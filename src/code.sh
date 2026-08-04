#!/bin/bash
# eggd_pyTMB, version 1.0.0

# The following line causes bash to exit at any point if there is any error
# and to output each line as it is executed -- useful for debugging
set -exo pipefail


# Writes a TMB report populated with NA values, used as a fallback when
# effective genome size or pyTMB itself fails to produce a valid result.
# Usage: _create_na_tmb_report <output_file> <input_vcf> <vaf> <maf> <min_depth> <min_alt_depth> [sample_name]
_create_na_tmb_report() {
    local output_file="$1"
    local input_vcf="$2"
    local vaf_val="$3"
    local maf_val="$4"
    local min_depth_val="$5"
    local min_alt_depth_val="$6"
    local sample_name="${7:-None}"
 
    if [[ -z "$output_file" || -z "$input_vcf" || -z "$vaf_val" || -z "$maf_val" || -z "$min_depth_val" || -z "$min_alt_depth_val" ]]; then
        echo "Usage: _create_na_tmb_report <output_file> <input_vcf> <vaf> <maf> <min_depth> <min_alt_depth> [sample_name]" >&2
        return 1
    fi
 
    cat > "$output_file" <<EOF
pyTMB version= 1.6.0
When= $(date +%Y-%m-%d)
 
Input= ${input_vcf}
Sample= ${sample_name}
 
Config caller= /home/dnanexus/pytmb/config/vcf.yml
Config databases= /home/dnanexus/pytmb/config/vep.yml
 
Filters:
-------
VAF= ${vaf_val}
MAF= ${maf_val}
minDepth= ${min_depth_val}
minAltDepth= ${min_alt_depth_val}
filterLowQual= True
filterIndels= False
filterCoding= False
filterNonCoding= True
filterSplice= True
filterSyn= True
filterNonSyn= False
filterCancerHotspot= False
filterPolym= True

Filter statistics:
------------------
  QUAL= NA
  VAF= NA
  DEPTH= NA
  SPLICING= NA
  NONCODING= NA
  SYN= NA

Total number of variants= NA
Non-informative variants= NA
Variants after filters= NA
Effective Genome Size= NA

TMB= NA
TMB_95CI_low= NA
TMB_95CI_high= NA
EOF

    echo "Written: ${output_file}"
}

# check_chr_consistency <bed_path> <bam_path> <gtf_path>
# Verifies that BED, BAM, and GTF files all use the same chromosome
# naming convention (either all "chr"-prefixed or all without prefix).
# Returns 0 if consistent, 1 otherwise (with details printed to stderr).
check_chr_consistency() {
    local bed_path="$1"
    local bam_path="$2"
    local gtf_path="$3"

    local bed_chr gtf_chr bam_chr
    local bed_style gtf_style bam_style

    # Extract a representative chromosome name from each file
    bed_chr=$(head -n1 "${bed_path}" | cut -f1)
    gtf_chr=$(zgrep -v -E '^#' "${gtf_path}" | head -n1 | cut -f1)
    bam_chr=$(samtools view -H "${bam_path}" | awk -F'\t' '/^@SQ/{for(i=1;i<=NF;i++) if ($i ~ /^SN:/) {print substr($i,4); exit}}')

    if [[ -z "${bed_chr}" || -z "${gtf_chr}" || -z "${bam_chr}" ]]; then
        echo "Error: could not determine chromosome naming from one or more files" >&2
        echo "  bed_chr='${bed_chr}' gtf_chr='${gtf_chr}' bam_chr='${bam_chr}'" >&2
        return 1
    fi

    # Classify each as "chr" or "nochr" style
    [[ "${bed_chr}" == chr* ]] && bed_style="chr" || bed_style="nochr"
    [[ "${gtf_chr}" == chr* ]] && gtf_style="chr" || gtf_style="nochr"
    [[ "${bam_chr}" == chr* ]] && bam_style="chr" || bam_style="nochr"

    if [[ "${bed_style}" != "${gtf_style}" || "${bed_style}" != "${bam_style}" ]]; then
        echo "Error: chromosome naming mismatch between BED, GTF, and BAM files" >&2
        echo "  BED (${bed_chr}) -> ${bed_style}" >&2
        echo "  GTF (${gtf_chr}) -> ${gtf_style}" >&2
        echo "  BAM (${bam_chr}) -> ${bam_style}" >&2
        return 1
    fi

    echo "Chromosome naming consistent across BED, GTF, and BAM (${bed_style})"
    return 0
}

main() {

    dx-download-all-inputs --parallel

    mkdir -p /home/dnanexus/out/tmb_report
    # Check chromosome naming consistency between BED, BAM, and GTF files.
    # Fail early if inconsistent to avoid downstream errors in pyEffGenomeSize and pyTMB.
    if ! check_chr_consistency "${bed_path}" "${bam_path}" "${gtf_path}"; then
        echo "Error: chromosome naming check failed — see output above" >&2
        exit 1
    fi

    # Derive sample name by stripping pipeline suffix from VCF filename
    # The VCF filename is stripped after the first underscore, which is assumed to be the sample name.
    # e.g if the VCF prefix is "AA-BB-CC_tnhaplotyper2", then SAMPLE="AA-BB-CC".
    SAMPLE="${vcf_prefix%%_*}"

    # Place BAM index at the exact path mosdepth/htslib derives from the BAM path
    BAI_PATH="${bam_path}.bai"
    mv "${bam_bai_path}" "${BAI_PATH}"
    if [[ ! -f "${BAI_PATH}" ]]; then
        echo "Error: BAM index not found at ${BAI_PATH}" >&2
        exit 1
    fi

    # Run pyEffGenomeSize and extract effective genome size
    if ! EFF_OUTPUT=$(python3 -m pytmb.cli.run_effgenomesize \
    --bed "${bed_path}" \
    --gtf "${gtf_path}" \
    --filterNonCoding \
    --bam "${bam_path}" \
    --minCoverage "${min_coverage}" \
    --minMapq "${min_mapquality}" \
    --oprefix "/home/dnanexus/out/${SAMPLE}" 2>&1); then
        echo "${EFF_OUTPUT}"
        echo "Error: pyEffGenomeSize failed — see output above" >&2
        exit 1
    fi
    echo "${EFF_OUTPUT}"
    EFF_GENOME_SIZE=$(echo "${EFF_OUTPUT}" | grep "Effective Genome Size" | awk '{print $(NF-1)}')

    # Extract effective genome size; fall back to a sentinel if the pipeline
    # finds no match, so set -e/pipefail doesn't kill the script here
    EFF_GENOME_SIZE=$(echo "${EFF_OUTPUT}" | grep "Effective Genome Size" | awk '{print $(NF-1)}') || EFF_GENOME_SIZE="NO_EFF_GENOME_SIZE"

    # Single check: catches the no-match sentinel, any non-numeric value, AND zero
    if [[ ! "${EFF_GENOME_SIZE}" =~ ^[1-9][0-9]*$ ]]; then
        echo "Error: Effective genome size ('${EFF_GENOME_SIZE}') is not a valid positive integer." >&2
        _create_na_tmb_report "/home/dnanexus/out/tmb_report/${SAMPLE}_TMB.txt" "${vcf_path}" "${vaf}" "${maf}" "${min_depth}" "${min_alt_depth}" "${SAMPLE}"
    else
        echo "Effective genome size is ${EFF_GENOME_SIZE} bp"  
        # run pyTMB, capture output to the report file, and capture exit status separately
        if python3 -m pytmb.cli.run_tmb \
            -i "${vcf_path}" \
            --varConfig "/home/dnanexus/pytmb/config/vcf.yml" \
            --dbConfig "/home/dnanexus/pytmb/config/vep.yml" \
            --effGenomeSize "${EFF_GENOME_SIZE}" \
            --vaf "${vaf}" --maf "${maf}" \
            --minDepth "${min_depth}" --minAltDepth "${min_alt_depth}" \
            --filterLowQual --filterNonCoding --filterSyn --filterSplice \
            --filterPolym --polymDb 1k,gnomad \
            > "/home/dnanexus/out/tmb_report/${SAMPLE}_TMB.txt" 2> "/home/dnanexus/out/tmb_report/${SAMPLE}_TMB.stderr.log"; then
            echo "pyTMB completed successfully"
        else
            echo "Error: run_tmb failed (exit code $?) — see stderr log" >&2
            cat "/home/dnanexus/out/tmb_report/${SAMPLE}_TMB.stderr.log" >&2
            _create_na_tmb_report "/home/dnanexus/out/tmb_report/${SAMPLE}_TMB.txt" "${vcf_path}" "${vaf}" "${maf}" "${min_depth}" "${min_alt_depth}" "${SAMPLE}"  
        fi
    fi

    dx-upload-all-outputs

}
