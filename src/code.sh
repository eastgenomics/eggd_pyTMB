#!/bin/bash
# eggd_pyTMB, version 1.0.0

# The following line causes bash to exit at any point if there is any error
# and to output each line as it is executed -- useful for debugging
set -exo pipefail

dx-download-all-inputs --parallel

mkdir -p /home/dnanexus/out/tmb_report

# Derive sample name by stripping pipeline suffix from VCF filename
SAMPLE="${vcf_prefix/_markdup_recalibrated_tnhaplotyper2_annotated/}"

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
  --minCoverage 100 \
  --minMapq 50 \
  --oprefix "/home/dnanexus/out/${SAMPLE}" 2>&1); then
    echo "${EFF_OUTPUT}"
    echo "Error: pyEffGenomeSize failed — see output above" >&2
    exit 1
fi
echo "${EFF_OUTPUT}"
EFF_GENOME_SIZE=$(echo "${EFF_OUTPUT}" | grep "Effective Genome Size" | awk '{print $(NF-1)}')

if [[ -z "${EFF_GENOME_SIZE}" ]]; then
    echo "Error: pyEffGenomeSize did not return an effective genome size. Check logs above." >&2
    exit 1
fi

if [[ "${EFF_GENOME_SIZE}" -eq 0 ]]; then
    echo "Error: Effective genome size is 0 bp" >&2
    exit 1
fi

echo "effGenomeSize for ${SAMPLE}: ${EFF_GENOME_SIZE}"

# Run pyTMB and capture report
python3 -m pytmb.cli.run_tmb \
  -i "${vcf_path}" \
  --varConfig "/home/dnanexus/pytmb/config/tnhaplotyper.yml" \
  --dbConfig "/home/dnanexus/pytmb/config/vep.yml" \
  --effGenomeSize "${EFF_GENOME_SIZE}" \
  --vaf "${vaf}" --maf "${maf}" \
  --minDepth "${min_depth}" --minAltDepth "${min_alt_depth}" \
  --filterLowQual --filterNonCoding --filterSyn --filterSplice \
  --filterPolym --polymDb 1k,gnomad \
  > "/home/dnanexus/out/tmb_report/${SAMPLE}_TMB.txt" 2>&1

if [[ ! -s "/home/dnanexus/out/tmb_report/${SAMPLE}_TMB.txt" ]]; then
    echo "Error: TMB report is empty or was not created for ${SAMPLE}" >&2
    exit 1
fi

dx-upload-all-outputs
