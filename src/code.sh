#!/bin/bash
# eggd_pyTMB, version 1.0.0

# The following line causes bash to exit at any point if there is any error
# and to output each line as it is executed -- useful for debugging
set -exo pipefail

dx-download-all-inputs

mkdir -p /home/dnanexus/out/tmb_report

# Add bundled binaries (mosdepth) to PATH
export PATH="/home/dnanexus/bin:${PATH}"

# Install all wheels (dependencies + pytmb package)
sudo -H python3 -m pip install --no-index --no-deps /home/dnanexus/packages/*

# Derive sample name from VCF filename
SAMPLE="${vcf_prefix}"

# Place BAM index at the exact path mosdepth/htslib derives from the BAM path
BAI_PATH="${bam_path}.bai"
mv "${bam_bai_path}" "${BAI_PATH}"
if [[ ! -f "${BAI_PATH}" ]]; then
    echo "Error: BAM index not found at ${BAI_PATH}" >&2
    exit 1
fi

# Run pyEffGenomeSize and extract effective genome size
EFF_OUTPUT=$(pyEffGenomeSize \
  --bed "${bed_path}" \
  --gtf "${gtf_path}" \
  --filterNonCoding \
  --bam "${bam_path}" \
  --minCoverage 100 \
  --minMapq 50 \
  --oprefix "/home/dnanexus/out/${SAMPLE}" 2>&1) || true
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
pyTMB \
  -i "${vcf_path}" \
  --varConfig "/home/dnanexus/pytmb/config/tnhaplotyper.yml" \
  --dbConfig "/home/dnanexus/pytmb/config/vep.yml" \
  --effGenomeSize "${EFF_GENOME_SIZE}" \
  --vaf 0.05 --maf 0.0001 \
  --minDepth 50 --minAltDepth 2 \
  --filterLowQual --filterNonCoding --filterSyn --filterSplice \
  --filterPolym --polymDb 1k,gnomad \
  > "/home/dnanexus/out/tmb_report/${SAMPLE}_TMB.txt" 2>&1

if [[ ! -s "/home/dnanexus/out/tmb_report/${SAMPLE}_TMB.txt" ]]; then
    echo "Error: TMB report is empty or was not created for ${SAMPLE}" >&2
    exit 1
fi

dx-upload-all-outputs
