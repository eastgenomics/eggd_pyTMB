# eggd_pyTMB (DNAnexus Platform App)

## What does this app do?
Calculates Tumour Mutational Burden (TMB) from an annotated VCF using [pyTMB](https://github.com/bioinfo-pf-curie/TMB).

The app first checks that the BED, BAM, and GTF inputs use a consistent chromosome naming convention (all `chr`-prefixed or all without prefix), failing early if they disagree. It then uses mosdepth to compute a sample-specific, coverage-based effective genome size from the input BAM and capture BED (subject to the `min_coverage` and `min_mapquality` thresholds), and runs pyTMB against the annotated VCF using the supplied VAF, population allele frequency, depth and alt-depth filters — with the variant-config YAML selected via `vcf_yaml` — to produce a TMB score, its 95% confidence interval, and filtering statistics.

If the effective genome size cannot be computed (or pyTMB produces no report data), the app writes a fallback TMB report populated with `NA` values instead of failing the job.

## What are the typical use cases for this app?
This app may be executed as a standalone app to calculate TMB for a tumour sample as part of a somatic variant calling workflow, where a VEP-annotated VCF and the corresponding BAM are already available.

## What are the inputs?
**Packages (DNAnexus assets)**

- `bedtools` (v2.30.0)
- `htslib` (v1.15.0)
- `mosdepth` (v0.3.3)
- `pyTMB` (v1.6.0)
- `samtools` (v1.19.2)

**Inputs (required)**:

- `vcf`: VEP-annotated VCF file (`.vcf`, `.vcf.gz`, `.bcf`, `.bcf.gz`).
- `bam`: BAM file used to compute the sample-specific, coverage-based effective genome size.
- `bam_bai`: BAM index file (`.bam.bai` / `.bai`) required by `mosdepth` for coverage computation.
- `bed`: Capture design BED file. Should be 0-based, sorted, and with no header.
- `gtf`: Genome annotation GTF file (`.gtf` / `.gtf.gz`).
- `vcf_yaml` (`string`): VCF YAML configuration, selecting the variant-config file used by pyTMB. Options are: `tso500` or `tnhaplotyper2`. Default: `tnhaplotyper2`.

**Other inputs (optional)**:

`vaf` (`float`): Minimum variant allele frequency. Variants with VAF below this threshold are excluded. Default: `0.05`.

`maf` (`float`): Maximum population allele frequency. Variants with population MAF above this threshold are excluded as germline. Default: `0.0001`.

`min_depth` (`int`): Minimum total read depth. Variants with total depth below this threshold are excluded. Default: `50`.

`min_alt_depth` (`int`): Minimum alt allele depth. Variants with fewer than this many alt reads are excluded. Default: `2`.

`min_coverage` (`int`): Minimum coverage per region of the BED file, used when computing the effective genome size. Default: `50`.

`min_mapquality` (`int`): Minimum mapping quality threshold; reads below this value are ignored when computing the effective genome size. Default: `50`.

## What are the outputs?

This app outputs a single plain-text TMB report (`<SAMPLE>_TMB.txt`) produced by `pyTMB`, containing the TMB score, its 95% confidence interval, and filtering statistics. `<SAMPLE>` is derived from the input VCF filename (the portion before the first underscore).

If the effective genome size could not be determined, or pyTMB produced no report data, the report is instead populated with `NA` values for the score, confidence interval, and filtering statistics.

## How to run this app from command line?
**Example**:
```
dx run eggd_pyTMB \
  -ivcf="sample.vep_annotated.vcf.gz" \
  -ibam="sample.bam" \
  -ibam_bai="sample.bam.bai" \
  -ibed="capture_design.bed" \
  -igtf="annotation.gtf" \
  -ivcf_yaml="tnhaplotyper2" \
  -ivaf=0.05 \
  -imaf=0.0001 \
  -imin_depth=50 \
  -imin_alt_depth=2 \
  -imin_coverage=50 \
  -imin_mapquality=50
```

The above will do the following:
- verify that the BED, BAM, and GTF files use a consistent chromosome naming convention
- compute a coverage-based effective genome size from `sample.bam` restricted to `capture_design.bed`, using the `min_coverage` and `min_mapquality` thresholds
- filter the annotated VCF (using the `tnhaplotyper2` variant config) to variants with VAF ≥ 0.05, population AF ≤ 0.0001, total depth ≥ 50 and alt depth ≥ 2
- calculate the TMB score and 95% CI over the resulting variant set

This is the source code for an app that runs on the DNAnexus Platform.
For more information about how to run or modify it, see https://documentation.dnanexus.com/.
