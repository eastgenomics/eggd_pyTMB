# eggd_pyTMB (DNAnexus Platform App)
This could be used as a starting point when developing new apps for DNAnexus

<!-- Insert a description of your app here -->
## What does this app do?
Calculates Tumour Mutational Burden (TMB) from an annotated VCF using [pyTMB](https://github.com/bioinfo-pf-curie/TMB).

The app uses mosdepth to compute a sample-specific, coverage-based effective genome size from the input BAM and capture BED, then runs pyTMB against the annotated VCF (restricted to the capture BED region) using the supplied VAF, population allele frequency, depth and alt-depth filters to produce a TMB score, its 95% confidence interval, and filtering statistics.

## What are the typical use cases for this app?

This app may be executed as a standalone app to calculate TMB for a tumour sample as part of a somatic variant calling workflow, where a VEP-annotated VCF and the corresponding BAM are already available.

## What are the inputs?

**Packages (DNAnexus assets)**

- `bedtools` (v2.30.0)
- `htslib` (v1.15.0)
- `mosdepth` (v0.3.3)
- `pyTMB` (v1.6.0)

**File inputs (required)**:

- `vcf`: VEP-annotated VCF file (`.vcf`, `.vcf.gz`, `.bcf`, `.bcf.gz`).
- `bam`: BAM file used to compute the sample-specific, coverage-based effective genome size.
- `bam_bai`: BAM index file (`.bam.bai` / `.bai`) required by `mosdepth` for coverage computation.

**Other inputs (optional)**:

`bed` (`file`): Capture design BED file.

`gtf` (`file`): Genome annotation GTF file (`.gtf` / `.gtf.gz`).

`vaf` (`float`): Minimum variant allele frequency. Variants with VAF below this threshold are excluded. Default: `0.05`.

`maf` (`float`): Maximum population allele frequency. Variants with population MAF above this threshold are excluded as germline. Default: `0.0001`.

`min_depth` (`int`): Minimum total read depth. Variants with total depth below this threshold are excluded. Default: `50`.

`min_alt_depth` (`int`): Minimum alt allele depth. Variants with fewer than this many alt reads are excluded. Default: `2`.

## What are the outputs?

This app outputs a single plain-text TMB report (`*.txt`) produced by `pyTMB`, containing the TMB score, its 95% confidence interval, and filtering statistics.

## How to run this app from command line?

**Example**:
```
dx run eggd_pyTMB \
  -ivcf="sample.vep_annotated.vcf.gz" \
  -ibam="sample.bam" \
  -ibam_bai="sample.bam.bai" \
  -ibed="capture_design.bed" \
  -ivaf=0.05 \
  -imaf=0.0001 \
  -imin_depth=50 \
  -imin_alt_depth=2
```

The above will do the following:
- compute a coverage-based effective genome size from `sample.bam` restricted to `capture_design.bed`
- filter the annotated VCF to variants with VAF ≥ 0.05, population AF ≤ 0.0001, total depth ≥ 50 and alt depth ≥ 2
- calculate the TMB score and 95% CI over the resulting variant set

This is the source code for an app that runs on the DNAnexus Platform.
For more information about how to run or modify it, see https://documentation.dnanexus.com/.
