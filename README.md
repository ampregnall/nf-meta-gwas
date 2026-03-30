# nf-meta-gwas

A Nextflow pipeline for multi-population GWAS summary statistics processing, meta-analysis, and downstream functional interpretation including fine-mapping, heritability estimation, and tissue/cell type enrichment analysis.

## Introduction

**nf-meta-gwas** takes raw GWAS summary statistics from one or more cohorts and populations, harmonises them to a common reference genome (GRCh38), performs LD-score regression (LDSC) intercept correction for genomic inflation, and runs fixed-effects inverse-variance weighted (IVW) meta-analysis within and across populations. Downstream modules perform lead variant extraction, SNP heritability estimation, Approximate Bayes Factor (ABF) credible set fine-mapping, and tissue/cell type enrichment analysis via LDSC-CTS.

The pipeline is written in [Nextflow](https://www.nextflow.io/) and uses [Apptainer](https://apptainer.org/) (formerly Singularity) containers to ensure reproducibility. All compute-intensive steps are designed to run on an HPC cluster via the LSF scheduler, though a local execution profile is also provided.

### Citation

nf-meta-gwas is currently in active development and has not yet been formally published. If you use this pipeline in your research, please cite the GitHub repository:

> Pregnall, A. (2025). *nf-meta-gwas*. GitHub. https://github.com/ampregnall/nf-meta-gwas

---

## Pipeline overview

```
Raw summary statistics (per cohort × population)
        │
        ▼
┌─────────────────────────────────────────────────────────────┐
│  MUNGE_SUMSTATS                                             │
│  • Build inference and liftover to GRCh38 (if needed)      │
│  • Basic QC: remove duplicates, normalize alleles           │
│  • Harmonize against reference genome, assign rsIDs         │
│  • Allele frequency inference from population VCF           │
│  • MAC filter                                               │
│  • DAF plot (PDF + PNG)                                     │
└─────────────────────────────────────────────────────────────┘
        │
        ▼
┌─────────────────────────────────────────────────────────────┐
│  LDSC_CORRECTION                                            │
│  • Filter to HapMap3 SNPs                                   │
│  • Estimate h² and LDSC intercept                           │
│  • Correct SE/Z/P if intercept > 1                          │
│  • Manhattan + QQ plot (PDF + PNG)                          │
│  • Export per-cohort h² results                             │
└─────────────────────────────────────────────────────────────┘
        │
        ▼
┌─────────────────────────────────────────────────────────────┐
│  META_ANALYZE (per chromosome)                              │
│  • Fixed-effects IVW meta-analysis                          │
│  • Within-population and across-population                  │
│  • Cochran's Q heterogeneity statistics                     │
└─────────────────────────────────────────────────────────────┘
        │
        ├──────────────────────────┬─────────────────────────┐
        ▼                          ▼                         ▼
┌───────────────────┐  ┌───────────────────────┐  ┌─────────────────────┐
│ EXTRACT_LEAD_     │  │ HERITABILITY          │  │ TISSUE_ENRICHMENT   │
│ VARIANTS          │  │ • LDSC h² on meta     │  │ • LDSC-CTS          │
│ • genome-wide     │  │   sumstats            │  │   partitioned h²    │
│   significant     │  │ • per-population      │  │ • tissue and cell   │
│   loci            │  │   only                │  │   type enrichment   │
│ • Manhattan plot  │  └───────────────────────┘  └─────────────────────┘
└───────────────────┘
        │
        ▼
┌─────────────────────────────────────────────────────────────┐
│  ABF_FINEMAPPING                                            │
│  • Extract variants within 500 kb of each lead locus        │
│  • Approximate Bayes Factor credible sets (99%)             │
└─────────────────────────────────────────────────────────────┘
```

---

## Requirements

| Tool | Minimum version |
|------|----------------|
| [Nextflow](https://www.nextflow.io/docs/latest/install.html) | 24.04 |
| [Apptainer](https://apptainer.org/docs/user/latest/quick_start.html) | 1.3 |

Container images are pulled automatically from the GitHub Container Registry (`ghcr.io/ampregnall/nf-meta-gwas`) and cached locally. No manual installation of R or Python packages is required.

---

## Quick start

```bash
# Clone the repository
git clone https://github.com/ampregnall/nf-meta-gwas.git
cd nf-meta-gwas

# Run on a local machine (small test)
nextflow run main.nf \
  --input assets/test_samplesheet.csv \
  -profile standard

# Run on the Voltron HPC cluster
nextflow run main.nf \
  --input assets/samplesheet.csv \
  -profile voltron
```

---

## Inputs

### Samplesheet

Pass the samplesheet path with `--input`. The file must be a comma-separated CSV with a header row and the following columns:

| Column | Type | Description |
|--------|------|-------------|
| `phenotype` | string | Phenotype identifier (e.g. `inguinal_hernia`) |
| `trait_type` | string | `binary` or `continuous` |
| `population` | string | Population code — one of `eur`, `afr`, `amr`, `eas`, `csa`, `mid` |
| `cohort` | string | Cohort identifier (e.g. `ukbb`) |

Example:

```csv
phenotype,trait_type,population,cohort
inguinal_hernia,binary,eur,ukbb
inguinal_hernia,binary,afr,ukbb
```

### Raw summary statistics

Raw sumstats files must be placed at:

```
data/raw/{phenotype}/{cohort}-{phenotype}-{population}.sumstats.txt.gz
```

Files are expected in the **gwaslab** input format (tab-delimited, gzip-compressed). Required columns depend on trait type:

| Column | Description |
|--------|-------------|
| `CHR` | Chromosome |
| `POS` | Base-pair position (any build — liftover to GRCh38 is performed automatically if hg19 is detected) |
| `RSID` | rsID (optional — assigned from dbSNP155 during harmonization if absent) |
| `EA` | Effect allele |
| `NEA` | Non-effect allele |
| `EAF` | Effect allele frequency |
| `BETA` | Effect size estimate |
| `SE` | Standard error |
| `P` | P-value |
| `N` | Total sample size |
| `N_CASE` | Case count (binary traits only) |
| `N_CONTROL` | Control count (binary traits only) |

### Reference files

The following reference files are required and configured in `nextflow.config`. Default paths are set for the Voltron HPC environment:

| Parameter | Description |
|-----------|-------------|
| `params.fasta` | GRCh38 reference genome FASTA |
| `params.dbsnp` | dbSNP VCF for rsID assignment (GCF_000001405.40) |
| `params.chain` | LiftOver chain file (hg19 → hg38) |
| `params.gtf` | Ensembl GTF for gene annotation |
| `params.population_vcf` | Per-population 1000 Genomes VCF (hg38, 30×) for allele frequency inference |
| `params.ldsc_reference` | Per-population Pan-UKBB LD score reference panels (used for LDSC correction, heritability, and enrichment) |
| `params.ldcts` | Cell-type LD score file for LDSC-CTS enrichment analysis |
| `params.ldcts_weights` | LD score weights for LDSC-CTS (typically HapMap3, no MHC) |
| `params.ldbaseline` | Baseline LD score annotations for LDSC-CTS |

---

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--input` | `null` | Path to samplesheet CSV (required) |
| `--output` | `results` | Output directory |
| `--mac` | `50` | Minimum minor allele count filter |
| `--fasta` | see config | GRCh38 reference FASTA |
| `--dbsnp` | see config | dbSNP VCF |
| `--chain` | see config | LiftOver chain file |
| `--gtf` | see config | Ensembl GTF file |
| `--ldcts` | see config | Cell-type LD scores for LDSC-CTS |
| `--ldcts_weights` | see config | LD score weights for LDSC-CTS |
| `--ldbaseline` | see config | Baseline LD score annotations for LDSC-CTS |

---

## Outputs

```
results/
├── daf/
│   └── {phenotype}/
│       └── {cohort}-{phenotype}-{population}-daf.{pdf,png}
├── manhattan-qq/
│   └── {phenotype}/
│       └── {cohort}-{phenotype}-{population}-manhattan-qq.{pdf,png}
├── ldsc/
│   └── {phenotype}/
│       └── {phenotype}.ldsc_h2.txt
├── lead-variants/
│   └── {phenotype}/
│       ├── {phenotype}-{population}-manhattan-qq.{pdf,png}
│       └── {phenotype}-{population}.lead.variants.txt
├── heritability/
│   └── {phenotype}/
│       └── {phenotype}-{population}.ldsc_h2.txt
├── finemapping/
│   └── {phenotype}/
│       └── {phenotype}-{population}.credset.txt
└── enrichment/
    └── {phenotype}/
        └── {phenotype}-{population}.enrichment.txt

data/
├── sumstats-processed/
│   └── {phenotype}/
│       └── {cohort}-{phenotype}-{population}.sumstats.processed.txt.gz
└── meta-analysis/
    └── {phenotype}/
        └── {phenotype}-{population}.txt.gz
```

**`ldsc_h2.txt`** columns: `phenotype`, `cohort`, `population`, `h2`, `h2_se`, `Intercept`, `Intercept_se`, `Lambda_GC`, `Mean_Chi2`, `Ratio`, `Ratio_se`

**Meta-analysis output** columns: `SNPID`, `CHR`, `POS`, `rsID`, `EA`, `NEA`, `BETA`, `SE`, `P`, `Z`, `EAF`, `N_CONTRIBUTIONS`, `N` (+ `N_CASE`, `N_CONTROL` for binary traits), `Q`, `Q_DF`, `Q_PVAL`, `I2`

**`lead.variants.txt`**: genome-wide significant loci annotated with nearest gene (gwaslab format)

**`credset.txt`** columns: `SNPID`, `rsID`, `LOCUS`, `CHR`, `POS`, `EA`, `NEA`, `BETA`, `SE`, `NEAREST_GENE`, `BF`, `BF_PIP`

**`enrichment.txt`**: LDSC-CTS output with per-tissue/cell-type enrichment statistics

> **Note:** Heritability results (`results/heritability/`) are only produced for per-population meta-analysis results. The across-population (`all`) meta uses EUR LD panels by default; update `params.ldsc_reference.all` in `nextflow.config` as appropriate for your study.

---

## Execution profiles

| Profile | Executor | Notes |
|---------|----------|-------|
| `standard` | local | For testing on a workstation |
| `voltron` | LSF (`voltron_normal` queue) | HPC execution with Apptainer |

---

## References

If you use nf-meta-gwas, please also cite the following tools:

- **GWASLab** — He, Y. et al. (2023). GWASLab: a Python package for processing and visualizing GWAS summary statistics. *Preprint*. https://doi.org/10.1101/2023.01.15.524141
- **LDSC** — Bulik-Sullivan, B. et al. (2015). LD Score regression distinguishes confounding from polygenicity in genome-wide association studies. *Nature Genetics*, 47, 291–295. https://doi.org/10.1038/ng.3211
- **HyPrColoc** — Foley, C.N. et al. (2021). A fast and efficient colocalization algorithm for identifying shared genetic risk factors across multiple traits. *Nature Communications*, 12, 764. https://doi.org/10.1038/s41467-020-20885-8

---

## To-do

- [x] Munge and harmonize summary statistics
- [x] LDSC intercept correction
- [x] Fixed-effects IVW meta-analysis (within- and across-population)
- [x] Lead variant extraction and meta-analysis Manhattan plots
- [x] SNP heritability estimation on meta-analysis summary statistics
- [x] ABF credible set fine-mapping
- [x] Tissue and cell type enrichment analysis with LDSC-CTS
- [ ] eQTL and pQTL colocalization with HyPrColoc
- [ ] Nextflow pipeline tests
- [ ] Summary reports of all results
