# nf-meta-gwas

A Nextflow pipeline for multi-population GWAS summary statistics processing, meta-analysis, and downstream functional interpretation including fine-mapping, gene prioritisation, heritability estimation, and tissue/cell type enrichment analysis.

## Introduction

**nf-meta-gwas** takes raw GWAS summary statistics from one or more cohorts and populations, harmonises them to a common reference genome (GRCh38), performs LD-score regression (LDSC) intercept correction for genomic inflation, and runs fixed-effects inverse-variance weighted (IVW) meta-analysis within and across populations.

Downstream modules cover the full post-GWAS functional interpretation workflow:

- **Visualisation** — Manhattan and QQ plots for every input GWAS and meta-analysis result, with ggplot-style styling, per-chromosome significance colouring, and nearest-gene annotation for genome-wide significant loci
- **Lead variant extraction** — genome-wide significant loci identified and aggregated across populations
- **ABF fine-mapping** — Approximate Bayes Factor 99% credible sets per locus
- **Gene prioritisation** — MAGMA gene-level analysis, MAGMA tissue expression analysis (GTEx v8), PoPS polygenic priority scores, and FLAMES integrative gene scoring
- **Heritability estimation** — SNP heritability on meta-analysis summary statistics
- **Tissue/cell type enrichment** — LDSC-CTS partitioned heritability

The pipeline is written in [Nextflow](https://www.nextflow.io/) and uses [Apptainer](https://apptainer.org/) containers to ensure reproducibility. All compute-intensive steps are designed to run on an HPC cluster via the LSF scheduler.

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
└─────────────────────────────────────────────────────────────┘
        │
        ├──────────────────────────────────────────────────────────►  PLOT_INPUT_GWAS
        │                                                              Manhattan + QQ per input GWAS
        ▼
┌─────────────────────────────────────────────────────────────┐
│  LDSC_CORRECTION                                            │
│  • Filter to HapMap3 SNPs                                   │
│  • Estimate h² and LDSC intercept                           │
│  • Correct SE/Z/P if intercept > 1                          │
│  • Export per-cohort h² results                             │
└─────────────────────────────────────────────────────────────┘
        │
        ▼
┌─────────────────────────────────────────────────────────────┐
│  META_ANALYZE (per chromosome)                              │
│  • Fixed-effects IVW meta-analysis                          │
│  • Within-population and across-population                  │
│  • Cochran's Q heterogeneity statistics                     │
│  • MLOG10P computed natively from Z (handles underflow)     │
└─────────────────────────────────────────────────────────────┘
        │
        ├──────────────┬──────────────┬──────────────┬──────────────────────────────┐
        ▼              ▼              ▼              ▼                              ▼
┌──────────────┐ ┌──────────┐ ┌──────────────┐ ┌─────────────────┐   ┌────────────────────┐
│ EXTRACT_LEAD │ │HERITABI- │ │ TISSUE_      │ │ PLOT_META_GWAS  │   │ PREPARE_MAGMA_INPUT│
│ _VARIANTS    │ │LITY      │ │ ENRICHMENT   │ │ Manhattan + QQ  │   │ → MAGMA_GENE       │
│ • GWS loci   │ │ • LDSC   │ │ • LDSC-CTS   │ │ per meta result │   │   ├─ MAGMA_TISSUE  │
│ • nearest    │ │   h² on  │ │ • tissue and │ │                 │   │   │  (GTEx v8)      │
│   gene       │ │   meta   │ │   cell type  │ └─────────────────┘   │   └─ POPS          │
└──────┬───────┘ └──────────┘ └──────────────┘                       └────────┬───────────┘
       │                                                                        │
       ▼                                                                        │
┌─────────────────────────────────────────────────────────────┐                │
│  ABF_FINEMAPPING  (GWS loci only)                           │                │
│  • 99% ABF credible sets per locus                          │                │
│  • SNPID format: CHR:POS:NEA:EA                             │                │
└──────────────────────────┬──────────────────────────────────┘                │
                           │                                                    │
                           ▼                                                    │
┌─────────────────────────────────────────────────────────────┐                │
│  PREPARE_FLAMES_INPUTS                                      │                │
│  • Per-locus credset files (SNPID flipped to CHR:POS:EA:NEA)│                │
│  • Genomic locus file (±500 kb windows)                     │                │
│  • FLAMES index file (relative paths)                       │                │
└──────────────────────────┬──────────────────────────────────┘                │
                           │◄──────────────────────────────────────────────────┘
                           │  (join: MAGMA genes.out, MAGMA gsa.out, PoPS preds)
                           ▼
┌─────────────────────────────────────────────────────────────┐
│  FLAMES                                                     │
│  • Annotate loci with SNP-to-gene evidence, MAGMA-Z, PoPS  │
│  • Score loci with XGBoost + PoPS model                     │
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

The following reference files are configured in `nextflow.config`. Voltron HPC paths are pre-filled where available. MAGMA, PoPS, and FLAMES params default to `null` — these analyses are **opt-in**: set a param to a non-null path to enable it, leave it `null` to skip. Harmonisation, meta-analysis, and fine-mapping reference files are required for all runs.

| Parameter | Description |
|-----------|-------------|
| `params.fasta` | GRCh38 reference genome FASTA |
| `params.dbsnp` | dbSNP VCF for rsID assignment (GCF_000001405.40) |
| `params.chain` | LiftOver chain file (hg19 → hg38) |
| `params.gtf` | Ensembl GTF for gene annotation (Manhattan plot labelling) |
| `params.population_vcf` | Per-population 1000 Genomes VCF (hg38, 30×) for allele frequency inference |
| `params.ldsc_reference` | Per-population Pan-UKBB LD score reference panels |
| `params.ldcts` | Cell-type LD score file for LDSC-CTS enrichment analysis |
| `params.ldcts_weights` | LD score weights for LDSC-CTS (HapMap3, no MHC) |
| `params.ldbaseline` | Baseline LD score annotations for LDSC-CTS |
| `params.magma_bfile` | MAGMA 1000G LD reference bfile prefix (`.bed/.bim/.fam`) |
| `params.magma_gene_annot` | MAGMA pre-built gene annotation file (GRCh38 Ensembl) |
| `params.magma_gtex_covar` | GTEx v8 tissue-average expression matrix for MAGMA tissue analysis |
| `params.pops_gene_annot` | PoPS gene annotation file |
| `params.pops_feature_prefix` | PoPS munged feature matrix prefix |
| `params.pops_num_feature_chunks` | Number of PoPS feature chunks (default: `116`) |
| `params.pops_control_features` | PoPS control features file |
| `params.flames_annotation_data` | FLAMES annotation data directory |

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
| `--magma_bfile` | `null` | MAGMA LD reference bfile prefix |
| `--magma_gene_annot` | `null` | MAGMA gene annotation file |
| `--magma_gtex_covar` | `null` | GTEx v8 covariate matrix for MAGMA tissue analysis |
| `--pops_gene_annot` | `null` | PoPS gene annotation |
| `--pops_feature_prefix` | `null` | PoPS feature matrix prefix |
| `--pops_num_feature_chunks` | `116` | Number of PoPS feature chunks |
| `--pops_control_features` | `null` | PoPS control features file |
| `--flames_annotation_data` | `null` | FLAMES annotation data directory |

---

## Outputs

```
results/
├── plots/
│   ├── manhattan/
│   │   └── {phenotype}/
│   │       └── {prefix}.manhattan.{png,pdf}
│   └── qq/
│       └── {phenotype}/
│           └── {prefix}.qq.{png,pdf}
├── ldsc/
│   └── {phenotype}/
│       └── {phenotype}.ldsc_h2.txt
├── lead-variants/
│   └── {phenotype}/
│       └── {phenotype}-{population}.lead.variants.txt
├── heritability/
│   └── {phenotype}/
│       └── {phenotype}-{population}.ldsc_h2.txt
├── finemapping/
│   └── {phenotype}/
│       └── {phenotype}-{population}.credset.txt
├── enrichment/
│   └── {phenotype}/
│       └── {phenotype}-{population}.enrichment.txt
├── magma/
│   └── {phenotype}/
│       ├── {phenotype}-{population}.genes.out
│       ├── {phenotype}-{population}.genes.raw
│       └── {phenotype}-{population}-gtex.gsa.out
├── pops/
│   └── {phenotype}/
│       └── {phenotype}-{population}-pops.preds
└── flames/
    └── {phenotype}/
        ├── {phenotype}-{population}.genomic-locus.txt
        ├── {phenotype}-{population}.index.txt
        ├── flames-annotate/
        │   └── FLAMES_annotated_{locus}.txt
        └── flames-scores/

data/
└── meta-analysis/
    └── {phenotype}/
        └── {phenotype}-{population}.txt.gz
```

### Key output formats

**Manhattan plots** — ggplot-style (`theme_classic`), 8×5 inches, 300 dpi, PNG and PDF. Odd chromosomes: `#045ea7`; even chromosomes: `#82afd3`; chromosomes with ≥1 genome-wide significant variant: `#990000`. The most significant GWS variant per chromosome is annotated with the nearest Ensembl gene name. No annotations are added when no variants reach the genome-wide significance threshold (5×10⁻⁸). `-log₁₀(p)` is computed natively from BETA/SE using the log-scale survival function to handle p-value underflow.

**QQ plots** — ggplot-style, observed vs expected `-log₁₀(p)`. Points below `-log₁₀(p) = 1` are downsampled to 5% for plotting efficiency while preserving correct expected quantiles.

**Meta-analysis output** (`data/meta-analysis/`) columns: `SNPID`, `CHR`, `POS`, `rsID`, `EA`, `NEA`, `BETA`, `SE`, `P`, `Z`, `MLOG10P`, `EAF`, `N_CONTRIBUTIONS`, `N` (+ `N_CASE`, `N_CONTROL` for binary traits), `Q`, `Q_DF`, `Q_PVAL`, `I2`, `MAF`. The `MLOG10P` column is computed natively from the Z-score using `pnorm(|Z|, lower.tail=FALSE, log.p=TRUE)` to avoid floating-point underflow for highly significant loci.

**`lead.variants.txt`** — genome-wide significant loci annotated with nearest gene (gwaslab format)

**`credset.txt`** columns: `SNPID` (CHR:POS:NEA:EA), `rsID`, `LOCUS`, `CHR`, `POS`, `EA`, `NEA`, `BETA`, `SE`, `NEAREST_GENE`, `BF`, `BF_PIP`

**`{phenotype}-{population}-pops.preds`** — per-gene PoPS scores

**`FLAMES_annotated_{locus}.txt`** — per-locus SNP-to-gene annotation table incorporating MAGMA-Z, PoPS scores, and functional evidence

**`flames-scores/`** — final FLAMES XGBoost gene prioritisation scores per locus

> **Note:** Heritability results are only produced for per-population meta-analysis results. The across-population (`all`) meta uses EUR LD panels by default; update `params.ldsc_reference.all` in `nextflow.config` as appropriate for your study.

> **Note — opt-in analyses:** MAGMA, PoPS, and FLAMES are opt-in and will only run when the required reference file params are non-null (see [Parameters](#parameters)). Set `params.magma_bfile` (and the other MAGMA/PoPS params) to enable MAGMA + PoPS. Set `params.flames_annotation_data` in addition to enable FLAMES. Providing MAGMA params but not `flames_annotation_data` runs MAGMA + PoPS only. Harmonisation, meta-analysis, and fine-mapping always run — this pipeline is purpose-built for meta-analysis and these core steps cannot be skipped.

> **Note — GWS gating:** FLAMES additionally requires genome-wide significant loci; phenotype-population combinations without GWS hits are automatically skipped by the upstream ABF fine-mapping filter.

---

## Execution profiles

| Profile | Executor | Notes |
|---------|----------|-------|
| `standard` | local | For testing on a workstation |
| `voltron` | LSF (`voltron_normal` queue) | HPC execution with Apptainer |

---

## Containers

All containers are built automatically via GitHub Actions on push to `main` or when a versioned tag (`{tool}-v*`) is pushed.

| Container | Image | Description |
|-----------|-------|-------------|
| gwaslab | `ghcr.io/ampregnall/nf-meta-gwas/gwaslab:latest` | GWASLab, bcftools, tabix, plotnine — munging, harmonization, plotting |
| meta-analysis | `ghcr.io/ampregnall/nf-meta-gwas/meta-analysis:latest` | R with tidyverse, arrow, box — meta-analysis, fine-mapping, FLAMES input prep |
| magma | `ghcr.io/ampregnall/nf-meta-gwas/magma:latest` | MAGMA v1.10 static binary + pandas |
| pops | `ghcr.io/ampregnall/nf-meta-gwas/pops:latest` | PoPS pinned to FinucaneLab/pops@76eb86c |
| flames | `ghcr.io/ampregnall/nf-meta-gwas/flames:latest` | FLAMES pinned to ampregnall/FLAMES@2008304 |

---

## References

If you use nf-meta-gwas, please also cite the following tools:

- **GWASLab** — He, Y. et al. (2023). GWASLab: a Python package for processing and visualizing GWAS summary statistics. *Preprint*. https://doi.org/10.1101/2023.01.15.524141
- **LDSC** — Bulik-Sullivan, B. et al. (2015). LD Score regression distinguishes confounding from polygenicity in genome-wide association studies. *Nature Genetics*, 47, 291–295. https://doi.org/10.1038/ng.3211
- **MAGMA** — de Leeuw, C.A. et al. (2015). MAGMA: Generalized gene-set analysis of GWAS data. *PLOS Computational Biology*, 11, e1004219. https://doi.org/10.1371/journal.pcbi.1004219
- **PoPS** — Weeks, E.M. et al. (2023). Leveraging polygenic enrichments of gene features to predict genes underlying complex traits and diseases. *Nature Genetics*, 55, 1267–1276. https://doi.org/10.1038/s41588-023-01443-6
- **FLAMES** — Ulirsch, J.C. et al. (2021). Interrogation of human hematopoiesis at single-cell and single-variant resolution. *Nature Genetics*, 51, 683–693. https://doi.org/10.1038/s41588-019-0362-6
- **HyPrColoc** — Foley, C.N. et al. (2021). A fast and efficient colocalization algorithm for identifying shared genetic risk factors across multiple traits. *Nature Communications*, 12, 764. https://doi.org/10.1038/s41467-020-20885-8

---

## To-do

- [x] Munge and harmonize summary statistics
- [x] LDSC intercept correction
- [x] Fixed-effects IVW meta-analysis (within- and across-population)
- [x] Lead variant extraction
- [x] Manhattan and QQ plots for all input and meta-analysis GWAS
- [x] SNP heritability estimation on meta-analysis summary statistics
- [x] ABF credible set fine-mapping
- [x] Tissue and cell type enrichment analysis with LDSC-CTS
- [x] MAGMA gene-level and tissue expression analysis
- [x] PoPS polygenic priority scores
- [x] FLAMES integrative gene prioritisation
- [ ] eQTL and pQTL colocalization with HyPrColoc
- [ ] Nextflow pipeline tests
- [ ] Summary reports of all results
