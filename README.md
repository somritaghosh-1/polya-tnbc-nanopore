# PolyA Tail Length Landscape in Triple-Negative Breast Cancer
 
<p align="center">
  <img src="https://img.shields.io/badge/Platform-Oxford%20Nanopore%20MinION-0084A9?style=flat-square" alt="Platform">
  <img src="https://img.shields.io/badge/Kit-SQK--RNA004-2E75B6?style=flat-square" alt="Kit">
  <img src="https://img.shields.io/badge/Reads-247%2C258-1F3864?style=flat-square" alt="Reads">
  <img src="https://img.shields.io/badge/Cancer%20Genes-96-C00000?style=flat-square" alt="Genes">
  <img src="https://img.shields.io/badge/Language-R-276DC3?style=flat-square&logo=r" alt="R">
  <img src="https://img.shields.io/badge/License-MIT-green?style=flat-square" alt="License">
</p>
 
<p align="center">
  <strong>Direct RNA Nanopore Sequencing | MDA-231 TNBC | 96 Cancer Genes</strong><br>
  <em>Somrita Ghosh — MSc Molecular Biology and Biotechnology, Queen's University Belfast (2024)</em><br>
  Patrick G Johnston Centre for Cancer Research | Institute for Global Food Security
</p>
 
<br>
 
> **Part 1 of 2** — This repository is the first of two complementary epitranscriptomic analyses. Part 2 covering the m6A modification landscape is available at [m6a-epitranscriptomics-tnbc](https://github.com/somritaghosh-1/m6a-epitranscriptomics-tnbc). Together they form the first combined polyA tail length and m6A modification profile of a clinically validated cancer gene panel in triple-negative breast cancer, generated from a single Oxford Nanopore direct RNA sequencing experiment.
 
<br>
 
## Table of Contents
 
- [Overview](#overview)
- [Key Findings](#key-findings)
- [Figures](#figures)
- [Repository Structure](#repository-structure)
- [Methods](#methods)
- [How to Reproduce](#how-to-reproduce)
- [Data Description](#data-description)
- [Biological Context](#biological-context)
- [Limitations](#limitations)
- [Future Directions](#future-directions)
- [Citation](#citation)
- [Contact](#contact)
 
<br>
 
## Overview
 
This repository contains the data, analysis code, and results from my MSc research project characterising polyA tail length regulation across 96 cancer-associated genes in MDA-231 triple-negative breast cancer (TNBC) cells using Oxford Nanopore direct RNA sequencing.
 
The polyA tail is a stretch of adenosine residues at the 3' end of eukaryotic mRNAs and is a key determinant of mRNA stability, translational efficiency, and degradation rate. Until third-generation direct RNA sequencing, genome-wide polyA tail profiling required laborious indirect methods. The Oxford Nanopore SQK-RNA004 kit enables simultaneous measurement of polyA tail length, RNA base modifications (m6A, pseudouridine), and transcript abundance from a single native RNA sequencing experiment without reverse transcription or amplification bias.
 
**Why TNBC:** Triple-negative breast cancer lacks oestrogen receptor, progesterone receptor, and HER2 amplification, making it unresponsive to hormone therapies and HER2-targeted drugs. It relies heavily on chemotherapy and carries poor prognosis. Understanding post-transcriptional gene regulation in TNBC may reveal new layers of therapeutic vulnerability not captured by standard RNA-seq approaches.
 
<br>
 
## Key Findings
 
| Finding | Detail |
|---|---|
| **Total reads** | 247,258 individual polyA tail length measurements |
| **Genes detected** | 88 of 96 (91.7%) of the cancer gene panel |
| **PolyA range** | 83.5 nt (RPL37A) to 167.2 nt (ITGAV), a 2-fold range |
| **Statistical test** | Kruskal-Wallis H = 22,984, p < 0.001 across 59 high-confidence genes |
| **Integrin stabilisation** | Five of six detected integrin subunits display above-average polyA tails |
| **TGF-beta1 paradox** | TGFB1 = 90.8 nt despite being a key EMT and immunosuppression driver |
| **Oncogene divergence** | AKT1 (89.8 nt) and CDK4 (90.2 nt) short-tailed; CDK2 long-tailed (127.4 nt) |
| **Technical validation** | YWHAZ sits precisely at dataset mean (109.8 nt) across 16,458 reads |
 
### Biological Highlights
 
**Integrin mRNA stabilisation:** Five of six detected integrin subunits (ITGB1 at 143.0 nt, ITGA1 at 134.7 nt, ITGA2 at 134.7 nt, ITGA3 at 122.0 nt, and ITGAV at 167.2 nt) display above-average polyA tail lengths, consistent with post-transcriptional mRNA stabilisation driving constitutive integrin expression and MDA-231 invasiveness. ITGB1 (mean 143.0 nt, n = 20,214 reads) is the most stably expressed gene in the dataset. ITGB5 (106.6 nt) sits at the dataset mean and does not show the same stabilisation pattern as the alpha subunits.
 
**TGF-beta1 short polyA paradox:** TGFB1 (mean 90.8 nt) displays one of the shortest polyA tails in the dataset despite being a key driver of epithelial-mesenchymal transition and immunosuppression in TNBC. This implies that TGF-beta1 activity is sustained by transcriptional upregulation or protein stability rather than mRNA stabilisation, a distinction with direct therapeutic relevance.
 
**Oncogene divergence:** AKT1 (89.8 nt) and CDK4 (90.2 nt) show below-average polyA tails despite being oncogenically active, suggesting their high protein output is achieved through transcriptional rather than post-transcriptional mechanisms. CDK2, by contrast, shows a long tail (127.4 nt), revealing differential post-transcriptional regulation within the oncogenic kinase landscape.
 
<br>
 
## Figures
 
### Figure 1 — Overall PolyA Tail Length Distribution
 
![Overall Distribution](https://raw.githubusercontent.com/somritaghosh-1/polya-tnbc-nanopore/main/Fig1_Overall_PolyA_Distribution.png)
 
*Distribution of 247,258 individual polyA tail length measurements across 88 cancer genes in MDA-231 TNBC cells. The distribution is unimodal and right-skewed. Red dashed line = dataset mean (108.3 nt); orange dotted line = dataset median (105.0 nt). The extended right tail reflects translationally active, stabilised transcripts.*
 
<br>
 
### Figure 2 — Mean PolyA Tail Length per Cancer Gene
 
![Mean PolyA Per Gene](https://raw.githubusercontent.com/somritaghosh-1/polya-tnbc-nanopore/main/Fig2_Mean_PolyA_Per_Gene.png)
 
*Mean polyA tail length for all genes with 100 or more reads (n = 72 genes), sorted by length. Blue bars indicate significantly long tails (>130 nt); red bars indicate significantly short tails (<95 nt); grey bars are intermediate. Error bars represent 95% CI. Read counts are shown for each gene. Dashed line = dataset mean (108.3 nt).*
 
<br>
 
### Figure 3 — PolyA Tail Length by Functional Gene Category
 
![Functional Groups](https://raw.githubusercontent.com/somritaghosh-1/polya-tnbc-nanopore/main/Fig3_PolyA_By_Functional_Group.png)
 
*PolyA tail length distributions grouped by biological function. Box = IQR; whiskers = 10th to 90th percentile; diamond = mean; individual data points are overlaid (subsampled for clarity). Cell adhesion genes (integrins) display the longest tails; housekeeping genes display the shortest. Dashed line = dataset mean (108.3 nt).*
 
<br>
 
### Figure 4 — Integrin Family: Coordinated mRNA Stabilisation
 
![Integrin Family](https://raw.githubusercontent.com/somritaghosh-1/polya-tnbc-nanopore/main/Fig4_Integrin_Family_PolyA.png)
 
*Mean polyA tail length for all detected integrin subunits with 100 or more reads. Gradient blue bars indicate relative tail length intensity. Values show mean tail length and deviation from dataset mean. Error bars = 95% CI. Red dashed line = dataset mean (108.3 nt).*
 
<br>
 
### Figure 5 — ITGB1 vs RPL37A: Biological Contrast
 
![ITGB1 vs RPL37A](https://raw.githubusercontent.com/somritaghosh-1/polya-tnbc-nanopore/main/Fig5_ITGB1_vs_RPL37A_Contrast.png)
 
*Per-read polyA tail length distributions for ITGB1 (blue, n = 20,214) vs RPL37A (red, n = 16,374). Mean difference = 59.5 nt (Mann-Whitney U, p < 0.001). ITGB1, the master regulator of cell-matrix adhesion in MDA-231, displays constitutive mRNA stabilisation consistent with its role driving tumour invasiveness. RPL37A undergoes rapid turnover consistent with tight translational autoregulation.*
 
<br>
 
## Repository Structure
 
```
polya-tnbc-nanopore/
├── combined_pti_values.csv    # Per-read polyA tail lengths (96 gene columns, 20,214 rows, 247,258 non-NA values)
├── polyA_analysis.R           # Full R analysis script
├── Fig1_Overall_PolyA_Distribution.png
├── Fig2_Mean_PolyA_Per_Gene.png
├── Fig3_PolyA_By_Functional_Group.png
├── Fig4_Integrin_Family_PolyA.png
├── Fig5_ITGB1_vs_RPL37A_Contrast.png
├── LICENSE
└── README.md
```
 
<br>
 
## Methods
 
| Parameter | Value |
|---|---|
| Cell Line | MDA-231 (Triple-Negative Breast Cancer) |
| Sequencing Platform | Oxford Nanopore MinION |
| Library Kit | SQK-RNA004 (Direct RNA Sequencing) |
| Basecaller | Dorado (modification-aware, `--emit-moves`) |
| Modifications Also Detected | N6-methyladenosine (m6A), Pseudouridine |
| Gene Panel | TaqMan Array Micro Fluidic Card 96a (96 cancer genes) |
| Genes Detected | 88/96 (91.7%) |
| Total Reads Analysed | 247,258 |
| Statistical Tests | Kruskal-Wallis H-test, Mann-Whitney U, Pearson correlation |
| Analysis Language | R (ggplot2, dplyr, tidyr, ggpubr, RColorBrewer) |
 
**Cell culture:** MDA-MB-231 cells were maintained in DMEM supplemented with 10% FBS and 1% penicillin/streptomycin at 37°C with 5% CO2. Cells were confirmed mycoplasma-negative prior to RNA extraction.
 
**Library preparation:** Total RNA was extracted using TRIzol. Polyadenylated RNA was enriched prior to library preparation using the SQK-RNA004 Direct RNA Sequencing Kit. Direct sequencing preserves native RNA without cDNA synthesis or amplification bias, enabling simultaneous detection of polyA tail length, m6A, and pseudouridine from a single experiment.
 
**PolyA tail length estimation:** PolyA tail lengths were estimated directly from native RNA signal using Dorado with the `--emit-moves` flag, which activates the Dorado poly(A) estimator. This measures tail length from raw ionic current at single-molecule resolution without reverse transcription or amplification bias. Per-read estimates (in nucleotides) were extracted from the `pt` tag in the Dorado output BAM and assigned to genes by alignment to the Ensembl GRCh38 transcriptome.
 
**Statistical analysis:** Per-gene summary statistics including mean, median, SD, and 95% CI were calculated for all genes with 10 or more reads. The Kruskal-Wallis H-test was applied across genes with 500 or more reads. Pearson correlation assessed the relationship between sequencing depth and mean polyA length. Functional group comparisons used pooled per-read values across gene members within each category.
 
<br>
 
## How to Reproduce
 
**Requirements**
 
```
R >= 4.0
ggplot2
dplyr
tidyr
ggpubr
RColorBrewer
```
 
**Install packages**
 
```r
install.packages(c("ggplot2", "dplyr", "tidyr", "ggpubr", "RColorBrewer"))
```
 
**Run analysis**
 
```r
setwd("path/to/polya-tnbc-nanopore")
source("polyA_analysis.R")
```
 
This will generate all five figures and a summary statistics CSV file in your working directory.
 
<br>
 
## Data Description
 
| File | Description | Dimensions |
|---|---|---|
| `combined_pti_values.csv` | Per-read polyA tail length estimates in nucleotides for each of the 96 genes in the cancer panel. Each column represents one gene and each row represents one sequencing read. NA values indicate reads not assigned to that gene. | 20,214 rows × 96 columns |
 
> **Note on read count:** The 20,214 rows equal the read depth of ITGB1, the most highly covered gene. The total reads in the analysis (247,258) is the **sum of non-NA values across all 96 columns**, not 20,214 × 96. To verify in R: `sum(!is.na(combined_pti_values))`
 
**Eight genes produced no detectable signal:** CDKN2A, FGFR2, GZMA, IFNA1, IFNB1, IGF1, TEK, and TWIST1. This is consistent with absent or sub-threshold expression in MDA-231 cells. CDKN2A is frequently deleted or epigenetically silenced in this cell line.
 
**Coverage tiers:**
 
| Tier | Coverage | Genes (n) | Confidence |
|---|---|---|---|
| 1 | 1,000 reads or more | 47 | High |
| 2 | 100 to 999 reads | 25 | Moderate |
| 3 | 10 to 99 reads | 12 | Low |
| 4 | Fewer than 10 reads | 4 | Unreliable |
 
<br>
 
## Biological Context
 
**Why polyA tail length matters in cancer**
 
The polyA tail is added co-transcriptionally to the 3' end of eukaryotic mRNAs by the cleavage and polyadenylation machinery. Tail length is regulated dynamically. Long tails (greater than 130 nt) are associated with mRNA stabilisation, sustained translational activity, and high protein output. Short tails (less than 95 nt) are associated with rapid mRNA turnover and reduced protein synthesis. Deadenylation by the CCR4-NOT complex is the rate-limiting step in mRNA decay and is recruited by m6A reader YTHDF2, connecting this dataset directly to [Project 2](https://github.com/somritaghosh-1/m6a-epitranscriptomics-tnbc).
 
Until Oxford Nanopore direct RNA sequencing, measuring polyA tail length genome-wide required laborious indirect methods such as PAT assays and TAIL-seq. The SQK-RNA004 kit enables direct, single-molecule measurement from native RNA.
 
**Key biological findings**
 
Five of six detected integrin subunits (ITGB1 at 143.0 nt, ITGA1 at 134.7 nt, ITGA2 at 134.7 nt, ITGA3 at 122.0 nt, and ITGAV at 167.2 nt) exceed the dataset mean. Integrins are the primary mediators of MDA-231 adhesion and invasion through extracellular matrix interactions. Constitutive post-transcriptional stabilisation of these integrin transcripts adds a new regulatory dimension to understanding MDA-231 invasiveness. ITGB5 (106.6 nt) sits at the dataset mean and does not display the same stabilisation pattern.
 
TGFB1 (90.8 nt) drives EMT and immunosuppression but is not sustained by mRNA stabilisation. Therapeutic strategies targeting TGF-beta1 output should focus on transcriptional or protein-level mechanisms rather than mRNA stability.
 
CDK2 (127.4 nt) vs AKT1 (89.8 nt) and CDK4 (90.2 nt) reveals that not all oncogenes exploit post-transcriptional stabilisation, a distinction relevant to understanding how TNBC maintains proliferative signalling.
 
Reference gene YWHAZ (14-3-3 zeta) sits precisely at the dataset mean (109.8 nt) across 16,458 reads, providing robust internal validation of technical reproducibility.
 
<br>
 
## Limitations
 
- **Single cell line, single run, no biological replicates.** All data derive from one sequencing experiment on MDA-231 cells. Findings reflect a snapshot of this cell line under standard culture conditions and cannot be generalised to primary TNBC tumours without further validation.
- **No normal comparator.** There is no matched normal breast epithelial dataset. It is not possible to determine which polyA tail length patterns are TNBC-specific versus general features of cultured epithelial cells.
- **96-gene panel design.** The TaqMan Array Micro Fluidic Card 96a was designed for expression profiling, not epitranscriptomics. Panel composition may introduce ascertainment bias toward known cancer-associated genes.
- **Low-coverage genes.** Twelve genes have 10–99 reads and four have fewer than 10 reads. Mean polyA estimates for these genes are unreliable and should not be interpreted independently.
- **Functional interpretation is correlative.** PolyA tail length is a proxy for mRNA stability. Direct functional validation (e.g. mRNA half-life assays, polyA tail shortening experiments) is required to confirm causal relationships.
 
<br>
 
## Future Directions
 
- [ ] Integration of m6A and pseudouridine modification data from the same MinION sequencing run — see [Project 2](https://github.com/somritaghosh-1/m6a-epitranscriptomics-tnbc)
- [ ] Comparison with PEO1 and PEO4 ovarian cancer cell lines (preliminary data collected during MSc)
- [ ] Comparison with normal breast epithelial cells (MCF-10A) to identify TNBC-specific stabilisation patterns
- [ ] Biological replicates (≥3) to confirm gene-level polyA patterns
- [ ] TCGA pan-cancer immune correlation analysis (CD8+ T cell infiltration across 40 cancer types)
- [ ] Validation of polyA-stability relationships via PAT assays
- [ ] Functional validation of integrin mRNA stabilisation via CCR4-NOT perturbation
 
<br>
 
## Citation
 
If you use this data or analysis in your research, please cite:
 
```
Ghosh, S. (2024). PolyA Tail Length Landscape in Triple-Negative Breast Cancer.
MSc Research Project, Queen's University Belfast.
GitHub: https://github.com/somritaghosh-1/polya-tnbc-nanopore
```
 
<br>
 
## Contact
 
**Somrita Ghosh**
MSc Molecular Biology and Biotechnology | Queen's University Belfast
 
[![Email](https://img.shields.io/badge/Email-somritaghosh56%40gmail.com-D14836?style=flat-square&logo=gmail)](mailto:somritaghosh56@gmail.com)
[![LinkedIn](https://img.shields.io/badge/LinkedIn-Somrita%20Ghosh-0077B5?style=flat-square&logo=linkedin)](https://www.linkedin.com/in/somrita-ghosh-014826211)
[![GitHub](https://img.shields.io/badge/GitHub-somritaghosh--1-181717?style=flat-square&logo=github)](https://github.com/somritaghosh-1)
 
<p align="center">
  <sub>Data generated at Queen's University Belfast | Patrick G Johnston Centre for Cancer Research | 2024</sub>
</p>
 
