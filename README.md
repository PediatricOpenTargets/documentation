# MTP Documentation 




#### Table of Contents
- [Datasets](#datasets)
- [Histology and Molecular Subtyping](#Histology)
- [RNA-seq](#RNA-seq)
- [Somatic Variants](#Somatic_Variants)

## Datasets

**Open Pediatric Brain Tumor Atlas (OpenPBTA):** In September of 2018, the Children's Brain Tumor Network (CBTN) released the Pediatric Brain Tumor Atlas (PBTA), a genomic dataset (whole genome sequencing, whole exome sequencing, RNA sequencing, proteomic, and clinical data) for nearly 1,000 tumors, available from the Gabriella Miller Kids First Portal. The Open Pediatric Brain Tumor Atlas (OpenPBTA) Project is a global open science initiative to comprehensively define the molecular landscape of tumors of 943 patients from the CBTN and the PNOC003 DIPG clinical trial from the Pediatric Pacific Neuro-oncology Consortium through real-time, collaborative analyses and collaborative manuscript writing on GitHub.

**Therapeutically Applicable Research to Generate Effective Treatments (TARGET):** The Therapeutically Applicable Research to Generate Effective Treatments (TARGET) Initiative is an NCI-funded collection of disease-specific projects that seeks to identify the genomic changes of pediatric cancers. 'The overall goal is to collect genomic data to accelerate the development of more effective therapies. OpenPedCan analyses include the seven diseases present in the TARGET dataset: Acute Lymphoblastic Leukemia (ALL), Acute Myeloid Leukemia (AML), Clear cell sarcoma of the kidney, Neuroblastoma, Osteosarcoma, Rhabdoid tumor, and Wilm’s Tumor.

**Gabriella Miller Kids First Neuroblastoma (Kids First) The Gabriella Miller Kids First Pediatric Research Program (Kids First):** is a large-scale effort to accelerate research and gene discovery in pediatric cancers and structural birth defects. The program includes whole genome sequencing (WGS) from patients with pediatric cancers and structural birth defects and their families. OpenPedCan analyses include Neuroblastoma data from the Kids First project.

**The Genotype-Tissue Expression (GTEx):** GTEx project is an ongoing effort to build a comprehensive public data resource and tissue bank to study tissue-specific gene expression, regulation and their relationship with genetic variants. Samples were collected from 54 non-diseased tissue sites across nearly 1000 individuals, primarily for molecular assays including WGS, WES, and RNA-Seq. OpenPedCan project includes 17382 GTEx RNA-Seq samples from GTEx v8 release, which span across 31 GTEx groups in the v10 release.

**The Cancer Genome Atlas Program (TCGA):** TCGA is a landmark cancer genomics program that molecularly characterized over 20,000 primary cancer and matched normal samples spanning 33 cancer types. It is a joint effort between NCI and the National Human Genome Research Institute. OpenPedCan project includes 10414 TCGA RNA-Seq samples (716 normal and 9698 tumor) from 33 cancer types in the v10 release.

---

## Histology

---

## RNA-seq

This document describes the data processing steps. For a description of the data see LINK TO ABOUT MOLECULAR TARGETS and to browse the Pediatric Cancer Data see [https://ppdc-otp-dev.bento-tools.org/pediatric-cancer-data-navigation](https://ppdc-otp-dev.bento-tools.org/pediatric-cancer-data-navigation)

### RNA-seq Data Processing

The RNA-seq Alignment Workflow begins by trimming adapters, only if adapters are provided, using Cutadapt. Reads were then aligned using STAR in two-pass mode to reference genome GRCh38. While all data is paired-end, methods are provided for single-end alignment if you are interested in processing your data in the same manner. Transcripts are quantified using RSEM  with the Gencode v27 annotation. Fusion calling is done using both Arriba and STAR-Fusion and then filtered for high confidence fusion calls using annoFuse. QC metrics for the alignment are summarized using RNA-seQC. If you would like to view the code in more detail, please see the GitHub repository at [https://github.com/kids-first/kf-rnaseq-workflow](https://github.com/kids-first/kf-rnaseq-workflow) and if you would like to run the pipeline, please see the CAVATICA App at [https://cavatica.sbgenomics.com/public/apps/cavatica/apps-publisher/kfdrc-rnaseq-workflow/6](https://cavatica.sbgenomics.com/public/apps/cavatica/apps-publisher/kfdrc-rnaseq-workflow/6) (for more documentation on CAVATICA, please click READ MORE).

### RNA-seq Data Table Generation

#### Fusions

Fusions are filtered using custom R scripts. Fusion calls are retained if they are called by both STAR-Fusion and Arriba and if the fusion was specific and present in 3 or more samples in a single disease. Fusions were then annotated with gene and fusion specific information as well as whether they are known cancer genes from OncoKB, TCGA, and COSMIC. Summary frequencies are calculated using R. See https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/fusion_filtering for specific code and further details.

| Annotation | Description | Values |
| --- | --- | --- |
| Fusion Name | Genes fused with the name of the genes fused separated by “--”. Gene order is the order they fused in with the dashes representing the breakpoint. |  |
| Fusion Type |  | frameshift, other |
| Gene Position | Whether the gene is to the left or right of the breakpoint. | Gene1A = gene to left of breakpoint, Gene1B = gene to right of breakpoint |
| Fusion Annotation |  |  |
| Breakpoint Location |  | genic |
| Annotations |  |  |
| Kinase Domain Retained Gene1A | Whether the kinase domain is retained; left blank if gene is not a kinase. | Yes/No |
| Kinase Domain Retained Gene1B | Whether the kinase domain is retained; left blank if gene is not a kinase. | Yes/No |
| Reciprocal exists either gene kinase | Whether or not the reciprocal fusion with the gene order around the breakpoint swapped exists | TRUE/FALSE |
| Gene1A/Gene1B/Gene2A/Gene2B Annotation | Left blank if no annotations | TumorSuppressorGene, Oncogene, TranscriptionFactor |
| Gene Ensembl ID | Official gene ID from Ensembl  |  |
| Dataset | See Pediatric Cancer Data Sources at https://ppdc-otp-dev.bento-tools.org/about | All Cohorts = all data sets combined, TARGET (Therapeutically Applicable Research to Generate Effective Treatments), GMFK (Gabriella Miller Kids First Neuroblastoma)  |
| Disease | See Histology sections. |  |
| Total alterations / Subjects in Dataset | Total number of samples with the given fusions over the total number of disease samples in the given Dataset |  |
| Frequency in overall dataset | Fraction of the samples for the given disease in the given dataset that have the fusion |  |
| Total primary tumors altered / Primary tumors in dataset | Same as Total alterations, but for primary tumors only |  |
| Frequency in primary tumors | Same as overall frequency, but for primary tumors only |  |
| Total relapse tumors altered / Relapse tumors in dataset | Same as Total alterations, but for relapse tumors only |  |
| Frequency in relapse tumors | Same as Total alterations, but for relapse tumors only |  |
| Gene full name | Standard full name of the gene (from HGNC I assume?) |  |
| PMTL | Whether the gene is a relevant target on the PMTL (Pediatric Molecular Target List); left blank if not |  |
| OncoKB Cancer Gene | Yes or no whether the gene is a annotated cancer gene listed in OncoKB https://www.oncokb.org/ | Y/N |
| OncoKB Oncogene | TSB | Whether the gene is an oncogene or tumor suppressor (TSG) | oncogene, TSG = tumor suppressor gene |

#### Transcript Expression

TPMs (transcripts per million reads) were calculated using RSEM and plotted using R. See CAVATICA public app for more details on RSEM. [https://cavatica.sbgenomics.com/public/apps/cavatica/apps-publisher/kfdrc-rnaseq-workflow/6](https://cavatica.sbgenomics.com/public/apps/cavatica/apps-publisher/kfdrc-rnaseq-workflow/6) 

---

## Somatic_Variants
