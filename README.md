# MTP Documentation 




#### Table of Contents
- [Datasets](#datasets)
- [Histology and Molecular Subtyping](#Histology)
- [RNA-seq](#RNA-seq)
- [Somatic Variants](#Somatic)

## Datasets

**Open Pediatric Brain Tumor Atlas (OpenPBTA):** In September of 2018, the Children's Brain Tumor Network (CBTN) released the Pediatric Brain Tumor Atlas (PBTA), a genomic dataset (whole genome sequencing, whole exome sequencing, RNA sequencing, proteomic, and clinical data) for nearly 1,000 tumors, available from the Gabriella Miller Kids First Portal. The Open Pediatric Brain Tumor Atlas (OpenPBTA) Project is a global open science initiative to comprehensively define the molecular landscape of tumors of 943 patients from the CBTN and the PNOC003 DIPG clinical trial from the Pediatric Pacific Neuro-oncology Consortium through real-time, collaborative analyses and collaborative manuscript writing on GitHub.

**Therapeutically Applicable Research to Generate Effective Treatments (TARGET):** The Therapeutically Applicable Research to Generate Effective Treatments (TARGET) Initiative is an NCI-funded collection of disease-specific projects that seeks to identify the genomic changes of pediatric cancers. 'The overall goal is to collect genomic data to accelerate the development of more effective therapies. OpenPedCan analyses include the seven diseases present in the TARGET dataset: Acute Lymphoblastic Leukemia (ALL), Acute Myeloid Leukemia (AML), Clear cell sarcoma of the kidney, Neuroblastoma, Osteosarcoma, Rhabdoid tumor, and Wilm’s Tumor.

**Gabriella Miller Kids First Neuroblastoma (Kids First) The Gabriella Miller Kids First Pediatric Research Program (Kids First):** is a large-scale effort to accelerate research and gene discovery in pediatric cancers and structural birth defects. The program includes whole genome sequencing (WGS) from patients with pediatric cancers and structural birth defects and their families. OpenPedCan analyses include Neuroblastoma data from the Kids First project.

**The Genotype-Tissue Expression (GTEx):** GTEx project is an ongoing effort to build a comprehensive public data resource and tissue bank to study tissue-specific gene expression, regulation and their relationship with genetic variants. Samples were collected from 54 non-diseased tissue sites across nearly 1000 individuals, primarily for molecular assays including WGS, WES, and RNA-Seq. OpenPedCan project includes 17382 GTEx RNA-Seq samples from GTEx v8 release, which span across 31 GTEx groups in the v10 release.

**The Cancer Genome Atlas Program (TCGA):** TCGA is a landmark cancer genomics program that molecularly characterized over 20,000 primary cancer and matched normal samples spanning 33 cancer types. It is a joint effort between NCI and the National Human Genome Research Institute. OpenPedCan project includes 10414 TCGA RNA-Seq samples (716 normal and 9698 tumor) from 33 cancer types in the v10 release.

---

## Histology and Molecular Subtyping

---

## RNA-seq

This document describes the data processing steps. For a description of the data see LINK TO ABOUT MOLECULAR TARGETS and to browse the Pediatric Cancer Data see [https://ppdc-otp-dev.bento-tools.org/pediatric-cancer-data-navigation](https://ppdc-otp-dev.bento-tools.org/pediatric-cancer-data-navigation)

### RNA-seq Data Processing

The RNA-seq Alignment Workflow begins by trimming adapters, only if adapters are provided, using Cutadapt. Trimmed reads were then aligned using STAR in two-pass mode to reference genome GRCh38. Transcripts are quantified using both RSEM and by pseudoalignment using Kallisto using the Gencode v27 annotation. Fusion calling is done using both Arriba and STAR-Fusion and then filtered for high confidence fusion calls using annoFuse. QC metrics for the alignment are summarized using RNA-seQC. If you would like to view the code in more detail, please see the GitHub repository at [https://github.com/kids-first/kf-rnaseq-workflow](https://github.com/kids-first/kf-rnaseq-workflow) and if you would like to run the pipeline, please see the CAVATICA App at [https://cavatica.sbgenomics.com/public/apps/cavatica/apps-publisher/kfdrc-rnaseq-workflow/6](https://cavatica.sbgenomics.com/public/apps/cavatica/apps-publisher/kfdrc-rnaseq-workflow/6) (for more documentation on CAVATICA, please click READ MORE).

### RNA-seq Data Table Generation

#### Fusions

Fusions are filtered using custom R scripts. Fusion calls are retained if they are called by both STAR-Fusion and Arriba and if the fusion was specific and present in 3 or more samples in a single broad_histology/disease. Fusions were then annotated with gene and fusion specific information as well as whether they are known cancer genes from OncoKB, TCGA, and COSMIC. Summary frequencies are calculated using R. See [https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/fusion_filtering](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/fusion_filtering) for specific code and further details.

#### Transcript Expression

TPMs (transcripts per million reads) were calculated using RSEM and plotted using R. See CAVATICA public app for more details on RSEM.

---

## Somatic Variants
