# Cell-Meta-Analysis-Bhaskar

**A Cell-Type-Resolved Meta-Analysis Reveals Glial DNA Methylation Changes Associated with Aging and Alzheimer's Disease**

Uchit Bhaskar, Mark Z. Kos*, Melanie A. Carless*

Corresponding authors: mark.kos@utrgv.edu OR melanie.carless@utsa.edu

**Summary:**
Epigenome-wide association studies implicate DNA methylation in the development and progression of Alzheimer's disease (AD). Although recent studies show that the epigenetics of non-neuronal cell types contribute to disease risk, the role of the methylome in individual glial cell types (i.e., astrocytes, oligodendrocytes) in biological aging and AD pathogenesis is unclear. In this study, we examined archived DNA methylation data across 13 cohorts and performed cell-type deconvolution in silico to identify novel epigenetic signatures associated with aging and AD in glial cells. We observed pronounced age-associated methylation in astrocytes within the prefrontal cortex, whereas oligodendrocytes of the entorhinal cortex show the most differential methylation with AD status. Astrocytes, along with neurons, within the prefrontal cortex, emerge as key players in Braak stage-associated methylation, exhibiting strong concordance with previously reported associations at the brain tissue level. Age-associated changes in oligodendrocytes exhibit strong directional correlation with, and amplification of age-related effects with AD that affect neurodevelopmental processes, while AD-related methylation changes at age-associated sites in astrocytes diverge from those representative of normative aging processes. Our study expands on previous findings and reveals glial-specific methylation patterns associated with epigenetic aging and AD.

Paper can be accessed on BioRxiv (https://www.biorxiv.org/content/10.64898/2026.05.04.722662v1)

This page contains the main analysis scripts conducted for this paper:

1. preprocessing.r: Script for all preprocessing, data quality and harmonization and includes
     - reading input files
     - Sample-level QC steps to filter for a) bisulfite conversion efficiency, b) abnormal background intensity of negative probes,  c) outlier mean methylated or          unmethylated intensities, d) discordance between predicted and reported sex, and e) genetic identity mismatches detected using SNP probes.
     - Probe-level QC steps including a) background correction and dye-bias normalization using NOOB, b) filtering probes with detection P values > 0.01, and c)            those with < 3 beadcounts
     - BMIQ-based normalization, and
     - EpiSCORE-based cell type proportion estimation
     
2. cell-specific-association-analysis.r: Script for cohort level cell-specific association analysis and includes
     - harmonization of probe-set across 450k/EPIC arrays
     - implementation of surrogate variables for model fit across different brain regions
     - CellDMC based association testing for one of three phenotypes (age, AD, or Braak stage) for each brain cell type 

3. main_meta_analysis.r: Script for overall meta-analysis, and plots generated for the paper, including
     - integration of CellDMC-based estimates across cohorts
     - BACON-based correction for bias and inflation
     - Fixed-effect inverse variance weighting-based meta-analysis per region, per cell type using the _metagen_ function
     - Codes for summary plots, manhattan plots and volcano plots
     - Gene onotology and KEGG pathway enrichment analysis
     - Effect-size rank enrichment analysis
     - Age-to-disease trajectory analysis

All analysis was run using R v4.5.1, Bioconductor v3.22 and Python v2.7. Programs and packages were installed via anaconda on a high-performance computing cluster system (Red Hat Enterprise Linux 8.10).  
