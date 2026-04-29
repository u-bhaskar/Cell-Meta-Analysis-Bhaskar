# Load required libraries ----------------------------------------------------------------
suppressPackageStartupMessages({
  library(minfi)
  library(wateRmelon)
  library(minfiData)
  library(ewastools)
  library(ChAMP)
  library(IlluminaHumanMethylation450kanno.ilmn12.hg19) #adjust if EPIC array
  library(IlluminaHumanMethylation450kmanifest)
  library(maxprobes)
  library(RnBeads)
  library(RnBeads.hg19)
  library(lumi)
  library(dplyr)
  library(data.table)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(cowplot)
  library(devtools)
  library(presto)
  library(EpiSCORE)
  })

# The following steps need to be performed once per cohort
# Load sample metadata file
  pheno <- read.csv(".csv")

# Load IDAT files as rgSet object
  raw_dir <- c("") #path to directory containing IDAT files
  rgSet <- read.metharray.exp(base = raw_dir, recursive = TRUE, 
                              extended = TRUE, verbose = TRUE, force = TRUE)

# Perform Sample QC-----------------------------------------------------------------------------
# Populate Basename in pheno data
  if(is.null(pheno$Basename)){
    pheno$Basename <- file.path(raw_dir, paste0(sampleNames(rgSet)))
  }
  
# Ensure pheno and rgSet files match
  pheno <- pheno[order(rownames(pheno)), ]
  rownames(pheno) <- colnames(rgSet)
  stopifnot(identical(colnames(rgSet), rownames(pheno))) 
  
# a) Bisulfite conversion efficiency: identify samples with bisulfite conversion efficiency < 0.8
  bsEff <- bscon(rgSet)
  lowBS <- rownames(bsEff)[bsEff < 0.8]

# b) Remove samples with extreme intensities using negative controls as references (following Smith et. al. 2021)
  #Control probe matrices
  contInfo <- minfi::getProbeInfo(rgSet, type = "Control")
  
  redCont <- minfi::getRed(rgSet)[contInfo$Address, , drop = FALSE]
  grnCont <- minfi::getGreen(rgSet)[contInfo$Address, , drop = FALSE]
  rownames(redCont) <- contInfo$ExtendedType
  rownames(grnCont) <- contInfo$ExtendedType
  
  # Identify negative control probes
  negIdx <- which(contInfo$Type == "NEGATIVE")
  
  negRedMean <- colMeans(redCont[negIdx, , drop = FALSE], na.rm = TRUE)
  negGrnMean <- colMeans(grnCont[negIdx, , drop = FALSE], na.rm = TRUE)
  redMean <- colMeans(minfi::getRed(rgSet), na.rm = TRUE)
  grnMean <- colMeans(minfi::getGreen(rgSet), na.rm = TRUE)
  
  extremeIntensity <- names(which(negRedMean > 1000 | negGrnMean > 1000 | redMean < 2000 | grnMean < 2000))
  
# c) Detect samples with mean M/U intensities +/- 3SD of the mean
  M <- minfi::getMeth(preprocessRaw(rgSet))
  U <- minfi::getUnmeth(preprocessRaw(rgSet))
  
  Mmean <- colMeans(M, na.rm = TRUE)
  Umean <- colMeans(U, na.rm = TRUE)
  
  signalOutliers <- names(which(
  Mmean > mean(Mmean) + 3*sd(Mmean) |
  Mmean < mean(Mmean) - 3*sd(Mmean) |
  Umean > mean(Umean) + 3*sd(Umean) |
  Umean < mean(Umean) - 3*sd(Umean)))
  
# d) Detect genetic correlation between samples based on SNP probes
  snpBetas <- getSnpBeta(rgSet)
  
  if(nrow(snpBetas) > 0){
  
  snpCor <- cor(snpBetas, use = "pairwise.complete.obs")
  diag(snpCor) <- NA
  
  snpFail <- rep(FALSE, ncol(snpCor))
  names(snpFail) <- colnames(snpCor)
  
  for(i in seq_len(ncol(snpCor))){
  
  highCorr <- names(which(snpCor[, i] > 0.65))
  
  if (length(highCorr) > 0){
  if (is.null(donor_id_col) ||
    length(unique(pheno[[donor_id_col]])) == nrow(pheno)) {
  
  message("All donors unique — skipping SNP relatedness check.")
  snpFailNames <- character(0)

} else {
  donor_i <- pheno[[donor_id_col]][i]
  donor_j <- pheno[[donor_id_col]][match(highCorr, rownames(pheno))]
  
  if(any(donor_j != donor_i)){
  snpFail[i] <- TRUE
  }
  }
  }
  }} else {
  message("No SNP probes detected; skipping SNP identity check.")
  snpFail <- rep(FALSE, ncol(rgSet))
  names(snpFail) <- colnames(rgSet)
  }
  
# e) Sex Prediction check
  mSetRaw <- preprocessRaw(rgSet)
  gmSet <- mapToGenome(mSetRaw)
  
  predSex <- getSex(gmSet)$predictedSex
  reportedSex <- pheno$sex
  
  sexMismatch <- names(which(predSex != reportedSex))

# f) Combine all flags and remove samples that failed QC
  failedSamples <- unique(c(lowBS, extremeIntensity, signalOutliers, names(snpFail)[snpFail], sexMismatch))
  
  rgSet <- rgSet[, !(colnames(rgSet) %in% failedSamples)]
  pheno <- pheno[!(rownames(pheno) %in% failedSamples), ] 

# Perform probe-level QC-----------------------------------------------------------------------------
# a) Generate detection P-values
  detPewas <- ewastools::detectionP.minfi(rgSet)  
  detPewas <- detPewas[, order(colnames(detPewas))]
  
# b) Generate beadcount data
  beadCount <- beadcount(rgSet)
  
# c) Background correction and dye bias normalization
  mSet <- preprocessNoob(rgSet, dyeCorr = TRUE, dyeMethod = "single") # generate methyl-set object
  rSet <- ratioConvert(mSet, what = "both", keepCN = TRUE) #Convert mSet to RatioSet
  # Convert Ratio Set → GenomicRatioSet
  grSet <- mapToGenome(rSet)

# d) Remove Sex-specific probes
  keepAutosomal <- !(seqnames(grSet) %in% c("chrX", "chrY"))
  grSet <- grSet[keepAutosomal, ]
  
# e) Remove probes with SNPs at CpG sites
  grSet <- mapToGenome(grSet)
  grSet <- dropLociWithSnps(grSet)
  
# f)Remove Cross-Reactive Probes
  # maxprobes package provides cross-reactive probe lists
  crossReact <- xreactive_probes("450K") # or "EPIC"
  keepNonReactive <- !(featureNames(grSet) %in% crossReact)
  grSet <- grSet[keepNonReactive, ]
  
# Generate Beta Values
  betaValues <- getBeta(grSet)
  
# Filter samples based on detP values
  detPewas <- detPewas[rownames(detPewas) %in% rownames(betaValues), colnames(detPewas) %in% colnames(betaValues)]
  detPewas <- detPewas[order(rownames(detPewas)), ]
  betaValues <- betaValues[order(rownames(betaValues)), ]

# subset bead count
  beadCount <- beadCount[rownames(betaValues), ]

# Use CHAMP to filter out bad probes
  champ <- champ.filter(beta = betaValues,
                        pd = pheno,
                        detP = detPewas,
                        beadcount = beadCount,
                        ProbeCutoff = 0.05, #remove probes only if at least 5% of samples fail detP
                        SampleCutoff = 0.1, #remove samples with >10% probes failing detP
                        detPcut = 0.01, #remove probes with detP > 0.01
                        filterDetP = TRUE,
                        beadCutoff = 0.05, # <3 beads in at least 5% samples
                        filterBeads = TRUE,
                        arraytype = "450K"
  )
  
  # Extract list of beta values from champ filtering
  betaFiltered <- champ$beta
  phenoUpdated <- pheno[match(colnames(betaFiltered), rownames(pheno)), , drop = FALSE]
  
# If IDAT file input is not available, but Signal intensity files are, use the following commands for input-----------------------------------------------------------------------------
# Then use mSet object and proceed to probe QC as before
  meth <- data.table::fread(meth_file, data.table = FALSE) # meth_file should be the path to methylation data file
  meth <- as.matrix(meth)

  unmeth <- data.table::fread(unmeth_file, data.table = FALSE) # unmeth_file should be the path to the unmethylated data file
  unmeth <- as.matrix(unmeth)

  detP <- data.table::fread(detP_file, data.table = FALSE)  #detP_file should be the path to the detection P-value data file
  detP <- as.matrix(detP)

# Build methyl-set object
  mSet <- minfi::MethylSet(
  Meth = meth,
  Unmeth = unmeth,
  colData = pheno,
  annotation = c(array = "IlluminaHumanMethylation450k", annotation = "ilmn12.hg19")
  )

# BMIQ-based Normalization of Beta Values-----------------------------------------------------------------------------
# Compute m-values from beta - Useful for statistics
  mvals <- lumi::beta2m(betaFiltered)
  
# Normalize (BMIQ method) using Champ
  betaBMIQ <- ChAMP::champ.norm(
    beta = betaFiltered,
    arraytype = "450K",
    method = "BMIQ",
    cores = 1
  )

# Ensure pheno data is aligned
  phenoUpdated <- phenoUpdated[match(colnames(betaBMIQ), rownames(phenoUpdated)), ]

# Estimate cell type proportions using EpiSCORE-----------------------------------------------------------------------------
# a) Read brain methylation reference file
  brainRef <- read.csv("Brain_Mref.csv", stringsAsFactors = FALSE) # can be downloaded from: https://www.biosino.org/episcore/file/Brain_Mref.csv
  rownames(brainRef) <- brainRef$Genes
  brainRef <- as.matrix(brainRef)
 
# b) Collapse CpGs to Promoter Levels, using Entrez IDs
  promoterBeta <- EpiSCORE::constAvBetaTSS(betaBMIQ, type = episcoreType)
  
# c) Estimate cell type proportions using weighted robust partial correlation method
  cellEstimates <- EpiSCORE::wRPC(promoterBeta, ref = brainRef, useW = TRUE, wth = 0.4, maxit = 300)
  cellFractions <- as.data.frame(cellEstimates$estF)

# d) Append cell fractions to pheno data file
  phenoUpdated <- bind_cols(phenoUpdated, cellFractions)

# Generate final pre-processed files-----------------------------------------------------------------------------
  commonSamples <- Reduce(intersect, list(colnames(betaBMIQ), rownames(phenoUpdated),
                                          rownames(cellFractions)))
  betaFinal <- betaBMIQ[, commonSamples, drop = FALSE]
  phenoFinal <- phenoUpdated[commonSamples, , drop = FALSE]
  cellFracFinal <- cellFractions[commonSamples, , drop = FALSE]
