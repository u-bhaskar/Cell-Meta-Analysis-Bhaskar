# Load required packages------------------------------------------------
suppressPackageStartupMessages({
    library(EpiDISH)
    library(sva)
    library(minfi)
    library(ggplot2)
    library(dplyr)
    library(tidyr)
    library(pheatmap)
    library(RColorBrewer)
    library(gridExtra)
    library(qqman)
    library(cowplot)
    library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
    library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)})

# Load annotations and generate common CpG probe list
  anno450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  annoEPIC <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
  commonCpGs <- intersect(rownames(anno450k), rownames(annoEPIC))

# Subset preprocessed beta file to include only common CpGs
  beta <- betaFinal[rownames(betaFinal) %in% commonCpGs, , drop = FALSE]
  pheno <- phenoFinal
  cellFrac <- cellFracFinal

# Filter cell fractions to remove low abundance cell types, combine oligo/OPC fractions
  cellFracFiltered <- cellFrac %>%
      dplyr::select(Neuron, Astro, Oligo, OPC, Endo, Microglia) %>%
      mutate(
        Oligo_OPC = Oligo + OPC
      ) %>%
      dplyr::select(-Oligo, -OPC)
  
# Identify low abundance cell types
  meanFrac <- colMeans(cellFracFiltered)
  keepCells <-meanFrac >= 0.005 # remove cell types with <0.5% mean fraction
  
  if(!all(keepCells)) {
    dropped <- names(meanFrac)[!keepCells]
    cellFracFiltered <- cellFracFiltered[, keepCells, drop = FALSE]
  }
  
# Generate final beta, cell frac and pheno input files for CellDMC analysis
  common <- Reduce(
    intersect, list(colnames(beta), rownames(pheno), rownames(cellFracFiltered))
  )
  
  beta <- beta[, common, drop = FALSE]
  pheno <- pheno[common, , drop = FALSE]
  cellFrac <- cellFracFiltered[common, , drop = FALSE]

# Specify models to run for CellDMC (i.e., age, AD, and Braak Stage Associations)-------------------------------------------------------------------------

age_by_region <- list(
    description = "Age effects stratified by region",
    outcome = "age",
    covariates = c("sex", "diagnosis"),
    cell_fractions = TRUE, 
    stratify_by = "region"
  )

AD_by_region <- list(
      description = "AD association stratified by region",
      outcome = "diagnosis",
      covariates = c("age", "sex"),
      cell_fractions = TRUE, 
      stratify_by = "region"
    )

braak <- list(
    description = "Braak (continuous) pathology burden, region-stratified",
    outcome = "braak",
    covariates = c("age", "sex"),
    cell_fractions = TRUE,
    stratify_by = "region"
    )

# Specify which model to run
model <- c("") # run each model once per cohort

# Check sample size
# Min. number of samples per cohort: 30, Min. samples per trait: 5, Min. samples per region: 5. 
  n_total <- nrow(pheno)
  if(n_total < 30){
    message("Too few samples, cannot run model")
  }

  outcome <- model$outcome
  trait <- pheno[[outcome]]
  if(is.factor(trait) && nlevels(trait) > 1){
    groupCounts <- table(trait)
    if(any(groupCounts < 5)){
      message("Outcome groups too small, cannot run model")
    }
  }

  strat <- model$stratify_by
  strat_counts <- table(pheno[[strat]])
  if(any(strat_counts < 5)){
      message("One or more strata too small for region stratified analysis, cannot run model")
    }

# Run region-stratified association analysis----------------------------------------------------------------------------
# Split data by regions
stratLevels <- unique(pheno[[strat]])

strats <- lapply(as.character(stratLevels), function (s) {
    keep <- pheno[[strat]] == s
    phenoStrat <- droplevels(pheno[keep, , drop = FALSE])
    cellFracStrat <- cellFrac[rownames(phenoStrat), , drop = FALSE]
    ctVar <- apply(cellFracStrat, 2, var)
    keepCellTypes <- names(ctVar)[ctVar > 0]
    cellFracStrat <- cellFracStrat[, keepCellTypes, drop = FALSE]
    
    list(
      betaMod = beta[, rownames(phenoStrat), drop = FALSE],
      phenoMod = phenoStrat,
      cellFracMod = cellFracStrat,
      strata = s
    )
  })

for(s in names(strats)){
  
  beta <- strats[[s]]$betaMod
  pheno <- strats[[s]]$phenoMod
  cellFrac <- strats[[s]]$cellFracMod
  phenoSVA <- cbind(pheno, cellFrac)
  
  #convert beta to m-vals
  mVals <- lumi::beta2m(beta)

  # Compute number of surrogate variables (SVs) to use -----------------------
   # Collect covariates and outcome
   outcome <- model$outcome  
   baseCovs <- intersect(model$covariates, colnames(phenoSVA))
   keepCovs <- unique(c(outcome, baseCovs, cfCols))
 
    if (is.factor(phenoSVA[[outcome]]) && nlevels(droplevels(phenoSVA[[outcome]])) < 2) {
    message("Outcome has <2 levels in this stratum - skipping SVA")
    }
  
    cc <- stats::complete.cases(phenoSVA[, keepCovs, drop = FALSE])
    phenoSVA <- phenoSVA[cc, , drop = FALSE]
    mVals <- mVals[, cc, drop = FALSE]
    var_post <- apply(mVals, 1, var)
    mVals <- mVals[var_post > 1e-6, , drop = FALSE]

    # Build SVA design, outcome as predictor
    # mod = ~ outcome + covars; mod0 = ~ covars

    if (length(keepCovs) > 1) {
    formFull <- as.formula(paste("~", paste(keepCovs, collapse = " + ")))
    formNull <- as.formula(paste("~", paste(setdiff(keepCovs, outcome), collapse = " + ")))
    } else {
    formFull <- as.formula(paste("~", outcome))
    formNull <- ~1
    }

   mod <- stats::model.matrix(formFull, data = phenoSVA)
   mod0 <- stats::model.matrix(formNull, data = phenoSVA)
   common_cols <- intersect(colnames(mod0), colnames(mod))
   mod0 <- mod0[, common_cols, drop = FALSE]
  
  # Compute residual df and k-cap
  n_samples <- ncol(mVals)
  df_resid <- n_samples - qr(mod)$rank
  
  # Estimate Number of SVs based on Lambda (similar to Smith et. al. 2021)
  lambdaTarget <- 1.20
  max_sv <- 20
  max_allowed <- max(0L, min(max_sv, floor(df_resid/4)))
  
  # Compute lambda for the outcome/trait
  outcomeCols <- grep(paste0("^", outcome, "$|^", outcome), colnames(mod), perl = TRUE)
  outcome_col <- outcomeCols[1]
 
  lambdaOut <- function(extra) {
  design_tmp <- if(is.null(extra)) mod else cbind(mod, extra)
  fit <- limma::lmFit(mVals, design_tmp)
  fit <- limma::eBayes(fit)
  p <- fit$p.value[, outcome_col]
  p <- p[is.finite(p) & p > 0 & p <= 1]
  if (length(p) < 100) return (NA_real_)
  median(stats::qchisq(1-p, df = 1), na.rm = TRUE)/stats::qchisq(0.5, df=1)
  }
 
  # Check lambda at baseline
  lam0 <- lambdaOut(NULL)
  message("k = 0 | lambda = ", round(lam0, 3))
 
  # If lambda < 1.2, skip SVA
  if (is.finite(lam0) && lam0 <= lambdaTarget) {
  message("No SVs required: baseline lambda <= target")
  return(NULL)
  }
 
  # Try increasing values of k
  chosen_sv <- NULL
  chosen_k <- 0L

  for (k in seq_len(max_allowed)) {
  svobj <- try(sva::sva(dat = mVals, mod = mod, mod0 = mod0, n.sv = k))
  SV <- as.data.frame(svobj$sv)
  colnames(SV) <- paste0("SV", seq_len(ncol(SV)))
  rownames(SV) <- rownames(phenoSVA)
  lamk <- lambdaOut(SV)
  chosen_sv <- SV
    }

  # Append SVs to pheno data  
  phenoSVA <- cbind(phenoSVA, chosen_sv)

# Run CellDMC ------------------------------
  #Build design matrix
  baseFormula <- paste(
    "~",
    paste(c(model$covariates), collapse = " + ")
  )
  modelMatrix <- model.matrix(as.formula(baseFormula), data = phenoSVA)

  # drop constant model matrix columns
  zeroVar <- apply(modelMatrix, 2, function(x) isTRUE(all(x == x[1])))
  modelMatrix <- modelMatrix[, !zeroVar, drop = FALSE]

  # Append SVs, if present, as covariates to design
  svCols <- grep("^SV[0-9]+$", colnames(phenoSVA), value = TRUE)
  if(length(svCols) > 0){
  modelMatrix <- cbind(modelMatrix, as.matrix(phenoSVA[, svCols, drop = FALSE]))}

  # Run CellDMC                 
  fit <- EpiDISH::CellDMC(
    beta.m = as.matrix(beta),
    pheno.v = phenoSVA[[outcome]],
    frac.m = as.matrix(cellFrac),
    cov.mod = modelMatrix,
    adjPMethod = "fdr",
    adjPThresh = 0.05,
    sort = TRUE
  )

  # Generate association result tables for each cell type
  cellTypes <- names(fit$coe)
  resList <- lapply(cellTypes, function (ct){
    
    df <- fit$coe[[ct]] %>%
      as.data.frame() %>%
      tibble::rownames_to_column("CpG") %>%
      mutate(
        CellType = ct,
        Direction = ifelse(Estimate > 0, "Hyper", "Hypo"),
        FDR_sig = adjP < 0.05,
        Model = model,
        Array = "450k", #change if EPIC
        Region = s
      )
    
    df
  })
  
  combinedResults <- bind_rows(resList)
  combinedResultsTable <- combinedResults %>%
    dplyr::select(CpG, CellType, Estimate, SE, t, p, adjP) %>%
    arrange(CellType, adjP)
  
  suffix <- if(!is.null(strat)) paste0("_", strat) else "" 
  
  write.csv(combinedResultsTable,
            file.path(c(".csv"))) # Save per cohort per model               
}
