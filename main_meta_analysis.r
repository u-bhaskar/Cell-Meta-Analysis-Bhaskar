#------------------------ Load required packages and inputs --------------------------
suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
  library(ggplot2)
  library(cowplot)
  library(ggrepel)
  library(meta)
  library(bacon)
  library(missMethyl)
  library(minfi)
  library(GO.db)
  library(AnnotationDbi)
  library(KEGGREST)
  library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  library(ggridges)
  library(viridis)
})

# Load input data
phenotype <- c("") # set phenotype to one of "Age", "AD" or "Braak"
input_dir <- c("") # set path to input directory containing files for meta-analysis, unique to each phenotype
files <- list.files(input_dir, pattern = "\\.csv$", full.names = TRUE)
file_list <- data.table(path = files, file = basename(files))
dt <- rbindlist(lapply(file_list$path, fread, use.names = TRUE, fill = TRUE))

# Load annotation data
annData <- minfi::getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
anno <- as.data.table(data.frame(
      CpG = rownames(annData),
      chr = as.character(annData$chr),
      bp  = annData$pos,
      UCSC_RefGene_Name  = annData$UCSC_RefGene_Name,
      UCSC_RefGene_Group = annData$UCSC_RefGene_Group,
      Relation_to_Island = annData$Relation_to_Island
    ), keep.rownames = FALSE)
dt <- dt[anno, on = .(CpG), nomatch = 0L]

# Correct for bias/inflation using BACON
keys <- unique(dt[, .(Cohort, Region, CellType)]) # get a unique list of cohort x region x celltype combinations in the data
out <- vector("list", nrow(keys))

for (i in seq_len(nrow(keys))) {
  cohort_i <- keys$Cohort[i]
  region_j <- keys$Region[i]
  celltype_k <- keys$CellType[i]
  
  sub <- dt[Cohort == cohort_i & Region == region_j & CellType == celltype_k & is.finite(Estimate) & is.finite(SE) & SE > 0]
  n_sub <- nrow(sub)

  # Calculate lambda prior to BACON correction
  # If lambda > 2.0, cohort will be omitted from analysis
  lambda_pre <- median(qchisq(1- sub$p, 1), na.rm = TRUE)/ qchisq(0.5, 1)
  if(is.finite(lambda_pre) && lambda_pre >= 2.0){next}

  # Perform BACON-based correction
  bc <- bacon::bacon(effectsizes = sub$Estimate, standarderrors = sub$SE, verbose = FALSE)
  # Extract corrected ES and SEs, and estimate bias and inflation metrics
  es_corr <- bacon::es(bc)
  se_corr <- bacon::se(bc)
  bc_bias <- as.numeric(bacon::bias(bc))
  bc_infl <- as.numeric(bacon::inflation(bc))

  # Overwrite in sub
  sub <- sub[, Estimate := as.numeric(es_corr)]
  sub <- sub[, SE := as.numeric(se_corr)]
  sub <- sub[, p := 2*pnorm(-abs(Estimate/SE))]

  # Calculate lambda post BACON correction
  lambda_post <- median(qchisq(1- sub$p, 1), na.rm = TRUE)/ qchisq(0.5, 1)

  # If lambda after BACON correction is > 1.5, exclude those cohorts
  if (is.finite(lambda_post) && lambda_post > 1.5){next}

  # Retain cohorts that pass filtering metrics
  out[[i]] <- sub
  }

# Return combinations that pass filtering
dt <- rbindlist(out, use.names = TRUE, fill = TRUE)

#------------------------ Fixed effect meta analyis --------------------------
# analysis is done at the phenotype x region x celltype level

setkeyv(dt, c("Region","CellType","CpG"))
keys <- unique(dt[, .(Region, CellType, CpG)]) # get a unique set of region x celltype x CpG combination per phenotype for meta-analysis
out <- vector("list", nrow(keys))

for (i in seq_len(nrow(keys))) {
  cpg_i <- keys$CpG[i]
  region_j <- keys$Region[i]
  celltype_k <- keys$CellType[i]
  
  sub <- dt[CpG == cpg_i & Region == region_j & CellType == celltype_k]
  sub <- sub[is.finite(Estimate) & is.finite(SE) & SE > 0]
  
  if (nrow(sub) < 2) next # if number of cohorts available for the particular combination is <2, do not perform meta-analysis

  # perform meta-analysis
  meta_res <- meta::metagen(TE = sub$Estimate, seTE = sub$SE, studlab = sub$Cohort, comb.fixed = TRUE, comb.random = FALSE, method.tau = "REML", hakn = FALSE, sm = "")
  
  # build final output data 
  cpg_meta <- data.table(
    n_studies = meta_res$k,
    CpG = sub$CpG[1],
    beta_hat = as.numeric(meta_res$TE.fixed),
    se_hat   = as.numeric(meta_res$seTE.fixed),
    pval     = as.numeric(meta_res$pval.fixed),
    Q        = as.numeric(meta_res$Q),
    I2       = as.numeric(meta_res$I2),
    H2       = as.numeric(meta_res$H^2),
    sign_agree = mean(sign(sub$Estimate) == sign(meta_res$TE.fixed))
  )

  if (!is.null(cpg_meta)){
    cpg_meta[, `:=`(CpG = cpg_i, Region = region_j, CellType = celltype_k)]
    out[[i]] <- cpg_meta
  }
}

res <- Filter(Negate(is.null), out)
fe <- rbindlist(res, use.names = TRUE, fill = TRUE)
fe <- fe[, `:=`(
      z = beta_hat / se_hat
      FDR = p.adjust(pval, method = "BH"),
      Direction = ifelse(beta_hat > 0, "Hyper", "Hypo"),
      Sig = FDR < 0.05,
      abs_beta_hat := abs(beta_hat)
    )]
setorder(fe, Region, CellType, FDR, -abs_beta_hat)

# Add annotation data to fe results
feA <- fe[anno, on = .(CpG), nomatch = 0L] 
fe_path <- c("") # set filepath to save results, file name should contain phenotype x region x celltype combination analyzed
fwrite(feA, fe_path, sep = "\t") 

#------------------------ Main plots --------------------------
# Note: Perform all meta-analysis before plot generation for a given phenotype
# Load meta-analysis data
meta_dir <- c("") # set path to directory containing meta analyzed tsv files for each phenotype x region x celltype combination
meta_files <- list.files(meta_dir, pattern = "\\.tsv\\.gz$", full.names = TRUE)
file_list <- data.table(path = meta_files, file = basename(meta_files))
meta <- rbindlist(lapply(file_list$path, fread, use.names = TRUE, fill = TRUE))

# We'll create functions to generate plots
# summary bubble/dot plots
summary_dotplot <- function(meta, phenotype) {
  # Aggregate counts of significant DMCs
  # We count hyper and hypo separately to show directionality via color
  summary_dt <- meta[Sig == TRUE, .(
    Total = .N,
    Hyper = sum(beta_hat > 0, na.rm = TRUE),
    Hypo  = sum(beta_hat < 0, na.rm = TRUE)
  ), by = .(Region, CellType)]
  summary_dt[, Direction := (Hyper - Hypo) / Total]
  
  p <- ggplot(summary_dt, aes(x = CellType, y = Region)) +
    geom_point(aes(size = Total, fill = Direction), shape = 21, stroke = 1) +
    scale_size_continuous(range = c(2, 18), name = "N Sig CpGs",
    breaks = seq(0, max(summary_dt$Total, na.rm = TRUE), length.out = 5) %>% round() %>% unique()) +
    scale_fill_gradient2(low = "#f2c45f", mid = "white", high = "#298c8c", 
                         midpoint = 0, name = "Proportion") +
    theme_cowplot() +
    labs(
      title = sprintf("Differential Methylation Burden: %s", phenotype),
      subtitle = sprintf("Comparison of %s-significant sites across cell types and regions", sig_type),
      x = "Cell Type", y = "Brain Region"
    )  

  filepath <- c("") # path to save the figure
  ggsave(filepath, p, width = 10, height = 8)
  }

# Manhattan plot
plot_manhattan <- function(meta, phenotype) {
  keys <- unique(meta[, .(Region, CellType)]) #generate list of region x cell type combinations for given phenotype
  for (i in seq_len(nrow(keys))) {
    region_i <- keys$Region[i] 
    celltype_j <- keys$CellType[i]
    sub <- meta[Region == region_i & CellType == celltype_j]
    if (!nrow(sub)) next

    anno <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
    anno <- as.data.table(anno, keep.rownames = "CpG")

    sub_plot <- copy(sub) # copy the data to generate a "sub" plot for each celltype x region combination

    sub_plot[anno, on = "CpG", `:=`(
      chr = as.character(i.chr),
      bp = i.pos,
      UCSC_RefGene_Name = i.UCSC_RefGene_Name
      )]
    chrLengths <- sub_plot[, .(max_bp = max(bp, na.rm = TRUE)), by = chr][order(chr)]
    chrLengths[, bp_add := data.table::shift(cumsum(as.numeric(max_bp)), fill = 0, type = "lag")]
    chrLengths[, chr_center := bp_add + max_bp / 2]
  
    sub_plot[chrLengths, on = "chr", bp_add := i.bp_add]
    sub_plot[, `:=`(
      bpCummulative = bp + bp_add,
      log10_p = -log10(pval)
      )]

    # Set significance threshold
    sig_thr <- !is.na(sub_plot[[FDR]]) & sub_plot[[FDR]] < 0.05
    thr_line <- if (any(sig_thr)) -log10(max(sub_plot$pval[sig_mask], na.rm = TRUE)) else NA_real_

    # Generate the plot
    p <- ggplot(sub_plot, aes(bpCummulative, log10_p)) +
    geom_point(aes(color = factor(chr %% 2)), size = 0.6, alpha = 0.75, show.legend = FALSE) +
    scale_color_manual(values = c("0" = "orange3", "1" = "blue4")) +
    theme_cowplot() +
    scale_x_continuous(
    breaks = chrLengths$chr_center, 
    labels = c(1:22, "X", "Y")[1:nrow(chrLengths)], 
    expand = c(0.01, 0)
    )

    # Add significance threshold line
    if (is.finite(thr_line)) {
    p <- p + geom_hline(yintercept = thr_line, color = "red", linetype = "solid", linewidth = 1.2)
    }
    
    filepath <- c("") # path to save the figure
    ggsave(filepath, p, width = 10, height = 8)
    }
}

# Volcano plot
plot_volcano <- function(meta, phenotype) {
  keys <- unique(meta[, .(Region, CellType)])
  for (i in seq_len(nrow(keys))) {
    region_i <- keys$Region[i] 
    celltype_j <- keys$CellType[i]
    sub <- meta[Region == region_i & CellType == celltype_j & is.finite(get(p_col)) & get(p_col) > 0]
    if (!nrow(sub)) next

    sub[, neglog10p := -log10(get(p_col))]
    sub[, sig := (get(p_col) < 0.05)]
    
    xmax <- max(abs(sub[["beta_hat"]]), na.rm = TRUE)
    xmax <- max(xmax, 0.005) * 1.1 # we use 0.005 as threshold for age association, change if necessary
    sub[, beta := get("beta_hat")]
    sub[, category := fifelse(
    get(p_col) < 0.05 & beta >= 0.005, "Hyper_sig",
    fifelse(get(p_col) < 0.05 & beta <= -0.005, "Hypo_sig",
    fifelse(get(p_col) < 0.05, "Sig_small",
    "NS")))
    ]

    p <- ggplot(sub, aes(x = beta, y = neglog10p)) + 
    geom_point(aes(color = category), alpha = 0.7, size = 1.2) +
    scale_color_manual(values = c(
    "Hyper_sig" = "#298c8c",   
    "Hypo_sig"  = "#f2c45f",   
    "Sig_small" = "grey70",    
    "NS"        = "grey40"     
    )) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +\
    geom_vline(xintercept = c(-0.005, 0.005), linetype = "dashed", color = "black") +
    scale_x_continuous(
    limits = c(-0.1, 0.1),
    breaks = seq(-0.1, 0.1, by = 0.02)
    ) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 6)) +
    theme_cowplot()

    filepath <- c("") # path to save the figure
    ggsave(filepath, p, width = 8, height = 8)
   }
}

# running the plotting functions:
summary_dotplot(meta, phenotype)
plot_manhattan(meta, phenotype)
plot_volcano(meta, phenotype)


# ----------------------------------Gene ontology and pathway analysis--------------------------------------------------
enrichment_analysis <- function(meta, phenotype, collection = "GO"){
  keys <- unique(dt[, .(Region, CellType)])

  for (i in seq_len(nrow(keys))) {
    region_i <- keys$Region[i] 
    celltype_j <- keys$CellType[i]
    sub <- dt[Region == region_i & CellType == celltype_j]
    if (!nrow(sub)) next
    sig_cpg <- unique(sub[Sig == TRUE, CpG])
    all_cpg <- unique(sub$CpG)

    go_meth <- missMethyl::gometh(
      sig.cpg = sig_cpg,
      all.cpg = all_cpg,
      collection = collection
    )
    go_meth <- as.data.table(go_meth, keep.rownames = "ID")
    go_meth$Phenotype <- phenotype

    filepath <- c("") #file name and path to save GO results
    write.csv(go_meth, filepath)

    # Generate plots
    # calculate fold enrichment
    total_de <- sum(go_meth$DE, na.rm = TRUE)
    total_n  <- sum(go_meth$N, na.rm = TRUE)
    global_rate <- total_de / total_n
    go_meth[, Fold_Enrichment := (DE / N) / global_rate] 

    pcol <- if ("FDR" %in% names(go_meth)) "FDR" else "P.DE"
    go_meth <- go_meth[is.finite(get(pcol)) & get(pcol) > 0]
    go_meth[, neglog10 := -log10(get(pcol))]
    setorder(go_meth, -neglog10)
    
    plot_data <- go_meth[1:min(15, .N)]
    p <- ggplot(plot_data, aes(x = Fold_Enrichment, y = term_plot)) + 
    geom_segment(aes(x=0, xend = Fold_Enrichment, y = term_plot, yend = term_plot, color = neglog10), linewidth = 1) +
    geom_point(aes(size = DE, color = neglog10), alpha = 0.8) +
    scale_color_gradient(low = "navy", high = "firebrick3", name = expression(-log[10](pValue))) +
    scale_size_continuous(name = "Count (DE)") +
    theme_cowplot()

    filepath <- c("") #path to save the figure
    ggsave(filepath, p, width = 20, height = 10)
  }
}

# run the function itiratively for GO and KEGG pathway enrichment analyses
go_results <- enrichment_analysis(meta, phenotype, "GO")
kegg_results <- enrichment_analysis(meta, phenotype, "KEGG")

#---------------------------- Effect-size Rank Enrichment Analysis-------------------------------------------
# related to Fig. 4
# Evaluating cell-type enrichment of effect size for bulk tissue-level Braak CpGs (identified in Smith et al. 2021)
# This analysis is limited to phenotype = "Braak"
meta_dir_braak <- c("") # set path to directory containing meta analyzed tsv files for each Braak x region x celltype combination
region <- c("") # specify brain region for analysis
smith_list_filepath <- c("") # path to file containing list of CpGs from Smith et al 2021, for the given brain region
smith_cpgs <- fread(smith_list_filepath)

rank_enrichment_results <- list()
celltypes <- c("")

for (ct in celltypes) {
meta_braak_file <- file.path(meta_dir_braak, paste0("", ct, ".tsv.gz")) # make sure file name contains cell type label
dt <- fread(meta_braak_file)
dt[, .(CpG, chr, bp, beta_hat, se_hat, z, pval)]

  # perform wilcox rank sum test
  dt[, is_smith := CpG %in% smith_cpgs]
  w <- wilcox.test(
  abs(dt$z[dt$is_smith == 1]),
  abs(dt$z[dt$is_smith == 0]),
  alternative = "greater"
  )

  w_data <- data.table(
    n_smith_overlap = sum(dt$is_smith),
    p_wilcox = w$p.value,
    median_smith = median(abs(dt[is_smith == TRUE, z])),
    median_background = median(abs(dt[is_smith == FALSE, z])),
    direction_agreement = mean(sign(dt[is_smith == TRUE]$z) == sign(mean(dt[is_smith == FALSE]$z)))
    )
    
  # perform AUC perm-test
  observed_auc <- mean(rank(abs(dt$z))[dt$is_smith]) / nrow(dt)
  perm_auc <- replicate(10000, {
    perm <- sample(dt$is_smith)
    mean(rank(dt$z)[perm]) / nrow(dt)
  })

  auc_data <- data.table(auc = obs_auc, p_perm = mean(perm_auc >= obs_auc))

  rank_enrichment_results[[paste(ct)]] <- cbind(data.table(CellType = ct), w_data, auc_data)
 }

rank_enrichment_compiled <- rbindlist(rank_enrichment_results, fill = TRUE)
fwrite(rank_enrichment_compiled, file.path("")) # save results

# generate ridge plots after running all cell types
plot_braak_ridges <- function(all_dt, smith_list, region) {
  all_dt[, Smith := CpG %in% smith_list]
  all_dt[, CellType := factor(CellType, levels = c("Astro", "Oligo_OPC", "Endo", "Neuron"))]
  
  ggplot(all_dt, aes(x = abs(z), y = CellType)) +
  geom_density_ridges(
    data = subset(all_dt, Smith == FALSE), fill = "grey80", alpha = 0.8, scale = 1.3) +
  geom_density_ridges(data = subset(all_dt, Smith == TRUE), aes(fill = CellType), alpha = 0.9, scale = 1.3) +
    scale_fill_manual(values=c(
    "Neuron" = "#023743",
    "Astro" = "#FED789",
    "Oligo_OPC" = "#72874E",
    "Endo" = "#476F84")) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    theme_classic()
}

# leading edge directional concordance analysis
# this was done only between neurons and astrocytes in the PFC, in association with Braak stage for Smith CpGs
neuron <- fread("") # path to meta-analyzed Braak association file for PFC neurons
astro <- fread("") # path to meta-analyzed Braak association file for PFC astrocytes
smith <- fread("") # path to file containing braak-association ES/SE data for the PFC from Smith et al. 2021

# calculate z-values
neuron[, z_cell := beta_hat / se_hat]
astro[, z_cell := beta_hat / se_hat]
smith[, z_smith := beta_hat / se_hat]
neuron[, CellType := "Neuron"]
astro[, CellType := "Astro"]

dt <- rbind(neuron, astro)
dt <- merge(dt, smith[, .(CpG, z_smith)], by = "CpG")
dt <- dt[is.finite(z_cell) & is.finite(z_smith)]

# leading edge function
leading_edge <- function(d, step = 5) {
  d <- d[order(-abs(z_smith))]
  n <- nrow(d)
  ks <- seq(step, n, by = step)
  
  data.table(
    frac = ks / n,
    concordance = sapply(ks, function(k) {
      mean(sign(d$z_cell[1:k]) == sign(d$z_smith[1:k]))
    })
  )
}

# compute leading edge for neuron and astrocytes
le_neuron  <- leading_edge(dt[CellType == "Neuron"])
le_astro <- leading_edge(dt[CellType == "Astro"])

# Plot leading edge concordance
p <- ggplot() + 
geom_smooth(data = le_astro, aes(x = frac, y = concordance), method = "loess", span = 0.2) +
geom_smooth(data = le_neuron, aes(x = frac, y = concordance), method = "loess", span = 0.2) +
geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey40") #reference line +
scale_x_continuous(
    labels = scales::percent_format(accuracy = 1),
    breaks = c(0.1, 0.25, 0.5, 0.75, 1)
  ) +
theme_classic()

filepath <- c("")
ggsave(filepath, p, width=10, height=6)

# comparison between neuron vs astrocyte concordance for top 50 ranked Smith-CpGs by |z|
astro_obs <- mean(sign(le_astro$z_cell[1:50]) == sign(le_astro$z_smith[1:50]))
neur_obs  <- mean(sign(le_neur$z_cell[1:50]) == sign(le_neur$z_smith[1:50]))
delta_obs <- astro_obs - neur_obs

perm_delta <- replicate(1000, {
    idx <- sample(seq_len(nrow(le_astro)))
    astro_perm <- mean(sign(le_astro$z_cell[idx][1:50]) == sign(le_astro$z_smith[1:50]))
    neur_perm <- mean(sign(le_neur$z_cell[idx][1:50]) == sign(le_neur$z_smith[1:50]))
    astro_perm - neur_perm
  })


#----------------------------------------- Age-to-disease trajectory analysis--------------------------------------------------
# related to Fig. 5
# to compute directional and ES magnitude-related change in age-associated DMPs between age- and AD-associations across cell types x regions
# run once per region x celltype combination

age <- fread("") # load age-association data for cell type x region
ad <- fread("") # load AD-association data for the same cell type x region

age_sig <- age[pval < 0.05]
dt <- merge(age_sig[, .(CpG, beta_age = beta_hat)], 
                ad[, .(CpG, beta_ad = beta_hat)], by = "CpG")
    
# Remove any non-finite values for regression
dt <- dt[is.finite(beta_age) & is.finite(beta_ad)]

# Perform linear regression analysis for AD ~ Age
# The slope (beta) measures ES acceleration/divergence
model <- lm(beta_ad ~ beta_age, data = dt)
slope_val <- coef(model)[2]
slope_p    <- summary(model)$coefficients["beta_age", "Pr(>|t|)"]

# compute correlation
cor_test   <- cor.test(dt$beta_age, dt$beta_ad)
r_val      <- cor_test$estimate
cor_p      <- cor_test$p.value

# scatter plot
p_scatter <- ggplot(dt, aes(x = beta_age, y = beta_ad)) +
      geom_point(size = 1) + 
      geom_smooth(method = "lm", color = "black", linetype = "solid") +
      theme_cowplot() 

filepath <- c("")
ggsave(filepath, p_scatter, width = 6, height = 6)
