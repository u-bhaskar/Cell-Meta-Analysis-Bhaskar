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
meta_files <- list.files(meta, pattern = "\\.tsv\\.gz$", full.names = TRUE)
file_list <- data.table(path = files, file = basename(files))
meta <- rbindlist(lapply(file_list$path, fread, use.names = TRUE, fill = TRUE))
















for (sig in c("FDR","Bonf","Nominal")) {
  logmsg("[%s] Summary bars (per region) ...", sig);     plot_summary_bars_per_region(meta, phenotype, sig = sig)
  logmsg("[%s] Summary DotPlot (FACET ALL) ...", sig); plot_glia_summary_dotplot(meta, phenotype, sig = sig)
  logmsg("[%s] Manhattan ...", sig);                      plot_manhattan_all(meta, phenotype, sig = sig)
  logmsg("[%s] Volcano ...", sig);                        plot_volcano_all(meta, phenotype, sig=sig)
  logmsg("[%s] Context bars (Island & Genic) ...", sig);  plot_context_bars_all(meta, phenotype, sig = sig)
  logmsg("[%s] CVI heatmap ...", sig);                    plot_cvi(meta, phenotype, sig = sig)
  logmsg("[%s | %s] Running gometh...", phenotype, sig)

go_res   <- run_gometh(meta, phenotype, sig, "GO")
kegg_res <- run_gometh(meta, phenotype, sig, "KEGG")
}

logmsg("Effect-size correlation across regions (Panel C) ...")
plot_region_cor_heatmaps(meta, phenotype)

logmsg("DONE plotting for %s", phenotype)

#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
  library(ggplot2)
  library(cowplot)
  library(ggrepel)
  library(minfi)
  library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  suppressWarnings(suppressMessages(requireNamespace("missMethyl", quietly = TRUE)))
  suppressWarnings(suppressMessages(requireNamespace("GO.db", quietly = TRUE)))
  suppressWarnings(suppressMessages(requireNamespace("AnnotationDbi", quietly = TRUE)))
  suppressWarnings(suppressMessages(requireNamespace("KEGGREST", quietly = TRUE)))
})

# -------------------- Config --------------------
cfg <- list(
  res_dir        = "03_MetaAnalysis/02_Results_Region_Stratified",
  fig_dir        = "03_MetaAnalysis/03_Figures_Region_Stratified",
  log_dir        = "03_MetaAnalysis/04_Logs",
  alpha          = 0.05,
  volcano_es_thr = 0.01,  # absolute beta cut
  volcano_p_thr  = 0.05, # p-value threshold line (you can set to Bonf/FDR converted to p)
  top_labels     = 20,       # label up to N points on Manhattan
  gometh_min_sig = 100,      # minimum sig CpGs for GO/KEGG enrichment
  gometh_collections = c("GO","KEGG"),
  seed           = 123
)
dir.create(cfg$fig_dir, recursive = TRUE, showWarnings = FALSE)

logmsg <- function(...) cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "|", sprintf(...), "\n")

# Put this below cfg and ggsave_tiff definitions
use_common_theme <- function(base_size = 18, base_family = "") {
  theme_set(
    theme_cowplot(font_size = base_size, font_family = base_family) +
      theme(
        text = element_text(family = "Helvetica"),
        plot.title   = element_text(hjust=0.5, face="bold", size=24),
        plot.subtitle= element_text(hjust=0.5, face="italic", size=20),
        axis.title   = element_text(face="bold", size=20),
        axis.text    = element_text(face="bold", size=20),
        legend.title = element_text(face="bold", size = 18),
        legend.text  = element_text(face="bold", size = 18),
        strip.text   = element_text(face="bold", size=20),
        panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA)
      )
  )
}
use_common_theme()

# TIFF saver (bold-friendly)
ggsave_tiff <- function(filename, plot, width=10, height=6, dpi=600) {
  ggsave(
    filename,
    plot = plot,
    width = width,
    height = height,
    device = "tiff",
    bg = "white",
    dpi = dpi,
    units = "in",
    compression = "lzw"
  )
}

# -------------------- IO helpers --------------------

# -------------------- Summary bars --------------------

plot_glia_summary_dotplot <- function(dt, phenotype, sig = c("FDR", "Bonf", "Nominal")) {
  sig_type <- match.arg(sig)
  flag <- switch(sig,
  "FDR" = "Sig_FDR",
  "Bonf" = "Sig_Bonf",
  "Nominal" = "Sig_Nominal"
)
  
  # Aggregate counts of significant DMCs
  # We count hyper and hypo separately to show directionality via color
  summary_dt <- dt[get(flag) == TRUE, .(
    Total = .N,
    Hyper = sum(beta_hat > 0, na.rm = TRUE),
    Hypo  = sum(beta_hat < 0, na.rm = TRUE)
  ), by = .(Region, CellType)]
  
  if (nrow(summary_dt) == 0) {
    logmsg("No significant results found for %s (%s). Skipping dotplot.", phenotype, sig_type)
    return(NULL)
  }
  
  # Calculate a 'Direction' metric: 1 if mostly hyper, -1 if mostly hypo
  summary_dt[, Direction := (Hyper - Hypo) / Total]
  
  p <- ggplot(summary_dt, aes(x = CellType, y = Region)) +
    # The 'size' of the dot represents the total number of significant sites
    # The 'fill' color represents the hyper/hypo balance
    geom_point(aes(size = Total, fill = Direction), shape = 21, stroke = 1) +
    scale_size_continuous(range = c(2, 18), name = "N Sig CpGs",
    breaks = seq(0, max(summary_dt$Total, na.rm = TRUE), length.out = 5) %>% round() %>% unique()) +
    scale_fill_gradient2(low = "#f2c45f", mid = "white", high = "#298c8c", 
                         midpoint = 0, name = "Proportion") +
    theme_minimal() +
    theme_cowplot() +
    labs(
      title = sprintf("Differential Methylation Burden: %s", phenotype),
      subtitle = sprintf("Comparison of %s-significant sites across cell types and regions", sig_type),
      x = "Cell Type", y = "Brain Region"
    ) +
    theme(plot.title = element_text(face="bold", size=24),
      axis.text.x = element_text(face = "bold", size = 20),
      axis.text.y = element_text(face = "bold", size = 20),
      axis.title  = element_text(face="bold", size=20, margin = margin(t = 10, unit = "pt")),
          legend.title= element_text(face="bold"),
          legend.text = element_text(face="bold"),
          panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
      panel.grid.major = element_line(color = "grey90")
    )
  
  fn <- file.path(cfg$fig_dir, sprintf("Summary_Dotplot_%s_%s.tiff", phenotype, sig_type))
  ggsave_tiff(fn, p, width = 10, height = 8)
  logmsg("Glia dominance dotplot saved to: %s", fn)
}


# -------------------- Manhattan (robust chr mapping + subtitle with Nsig) --------------------
manhattan_one <- function(sub, phenotype, sig = c("FDR","Bonf", "Nominal"), region, celltype) {
  sub_plot <- copy(sub)
  
  sig <- match.arg(sig)
  sig_col <- switch(sig,
  "FDR" = "FDR",
  "Bonf" = "Bonf",
  "Nominal" = "pval"
)
  anno <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)


  anno <- as.data.table(anno, keep.rownames = "CpG")



  # Filter valid genomic rows
  sub_plot <- sub_plot[!is.na(chr) & !is.na(bp) & !is.na(pval)]

sub_plot[anno, on = "CpG", `:=`(
chr = as.character(i.chr),
bp = i.pos,
UCSC_RefGene_Name = i.UCSC_RefGene_Name
)]
  

sub_plot[, chr := gsub("chr", "", chr)]
  sub_plot[, chr := ifelse(chr == "X", "23", ifelse(chr == "Y", "24", chr))]
  sub_plot[, chr := as.numeric(chr)]
  sub_plot <- sub_plot[!is.na(chr)] 
  
  chrLengths <- sub_plot[, .(max_bp = max(bp, na.rm = TRUE)), by = chr][order(chr)]
  
  chrLengths[, bp_add := data.table::shift(cumsum(as.numeric(max_bp)), fill = 0, type = "lag")]
  
  chrLengths[, chr_center := bp_add + max_bp / 2]
  
  sub_plot[chrLengths, on = "chr", bp_add := i.bp_add]
  sub_plot[, `:=`(
  bpCummulative = bp + bp_add,
  log10_p = -log10(pval)
  )]
  

  
  # 4. Significance Threshold
  sig_mask <- !is.na(sub_plot[[sig_col]]) & sub_plot[[sig_col]] < 0.05
  thr_line <- if (any(sig_mask)) -log10(max(sub_plot$pval[sig_mask], na.rm = TRUE)) else NA_real_
  
  # 5. Top Labels
  lab_dt <- sub_plot[sig_mask][order(get(sig_col), pval)]
  if (nrow(lab_dt) > 0) {
    n_labels <- if(exists("cfg")) cfg$top_labels else 15
    lab_dt <- lab_dt[1:min(n_labels, .N)]
lab_dt[, Gene := sapply(strsplit(as.character(UCSC_RefGene_Name), ";"),
function(x) if(length(x) > 0) x[1] else "")]
lab_dt[, PlotLabel := ifelse(Gene != "", 
sprintf("%s\n(%s)", CpG, Gene), CpG)]
}
 
  
  # 6. Plotting
  p <- ggplot(sub_plot, aes(bpCummulative, log10_p)) +
    geom_point(aes(color = factor(chr %% 2)), size = 0.6, alpha = 0.75, show.legend = FALSE) +
    scale_color_manual(values = c("0" = "orange3", "1" = "blue4")) +
    theme_cowplot() +
    labs(
      title = sprintf("%s (%s) | %s | %s", phenotype, sig, region, celltype),
      subtitle = sprintf("Significant (%s < %.2g): %d CpGs", sig, 0.05, sum(sig_mask)),
      x = "Chromosome",
      y = expression(bold(-log[10](p)))
    ) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 24),
          plot.subtitle = element_text(hjust = 0.5, face = "italic", size = 20),
          panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA))
  
  # X-axis scale using numeric midpoints
  p <- p + scale_x_continuous(
    breaks = chrLengths$chr_center, 
    labels = c(1:22, "X", "Y")[1:nrow(chrLengths)], 
    expand = c(0.01, 0)
  )
  
  
if (is.finite(thr_line)) {
  p <- p + geom_hline(yintercept = thr_line, color = "red", linetype = "solid", linewidth = 1.2)
}

  
  if (nrow(lab_dt)) {
    p <- p + ggrepel::geom_label_repel(
      data = lab_dt, aes(label = PlotLabel), 
      size = 3, box.padding = 0.5, max.overlaps = 20, fontface = "bold"
    )
  }
  
  # 7. Safe Saving
  if(!dir.exists(cfg$fig_dir)) dir.create(cfg$fig_dir, recursive = TRUE)
  fn <- file.path(cfg$fig_dir, sprintf("Manhattan_%s_%s_%s_%s.tiff", phenotype, sig, region, celltype))
  
  # Use standard ggsave if ggsave_tiff isn't defined
  ggsave(fn, p, width = 12, height = 7, compression = "lzw", dpi = 300)
}

plot_manhattan_all <- function(dt, phenotype, sig = c("FDR","Bonf", "Nominal")) {
  sig <- match.arg(sig)
  keys <- unique(dt[, .(Region, CellType)])
  for (i in seq_len(nrow(keys))) {
    rr <- keys$Region[i]; ct <- keys$CellType[i]
    sub <- dt[Region == rr & CellType == ct]
    if (!nrow(sub)) next
    manhattan_one(sub, phenotype, sig = sig, region = rr, celltype = ct)
  }
}

# -------------------- Volcano --------------------
# -------- Volcano with thresholds & axis breaks --------
# dt: annotated meta (feA) for a single Region×CellType (or whole set with filtering)
plot_volcano_one <- function(dt, phenotype, region, celltype,
                             sig = c("FDR","Bonf", "Nominal"), beta_col = "beta_hat",
                             p_thresh = 0.05, effect_thresh = 0.005) {

  sig <- match.arg(sig)
  p_col <- switch(sig,
  "FDR" = "FDR",
  "Bonf" = "Bonf",
  "Nominal" = "pval"
)

  sub <- dt[Region == region & CellType == celltype & is.finite(get(p_col)) & get(p_col) > 0]
  if (!nrow(sub)) return(invisible(NULL))

  sub[, neglog10p := -log10(get(p_col))]
  sub[, sig := (get(p_col) < p_thresh)]
  xmax <- max(abs(sub[[beta_col]]), na.rm = TRUE)
  xmax <- max(xmax, effect_thresh) * 1.1
  sub[, beta := get(beta_col)]
  sub[, beta := get(beta_col)]
sub[, category := fifelse(
  get(p_col) < p_thresh & beta >= effect_thresh, "Hyper_sig",
  fifelse(get(p_col) < p_thresh & beta <= -effect_thresh, "Hypo_sig",
  fifelse(get(p_col) < p_thresh, "Sig_small",
  "NS")))
]
  
  p <- ggplot(sub, aes(x = beta, y = neglog10p)) +
  geom_point(aes(color = category), alpha = 0.7, size = 1.2) +

  scale_color_manual(values = c(
    "Hyper_sig" = "#298c8c",   # your hyper color
    "Hypo_sig"  = "#f2c45f",   # your hypo color
    "Sig_small" = "grey70",    # sig but small effect
    "NS"        = "grey40"     # non-sig
  )) +

  geom_hline(yintercept = -log10(p_thresh), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-effect_thresh, effect_thresh), linetype = "dashed", color = "black") +

  scale_x_continuous(
    limits = c(-0.1, 0.1),
    breaks = seq(-0.1, 0.1, by = 0.02)
  ) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 6)) +
    labs(
      title = sprintf("Volcano — %s : %s / %s\n %s<%.2g, |Δβ|≥%.3f", phenotype, region, celltype, sig, p_thresh, effect_thresh),
      x = expression(beta),
      y = expression(-log[10](p))
    ) +
    theme_cowplot() +
    theme(
      plot.title  = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.title  = element_text(face = "bold"),
      axis.text   = element_text(face = "bold"),
      panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA)
    )

  fn <- file.path(cfg$fig_dir, sprintf("Volcano_%s_%s_%s_%s.tiff", phenotype, sig, region, celltype))
  ggsave_tiff(fn, p, width = 6.5, height = 6.0)
}

# helper to loop all Region×CellType combos
plot_volcano_all <- function(dt, phenotype, sig, p_thresh = 0.05, effect_thresh = 0.005) {
  keys <- unique(dt[, .(Region, CellType)])
  for (i in seq_len(nrow(keys))) {
    rr <- keys$Region[i]; ct <- keys$CellType[i]
    plot_volcano_one(dt, phenotype, rr, ct, sig, p_thresh = p_thresh, effect_thresh = effect_thresh)
  }
}


# -------- GO/KEGG term-name plotting --------
# df: gometh() result; must contain columns 'P.DE' (raw p) and either 'Term' or 'Pathway'
# db_label: "GO" or "KEGG" for title
plot_gometh <- function(df, phenotype, sig = c("FDR","Bonf", "Nominal"), region, celltype, db_label = "GO", top_n = 15L) {
  if (is.null(df) || !nrow(df)) return(invisible(NULL))

dd <- as.data.table(df)
  # detect term text column
 
  if (db_label == "GO") {
  if (!all(c("ID", "TERM") %in% names(dd))) {
    stop("GO results must contain 'Term' (ID) and 'TERM' (name)")
  }
  dd[, TermLabel := TERM]   # <-- THIS is what you want plotted

} else if (db_label == "KEGG") {
  if (!all(c("ID", "Description") %in% names(dd))) {
    stop("KEGG results must contain 'Term' (ID) and 'Description' (name)")
  }
  dd[, TermLabel := Description]

} else {
  stop("db_label must be 'GO' or 'KEGG'")
}

# 2. Calculate Fold Enrichment (The "Enrichment Score")
  # We use the full results 'df' to get global totals for the background rate
  total_de <- sum(dd$DE, na.rm = TRUE)
  total_n  <- sum(dd$N, na.rm = TRUE)
  global_rate <- total_de / total_n
  dd[, Fold_Enrichment := (DE / N) / global_rate]
 
 pcol <- if ("P.DE" %in% names(dd)) "P.DE" else if ("FDR" %in% names(dd)) "FDR" else NULL

if (is.null(pcol)) stop("No p-value column found")

dd <- dd[is.finite(get(pcol)) & get(pcol) > 0]
dd[, neglog10 := -log10(get(pcol))]

  # allow both raw P and adjusted columns in gometh output
  setorder(dd, -neglog10)
  plot_data <- dd[1:min(top_n, .N)]

  # wrap long labels
  plot_data[, term_plot := stringr::str_trunc(TermLabel, 60)]
  plot_data[, term_plot := factor(term_plot, levels = rev(term_plot))]

  title_txt <- sprintf("%s Enrichment — %s (%s) | %s / %s",
                     db_label, phenotype, sig, region, celltype)
  p <- ggplot(plot_data, aes(x = Fold_Enrichment, y = term_plot)) + 
  geom_segment(aes(x=0, xend = Fold_Enrichment, y = term_plot, yend = term_plot, color = neglog10), linewidth = 1) +
  geom_point(aes(size = DE, color = neglog10), alpha = 0.8) +
  scale_color_gradient(low = "navy", high = "firebrick3", name = expression(-log[10](pValue))) +
  scale_size_continuous(name = "Count (DE)") +
    #coord_flip() +
    labs(title = title_txt, x = "Fold Enrichment", y = NULL) +
    theme_cowplot() +
    theme(
      plot.title  = element_text(hjust = 0.5, face = "bold", size = 24),
      axis.title  = element_text(face = "bold", size = 20),
      axis.text   = element_text(face = "bold", size = 20),
      legend.text = element_text(face = "bold", size = 18),
      #panel.grid.major = element_line(color = "grey90"),
      panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA)
    )

  fn <- file.path(cfg$fig_dir, sprintf("Enrichment_%s_%s_%s_%s_%s.tiff", db_label, phenotype, sig, region, celltype))
  ggsave_tiff(fn, p, width = 20, height = 8.5)
}

run_gometh <- function(dt, phenotype, sig, collection = "GO") {

  require(missMethyl)

  flag <- switch(sig,
  "FDR"     = "Sig_FDR",
  "Bonf"    = "Sig_Bonf",
  "Nominal" = "Sig_Nominal"
)

  keys <- unique(dt[, .(Region, CellType)])
  for (i in seq_len(nrow(keys))) {
    rr <- keys$Region[i]; ct <- keys$CellType[i]
    sub <- dt[Region == rr & CellType == ct]
    if (!nrow(sub)) next
    sig_cpg <- unique(sub[get(flag) == TRUE, CpG])
  all_cpg <- unique(sub$CpG)

  if (length(sig_cpg) < 10) {
    message(sprintf("[%s | %s] Not enough CpGs for %s", phenotype, sig, collection))
    return(NULL)
  }

  gom <- missMethyl::gometh(
    sig.cpg = sig_cpg,
    all.cpg = all_cpg,
    collection = collection
  )
  gom <- as.data.table(gom, keep.rownames = "ID")
  gom$Phenotype <- phenotype
  gom$SigType <- sig

  gom <- as.data.frame(gom)
  #gom <- add_ids_if_possible(gom, collection)
  fn <- file.path(cfg$res_dir, sprintf("Enrichment_%s_%s_%s_%s_%s.csv", collection, phenotype, sig, rr, ct))
  write.csv(gom, fn)
  plot_gometh(gom, phenotype, sig, rr, ct, collection)      
  }

}




