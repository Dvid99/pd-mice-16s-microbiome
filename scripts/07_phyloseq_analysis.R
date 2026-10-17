#!/usr/bin/env Rscript
# Step 7 — Community analysis in R with phyloseq / vegan / lme4:
#   1. build a phyloseq object from the exported QIIME 2 results
#   2. subset Clostridiales and plot family abundances per sample
#   3. Shannon diversity by donor status, tested correctly for repeated measures
#   4. sample with the lowest Shannon index
#   5. phylogenetic tree of Betaproteobacteria genera, coloured by group
#   6. PCoA of unweighted UniFrac + PERMANOVA on independent samples
#
# Usage: Rscript scripts/07_phyloseq_analysis.R

suppressPackageStartupMessages({
  library(phyloseq)
  library(vegan)
  library(ape)
  library(ggplot2)
  library(lmerTest)
})

in_dir <- "results/06_export"
out_dir <- "results/07_phyloseq"
fig_dir <- file.path(out_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(42)

# ---- 1. phyloseq object ------------------------------------------------------
otu <- read.delim(file.path(in_dir, "feature_table.tsv"), skip = 1, row.names = 1,
                  check.names = FALSE)
otu <- as.matrix(otu)

tax_raw <- read.delim(file.path(in_dir, "taxonomy.tsv"), row.names = 1)
ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
tax <- t(sapply(strsplit(tax_raw$Taxon, ";\\s*"), function(x) {
  x <- sub("^[a-z]__", "", x)                 # "k__Bacteria" -> "Bacteria"
  x[x == ""] <- NA
  length(x) <- length(ranks)
  x
}))
dimnames(tax) <- list(rownames(tax_raw), ranks)

meta <- read.delim(file.path(in_dir, "metadata.tsv"), check.names = FALSE,
                   comment.char = "", stringsAsFactors = FALSE)
meta <- meta[!grepl("^#q2:", meta[[1]]), ]   # drop the QIIME 2 type row if present
rownames(meta) <- meta[[1]]
meta$days_post_transplant <- as.numeric(meta$days_post_transplant)

ps <- phyloseq(otu_table(otu, taxa_are_rows = TRUE), tax_table(tax),
               sample_data(meta), read_tree(file.path(in_dir, "tree.nwk")))
ps <- prune_taxa(taxa_sums(ps) > 0, ps)
cat("phyloseq object:", ntaxa(ps), "ASVs x", nsamples(ps), "samples\n")
cat("Reads per sample: min", min(sample_sums(ps)), "| median", median(sample_sums(ps)),
    "| max", max(sample_sums(ps)), "\n")

# Rarefy once to an even depth (same as QIIME 2 core metrics) for diversity
depth <- 2000
ps_rare <- rarefy_even_depth(ps, sample.size = depth, rngseed = 42, verbose = FALSE)
dropped <- setdiff(sample_names(ps), sample_names(ps_rare))
cat("Rarefied to", depth, "reads;", length(dropped), "sample(s) dropped:",
    paste(dropped, collapse = ", "), "\n\n")

# ---- 2. Clostridiales: families per sample -----------------------------------
ps_clos <- subset_taxa(ps, Order == "Clostridiales")
cat("Clostridiales: ", ntaxa(ps_clos), "of", ntaxa(ps), "ASVs\n")
ps_fam <- tax_glom(ps_clos, taxrank = "Family", NArm = FALSE)
ps_fam_rel <- transform_sample_counts(ps_fam, function(x) 100 * x / sum(x))
p <- plot_bar(ps_fam_rel, fill = "Family") +
  geom_bar(stat = "identity", colour = NA) +
  labs(y = "Relative abundance within Clostridiales (%)", x = NULL,
       title = "Clostridiales families per sample") +
  facet_grid(~ donor_status, scales = "free_x", space = "free_x") +
  theme_minimal(base_size = 9) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 5), legend.position = "bottom")
ggsave(file.path(fig_dir, "clostridiales_families.png"), p, width = 11, height = 6, dpi = 200, bg = "white")

# ---- 3. Shannon diversity by donor status ------------------------------------
alpha <- data.frame(sample_data(ps_rare))
alpha$shannon <- estimate_richness(ps_rare, measures = "Shannon")$Shannon   # natural log, as microbiome::alpha
write.table(alpha[, c("mouse_id", "cage_id", "donor_status", "genotype",
                      "days_post_transplant", "shannon")],
            file.path(out_dir, "shannon_by_sample.tsv"), sep = "\t", quote = FALSE, col.names = NA)

means <- tapply(alpha$shannon, alpha$donor_status, mean)
cat("\nMean Shannon by donor status:\n"); print(round(means, 3))

# (a) Naive test: every sample treated as independent. 48 samples come from only
#     12 mice (4 time points each) housed in 6 cages -> pseudoreplication.
naive <- wilcox.test(shannon ~ donor_status, data = alpha, exact = FALSE)
# (b) One value per mouse (mean over time points), n = 12 independent units
per_mouse <- aggregate(shannon ~ mouse_id + donor_status + cage_id, data = alpha, FUN = mean)
mouse_test <- wilcox.test(shannon ~ donor_status, data = per_mouse, exact = FALSE)
# (c) Linear mixed model: donor status + time, random intercept per mouse nested in cage
lmm <- lmer(shannon ~ donor_status + days_post_transplant + (1 | cage_id / mouse_id),
            data = alpha)
lmm_tab <- summary(lmm)$coefficients
lmm_p <- lmm_tab[grep("^donor_status", rownames(lmm_tab)), "Pr(>|t|)"]

tests <- data.frame(
  test = c("Wilcoxon, all samples (naive)", "Wilcoxon, per-mouse means",
           "LMM: shannon ~ donor + time + (1|cage/mouse)"),
  n_units = c(nrow(alpha), nrow(per_mouse), nrow(alpha)),
  p_value = signif(c(naive$p.value, mouse_test$p.value, lmm_p), 3))
write.table(tests, file.path(out_dir, "shannon_donor_tests.tsv"), sep = "\t", quote = FALSE,
            row.names = FALSE)
cat("\n"); print(tests, row.names = FALSE)

p <- ggplot(alpha, aes(donor_status, shannon)) +
  geom_boxplot(outlier.shape = NA, width = 0.5, fill = "grey92") +
  geom_jitter(aes(colour = factor(days_post_transplant)), width = 0.12, size = 1.8) +
  geom_point(data = per_mouse, shape = 23, size = 3, fill = "black") +
  scale_colour_viridis_d(name = "Day") +
  labs(x = "Donor", y = "Shannon index",
       title = "Shannon diversity by donor status",
       subtitle = sprintf("Per-mouse means (diamonds), n = %d mice | LMM p = %.3g",
                          nrow(per_mouse), lmm_p)) +
  theme_classic()
ggsave(file.path(fig_dir, "shannon_by_donor.png"), p, width = 6, height = 5, dpi = 200, bg = "white")

# ---- 4. Lowest Shannon index -------------------------------------------------
lo <- alpha[which.min(alpha$shannon), ]
cat(sprintf("\nLowest Shannon: %s (%.3f) — %s, %s donor, day %d\n",
            rownames(lo), lo$shannon, lo$genotype, lo$donor_status, lo$days_post_transplant))

# ---- 5. Betaproteobacteria tree ----------------------------------------------
ps_beta <- subset_taxa(ps, Class == "Betaproteobacteria")
cat("\nBetaproteobacteria:", ntaxa(ps_beta), "ASVs in families:",
    paste(sort(unique(na.omit(tax_table(ps_beta)[, "Family"]))), collapse = ", "), "\n")
if (ntaxa(ps_beta) > 1) {
  p <- plot_tree(ps_beta, ladderize = "left", color = "genotype_and_donor_status",
                 shape = "Family", label.tips = "Genus", size = "abundance",
                 text.size = 3.5, base.spacing = 0.03) +
    labs(title = "Betaproteobacteria: genera and samples by group") +
    theme(legend.position = "right")
  ggsave(file.path(fig_dir, "betaproteobacteria_tree.png"), p, width = 10, height = 6, dpi = 200, bg = "white")
}

# ---- 6. Beta diversity -------------------------------------------------------
uu <- as.matrix(read.delim(file.path(in_dir, "unweighted_unifrac.tsv"), row.names = 1,
                           check.names = FALSE))
keep <- intersect(rownames(uu), rownames(meta))
uu <- as.dist(uu[keep, keep])
pc <- cmdscale(uu, k = 2, eig = TRUE)
var_expl <- round(100 * pc$eig[1:2] / sum(pc$eig[pc$eig > 0]), 1)
pcoa <- data.frame(PC1 = pc$points[, 1], PC2 = pc$points[, 2], meta[keep, ])
p <- ggplot(pcoa, aes(PC1, PC2, colour = donor_status, shape = genotype)) +
  geom_point(size = 2.5, alpha = 0.85) +
  labs(x = sprintf("PCoA 1 (%.1f%%)", var_expl[1]), y = sprintf("PCoA 2 (%.1f%%)", var_expl[2]),
       title = "Unweighted UniFrac", colour = "Donor", shape = "Genotype") +
  theme_classic()
ggsave(file.path(fig_dir, "pcoa_unweighted_unifrac.png"), p, width = 6.5, height = 5, dpi = 200, bg = "white")

# PERMANOVA on one sample per mouse (last time point) to avoid pseudoreplication
last <- do.call(rbind, lapply(split(meta[keep, ], meta[keep, "mouse_id"]), function(d)
  d[which.max(d$days_post_transplant), ]))
rownames(last) <- last[[1]]                 # rbind renames rows; restore sample IDs
d_last <- as.dist(as.matrix(uu)[rownames(last), rownames(last)])
perm <- adonis2(d_last ~ donor_status + genotype, data = last, permutations = 9999, by = "margin")
cat("\nPERMANOVA (unweighted UniFrac, last time point per mouse, n =", nrow(last), "):\n")
print(perm)
capture.output(perm, file = file.path(out_dir, "permanova_last_timepoint.txt"))

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
cat("\nDone. Results in", out_dir, "\n")
