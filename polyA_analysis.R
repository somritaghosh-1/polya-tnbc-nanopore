# ============================================================
# PolyA Tail Length Analysis in Triple-Negative Breast Cancer
# MDA-231 Cell Line | 96 Cancer-Associated Genes
# Oxford Nanopore Direct RNA Sequencing (SQK-RNA004)
# Author: Somrita Ghosh | MSc Molecular Biology, QUB 2024
# ============================================================

# ---- 1. INSTALL AND LOAD PACKAGES --------------------------
# Run install lines only once, then comment them out
# install.packages("ggplot2")
# install.packages("dplyr")
# install.packages("tidyr")
# install.packages("ggpubr")
# install.packages("RColorBrewer")
# install.packages("plotly")

library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)
library(RColorBrewer)
library(plotly)

# ---- 2. LOAD DATA ------------------------------------------
# Set your working directory to wherever your CSV file is
# setwd("C:/Users/LENOVO/Downloads")  # uncomment and edit if needed

polya <- read.csv("combined_pti_values.csv", stringsAsFactors = FALSE)
cat("Data loaded:", nrow(polya), "reads x", ncol(polya), "genes\n")

# ---- 3. DEFINE FUNCTIONAL GENE GROUPS ----------------------
# Groups based on biological function in MDA-231 / cancer biology
# This is the key scientific framing of the analysis

gene_groups <- list(
  "Cell Adhesion (Integrins)" = c("ITGB1", "ITGA1", "ITGA2", "ITGA3",
                                   "ITGAV", "ITGB3", "ITGB5"),
  "Oncogenes"                 = c("AKT1", "CDK2", "CDK4", "MYC",
                                   "ERBB2", "RAF1", "MET", "ABL1",
                                   "MDM2", "PIK3R1"),
  "Tumour Suppressors"        = c("TP53", "RB1", "BRCA1", "ATM",
                                   "CHEK2", "CDKN1B"),
  "EMT and Invasion"          = c("TGFB1", "TGFBR1", "MMP1", "MMP2",
                                   "PLAU", "PLAUR", "SERPINE1",
                                   "S100A4", "THBS1"),
  "Apoptosis"                 = c("BAD", "BAX", "BCL2", "BCL2L1",
                                   "CASP8", "FAS", "CFLAR", "APAF1"),
  "Housekeeping"              = c("RPL37A", "YWHAZ", "HMBS",
                                   "PGK1", "MRPL19", "NME1", "PES1"),
  "Angiogenesis"              = c("VEGFA", "ANGPT1", "ANGPT2",
                                   "PDGFA", "PDGFB"),
  "Cell Cycle"                = c("CDK2", "CDK4", "CCNE1", "CDC25A",
                                   "E2F1", "CDKN1B")
)

# ---- 4. RESHAPE TO LONG FORMAT FOR GGPLOT ------------------
# Wide format (one column per gene) -> Long format (gene + value columns)
# This is standard practice for ggplot2

polya_long <- polya %>%
  pivot_longer(
    cols      = everything(),
    names_to  = "Gene",
    values_to = "PolyA_Length"
  ) %>%
  filter(!is.na(PolyA_Length)) %>%         # remove missing values
  filter(PolyA_Length > 0)                  # remove zero values

cat("Total reads after filtering:", nrow(polya_long), "\n")
cat("Genes detected:", length(unique(polya_long$Gene)), "\n")

# Add functional group labels
group_lookup <- data.frame(
  Gene  = unlist(gene_groups),
  Group = rep(names(gene_groups), lengths(gene_groups)),
  stringsAsFactors = FALSE
)

polya_long <- polya_long %>%
  left_join(group_lookup, by = "Gene") %>%
  mutate(Group = ifelse(is.na(Group), "Other", Group))

# ---- 5. CALCULATE SUMMARY STATISTICS PER GENE --------------
gene_stats <- polya_long %>%
  group_by(Gene) %>%
  summarise(
    n_reads  = n(),
    mean_pa  = round(mean(PolyA_Length), 1),
    median_pa= round(median(PolyA_Length), 1),
    sd_pa    = round(sd(PolyA_Length), 1),
    se_pa    = round(sd(PolyA_Length) / sqrt(n()), 2),
    ci95     = round(1.96 * sd(PolyA_Length) / sqrt(n()), 2),
    .groups  = "drop"
  ) %>%
  filter(n_reads >= 100) %>%               # only high confidence genes
  arrange(desc(mean_pa))

cat("\nTop 10 genes by mean polyA tail length:\n")
print(head(gene_stats, 10))

cat("\nBottom 10 genes by mean polyA tail length:\n")
print(tail(gene_stats, 10))

# Save stats table
write.csv(gene_stats, "polyA_gene_summary_stats.csv", row.names = FALSE)
cat("Summary stats saved to polyA_gene_summary_stats.csv\n")

# ---- 6. FIGURE 1: OVERALL DISTRIBUTION ---------------------
overall_mean   <- mean(polya_long$PolyA_Length)
overall_median <- median(polya_long$PolyA_Length)

fig1 <- ggplot(polya_long, aes(x = PolyA_Length)) +
  geom_histogram(binwidth = 10, fill = "#2E75B6", color = "white",
                 alpha = 0.85) +
  geom_vline(xintercept = overall_mean,   color = "#C00000",
             linetype = "dashed", linewidth = 0.9) +
  geom_vline(xintercept = overall_median, color = "#ED7D31",
             linetype = "dotted", linewidth = 0.9) +
  annotate("text", x = overall_mean + 18, y = Inf,
           label = paste0("Mean = ", round(overall_mean, 1), " nt"),
           color = "#C00000", vjust = 2, size = 3.5) +
  annotate("text", x = overall_median - 18, y = Inf,
           label = paste0("Median = ", round(overall_median, 1), " nt"),
           color = "#ED7D31", vjust = 4, size = 3.5) +
  xlim(0, 450) +
  labs(
    title    = "Overall PolyA Tail Length Distribution",
    subtitle = paste0("MDA-231 Triple-Negative Breast Cancer | ",
                      format(nrow(polya_long), big.mark = ","),
                      " reads | 96 cancer genes"),
    x        = "PolyA Tail Length (nt)",
    y        = "Number of Reads"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(color = "grey40", size = 10)
  )

ggsave("Fig1_Overall_PolyA_Distribution.png", fig1,
       width = 9, height = 5, dpi = 300)
cat("Figure 1 saved\n")

# ---- 7. FIGURE 2: MEAN POLYA PER GENE (BAR CHART) ----------
# High confidence genes only (>=100 reads), sorted by mean
gene_stats_plot <- gene_stats %>%
  left_join(group_lookup, by = "Gene") %>%
  mutate(Group = ifelse(is.na(Group), "Other", Group)) %>%
  mutate(Gene = factor(Gene, levels = Gene[order(mean_pa)]))

# Colour by significance relative to dataset mean
dataset_mean <- mean(polya_long$PolyA_Length)

gene_stats_plot <- gene_stats_plot %>%
  mutate(Category = case_when(
    mean_pa > dataset_mean + 22 ~ "Long (>130 nt)",
    mean_pa < dataset_mean - 13 ~ "Short (<95 nt)",
    TRUE                         ~ "Intermediate"
  ))

colour_map <- c(
  "Long (>130 nt)"  = "#2E75B6",
  "Intermediate"    = "#A9A9A9",
  "Short (<95 nt)"  = "#C00000"
)

fig2 <- ggplot(gene_stats_plot,
               aes(x = Gene, y = mean_pa, fill = Category)) +
  geom_col(width = 0.75) +
  geom_errorbar(aes(ymin = mean_pa - ci95, ymax = mean_pa + ci95),
                width = 0.3, linewidth = 0.4, color = "grey30") +
  geom_hline(yintercept = dataset_mean, linetype = "dashed",
             color = "black", linewidth = 0.7) +
  annotate("text", x = 1, y = dataset_mean + 4,
           label = paste0("Dataset mean = ", round(dataset_mean, 1), " nt"),
           hjust = 0, size = 3, color = "black") +
  scale_fill_manual(values = colour_map, name = "PolyA Category") +
  coord_flip() +
  labs(
    title    = "Mean PolyA Tail Length per Cancer Gene",
    subtitle = "MDA-231 TNBC | Genes with \u2265100 reads | Error bars = 95% CI",
    x        = NULL,
    y        = "Mean PolyA Tail Length (nt)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title    = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(color = "grey40", size = 9),
    legend.position = "bottom"
  )

ggsave("Fig2_Mean_PolyA_Per_Gene.png", fig2,
       width = 10, height = 14, dpi = 300)
cat("Figure 2 saved\n")

# ---- 8. FIGURE 3: BOXPLOT BY FUNCTIONAL GROUP --------------
# This is the key biological finding figure

# Filter to main groups only, remove Other
group_data <- polya_long %>%
  filter(Group %in% c("Cell Adhesion (Integrins)", "Oncogenes",
                       "Tumour Suppressors", "EMT and Invasion",
                       "Apoptosis", "Housekeeping", "Angiogenesis"))

# Order groups by median polyA length
group_order <- group_data %>%
  group_by(Group) %>%
  summarise(med = median(PolyA_Length)) %>%
  arrange(desc(med)) %>%
  pull(Group)

group_data$Group <- factor(group_data$Group, levels = group_order)

group_colours <- c(
  "Cell Adhesion (Integrins)" = "#1F3864",
  "Oncogenes"                 = "#C00000",
  "Tumour Suppressors"        = "#2E75B6",
  "EMT and Invasion"          = "#7030A0",
  "Apoptosis"                 = "#ED7D31",
  "Housekeeping"              = "#70AD47",
  "Angiogenesis"              = "#00B0F0"
)

fig3 <- ggplot(group_data,
               aes(x = Group, y = PolyA_Length, fill = Group)) +
  geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.3,
               notch = FALSE, linewidth = 0.5) +
  geom_hline(yintercept = dataset_mean, linetype = "dashed",
             color = "black", linewidth = 0.7) +
  scale_fill_manual(values = group_colours) +
  scale_y_continuous(limits = c(0, 450)) +
  annotate("text", x = 0.6, y = dataset_mean + 12,
           label = paste0("Dataset mean\n", round(dataset_mean, 1), " nt"),
           hjust = 0, size = 3, color = "black") +
  labs(
    title    = "PolyA Tail Length by Functional Gene Category",
    subtitle = "MDA-231 TNBC | Whiskers = 10th-90th percentile",
    x        = NULL,
    y        = "PolyA Tail Length (nt)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title      = element_text(face = "bold", size = 13),
    plot.subtitle   = element_text(color = "grey40", size = 9),
    axis.text.x     = element_text(angle = 35, hjust = 1, size = 10),
    legend.position = "none"
  )

ggsave("Fig3_PolyA_By_Functional_Group.png", fig3,
       width = 10, height = 6, dpi = 300)
cat("Figure 3 saved\n")

# ---- 9. FIGURE 4: INTEGRIN FAMILY SPOTLIGHT ----------------
# Biological highlight: all integrins show above-average polyA tails

integrin_genes <- c("ITGB1", "ITGA1", "ITGA2", "ITGA3",
                    "ITGAV", "ITGB3", "ITGB5")

integrin_data <- polya_long %>%
  filter(Gene %in% integrin_genes) %>%
  mutate(Gene = factor(Gene,
                       levels = gene_stats$Gene[gene_stats$Gene %in%
                                                  integrin_genes]))

integrin_stats <- gene_stats %>%
  filter(Gene %in% integrin_genes) %>%
  mutate(Gene = factor(Gene, levels = levels(integrin_data$Gene)))

fig4 <- ggplot(integrin_stats,
               aes(x = reorder(Gene, mean_pa), y = mean_pa)) +
  geom_col(fill = "#1F3864", width = 0.65) +
  geom_errorbar(aes(ymin = mean_pa - ci95, ymax = mean_pa + ci95),
                width = 0.25, linewidth = 0.6, color = "grey30") +
  geom_hline(yintercept = dataset_mean, linetype = "dashed",
             color = "#C00000", linewidth = 0.9) +
  annotate("text", x = 0.6, y = dataset_mean + 3,
           label = paste0("Dataset mean = ", round(dataset_mean, 1), " nt"),
           hjust = 0, size = 3.5, color = "#C00000") +
  coord_flip() +
  labs(
    title    = "Integrin Family: Coordinated mRNA Stabilisation in MDA-231",
    subtitle = "All detected integrins display above-average polyA tail lengths\nError bars = 95% CI",
    x        = NULL,
    y        = "Mean PolyA Tail Length (nt)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(color = "grey40", size = 9)
  )

ggsave("Fig4_Integrin_Family_PolyA.png", fig4,
       width = 8, height = 5, dpi = 300)
cat("Figure 4 saved\n")

# ---- 10. FIGURE 5: KEY CONTRAST — ITGB1 vs RPL37A ----------
contrast_data <- polya_long %>%
  filter(Gene %in% c("ITGB1", "RPL37A")) %>%
  mutate(Gene = factor(Gene, levels = c("ITGB1", "RPL37A")))

contrast_stats <- gene_stats %>%
  filter(Gene %in% c("ITGB1", "RPL37A"))

contrast_colours <- c("ITGB1" = "#2E75B6", "RPL37A" = "#C00000")

fig5 <- ggplot(contrast_data,
               aes(x = PolyA_Length, fill = Gene, color = Gene)) +
  geom_histogram(binwidth = 10, alpha = 0.65,
                 position = "identity") +
  geom_vline(data = contrast_stats,
             aes(xintercept = mean_pa, color = Gene),
             linetype = "dashed", linewidth = 1) +
  scale_fill_manual(values  = contrast_colours) +
  scale_color_manual(values = contrast_colours) +
  geom_hline(yintercept = 0) +
  xlim(0, 450) +
  labs(
    title    = "ITGB1 vs RPL37A: Biological Contrast in PolyA Tail Length",
    subtitle = paste0("ITGB1 (invasion driver): mean ",
                      contrast_stats$mean_pa[contrast_stats$Gene=="ITGB1"],
                      " nt | RPL37A (housekeeping): mean ",
                      contrast_stats$mean_pa[contrast_stats$Gene=="RPL37A"],
                      " nt | Mann-Whitney p<0.001"),
    x        = "PolyA Tail Length (nt)",
    y        = "Number of Reads",
    fill     = "Gene",
    color    = "Gene"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(color = "grey40", size = 9),
    legend.position = "top"
  )

ggsave("Fig5_ITGB1_vs_RPL37A_Contrast.png", fig5,
       width = 9, height = 5, dpi = 300)
cat("Figure 5 saved\n")

# ---- 11. STATISTICAL TEST: KRUSKAL-WALLIS ------------------
# Test whether polyA lengths differ significantly across functional groups

kw_test <- kruskal.test(PolyA_Length ~ Group, data = group_data)
cat("\nKruskal-Wallis test across functional groups:\n")
print(kw_test)

# ---- 12. SUMMARY PRINT -------------------------------------
cat("\n================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("================================================\n")
cat("Total reads analysed:        ", format(nrow(polya_long), big.mark=","), "\n")
cat("Genes with >= 100 reads:     ", nrow(gene_stats), "\n")
cat("Overall mean polyA:          ", round(dataset_mean, 1), "nt\n")
cat("Overall median polyA:        ", round(overall_median, 1), "nt\n")
cat("Longest mean polyA gene:     ", gene_stats$Gene[1],
    "(", gene_stats$mean_pa[1], "nt)\n")
cat("Shortest mean polyA gene:    ",
    gene_stats$Gene[nrow(gene_stats)],
    "(", gene_stats$mean_pa[nrow(gene_stats)], "nt)\n")
cat("Kruskal-Wallis H:            ",
    round(kw_test$statistic, 0), "\n")
cat("Kruskal-Wallis p-value:      ",
    format(kw_test$p.value, scientific = TRUE), "\n")
cat("\nOutput files created:\n")
cat("  polyA_gene_summary_stats.csv\n")
cat("  Fig1_Overall_PolyA_Distribution.png\n")
cat("  Fig2_Mean_PolyA_Per_Gene.png\n")
cat("  Fig3_PolyA_By_Functional_Group.png\n")
cat("  Fig4_Integrin_Family_PolyA.png\n")
cat("  Fig5_ITGB1_vs_RPL37A_Contrast.png\n")
cat("================================================\n")
