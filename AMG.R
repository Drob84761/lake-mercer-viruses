##### VIRAL DIVERSITY ANALYSIS #######

library(phyloseq)
library(ggplot2)
library(gridExtra)
library(vegan)
library(pheatmap)
library(DESeq)
library(DESeq2)
library(dplyr)
library(readxl)
library(tidyverse)
library(grid)
library(reshape2)
library(RColorBrewer)
#install.packages("wesanderson")
library(wesanderson)
library("knitr")
library("BiocStyle")
library(devtools)
#install_github("Russel88/MicEco")
library(MicEco)
if (!requireNamespace("remotes", quietly=TRUE))
  #install.packages("remotes")
  #remotes::install_github("YuLab-SMU/MicrobiotaProcess")
  library(MicrobiotaProcess)
#install.packages("ggthemes")

theme_set(theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank(), axis.line = element_line(colour = NULL)))                          ## I also like theme_calc, theme_classic, theme_gray


#BiocManager::install("vmikk/metagMisc")
getwd()
setwd("/users/davidrobinson/Desktop/FINAL_PNAS_FILES")
otu_mat <- read.delim("all_samples_amg_2.txt", header=T, row.names="Contig")
tax_mat <- read.delim("amg_summary.txt", header=T, row.names="Contig")
samples_df <- read.delim("map_salsa.txt", header=T, row.names="samples")



tax_mat <- as.matrix(tax_mat)
otu_mat <- as.matrix(otu_mat)


OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
#OTU2 = transform_sample_counts(OTU, function(x) x / sum(x, na.rm = TRUE))
#OTU2 = transform_sample_counts(otu_mat, as.integer)
samples = sample_data(samples_df)
TAX = tax_table(tax_mat)
TAX

viral_all<- phyloseq(OTU, TAX, samples)
viral_all
rank_names(viral_all)


# Define a category for each AMG based on annotation
tax_df <- data.frame(as(tax_table(viral_all), "matrix"))  # extract taxonomy data frame
tax_df$Category <- "Other"  # default category
# Assign categories based on keywords in 'sheet' or 'header'
tax_df$Category[grepl("carbon", tax_df$Category, ignore.case=TRUE) | 
                  grepl("carbon", tax_df$header, ignore.case=TRUE) | 
                  grepl("CAZY", tax_df$header, ignore.case=TRUE)] <- "Carbon"
tax_df$Category[grepl("Nitrogen", tax_df$Category, ignore.case=TRUE) | 
                  grepl("Amino Acid", tax_df$header, ignore.case=TRUE) | 
                  grepl("Peptidase", tax_df$header, ignore.case=TRUE)] <- "Nitrogen"
tax_df$Category[grepl("Nucleotide", tax_df$header, ignore.case=TRUE) | 
                  grepl("Nucleotide", tax_df$subheader, ignore.case=TRUE) | 
                  grepl("ribonucleotide|thymidylate|adenylate|ribosomal protein", tax_df$gene_description, ignore.case=TRUE)] <- "Nucleotide/Cofactor"

# Summarize the viral AMG abundance by category in each sample
otu_df <- as.data.frame(otu_table(viral_all))
otu_df$Category <- tax_df$Category  # add category info to each gene row
library(dplyr)
agg_df <- otu_df %>%
  group_by(Category) %>%
  summarize(across(where(is.numeric), sum, na.rm=TRUE)) %>%
  tidyr::pivot_longer(cols = -Category, names_to = "Sample", values_to = "Abundance") %>%
  group_by(Sample) %>%
  mutate(RelativeAbundance = Abundance / sum(Abundance))  # compute fraction per sample

# Plot stacked bar chart of category composition per sample
library(ggplot2)
ggplot(agg_df, aes(x = Sample, y = RelativeAbundance, fill = Category)) +
  geom_bar(stat="identity", position="stack", color="black") +
  scale_y_continuous(labels=scales::percent) +
  labs(y = "Relative Abundance of AMGs", x = "Sample", fill = "AMG Category") +
  ggtitle("Viral AMG Composition by Metabolic Function") +
  theme_minimal()




library(microeco)
# Create a presence/absence matrix of AMGs (features as rows, samples as columns)
otu_mat <- as(otu_table(viral_all), "matrix")
otu_mat[otu_mat > 0] <- 1  # convert counts to 1/0 presence-absence
# Initialize MicEco trans_venn object
venn_obj <- trans_venn$new(dataset = otu_mat, ratio = "numratio")
# Plot Venn diagram
venn_obj$plot_venn(fill_color = TRUE, petal_plot = FALSE)

# ordinate the phyloseq object (Bray-Curtis distance on raw abundances)
ord <- ordinate(viral_all, method = "PCoA", distance = "bray")
# Plot the ordination with samples colored by lake
p <- plot_ordination(viral_all, ord, type="samples", color="lake", label="name") +
  geom_point(size=4) +
  labs(title="PCoA of Viral AMG Profiles", color="Lake") +
  theme_bw()
print(p)


# Use your VC color palette
VC_colors <- c(
  "#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e",
  "#e6ab02", "#a6761d", "#666666", "#1f78b4", "#b2df8a",
  "#fb9a99", "#fdbf6f", "#cab2d6", "#ffff99", "#6a3d9a",
  "#ff7f00", "#b15928", "#8dd3c7", "#fef0d9", "#bebada",
  "#fb8072", "#80b1d3", "#fdb462", "#b3de69", "#fccde5",
  "#d9d9d9", "#bc80bd", "#ccebc5", "#ffed6f", "grey"
)

# Agglomerate and transform to relative abundance per sample
viral_category <- viral_all %>%
  tax_glom(taxrank = "sheet") %>%
  transform_sample_counts(function(x) x / sum(x)) %>%
  psmelt() %>%
  arrange(sheet)

# Recode sample names
viral_category$Sample <- recode(viral_category$Sample,
                                "SLMCELL" = "Mercer Cellular",
                                "SLMVIRAL" = "Mercer Viral",
                                "SLW" = "Whillans Cellular")

# Assign color palette
unique_categories <- sort(unique(viral_category$Category))
n_cat <- length(unique_categories)

if (length(VC_colors) < n_cat) {
  VC_colors <- rep(VC_colors, length.out = n_cat)
}
category_colors <- setNames(VC_colors[1:n_cat], unique_categories)

p_category <- ggplot(viral_category, aes(x = Sample, y = Abundance * 100, fill = sheet)) +
  geom_bar(stat = "identity", position = "stack", color = "NA", linewidth = 0.2) +
  scale_fill_manual(values = VC_colors) +
  theme_classic(base_size = 20) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank()
  ) +
  labs(
    y = "Relative Abundance (%)",
    fill = "AMG Category"
  ) 
  #ggtitle("Viral AMG Composition by Metabolic Function")

print(p_category)


# Summarize abundance by Sample and Category
viral_summary <- viral_category %>%
  group_by(Sample, Category) %>%
  summarise(RelAbund = sum(Abundance) * 100, .groups = "drop")

# Plot clean summary (fewer outlines)
p_category <- ggplot(viral_summary, aes(x = Sample, y = RelAbund, fill = Category)) +
  geom_bar(stat = "identity", position = "stack", color = "black", size = 0.2) +
  scale_fill_manual(values = category_colors) +
  theme_classic(base_size = 20) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank()
  ) +
  labs(
    y = "Relative Abundance (%)",
    fill = "AMG Category"
  ) 
  #ggtitle("Viral AMG Composition by Metabolic Function")

print(p_category)




######BRING THE HEAT

# ======================================================
# Contig-aligned AMG heatmap (TPM + log1p)
#   - Heatmap rows = gene_description × Lifestyle (so mixed-lifestyle genes show BOTH rows)
#   - Heatmap annotation: Lifestyle ONLY (no sheet)
#   - Bars: TOP 15 genes by total TPM (Lytic ↑, Temperate ↓; Unknown optional)
# Files:
#   - amg_summary.txt        (Contig, Lifestyle, sheet, gene_description)
#   - all_samples_amg_2.txt  (Contig + TPM columns per sample)
# ======================================================

library(tidyverse)
library(janitor)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(gridExtra)
library(grid)
library(scales)

# ---- knobs ----
MAX_ROWS_HEATMAP <- NA       # NA = show all gene×lifestyle rows; or set (e.g., 120) to cap rows
LOG_TRANSFORM     <- TRUE    # log10(TPM + 1) for heatmap color dynamic range (no z-scaling)
SHOW_UNKNOWN_BARS <- FALSE   # include Unknown lifestyle in bars (grey, upward)?
BORDER_COL        <- "grey85"
LYTIC_COL         <- "#1CA3A3"
TEMP_COL          <- "#E2B23F"
UNK_COL           <- "grey75"

# ---- 1) Load files and align by Contig ----
amg <- read.delim("amg_summary.txt", check.names = FALSE) |>
  clean_names()
abund <- read.delim("all_samples_amg_2.txt", check.names = FALSE) |>
  clean_names()

stopifnot("contig" %in% names(amg), "contig" %in% names(abund))

amg <- amg |>
  mutate(
    contig = as.integer(contig),
    lifestyle = case_when(
      str_detect(tolower(lifestyle %||% ""), "lytic|virulent") ~ "Lytic",
      str_detect(tolower(lifestyle %||% ""), "temperate|lyso") ~ "Temperate",
      TRUE ~ "Unknown"
    ),
    # sheet can exist but isn't used in the heatmap now
    gene_description = if_else(is.na(gene_description) | gene_description == "",
                               "unannotated_AMG", gene_description)
  )

abund <- abund |> mutate(contig = as.integer(contig))

# Sample columns & fixed order (as in the file)
sample_cols  <- setdiff(names(abund), "contig")
SAMPLE_ORDER <- sample_cols

# ---- 2) Join and tidy (on Contig) ----
dat <- amg |>
  select(contig, gene_description, lifestyle) |>
  inner_join(abund |> select(contig, all_of(SAMPLE_ORDER)), by = "contig")

dat_long <- dat |>
  pivot_longer(cols = all_of(SAMPLE_ORDER), names_to = "sample", values_to = "tpm") |>
  mutate(sample = factor(sample, levels = SAMPLE_ORDER))

# ======================================================
# 3) HEATMAP: rows = gene_description × Lifestyle
#    - sum TPM across contigs within each gene×lifestyle×sample
#    - optional cap on total rows
# ======================================================

# Aggregate per (gene_description × Lifestyle × sample)
gs_life_sample <- dat_long |>
  group_by(gene_description, lifestyle, sample) |>
  summarise(tpm = sum(tpm, na.rm = TRUE), .groups = "drop")

# Row key = "gene | Lifestyle", so mixed-lifestyle genes appear twice
gs_life_sample <- gs_life_sample |>
  mutate(row_key = paste0(gene_description, " | ", lifestyle))

# Optional cap on number of heatmap rows (choose top by total TPM across samples)
if (!is.na(MAX_ROWS_HEATMAP)) {
  keep_rows <- gs_life_sample |>
    group_by(row_key) |>
    summarise(total_tpm = sum(tpm), .groups = "drop") |>
    arrange(desc(total_tpm)) |>
    slice_head(n = MAX_ROWS_HEATMAP) |>
    pull(row_key)
  gs_life_sample <- gs_life_sample |> filter(row_key %in% keep_rows)
}

# Row order: by total TPM across samples (desc)
row_order_tbl <- gs_life_sample |>
  group_by(row_key) |>
  summarise(total = sum(tpm), .groups = "drop") |>
  arrange(desc(total)) |>
  mutate(row_ord = row_number())

# Collapse duplicates (safety) and pivot to matrix
gs_collapsed <- gs_life_sample |>
  group_by(row_key, sample) |>
  summarise(tpm = sum(tpm, na.rm = TRUE), .groups = "drop") |>
  mutate(sample = as.character(sample))

mat_df <- gs_collapsed |>
  left_join(row_order_tbl, by = "row_key") |>
  arrange(row_ord, sample) |>
  mutate(sample = factor(sample, levels = SAMPLE_ORDER)) |>
  arrange(row_ord, sample) |>
  select(row_key, sample, tpm) |>
  pivot_wider(
    names_from  = sample,
    values_from = tpm,
    values_fill = list(tpm = 0),
    values_fn   = list(tpm = sum)
  )

mat <- mat_df |>
  column_to_rownames("row_key") |>
  as.data.frame() |>
  mutate(across(everything(), as.numeric)) |>
  as.matrix()

plot_mat <- if (LOG_TRANSFORM) log10(mat + 1) else mat

# Row annotation: Lifestyle ONLY (parsed back from row_key)
row_anno <- tibble(row_key = rownames(plot_mat)) |>
  separate_wider_delim(row_key, delim = " | ", names = c("gene_description","Lifestyle"),
                       too_few = "align_start", cols_remove = FALSE) |>
  tibble::column_to_rownames("row_key")

# Colors
pal <- colorRampPalette(c("#084594","#BDD7E7","#F7F7F7","#FEE0D2","#A50F15"))(101)
ann_colors <- list(
  Lifestyle = c(Lytic = LYTIC_COL, Temperate = TEMP_COL, Unknown = UNK_COL)
)

# Heatmap (no clustering, gridlines, left annotation = Lifestyle)
hm <- pheatmap(
  plot_mat,
  color = pal,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  border_color = BORDER_COL,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 7,
  fontsize_col = 10,
  annotation_row = row_anno[, "Lifestyle", drop = FALSE],
  annotation_colors = ann_colors
)

# Print to RStudio
print(hm)

# ======================================================
# 4) BARS: TOP 15 genes by total TPM
#    - mean TPM per gene × Lifestyle (across samples)
#    - Lytic up, Temperate down (Unknown optional)
# ======================================================

# Total TPM per gene (to choose top 15)
gene_totals <- dat_long |>
  group_by(gene_description) |>
  summarise(total_tpm = sum(tpm, na.rm = TRUE), .groups = "drop") |>
  arrange(desc(total_tpm))

top15_genes <- gene_totals |>
  slice_head(n = 15) |>
  pull(gene_description)

# Mean TPM per Lifestyle for each top-15 gene
gene_life_mean <- dat_long |>
  filter(gene_description %in% top15_genes) |>
  group_by(gene_description, lifestyle, sample) |>
  summarise(tpm = sum(tpm, na.rm = TRUE), .groups = "drop") |>
  group_by(gene_description, lifestyle) |>
  summarise(mean_tpm = mean(tpm, na.rm = TRUE), .groups = "drop")

# Stable x order: top 15 by total TPM (high → low)
gene_order <- gene_totals |>
  filter(gene_description %in% top15_genes) |>
  pull(gene_description)

bars_df <- gene_life_mean |>
  mutate(
    gene_description = factor(gene_description, levels = gene_order),
    lifestyle = case_when(
      lifestyle %in% c("Lytic","Temperate","Unknown") ~ lifestyle,
      TRUE ~ "Unknown"
    ),
    sign_tpm = case_when(
      lifestyle == "Lytic" ~  mean_tpm,
      lifestyle == "Temperate" ~ -mean_tpm,
      TRUE ~ if (SHOW_UNKNOWN_BARS) mean_tpm else 0
    ),
    plotable = case_when(
      lifestyle %in% c("Lytic","Temperate") ~ TRUE,
      lifestyle == "Unknown" & SHOW_UNKNOWN_BARS ~ TRUE,
      TRUE ~ FALSE
    )
  ) |>
  filter(plotable)

bar_cols <- c(Lytic = LYTIC_COL, Temperate = TEMP_COL, Unknown = UNK_COL)

p_bars <- ggplot(bars_df, aes(x = gene_description, y = sign_tpm, fill = lifestyle)) +
  geom_hline(yintercept = 0, color = "grey50", linewidth = 0.4) +
  geom_col(width = 0.8, color = "black", linewidth = 0.2) +
  scale_fill_manual(values = bar_cols, drop = FALSE) +
  scale_y_continuous(labels = label_number(scale_cut = cut_short_scale())) +
  labs(x = NULL, y = "Mean TPM", fill = "Lifestyle") +
  theme_classic(base_size = 12) +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 35, hjust = 1))

# Print to RStudio
print(p_bars)

# ---- Save outputs ----
dir.create("amg_outputs", showWarnings = FALSE)

#ggsave("amg_outputs/fig_lifestyle_vertical_bars_top15_genes.png", p_bars,
       width = 9.5, height = 4.8, dpi = 300)

#png("amg_outputs/fig_heatmap_gene_description_by_lifestyle.png", width = 9.5, height = 10, units = "in", res = 300)
grid.draw(hm$gtable)
dev.off()

# Optional combined image (heatmap over bars)
png("amg_outputs/fig_combined_heatmap_plus_bars_top15.png", width = 9.5, height = 13, units = "in", res = 300)
grid.arrange(hm$gtable, ggplotGrob(p_bars), heights = c(3, 1.2))
dev.off()



# ================================
# A) COLLAPSED HEATMAP (Top 20 genes) + B) Lifestyle bars (Top 10 genes)
#     - Heatmap rows: one per gene (TPM summed across contigs & lifestyles)
#     - Row annotation: TPM-weighted dominant Lifestyle (not bold)
# ================================

library(dplyr)
library(tidyr)
library(pheatmap)
library(ggplot2)
library(ggplotify)   # as.ggplot()
library(patchwork)
library(scales)

# ----- Fix sample names globally so they propagate everywhere -----
dat_long <- dat_long %>%
  mutate(
    sample = case_when(
      sample == "SLMCELL"  ~ "Mercer Cellular",
      sample == "SLMVIRAL" ~ "Mercer Viral",
      sample == "SLW"      ~ "Whillans Cellular",
      TRUE ~ as.character(sample)
    )
  )

# Sample order should now use the new names
SAMPLE_ORDER <- c("Mercer Cellular", "Mercer Viral", "Whillans Cellular")

# ----- Colors for lifestyles -----
bar_cols <- c(Lytic = "#1CA3A3", Temperate = "#E2B23F", Unknown = "grey75")

# ----- 1) Collapse to gene × sample TPM -----
gene_sample_collapsed <- dat_long %>%
  group_by(gene_description, sample) %>%
  summarise(tpm = sum(tpm, na.rm = TRUE), .groups = "drop")

# ----- 2) Top 20 genes by total TPM -----
top20_genes <- gene_sample_collapsed %>%
  group_by(gene_description) %>%
  summarise(total_tpm = sum(tpm, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(total_tpm)) %>%
  slice_head(n = 20) %>%
  pull(gene_description)

gene_sample_top20 <- gene_sample_collapsed %>%
  filter(gene_description %in% top20_genes)

# ----- 3) Pivot to matrix (3 samples; then rename + order) -----
mat_top_df <- gene_sample_top20 %>%
  arrange(gene_description, sample) %>%
  pivot_wider(
    names_from  = sample,
    values_from = tpm,
    values_fill = list(tpm = 0),
    values_fn   = list(tpm = sum)
  ) %>%
  # robust renaming: works whether your columns are SLMCELL/slmcell/etc.
  rename_with(~"Mercer Cellular",    any_of(c("SLMCELL",  "slmcell"))) %>%
  rename_with(~"Mercer Viral",       any_of(c("SLMVIRAL", "slmviral"))) %>%
  rename_with(~"Whillans Cellular",  any_of(c("SLW",      "slw")))

# Convert to matrix
mat_top <- mat_top_df %>%
  tibble::column_to_rownames("gene_description") %>%
  as.matrix()

# ensure columns are in the desired order and only those 3
SAMPLE_ORDER <- c("Mercer Cellular", "Mercer Viral", "Whillans Cellular")
mat_top <- mat_top[, intersect(SAMPLE_ORDER, colnames(mat_top)), drop = FALSE]

# log10-transform
plot_mat_top <- log10(mat_top + 1)

# ----- 4) Dominant Lifestyle per gene (for row annotations) -----
row_anno_top <- dat_long %>%
  filter(gene_description %in% rownames(plot_mat_top)) %>%
  group_by(gene_description, lifestyle) %>%
  summarise(tpm = sum(tpm, na.rm = TRUE), .groups = "drop_last") %>%
  slice_max(tpm, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(gene_description, Lifestyle = lifestyle) %>%
  arrange(match(gene_description, rownames(plot_mat_top))) %>%
  tibble::column_to_rownames("gene_description")

# ----- 5) Heatmap -----
hm_top <- pheatmap::pheatmap(
  plot_mat_top,
  color = colorRampPalette(c("#084594","#BDD7E7","#F7F7F7","#FEE0D2","#A50F15"))(101),
  cluster_rows = FALSE, cluster_cols = FALSE,
  border_color = "grey85",
  show_rownames = TRUE, show_colnames = TRUE,
  labels_row = rownames(plot_mat_top),
  fontsize_row = 9, fontsize_col = 10,
  cellheight = 14,
  annotation_row = row_anno_top[, "Lifestyle", drop = FALSE],
  annotation_colors = list(Lifestyle = bar_cols),
  annotation_names_row = FALSE,
  silent = TRUE
)

hm_plot <- ggplotify::as.ggplot(hm_top$gtable) +
  theme(plot.margin = margin(8, 8, 4, 18))


# ----- 6) Lifestyle bars (Top 10 genes), horizontal, with Temperate labels -----
top10_genes <- dat_long %>%
  group_by(gene_description) %>%
  summarise(total_tpm = sum(tpm, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(total_tpm)) %>%
  slice_head(n = 10) %>%
  pull(gene_description)

bars_top10 <- dat_long %>%
  filter(gene_description %in% top10_genes) %>%
  group_by(gene_description, lifestyle, sample) %>%
  summarise(tpm = sum(tpm, na.rm = TRUE), .groups = "drop") %>%
  group_by(gene_description, lifestyle) %>%
  summarise(mean_tpm = mean(tpm, na.rm = TRUE), .groups = "drop")

# order so coord_flip shows largest at top
gene_order10 <- dat_long %>%
  group_by(gene_description) %>%
  summarise(total_tpm = sum(tpm, na.rm = TRUE), .groups = "drop") %>%
  filter(gene_description %in% top10_genes) %>%
  arrange(total_tpm) %>%
  pull(gene_description)

bars_10_df <- bars_top10 %>%
  mutate(
    gene_description = factor(gene_description, levels = gene_order10),
    lifestyle = ifelse(lifestyle %in% c("Lytic","Temperate","Unknown"), lifestyle, "Unknown"),
    y_val = ifelse(lifestyle == "Temperate", -mean_tpm, mean_tpm)
  )

lab_fun <- label_number(accuracy = 1, scale_cut = cut_short_scale())
labels_temp <- bars_10_df %>%
  filter(lifestyle == "Temperate") %>%
  mutate(label_txt = lab_fun(abs(y_val)))

p_bars_10 <- ggplot(bars_10_df, aes(x = gene_description, y = y_val, fill = lifestyle)) +
  geom_hline(yintercept = 0, color = "grey60", linewidth = 0.35) +
  geom_col(width = 0.8, color = "black", linewidth = 0.2) +
  geom_text(data = labels_temp,
            aes(label = label_txt),
            hjust = 1.05, vjust = 0.5, size = 3.2) +
  coord_flip() +
  scale_fill_manual(values = bar_cols, drop = FALSE, name = "Lifestyle") +
  scale_y_continuous(labels = lab_fun, expand = expansion(mult = c(0.06, 0.08))) +
  labs(x = NULL, y = "Mean TPM") +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "right",
    axis.text.y = element_text(size = 9),
    axis.text.x = element_text(size = 9),
    plot.margin = margin(4, 8, 8, 18)
  )

# ----- 7) Final stacked figure (A on top, B bottom; ~75/25 split) -----
final_tagged <- (hm_plot) / (p_bars_10) +
  plot_layout(heights = c(3, 1), guides = "keep") +
  plot_annotation(tag_levels = "A", tag_suffix = ".")

print(final_tagged)



