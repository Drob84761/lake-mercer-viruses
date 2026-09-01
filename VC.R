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
install.packages("ggthemes")

theme_set(theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank(), axis.line = element_line(colour = NULL)))                          ## I also like theme_calc, theme_classic, theme_gray


#BiocManager::install("vmikk/metagMisc")
getwd()
setwd("/users/davidrobinson/Desktop/FINAL_PNAS_FILES")
otu_mat <- read.delim("all_samples.txt", header=T, row.names="Contig")
tax_mat <- read.delim("genome_by_genome_overview.txt", header=T, row.names="Contig")
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


##############################################
## FINAL SCRIPT: Three Plot Versions Matching Paper Filtering
##############################################

# Load required libraries
library(phyloseq)
library(ggplot2)
library(dplyr)

##############################################
## COLOR VECTOR
##############################################

VC_colors <- c(
  "#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e",
  "#e6ab02", "#a6761d", "#666666", "#1f78b4", "#b2df8a",
  "#fb9a99", "#fdbf6f", "#cab2d6", "#ffff99", "#6a3d9a",
  "#ff7f00", "#b15928", "#8dd3c7", "#fef0d9", "#bebada",
  "#fb8072", "#80b1d3", "#fdb462", "#b3de69", "#fccde5",
  "#d9d9d9", "#bc80bd", "#ccebc5", "#ffed6f", "grey"
)

##############################################
## DATA PREPARATION
##############################################

vc_agg <- tax_glom(viral_all, taxrank = "VC")
vc_rel <- transform_sample_counts(vc_agg, function(x) x / sum(x))
df <- psmelt(vc_rel)
df$VC <- as.character(df$VC)

# Keep VCs if ANY sample >3%
df_max <- df %>%
  group_by(VC) %>%
  summarise(MaxAbundance = max(Abundance))

top_vcs_any_sample <- df_max$VC[df_max$MaxAbundance >= 0.03]

df$VC[!(df$VC %in% top_vcs_any_sample)] <- "<3% Abundance"

# Remove _0 suffix
df$VC <- gsub("_0$", "", df$VC)

# Summarize
df_grouped <- df %>%
  group_by(Sample, VC) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop")

df_grouped$Lake <- sample_data(viral_all)$lake[match(df_grouped$Sample, rownames(sample_data(viral_all)))]

unique_vcs <- unique(df_grouped$VC)
if (length(VC_colors) < length(unique_vcs)) {
  VC_colors <- c(VC_colors, rep("grey", length(unique_vcs) - length(VC_colors)))
}
VC_colors_named <- setNames(VC_colors[1:length(unique_vcs)], sort(unique_vcs))

##############################################
## VERSION 1: Minimalist Professional
##############################################

p1 <- ggplot(df_grouped, aes(x = Sample, y = Abundance * 100, fill = VC)) +
  geom_bar(stat = "identity", position = "stack", color = "black", size = 0.2) +
  scale_fill_manual(values = VC_colors_named) +
  theme_bw(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank()
  ) +
  labs(y = "Relative Abundance (%)", fill = "Viral Cluster") +
  ggtitle("Relative Abundance of Viral Clusters")

print(p1)

##############################################
## VERSION 2: Faceted by Lake
##############################################

p2 <- ggplot(df_grouped, aes(x = Sample, y = Abundance * 100, fill = VC)) +
  geom_bar(stat = "identity", position = "stack", color = "white", size = 0.2) +
  scale_fill_manual(values = VC_colors_named) +
  facet_wrap(~Lake, scales = "free_x") +
  theme_bw(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank()
  ) +
  labs(y = "Relative Abundance (%)", fill = "Viral Cluster") +
  ggtitle("Relative Abundance of Viral Clusters by Lake")

print(p2)

##############################################
## VERSION 3: Classic Theme
##############################################

# Recode Sample names to full labels
df_grouped$Sample <- recode(df_grouped$Sample,
                            "SLMCELL" = "Mercer Cellular",
                            "SLMVIRAL" = "Mercer Viral",
                            "SLW" = "Whillans Cellular")

# Classic Style Plot with Updated Labels
p3 <- ggplot(df_grouped, aes(x = Sample, y = Abundance * 100, fill = VC)) +
  geom_bar(stat = "identity", position = "stack", color = "black", linewidth = 0.2) +
  scale_fill_manual(values = VC_colors) +
  theme_classic(base_size = 20) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank()
  ) +
  labs(y = "Relative Abundance (%)", fill = "Viral Cluster") 
  #ggtitle("Relative Abundance of Viral Clusters")

print(p3)




##############################################
## PRINT
##############################################

print(p1)
print(p2)
print(p3)

# To save high-resolution:
# ggsave("VC_plot_version1.png", plot = p1, dpi=300, width=10, height=6)
# ggsave("VC_plot_version2.png", plot = p2, dpi=300, width=12, height=6)
# ggsave("VC_plot_version3.png", plot = p3, dpi=300, width=10, height=6)

# Load libraries
library(phyloseq)
library(ggplot2)
library(dplyr)

# Agglomerate by VC
vc_agg <- tax_glom(viral_all, taxrank = "VC")

# Transform to relative abundance
vc_rel <- transform_sample_counts(vc_agg, function(x) x / sum(x))

# Melt to long format
vc_df <- psmelt(vc_rel)

# Convert to character for relabeling
vc_df$VC <- as.character(vc_df$VC)

# Relabel <3% abundance
vc_df$VC[vc_df$Abundance < 0.03] <- "<3% Abundance"

# Group and sum abundances
vc_grouped <- vc_df %>%
  group_by(Sample, VC) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop")

# Plot
ggplot(vc_grouped, aes(x = Sample, y = Abundance * 100, fill = VC)) +
  geom_bar(stat = "identity", position = "stack", color = "black", size = 0.1) +
  scale_fill_manual(values = VC_colors) +
  theme_bw(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank()
  ) +
  labs(
    y = "Relative Abundance (%)",
    fill = "Viral Cluster"
  ) +
  ggtitle("Relative Abundance of Viral Clusters Across Samples")


##############################################
## COMPLEX UPSET PLOT OF VIRAL CONTIGS with Labels and Clear Axes
##############################################

# Remove the NA-named column
vc_pa_clean <- vc_pa[ , !is.na(colnames(vc_pa))]

# Load ComplexUpset
library(ComplexUpset)

# Create the plot
upset(
  vc_pa_clean,
  intersect = c("Mercer_Cell", "Mercer_Viral", "Whillans"),
  name = "Viral Contigs",
  base_annotations = list(
    'Intersection size' = intersection_size(
      counts = TRUE,
      text = list(
        vjust = -0.3,
        size = 3
      )
    )
  ),
  width_ratio = 0.2
) +
  #ggtitle("Intersection of Viral Contigs Across Samples") +
  labs(
    x = "Sample Intersections",
    y = "Number of Viral Contigs"
  ) +
  theme(
    plot.title = element_text(size=16, face="bold"),
    axis.text.x = element_text(size=12),
    axis.text.y = element_text(size=12)
  )


##################PCOA

# Recode lake names and ensure correct formatting
sample_data(viral_all)$Lake <- recode(sample_data(viral_all)$lake,
                                      "Mercer" = "Mercer",
                                      "Wissard" = "Whillans")

# Also rename shape variable if needed
sample_data(viral_all)$Type <- sample_data(viral_all)$type

# Re-run the ordination plot
plot_ordination(
  physeq = viral_all,
  ordination = viral_pcoa,
  color = "Lake"
) +
  geom_point(size = 6) +
  ggtitle("Viral Communities PCoA by Lake and Sample Type") +
  ylim(-0.6, 0.6) + 
  xlim(-0.6, 0.6) +
  scale_shape_manual(values = c(15, 16, 17)) +  # customize if you want different symbols
  guides(color = guide_legend(override.aes = list(shape = 16))) +
  labs(
    x = "P1-Percent Variation Explained 46.8%",
    y = "P2-Percent Variation Explained 52.2%"
  ) +
  theme_bw(base_size = 20) +
  theme(panel.grid = element_blank()) +
  geom_hline(yintercept = 0, color = "grey", linetype = "dashed") + 
  geom_vline(xintercept = 0, color = "grey", linetype = "dashed")



######SHARED CONTIG TABLE


#remotes::install_github("Russel88/MicEco")
library(MicEco)

venn <- ps_venn(depth, "name", fraction = 0, weight = FALSE, type = "percent",
                relative = FALSE, plot = TRUE) 
ggtitle("Shared Viral Contigs")


venn

library(dplyr)

# Read the file
depth <- read.delim("all_samples.txt")

# Convert abundances to presence/absence and create a pattern string for each contig
presence_patterns <- depth %>%
  mutate(across(c(SLMCELL, SLMVIRAL, SLW), ~ . > 0)) %>%
  rowwise() %>%
  mutate(PresentIn = paste(
    c("SLMCELL", "SLMVIRAL", "SLW")[c_across(c(SLMCELL, SLMVIRAL, SLW))],
    collapse = " + "
  )) %>%
  ungroup()

# Summarize the number of contigs in each combination
shared_by_sample <- presence_patterns %>%
  count(PresentIn, name = "Number_of_Contigs") %>%
  arrange(desc(Number_of_Contigs))

# Print summary table
print(shared_by_sample)


# Use the relative-abundance long table *before* you relabel "<3% Abundance"
df_raw <- psmelt(vc_rel) %>% mutate(VC = as.character(VC))

# If you recoded the sample labels, SLW might be "Whillans Cellular"; adjust as needed:
slw_id <- "SLW"                  # or "Whillans Cellular" if using recoded names

# Distinct VCs present (>0) in SLW
vc_present_slw <- df_raw %>%
  filter(Sample == slw_id, Abundance > 0) %>%
  distinct(VC)

nrow(vc_present_slw)   # total VCs detected in SLW (any abundance > 0)




# install.packages("ggVennDiagram")  # run once if needed
library(tidyverse)
library(ggVennDiagram)

# 1) Load abundance table (contigs as rownames, samples as columns)
depth <- read.delim("all_samples.txt", header = TRUE, row.names = "Contig", check.names = FALSE)

# Keep only the three samples we care about (adjust names if yours differ)
samples <- c("SLMCELL", "SLMVIRAL", "SLW")
stopifnot(all(samples %in% colnames(depth)))

# Drop contigs that are zero across all three samples
depth3 <- depth[, samples, drop = FALSE]
depth3 <- depth3[rowSums(depth3) > 0, , drop = FALSE]

# 2) Presence/absence by contig
pa <- as.matrix(depth3 > 0)

# 3) Build sets of contig IDs for the Venn
contig_sets <- list(
  "Mercer Cellular"  = rownames(depth3)[pa[, "SLMCELL"]],
  "Mercer Viral"     = rownames(depth3)[pa[, "SLMVIRAL"]],
  "Whillans Cellular"= rownames(depth3)[pa[, "SLW"]]
)

# 4) Plot Venn
ggVennDiagram(contig_sets, label_alpha = 0) +
  scale_fill_gradient(low = "white", high = "steelblue")

# 5) Also print a counts table for each intersection
venn_tbl <- tibble(
  Contig    = rownames(depth3),
  SLMCELL   = pa[, "SLMCELL"],
  SLMVIRAL  = pa[, "SLMVIRAL"],
  SLW       = pa[, "SLW"]
) %>%
  mutate(PresentIn = case_when(
    SLMCELL & SLMVIRAL & SLW ~ "All three",
    SLMCELL & SLMVIRAL & !SLW ~ "Mercer Cellular + Mercer Viral",
    SLMCELL & !SLMVIRAL & SLW ~ "Mercer Cellular + Whillans",
    !SLMCELL & SLMVIRAL & SLW ~ "Mercer Viral + Whillans",
    SLMCELL & !SLMVIRAL & !SLW ~ "Mercer Cellular only",
    !SLMCELL & SLMVIRAL & !SLW ~ "Mercer Viral only",
    !SLMCELL & !SLMVIRAL & SLW ~ "Whillans only",
    TRUE ~ "None"
  )) %>%
  count(PresentIn, name = "Number_of_Contigs") %>%
  arrange(desc(Number_of_Contigs))

print(venn_tbl)

# Optional: save figure
# ggsave("venn_contigs_SLW_Mercer.png", width = 6, height = 5, dpi = 300)


# Create Table 1 directly
table1 <- data.frame(
  Cast = 1:6,
  VLP_ml = c(14485.1, 3049.5, 12037.5, 19300.1, 14806.1, 17133.4),
  Cells_ml = c(23914.5, 13321.5, 25399.1, 23914.5, 20330.0, 22844.5)
)

# Calculate Viral-to-Bacterial Ratio (VBR)
table1$VBR <- round(table1$VLP_ml / table1$Cells_ml, 3)

# Print the table
print(table1)

# Optional: write to CSV
# write.csv(table1, "Table1_VLP_Cells_VBR.csv", row.names = FALSE)



#################PANEL

# === EXACT-LIKE-SCREENSHOT LAYOUT (wide, short) ===
# Tweak only these two to scale up/down without changing the look
FIG_W <- 10.0   # inches  (use 7.5 for PNAS double-column; look stays identical)
FIG_H <- 4.5    # inches
DPI   <- 300

library(tidyverse)
library(cowplot)
library(patchwork)
library(ggVennDiagram)

# ------- A) TAXONOMY -------
p3_base <- ggplot(df_grouped, aes(x = Sample, y = Abundance * 100, fill = VC)) +
  geom_bar(stat = "identity", position = "stack", color = "black", linewidth = 0.25) +
  scale_fill_manual(values = VC_colors_named) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x  = element_text(angle = 35, hjust = 1, vjust = 1, size = 11.5),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 13, margin = margin(r = 8)),
    plot.margin  = margin(6, 10, 2, 8),
    legend.position = "bottom",
    legend.justification = c(0,0),
    legend.title = element_text(size = 12.5),
    legend.text  = element_text(size = 11),
    legend.key.height = unit(5, "mm"),
    legend.key.width  = unit(7, "mm")
  ) +
  guides(fill = guide_legend(ncol = 9, title.position = "top")) +  # 2 rows (9 x 2)
  labs(y = "Relative Abundance (%)", fill = "Viral Cluster")

leg <- cowplot::get_legend(p3_base)

# add bottom padding so x labels don't sit under legend
p3_noleg_pad <- p3_base +
  theme(legend.position = "none",
        plot.margin = margin(6, 8, 24, 8))

left_col <- cowplot::plot_grid(
  p3_noleg_pad,
  patchwork::plot_spacer(),
  leg,
  ncol = 1,
  rel_heights = c(1, 0.04, 0.18)   # spacer + legend proportions tuned for this look
)

# ------- B) VENN (contigs only) -------
depth  <- read.delim("all_samples.txt", header = TRUE, row.names = "Contig", check.names = FALSE)
samples <- c("SLMCELL","SLMVIRAL","SLW")
depth3 <- depth[, samples, drop = FALSE]
depth3 <- depth3[rowSums(depth3) > 0, , drop = FALSE]
pa <- as.matrix(depth3 > 0)

contig_sets <- list(
  "Mercer Cellular"   = rownames(depth3)[pa[, "SLMCELL"]],
  "Mercer Viral"      = rownames(depth3)[pa[, "SLMVIRAL"]],
  "Whillans Cellular" = rownames(depth3)[pa[, "SLW"]]
)

venn_B <- ggVennDiagram(contig_sets, label = "count", label_alpha = 0) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  theme_void(base_size = 13.5) +
  theme(
    legend.position = "none",
    plot.margin = margin(8, 10, 8, 8),
    text = element_text(size = 13.5)
  ) +
  coord_cartesian(clip = "off")

# ------- Combine: taxonomy ~70%, venn ~30% (matches screenshot balance) -------
panel <- left_col | venn_B + plot_layout(widths = c(7, 3), guides = "keep")
panel <- panel + plot_annotation(tag_levels = "A",
                                 theme = theme(plot.tag = element_text(size = 15, face = "bold")))
print(panel)
# ------- Save -------
ggsave("figure_taxonomy_venn_wide.png", panel, width = FIG_W, height = FIG_H, dpi = DPI)
ggsave("figure_taxonomy_venn_wide.pdf",  panel, width = FIG_W, height = FIG_H)

