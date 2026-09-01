#####IPHOP

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


getwd()
setwd("/users/davidrobinson/Desktop/FINAL_PNAS_FILES")
library(dplyr)

# Load IPHOP host predictions
iphop <- read.csv("Host_prediction_to_genome_m90.csv")

# Keep the highest-confidence prediction per viral contig
iphop_top <- iphop %>%
  group_by(Virus) %>%
  slice_max(order_by = Confidence.score, with_ties = FALSE) %>%
  ungroup()

# Optional: write to file
write.csv(iphop_top, "iphop_top_prediction_per_contig.csv", row.names = FALSE)


# Load required package
library(dplyr)

# Read the IPHOP prediction file
iphop <- read.csv("iphop_top_prediction_per_contig.csv", header = TRUE)

# Summarize counts per Family
class_summary <- iphop %>%
  group_by(Class) %>%
  summarise(Count = n()) %>%
  arrange(desc(Count))

# View the summary table
print(class_summary, n = Inf)



#BiocManager::install("vmikk/metagMisc")
getwd()
setwd("/users/davidrobinson/Desktop/FINAL_PNAS_FILES")
otu_mat <- read.delim("all_samples.txt", header=T, row.names="Contig")
tax_mat <- read.delim("iphop_top_prediction_per_contig.txt", header=T, row.names="Virus")
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


VC_colors <- c(
  "#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e",
  "#e6ab02", "#a6761d", "#666666", "#1f78b4", "#b2df8a",
  "#fb9a99", "#fdbf6f", "#cab2d6", "#ffff99", "#6a3d9a",
  "#ff7f00", "#b15928", "#8dd3c7", "#fef0d9", "#bebada",
  "#fb8072", "#80b1d3", "#fdb462", "#b3de69", "#fccde5",
  "#d9d9d9", "#bc80bd", "#ccebc5", "#ffed6f", "grey"
)
viral_VC <- viral_all %>%
  tax_glom(taxrank = "Class") %>%
  transform_sample_counts(function(x) {x / sum(x)} ) %>%
  psmelt() %>%
  mutate(
    Class = ifelse(Abundance < 0.01, "Other (<2%)", as.character(Class))
  ) %>%
  group_by(Sample, Class) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop") %>%
  arrange(Class)

viral_VC$Sample <- recode(viral_VC$Sample,
                                "SLMCELL" = "Mercer Cellular",
                                "SLMVIRAL" = "Mercer Viral",
                                "SLW" = "Whillans Cellular")

p3 <- ggplot(viral_VC, aes(x = Sample, y = Abundance * 100, fill = Class)) +
  geom_bar(stat = "identity", position = "stack", color = "black", linewidth = 0.2) +
  scale_fill_manual(values = VC_colors) +
  theme_classic(base_size = 20) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank()
  ) +
  labs(y = "Relative Abundance (%)", fill = "Host Prediction") 
  #ggtitle("Relative Abundance of Viral Host Predictions")

print(p3)

########
# ---- Side-by-side panel: p3 (A) + PDF image (B) ----
# Packages (install if needed)
# install.packages(c("magick","cowplot","patchwork","ggplot2","grid"))
library(magick)
library(cowplot)
library(patchwork)
library(ggplot2)
library(grid)

pdf_path <- "class_family_lifestyle_squares_aligned_tighter.pdf"

# 1) Read first page of PDF at high DPI, trim margins, and convert to PNG raster
img_pdf <- image_read_pdf(pdf_path, density = 600)[1] |>    # higher DPI = crisper
  image_trim() |>                                           # remove surrounding white space
  image_background("white") |>                              # flatten transparency if present
  image_convert(format = "png")

# Optional: force a maximum height in pixels (keeps consistent sizing)
# Adjust "2400" to match your target figure height
img_pdf <- image_scale(img_pdf, "x2400")

# 2) Build the right panel (B) as a ggdraw image
pB <- ggdraw() + draw_image(img_pdf, x = 0, y = 1, hjust = 0, vjust = 1, scale = 1)

# 3) Make a compact legend from p3 and hide p3’s in-panel legend
leg <- cowplot::get_legend(
  p3 +
    theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.title = element_text(size = 10),
      legend.text  = element_text(size = 9),
      legend.key.size = unit(6, "pt"),
      legend.spacing.x = unit(6, "pt")
    )
)
p3_compact <- p3 + theme(legend.position = "none",
                         plot.margin = margin(5, 5, 5, 5))

# 4) Assemble A|B, with sensible widths and aligned vertical centers
AB <- plot_grid(
  p3_compact,
  pB,
  ncol = 2,
  labels = c("A", "B"),
  label_size = 16,
  rel_widths = c(1.1, 1.3),
  align = "h",
  axis = "tb"
)

# 5) Add the legend underneath the whole panel
panel <- plot_grid(
  AB,
  leg,
  ncol = 1,
  rel_heights = c(1, 0.12)
)

print(panel)

# 6) Save final outputs
ggsave("FIG_host_predictions_panel.pdf", panel,
       width = 12, height = 7.2, dpi = 300, useDingbats = FALSE)
ggsave("FIG_host_predictions_panel.png", panel,
       width = 12, height = 7.2, dpi = 300)
