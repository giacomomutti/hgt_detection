library("optparse")

option_list <- list(
  make_option(c("-i", "--input"),
              type = "character", default = NULL,
              help = "hgtector scores", dest = "hgtectordf"
  ),
  make_option(c("-b", "--blast"),
              type = "character", default = NULL,
              help = "blast filtered results", dest = "blast"
  ),
  make_option(c("-a", "--ai"),
              type = "character", default = NULL,
              help = "alien index scores", dest = "aidf"
  ),
  make_option(c("--both"),
              type = "character", default = NULL,
              help = "both ids", dest = "both"
  ),
  make_option(c("--onlyai"),
              type = "character", default = NULL,
              help = "alien index seeds", dest = "ais"
  ),
  make_option(c("--onlyhgtector"),
              type = "character", default = NULL,
              help = "hgtector seeds", dest = "hgtector"
  ),
  make_option(c("-o", "--output"),
              type = "character", default = "",
              help = "output plots directory", metavar = "character", dest = "outfile"
  )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# opt <- NULL
# opt$hgtectordf <- "output_hgt/Adeanei/hgtector/Adeanei_analysis/scores.tsv"
# opt$blast <- "output_hgt/Adeanei/blast/Adeanei_filtered.blast"
# opt$aidf <- "output_hgt/Adeanei/blast/Adeanei_filtered_AI.tsv"
# opt$both <- "output_hgt/Adeanei/ids/Adeanei_both.txt"
# opt$ais <- "output_hgt/Adeanei/ids/Adeanei_AIs.txt"
# opt$hgtector <- "output_hgt/Adeanei/ids/Adeanei_hgtector.txt"
# opt$outfile <- "output_hgt/Adeanei/plots/"

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(patchwork))


colors <- c("black", "#3B9AB2", "#78B7C5", "#EBCC2A", "#E1AF00", "#F21A00", "purple")
breaks <- c(seq(0, 0.99, length.out=6),1)
col_fun = circlize::colorRamp2(breaks, colors)
col_king <- c("Archaea"="#F2B705", "Eukaryota"="#D9043D", "Bacteria"="#033E8C", "Viruses" = "#FC8DCA")
col_hgt <- c("AIs" = "navyblue", "both" = "forestgreen", "hgtector" = "firebrick", "none" = "grey80")

hgtec <- read_delim(opt$hgtectordf, col_types = cols())
ais <- read_delim(opt$aidf, col_types = cols())

both_seeds <- readLines(opt$both)
names(both_seeds) <- rep("both", length(both_seeds))
ai_seeds <- readLines(opt$ais)
names(ai_seeds) <- rep("AIs", length(ai_seeds))
hgtector_seeds <- readLines(opt$hgtector)
names(hgtector_seeds) <- rep("hgtector", length(hgtector_seeds))
none_seeds <- setdiff(hgtec$protein, Reduce(base::union, list(both_seeds, ai_seeds, hgtector_seeds)))
names(none_seeds) <- rep("none", length(none_seeds))

mode <- enframe(c(both_seeds, ai_seeds, hgtector_seeds, none_seeds), name = "label", value="protein")

hgtec <- left_join(hgtec, ais, by=c("protein"="qaccver")) %>% 
  left_join(mode)

id_plot <- ggplot(hgtec, aes(close, distal, color=label)) +
  geom_point(alpha=.5, size=1) +
  scale_color_manual(values = col_hgt) +
  guides(color = guide_legend(override.aes = list(size=5, alpha=1, pch=15))) +
  theme_bw() # base_family = "Helvetica"

scatter <- paste0(opt$outfile, "/scatter_HGT.png")

ggsave(scatter, id_plot, dpi = 300, width = 5, height = 5)

####

blast_results <- read_delim(opt$blast, delim = "\t", col_types = cols())
  # filter(qaccver %in% mode$protein)

hits_protein <- pull(mode[mode$label!="none", "protein"])

blast_results_sub <- blast_results %>%
  filter(qaccver %in% hits_protein) %>% 
  group_by(qaccver, p) %>% 
  summarise(score=max(bit_self)) %>%
  pivot_wider(names_from = p, values_from = score, values_fill = 0) %>% 
  column_to_rownames("qaccver")

bestb <- t(as.matrix(blast_results_sub))

kingdoms <- blast_results %>% 
  filter(qaccver %in% hits_protein) %>% 
  select(k, p) %>% 
  distinct() %>% 
  filter(p %in% rownames(bestb))

bestb <- bestb[as.character(kingdoms$p), hits_protein]

dend_cols <- cluster_within_group(bestb, deframe(mode[mode$label!="none",2:1]))
dend_rows <- cluster_within_group(t(bestb), deframe(kingdoms[,2:1]))

ht <- Heatmap(bestb, "hgt_blast", col = col_fun, 
              row_split = 4,
              column_split = 3,
              column_labels = rep("", ncol(bestb)),
              cluster_columns = dend_cols,
              cluster_rows = dend_rows,
              rect_gp = gpar(col = NA, lwd = 0),
              use_raster = TRUE, raster_quality = 5,
              right_annotation = rowAnnotation(kingdom=deframe(kingdoms[,2:1]), col = list(kingdom=col_king),
                                               annotation_legend_param = list(kingdom = list(direction = "horizontal"))),
              bottom_annotation = columnAnnotation(mode=deframe(mode[mode$label!="none",2:1]), 
                                                   col = list(mode=col_hgt),
                                                   annotation_name_side = "left"),
              row_names_gp = gpar(fontsize = 5))

blast_hm <- paste0(opt$outfile, "/bitscore_tax_heatmap.pdf")
pdf(file=blast_hm, width = 20, height = 20)
draw(ht)
dev.off()
