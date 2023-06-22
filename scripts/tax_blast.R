library("optparse")

option_list <- list(
  make_option(c("-i", "--input"),
              type = "character", default = NULL,
              help = "input filtered blast", dest = "input"
  ),
  make_option(c("--ids"),
              type = "character", default = NULL,
              help = "directory of ids", dest = "ids"
  ),
  make_option(c("-o", "--output"),
              type = "character", default = NULL,
              help = "output filename", metavar = "character", dest = "outfile"
  )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)


suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ComplexHeatmap))

# opt <- NULL
# opt$input <- "output_hgt/Adeanei/blast/Adeanei_filtered.blast"
# opt$self <- "output_hgt/Adeanei/selfblast/Adeanei_self.blast"
# opt$seed <- 28296

both <- readLines("output_hgt/Adeanei/ids/Adeanei_both.txt")
names(both) <- rep("both", length(both))
AIs <- readLines("output_hgt/Adeanei/ids/Adeanei_AIs.txt")
names(AIs) <- rep("AIs", length(AIs))
hgtector <- readLines("output_hgt/Adeanei/ids/Adeanei_hgtector.txt")
names(hgtector) <- rep("hgtector", length(hgtector))

mode <- enframe(c(both, AIs, hgtector))

blast_results <- read_delim(opt$input, delim = "\t") %>% 
  filter(qaccver %in% mode$value)

blast_results_sub <- blast_results %>%
  group_by(qaccver, p) %>% 
  summarise(score=max(bit_self)) %>%
  pivot_wider(names_from = p, values_from = score, values_fill = 0) %>% 
  column_to_rownames("qaccver")

blast_results %>%
  filter(staxids == opt$seed) %>% 
  group_by(qaccver, p) %>% 
  summarise(score=max(bit_self)) %>% 
  filter(score<.5)

# frm <- ~k/p
# tr <- ape::as.phylo(frm, data = distinct(select(tax, k, p, c)), collapse=FALSE)
# tr$edge.length <- rep(1, nrow(tr$edge))

colors <- c("black", "#3B9AB2", "#78B7C5", "#EBCC2A", "#E1AF00", "#F21A00", "purple")
breaks <- c(seq(0, 0.99, length.out=6),1)
col_fun = circlize::colorRamp2(breaks, colors)
col_king <- c("Archaea"="#F2B705", "Eukaryota"="#D9043D", "Bacteria"="#033E8C", "Viruses" = "#FC8DCA")
col_hgt <- c("AIs" = "navyblue", "both" = "forestgreen", "hgtector" = "firebrick")


bestb <- t(as.matrix(blast_results_sub))

kingdoms <- blast_results %>% 
  select(k, p) %>% 
  distinct() %>% 
  filter(p %in% rownames(bestb))

bestb <- bestb[as.character(kingdoms$p), mode$value]

dend_cols <- cluster_within_group(bestb, deframe(mode[,2:1]))
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
              bottom_annotation = columnAnnotation(mode=deframe(mode[,2:1]), col = list(mode=col_hgt)),
              row_names_gp = gpar(fontsize = 5))

pdf(file=opt$outfile, width = 20, height = 20)
draw(ht)
dev.off()

hgtector_results <- read_tsv("output/hgtector/Adeanei_analysis/scores.tsv") %>% 
  filter(hits!=0)
hgtector_hits <- read_tsv("output/hgtector/Adeanei_analysis/hgts/Adeanei.txt", 
                          col_names = c("qid", "score", "staxids"))

##### Analyze HGT metrics
AI_df <- read_delim("output/blast/Adeanei_parsed_filtered_AI.tsv")

AI_df <- AI_df %>%
  left_join(hgtector_results, by=c("qaccver"="protein")) %>% 
  left_join(hgtector_hits, by=c("qaccver"="qid")) %>% 
  mutate(hit = ifelse(is.na(score), FALSE, TRUE)) 

ggplot(AI_df, aes(AI, AI2, color=out_pct_O)) + 
  geom_point() + 
  scale_color_viridis_c() +
  facet_grid(~hit)

ggplot(AI_df, aes(AHS, AI2, color=out_pct_O)) + 
  geom_point() + 
  scale_color_viridis_c() +
  facet_grid(~hit)

ggplot(AI_df, aes(close, AI, color=out_pct_O)) + 
  geom_point() + 
  scale_color_viridis_c() +
  facet_grid(~hit)

ggplot(hgtector_results, aes(distal)) +
  geom_histogram(binwidth = 20)

fig <- plot_ly(AI_df, x = ~AI, y = ~AI2, z = ~AHS, 
               marker = list(color = ~out_pct_O, 
                             colorscale = c('#FFE1A1', '#683531'), showscale = TRUE))
fig

### LOCATIONS
locations_file <- "data/assemblies/Angomonas_deanei/GCA_903995115.1/proteins.gff" 
genome_index <- "data/assemblies/Angomonas_deanei/GCA_903995115.1/GCA_903995115.1_Adeanei_nanopore_chromosomes_genomic.fna.fai"

locations <- rtracklayer::readGFF(locations_file)

chromosomes <- read_delim(genome_index, col_names = c("contig", "end")) %>% 
  mutate(start=1) %>% 
  select(contig, start, end)

results <- left_join(AI_df, y = as_tibble(locations), by = c("qaccver"="protein_id"))

ggplot() +
  geom_segment(aes(x=start, xend=end, y=contig, yend=contig), data = chromosomes,
               size=3, color="grey") +
  geom_segment(aes(x=start, xend=end, y=seqid, yend=seqid, color=AI2, size=hit), 
               data = filter(results, AI2>0)) +
  scale_color_viridis_c(option = "A") +
  theme(legend.position = "bottom")


