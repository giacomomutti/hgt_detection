library(tidyverse)
library(patchwork)
library(scales)
library(ComplexHeatmap)

theme_set(theme_bw(base_family = "Helvetica"))

### LOCATIONS
locations_file <- "data/genomes/Adeanei.gff" 
genome_index <- "data/genomes/Adeanei.fna.fai"

locations <- rtracklayer::readGFF(locations_file)

chromosomes <- read_delim(genome_index, col_names = c("contig", "end")) %>% 
  mutate(start=1) %>% 
  select(contig, start, end)

# blast_results <- read_delim("Adeanei_tax.blast", col_names = col_names)# %>% filter(pident>40, qcovhsp>40)
results <- left_join(hgtec, y = as_tibble(locations), by = c("protein"="protein_id"))

plot_chroms <- ggplot() +
  geom_segment(aes(x=start, xend=end, y=contig, yend=contig), data = chromosomes,
               size=3, color="grey") +
  geom_segment(aes(x=start, xend=end, y=seqid, yend=seqid, color=label), 
               data = results, size=3) +
  # facet_grid(vars(hit)) +
  # coord_flip() +
  scale_color_manual(values = col_hgt) +
  theme(legend.position = "bottom")

plot_chroms


results %>% 
  mutate(out_pct_O = ifelse(is.na(out_pct_O), 0, out_pct_O),
         out_pct_G = ifelse(is.na(out_pct_G), 100, out_pct_G)) %>% 
  ggplot() +
  # geom_rect(aes(xmin=start, xmax=end, ymin=0, ymax=1, fill=AI)) +
  geom_line(aes((start+end)/2, y=out_pct_G/100), color="black") +
  facet_grid(vars(seqid), scales = "free") +
  # scale_fill_manual(values = col_hgt) + 
  scale_y_continuous(limits = c(-.1, 1.1)) +
  theme(strip.text.y = element_text(angle = 0))


# tax_lca <- results %>%
#   mutate(kingdom = ifelse(taxid==131567, "cellular organism", kingdom)) %>% 
#   mutate(LCA = coalesce(species, genus, family, order, class, phylum, kingdom)) %>% 
#   mutate(phylum = ifelse(is.na(phylum), kingdom, phylum)) %>% 
#   mutate_if(sapply(tax_lca, is.character), as.factor)
# 
# lca_results <- results %>% 
#   filter(match!=0) %>% 
#   left_join(tax_lca, by = c("match" = "taxid"))
# 
# lca_labels <- lca_results %>% 
#   count(kingdom, phylum) %>% 
#   mutate(new_lca = ifelse(n<10, paste(kingdom), as.character(phylum))) %>% 
#   select(phylum, new_lca)

# tax <- read_delim("data/taxdumps/Adeanei/taxdump/names.dmp", delim = "\t", col_names = c("staxids", "X2", "name")) %>% 
#   select(staxids, name) %>%
#   mutate_if(is.character, as.factor)

plot_chroms_divided <- results %>% 
  left_join(tax, by = c("match"="staxids")) %>% 
  ggplot() +
  geom_segment(aes(x=start, xend=end, y=contig, yend=contig), data = chromosomes,
               size=.5, color="grey") +
  geom_segment(aes(x=start, xend=end, y=seqid, yend=seqid, color=new_lca, size=hit))
  facet_wrap(vars(kingdom), nrow = 3) +
  scale_color_manual(values = c("#E69F00","#56B4E9","#000000",
                                "#009E73","#0072B2","#CC79A7",
                                "#F0E442","#c7519c"))+
  # "#1f83b4","#12a2a8","black",
  #                             "#d63a3a", "#ffbf50", "#ff7f0e","#bcbd22",
  #                             "#6f63bb",
  #                             "#c7519c","#ba43b4","#8a60b0")) +
  # ggthemes::scale_color_tableau("Classic Cyclic") +
  guides(colour = guide_legend(title = "LCA", override.aes = list(linewidth=5))) +
  theme(axis.text.y = element_blank(), legend.position = "bottom")
plot_chroms_divided

plot_genome <- (plot_chroms | plot_chroms_divided)
ggsave("hgt_search/plots/angomonas_hgt_genome.pdf", plot_genome, width = 15, height = 8)

# library(gggenomes)
# 
# test <- read_fai("data/genomes/Adeanei.fna.fai") %>% 
#   mutate(seq_desc = seq_id)
# 
# test <- read_seqs("data/genomes/Adeanei.fna")
# 
# results$seq_id <- results$seqid
# 
# p1 <- gggenomes(seqs = test, genes = filter(results, label!="none")) + 
#   geom_seq() +
#   geom_bin_label() +
#   geom_gene(aes(color=label)) +
#   scale_color_manual(values = col_hgt) +
#   geom_ribbon(aes(x=(x+xend)/2, ymax=y+.24, ymin=y+.38-(.4*score),
#                   group=seq_id, linetype="GC-content"), use_features(emale_gc),
#               fill="blue", alpha=.5)
# 
# ggsave("test.pdf", p1, width = 10, height = 10)
