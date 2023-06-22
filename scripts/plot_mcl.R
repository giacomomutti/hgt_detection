#!/usr/bin/env Rscript
suppressMessages(library(tidyverse))
suppressMessages(library(igraph))
suppressMessages(library("optparse"))
suppressMessages(library(patchwork))

opt <- NULL
opt$dump <- "output_hgt/PKR01/selfblast/PKR01_self.dump"
opt$blast <- "output_hgt/PKR01/selfblast/PKR01_self.blast"

mcl_groups <- readLines(opt$dump)
# remomve singletons
mcl_groups <- mcl_groups[which(grepl("\t", mcl_groups))]

ortho_fams <- tibble()

idx <- 1

for (el in mcl_groups[1]){
  genes <- as.character(str_split(el, "\t", simplify = T))
  tmp <- tibble(genes=genes, cluster=rep(idx, length(genes)))
  ortho_fams <- rbind(ortho_fams, tmp)
  idx <- idx + 1
}

blast_cols <- c(
  "qseqid", "sseqid", "pident",
  "length", "mismatch", "gapopen",
  "qstart", "qend", "sstart", "send", "evalue", "bitscore",
  "qcovhsp", "qlen", "slen"
)

df <- read_tsv(opt$blast, col_names = blast_cols) 
df <- df %>% filter(qseqid %in% ortho_fams$genes)

df <- df %>% 
  left_join(ortho_fams, by=c("qseqid" = "genes")) %>% 
  group_by(cluster)

g <- graph_from_data_frame(df, directed = F)
E(g)$weight <- df$evalue
c_scale <- colorRamp(as.character(wesanderson::wes_palette("Zissou1")))

g <- simplify(g, edge.attr.comb = list(weight="mean", "ignore"))

E(g)$color = apply(c_scale(range01(E(g)$weight)), 1, function(x) rgb(x[1]/255, x[2]/255, x[3]/255))

layout_g <- layout_with_graphopt(g, niter = 1000, mass = 70) # , charge = 1e-6, mass = 70)

plot(g,
     vertex.size = 0, vertex.label = NA,
     edge.arrow.size = .00001, layout = layout_g
)
