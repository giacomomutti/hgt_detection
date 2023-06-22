#!/usr/bin/env Rscript
suppressMessages(library(tidyverse))
suppressMessages(library(igraph))
suppressMessages(library("optparse"))
suppressMessages(library(patchwork))

option_list <- list(
  make_option(c("-i", "--input"),
    type = "character", default = NULL,
    help = "blast output, see README for blast command", dest = "blast"
  ),
  make_option(c("-d", "--desc"),
    type = "character", default = NA,
    help = "description dataset for protein: id\tdescription", dest = "desc"
  ),
  make_option(c("-p", "--pident"),
    type = "integer", default = 80,
    help = "minimum percantage of identity. Default 80", dest = "pident"
  ),
  make_option(c("-c", "--cov"),
    type = "double", default = 80,
    help = "minimum percantage of coverage. Default 80", dest = "cov"
  ),
  make_option(c("-e", "--evalue"),
    type = "double", default = 1e-100,
    help = "minimum e value for which an edge is formed. default 1e-100", dest = "minevalue"
  ),
  make_option(c("--byevalue"),
    type = "logical", default = FALSE, action = "store_true",
    help = "if active get network by evalue otherwise by identity and coverage", dest = "byevalue"
  ),
  make_option(c("-m", "--mingenes"),
    type = "integer", default = 8,
    help = "minimum graph size for detailed analsys. Default 8", dest = "mingenes"
  ),
  make_option(c("-o", "--output"),
    type = "character", default = "self_blast",
    help = "output prefix", dest = "outfile"
  )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# testing
opt <- NULL
opt$blast <- "output_hgt/Adeanei/selfblast/Adeanei_self.blast"
# opt$desc <- "data/assemblies/Angomonas_deanei/GCA_903995115.1/protein_desc.txt"
opt$cov <- 0.8
opt$pident <- 80
# opt$outfile <- "self_clustering/plots/Angomonas"
# opt$mingenes <- 8
# opt$minevalue <- 1e-100
# opt$byevalue <- FALSE
args <- commandArgs(trailingOnly = F)
scriptPath <- normalizePath(dirname(sub("^--file=", "", args[grep("^--file=", args)])))

source(paste0(scriptPath, "/functions.R"))
theme_set(theme_bw(base_family = "Helvetica"))


blast_cols <- c(
  "qseqid", "sseqid", "pident",
  "length", "mismatch", "gapopen",
  "qstart", "qend", "sstart", "send", "evalue", "bitscore",
  "qcovhsp", "qlen", "slen"
)


df <- read_tsv(opt$blast, col_names = blast_cols)
norm_df <- df %>%
  filter(qseqid == sseqid) %>%
  mutate(self_bit = bitscore) %>%
  select(qseqid, self_bit)

## you can also use normalized bitscore!

if ("qcovhsp" %in% colnames(df)) {
  # in this case you can compute max_cov
  # remove self hits
  df <- df %>%
    filter(qseqid != sseqid) %>%
    left_join(norm_df) %>%
    # mutate(df, max_cov = length / pmin(qlen, slen)) %>%
    # do this to remove bidirectional hits
    # mutate(qseqid=pmin(qseqid, sseqid),
    #        sseqid=pmax(qseqid, sseqid)) %>%
    group_by(qseqid, sseqid) %>%
    summarise(
      evalue = min(evalue), bit_self = bitscore / self_bit,
      pident = mean(pident), max_cov = mean(qcovhsp)
    ) %>%
    ungroup()
} else {
  df <- df %>%
    # mutate(df, max_cov = length / pmin(qlen, slen)) %>%
    group_by(qseqid, sseqid) %>%
    left_join(norm_df) %>%
    summarise(
      evalue = min(evalue), bit_self = bitscore / self_bit,
      pident = mean(pident), max_cov = 100
    ) %>%
    filter(qseqid != sseqid) %>%
    ungroup()
}

# read_tsv(opt$blast, col_names = blast_cols) %>%
#   filter(qseqid!=sseqid) %>%
#   filter(pident==100, qcovs==100, qlen==slen) %>%
#   arrange(qseqid)

if (!is.na(opt$desc)) {
  # prot <- dendsort::dendsort(fastcluster::hclust(dist(t(matrix_presence), method = "man")))
  # } else {
  prot_desc <- read_tsv(opt$desc, col_names = c("id", "desc"))
  # mutate(desc = gsub("hypothetical protein, conserved", NA, desc),
  #        desc = gsub("hypothetical protein ADEAN_*", NA, desc),
  #        desc = gsub(", putative", "", desc))
}

## how many are only with itself?
ngenes <- length(unique(df$qseqid))
unique_hits <- round(nrow(filter(count(df, qseqid), n == 1)) / ngenes * 100, 2)

print(paste("there are", unique_hits, "% unique hits"))

# df <- df %>% mutate(case = case_when(pident < opt$pident & max_cov < opt$cov ~ "both",
#                                pident < opt$pident & max_cov >= opt$cov ~ "dissimilar",
#                                max_cov < opt$cov & pident >= opt$pident ~ "short",
#                                TRUE ~ "good"))

bn <- gsub("\\..*", "", basename(opt$blast))

if (opt$byevalue) {
  filtered_df <- filter_df(df,
    min_pident = opt$pident,
    min_cov = opt$cov, minevalue = opt$minevalue,
    byevalue = opt$byevalue
  )
  prefix <- paste0(opt$outfile, "/", bn, "_eval", as.character(opt$minevalue))
} else {
  filtered_df <- filter_df(df,
    min_pident = opt$pident,
    min_cov = opt$cov, minevalue = opt$minevalue,
    byevalue = opt$byevalue
  )
  prefix <- paste0(opt$outfile, "/", bn, "_id", as.character(opt$pident), "_cov", opt$cov)
}

df <- left_join(df, filtered_df) %>%
  mutate(case = ifelse(!is.na(weight), TRUE, FALSE))

# what would be the best identity percetage to cut?
x_d <- ggplot(df, aes(pident)) +
  geom_density() +
  # geom_vline(xintercept = opt$pident, color="firebrick", lty=2) +
  scale_x_continuous(breaks = seq(0, 100, 10), minor_breaks = seq(0, 100, 5))

y_d <- ggplot(df, aes(y = max_cov)) +
  geom_density() +
  # geom_hline(yintercept = opt$cov, color="firebrick", lty=2) +
  scale_y_continuous(breaks = seq(0, 100, 10), minor_breaks = seq(0, 100, 5))

scatter <- ggplot(df, aes(pident, max_cov, color = case)) +
  geom_point(alpha = .1) +
  # geom_hline(yintercept = opt$cov, color="firebrick", lty=2) +
  # geom_vline(xintercept = opt$pident, color="firebrick", lty=2) +
  scale_color_manual(values = c("#F23D4C", "#008978")) +
  theme(legend.position = "top")

layout <- "AAAAB
           AAAAB
           AAAAB
           CCCC#"

threshold_plot <- (scatter + y_d + x_d) + plot_layout(design = layout)

ggsave(paste0(prefix, "_threshold_plots.png"), threshold_plot,
  dpi = 300, width = 8, height = 7
)

# we can try with both 80 and 95

blast_graph <- get_network(filtered_df)

# plot(g, vertex.size= 1, edge.arrow.size=0.01, vertex.label=NA)
V(blast_graph)$color <- "grey90"
V(blast_graph)$frame.color <- NA
V(blast_graph)$size <- 0

layout_g <- layout_with_graphopt(blast_graph, niter = 1000, mass = 70) # , charge = 1e-6, mass = 70)

# paste0(opt$outfile, "_", as.character(opt$pident), "_", as.character(opt$cov*100), "_graphopt.pdf")
pdf(paste0(prefix, "_graphopt.pdf"),
  width = 15, height = 15
)
plot(blast_graph,
  vertex.size = 0, vertex.label = NA,
  edge.arrow.size = .00001, layout = layout_g
)
dev.off()

# layout_mds = layout_with_mds(blast_graph)
# pdf("self_clustering/plots/Angomonas_self_mds.pdf", width = 15, height = 15)
# plot(blast_graph, vertex.size=0, vertex.label=NA, edge.arrow.size=.01, layout=layout_mds)
# dev.off()

# you can visualize big families in more detail
split_graph <- decompose.graph(blast_graph)
lengths_graphs <- sapply(split_graph, length)

# ggplot(enframe(lengths_graphs), aes(value)) + geom_histogram(binwidth = 1) +
#   scale_x_continuous(breaks = seq(0,100,1))

V(blast_graph)$size <- 1
split_graph <- decompose.graph(blast_graph, min.vertices = opt$mingenes)
blast_subgraphs <- disjoint_union(split_graph)

# set names so that each cluster is represented by unique names, to avoid overlapping!
clust <- 1
clust_ids <- NULL
for (graph in split_graph) {
  clust_ids <- rbind(clust_ids, data.frame(id = V(graph)$name, cluster = rep(clust, vcount(graph))))
  clust <- clust + 1
}

if (!is.na(opt$desc)) {
  names_clusters <- clust_ids %>%
    left_join(prot_desc) %>%
    filter(!is.na(desc)) %>%
    group_by(cluster, desc) %>%
    summarise(desc = str_wrap(min(desc), width = 50), id = min(id))

  V(blast_subgraphs)$name <- names_clusters[match(V(blast_subgraphs)$name, names_clusters$id), ]$desc
}

l <- layout_with_graphopt(blast_subgraphs, mass = 40)

pdf(paste0(prefix, "_gns", opt$mingenes, "_description_graphopt.pdf"),
  width = 10, height = 10
)
plot(blast_subgraphs,
  layout = l,
  vertex.label.family = 'Helvetica',
  vertex.label.color = 'black',
  vertex.label.cex = .3
)
dev.off()

# find a way to plot the legend!
# image(, 1, as.matrix(seq_along(svals)), col=colors, axes=FALSE, xlab="", ylab="")
