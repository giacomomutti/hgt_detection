library("optparse")

option_list <- list(
  make_option(c("-i", "--input"),
              type = "character", default = NULL,
              help = "input tree", dest = "input"
  ),
  make_option(c("-a", "--aln"),
              type = "character", default = NULL,
              help = "input alignment", dest = "aln"
  ),
  make_option(c("-t", "--tax"),
              type = "character", default = "data/taxonomy/lng_db_fx.tsv",
              help = "Taxonomy of db", dest = "tax"
  ),
  make_option(c("-o", "--output"),
              type = "character", default = "",
              help = "output filename", metavar = "character", dest = "outfile"
  ),
  make_option(c("-s", "--self_hits"),
              type = "character", default = "",
              help = "self hits ids", metavar = "character", dest = "self"
  ),
  make_option(c("--selfid"),
              type = "character", default = "",
              help = "self taxid", metavar = "integer", dest = "selfid"
  )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

suppressPackageStartupMessages(library(ggtree))
suppressPackageStartupMessages(library(ggtreeExtra))
suppressPackageStartupMessages(library(tidyverse))

color = grDevices::colors()[grep('gr(a|e)y|black', grDevices::colors(), invert = T)]

# opt <- NULL
# opt$input <- "hgt_phylo/Adeanei/tree/CAD2222087.1.nwk"
# opt$aln <- "hgt_phylo/Adeanei/aligned/CAD2222087.1/CAD2222087.1_consensus.fa"
# opt$tax <- "output_hgt/Adeanei/taxonomy/Adeanei_tax.txt"
# opt$outfile <- "test.png"
# opt$self <- "hgt_phylo/Adeanei/ids/CAD2222087.1_self.txt"
# opt$selfid <- "28296"

seed <- gsub(".nwk", "", (basename(opt$input)))
tree <- read.tree(opt$input) #paste0("phylogeny/tree/",seed,".nwk")


if (length(tree$tip.label) < 2){
  print("Less than 2 sequences found, exiting")
  # this is shit but could not find a way to make snakemake understand
  file.create(opt$outfile)
  quit(status=0)
}
  
aln <- ape::as.AAbin(Biostrings::readAAMultipleAlignment(opt$aln))
attr(aln, "dimnames")[[1]] <- gsub(" .*", "", labels(aln))

# self_hits
self_hits <- intersect(readLines(opt$self), tree$tip.label)

# in case some sequences are lost due to trimming
aln <- aln[which(labels(aln) %in% tree$tip.label), ]

tree <- phytools::midpoint.root(tree)
# copied from rrphylo rescaleRR
tree$edge.length <- (tree$edge.length/max(diag(ape::vcv(tree))))

taxonomy <- readr::read_delim(opt$tax, col_names = c("id",  "k", "p", "c", "o", "f", "g", "s"), 
                              col_types = cols())
# map <- read_delim("data/taxonomy/taxid.map", col_names = c("id", "taxid"))
# map$taxid <- as.character(map$taxid)

df <- as_tibble(tree) %>% 
  mutate(is_self=ifelse(label %in% self_hits, TRUE, FALSE)) %>% 
  rowwise() %>% 
  mutate(id=case_when(node<=length(tree$tip.label) & !is_self ~ str_split(label, pattern = "_", simplify = T)[1],
                      is_self ~ opt$selfid),
         id=as.numeric(id)) %>% 
  # left_join(map)
  left_join(taxonomy) %>% 
  mutate(k2=ifelse(label==seed, "seed", k))

color_kingdoms <- c("forestgreen", "#D9043D", "#F2B705", "#033E8C", "#FC8DCA")#, "pink2", "black", "grey")
names(color_kingdoms) <- c("seed", "Eukaryota", "Archaea", "Bacteria", "Viruses")
# names(color_kingdoms) <- c("seed", unique(taxonomy$k))

color_kingdoms <- color_kingdoms[names(color_kingdoms) %in% unique(df$k2)]

p1 <- ggtree(tree) %<+% df + 
  # geom_tiplab(size=2, hjust = -.05) + 
  geom_tippoint(aes(color=k2, shape=is_self)) +
  scale_color_manual(values = color_kingdoms) +
  scale_shape_manual(values = c(19,1)) +
  theme(legend.position = "bottom")

df_text <- df %>% 
  select(label, p, c, o, f, g) %>% 
  pivot_longer(!label)

df_text$name <- factor(df_text$name, levels = c("p", "c", "o", "f", "g"), ordered = TRUE)

p <- p1 + 
  geom_fruit(data=df_text, geom=geom_tile,
             mapping=aes(y=label, x=name, fill=value),
             offset = .15, pwidth = 1.5) +
  scale_fill_manual(values = sample(color)) +
  geom_fruit(data=df_text, geom=geom_text, 
             mapping=aes(y=label, x=name, label=value),
             color = "black", offset = -.75, size=1.5, pwidth = 1.5) + 
  guides(fill="none") +
  ggnewscale::new_scale_fill()

p2 <- msaplot(p, aln, offset = 2.1, width = 3) +
  guides(fill="none")

n_seqs_coeff <- length(tree$tip.label)/13
height_plot <- ifelse(n_seqs_coeff<7,7,n_seqs_coeff)

ggsave(opt$outfile, p2, width = 12, height = height_plot, device = "png", dpi = 300)


# TODO
# seed taxonomy issue if not in db
# you can solve it by adding it 
# echo "309990" | taxonkit reformat --data-dir data/taxonomy -I1  -f '{k}\t{p}\t{c}\t{o}\t{f}\t{g}\t{s}' >> data/taxonomy/lng_db_fx.tsv
# support values
# add AI or hgtector results in title!
# add some blast info?
