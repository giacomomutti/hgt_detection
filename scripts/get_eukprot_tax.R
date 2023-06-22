suppressPackageStartupMessages(library(tidyverse))

library("optparse")

option_list <- list(
  make_option(c("-i", "--input"),
              type = "character", default = NULL,
              help = "input ncbi taxonomy of new proteomes, taxonkit lineage -i2 
              where the first column is the mnemo or genome id", dest = "input"
  ),
  make_option(c("-o", "--output"),
              type = "character", default = "",
              help = "output filename", metavar = "character", dest = "outfile"
  )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# unieuk is the most accurate eukaryotic taxonomy currently avaiable
unieuk <- as_tibble(read.table("data/taxonomy/unieuk_all.tsv", col.names = letters[1:14], sep = "\t", fill = TRUE)) %>% 
  mutate_all(na_if,"")

last_clade <- apply(unieuk, 1, function(x) tail(na.omit(x), 1))

# this is the taxonkit results based on ncbi taxonomy, more complete but less accurate
toadd <- read_delim(opt$input, col_names = F, col_types = cols())

# eukprot has already done some more or less arbitrary division in groups
eukprot_ids <- list.files(path = "data/taxonomy", pattern = "Euk*", full.names = T) %>% 
  map_df(~read_delim(., col_types = cols())) %>% 
  select(Supergroup_UniEuk, Taxogroup1_UniEuk, Taxogroup2_UniEuk) %>% 
  distinct()

# domain phylum class order  family genus species
# new_taxonomy <- tibble(sp=as.character(), k="Eukaryota", p=as.character(),
#                        c=as.character(), o=as.character(), f=as.character(), 
#                        g=as.character(), s=as.character())

new_taxonomy <- tibble(sp=as.character(), tax=as.character())

for (species in 1:nrow(toadd)) {
  if (species%%250==0) {
    print(species/nrow(toadd))
  }
  mnemo <- pull(toadd[species, "X1"])
  taxid <- pull(toadd[species, "X2"])

  els <- str_split(pull(toadd[species, 3]), pattern = ";", simplify = T)
  last_tax <- els[length(els)]
  els <- els[3:length(els)]
  found <- FALSE
  
  # not sure if needed but just in case
  new_class <- NA
  new_order <- NA
  new_family <- NA
  new_sk <- NA
  
  for (taxon in rev(els)){
    # in some cases the match was not exact!
    sub_df <- unieuk[which(grepl(paste0("^",taxon,"$"), last_clade)),]
    if (nrow(sub_df) > 0) {
      found <- TRUE

      new_class <- sub_df[which(sub_df %in% eukprot_ids$Supergroup_UniEuk)]
      new_class <- ifelse(ncol(new_class)>0, pull(new_class), NA)
      new_order <- sub_df[which(sub_df %in% eukprot_ids$Taxogroup1_UniEuk)]
      new_order <- ifelse(ncol(new_order)>0, pull(new_order), NA)
      new_family <- sub_df[which(sub_df %in% eukprot_ids$Taxogroup2_UniEuk)]
      new_family <- ifelse(ncol(new_family)>0, pull(new_family), NA)
      
      new_sk <- pull(sub_df[, 'b'])
      new_taxonomy <- add_row(new_taxonomy, sp=mnemo, 
                              tax=paste0("k__Eukaryota;p__",new_sk,";c__",new_class,
                                         ";f__", new_family, ";o__", new_order,
                                         ";g__", str_split(last_tax, " ")[[1]][1],";s__", last_tax))
      break
      }
  }
}

# new_taxonomy <- arrange(new_taxonomy, p,c,f,o,g)

write_delim(new_taxonomy, file = opt$outfile, delim = "\t", col_names = F)
