library("optparse")

option_list <- list(
  make_option(c("-i", "--input"),
    type = "character", default = NULL,
    help = "blast file. REQUIRES: -outfmt '6 std qcovhsp qlen slen staxids'", dest = "input"
  ),
  make_option(c("-t", "--taxonomy"),
    type = "character", default = NA,
    help = 'taxonomy file, you can get it with cut -f2 blast.txt  | cut -f1 -d\'_\' | sort | uniq | taxonkit reformat -I1  -f "{k}\\t{p}\\t{c}\\t{o}\\t{f}\\t{g}\\t{s}" > tax.txt', dest = "tax"
  ),
  make_option(c("-a", "--auto"),
    type = "character", default = NA,
    help = "self blast file to get self bitscore", dest = "self"
  ),
  make_option(c("-s", "--seed_id"),
    type = "character", default = NULL,
    help = "taxid of the seed species", dest = "seed"
  ),
  make_option(c("-o", "--output"),
    type = "character", default = "",
    help = "output filename", metavar = "character", dest = "outfile"
  ),
  make_option(c("-e", "--evalue"),
              type = "numeric", default = 1e-5,
              help = "max evalue", dest = "evalue"
  ),
  make_option(c("-l", "--len_ratio"),
              type = "numeric", default = 10,
              help = "max lenght ratio", dest = "len_ratio"
  ),
  make_option(c("-c", "--cov"),
              type = "numeric", default = 25,
              help = "coverage", dest = "cov"
  )
)

desc <- "Computes new variables from BLAST output in order to better filter the results"

opt_parser <- OptionParser(option_list = option_list, description = desc)
opt <- parse_args(opt_parser)

suppressMessages(library(tidyverse))

# opt <- NULL
# opt$input <- "data/blast/Adeanei.blast"
# opt$tax <- "output_hgt/Adeanei/taxonomy/Adeanei_tax.txt"
# opt$king <- "Eukaryota"
# opt$seed <- 28296
# opt$self <- "output_hgt/Adeanei/selfblast/Adeanei_self.blast"
# opt$cov <- 33
# opt$evalue <- 1e-20
# opt$l <- 10

col_names <- c(
  "qaccver", "saccver", "pident", "length", "mismatch",
  "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore",
  "qocvs", "qcovhsp", "qlen", "slen", "staxids"
)

dmnd_cols <- c(
  "qaccver", "saccver", "pident",
  "length", "mismatch", "gapopen",
  "qstart", "qend", "sstart", "send", "evalue", "bitscore",
  "qcovhsp", "qlen", "slen"
)

# taxonomic information to annotate blast results
tax <- read_delim(opt$tax, col_names = c("staxids", "k", "p", "c", "o", "f", "g", "s", "species"),
                  col_types = cols()) %>%
  mutate_if(is.character, as.factor)

# Get self bit score
selfblast <- read_delim(opt$self, col_names = col_names, col_types = cols()) %>%
  filter(qaccver == saccver) %>%
  mutate(self_bit = bitscore) %>%
  select(qaccver, self_bit)

all_seed <- unique(selfblast$qaccver)

# read original blast file
#                        0-11 12      13   14   15
# It requires: -outfmt '6 std qcovhsp qlen slen staxids' 
blast_results <- read_delim(opt$input, col_names = col_names, delim = "\t", col_types = cols()) %>% 
  # need to remove self hits if sequence is in DB!
  filter((staxids == opt$seed & pident<100 & qlen!=slen) | staxids!=opt$seed)

before <- nrow(blast_results)
before_seeds <- unique(blast_results$qaccver)

# some orphans might be because no hits in blast
cat(paste0(length(setdiff(all_seed, before_seeds)), " are orphans!\n"))

orphans_df <- tibble(id=setdiff(all_seed, before_seeds),
                     reason="no_hits")

# First, filter given evalue, coverage and length ratio
blast_results <- blast_results %>% 
  mutate(len_ratio = qlen / slen) %>% 
  filter(evalue<opt$evalue, qcovhsp>opt$cov, 
         len_ratio<opt$l, len_ratio>1/opt$l)

after <- nrow(blast_results)
after_seeds <- unique(blast_results$qaccver)

# This will remove X rows and some seed genes may lose all hits (orphans)
cat(paste0(after, " out of ", before, " hits (", round(after/before*100),"%) were kept after filtering!\n"))
cat(paste0(length(before_seeds)-length(after_seeds), " became orphans due to filtering\n"))

orphans_df <- rbind(orphans_df,
                    tibble(id=setdiff(before_seeds, after_seeds), reason="filtering"))

# Some query subject pairs might have multiple hits
# multiple_hits <- blast_results %>%
#   group_by(qaccver, saccver) %>% 
#   count() %>% 
#   filter(n>1)

# cat(paste0(length(unique(multiple_hits$qaccver)), " out of ", length(unique(blast_results$qaccver)), " have multiple hits with the same subject\n"))  
# cat("In these cases only the best query subject pair by btiscore will be considered\n")

out_blast_df <- blast_results %>%
  # how to handle same query subject pairs?
  # weighted mean based on coverage
  # except for region (take the widest)
  # group_by(qaccver, saccver) %>%
  # summarise(qaccver=min(qaccver), 
  #           saccver=min(saccver),
  #           pident=weighted.mean(pident, qcovhsp),
  #           length=weighted.mean(length, qcovhsp), 
  #           mismatch=weighted.mean(mismatch, qcovhsp),
  #           gapopen=weighted.mean(gapopen, qcovhsp),
  #           qstart=min(qstart), 
  #           qend=max(qend), 
  #           sstart=min(sstart), 
  #           send=max(send),
  #           evalue=weighted.mean(evalue, qcovhsp),
  #           bitscore=weighted.mean(bitscore, qcovhsp),
  #           qcovhsp=weighted.mean(qcovhsp, qcovhsp),
  #           qlen=min(qlen), 
  #           slen=min(slen),
  #           staxids=min(staxids)) %>% 
  # group_by(qaccver, saccver) %>% 
  # arrange(bitscore) %>% 
  # slice(1) %>% 
  # ungroup(saccver) %>% 
  group_by(qaccver) %>% 
  left_join(selfblast) %>%
  # Compute new variables that may be useful in exploring results
  mutate(
    n_hits = n(),
    evalue = ifelse(evalue == 0, 1e-180, evalue),
    trans_eval = -log10(evalue),
    norm_pident = pident * qcovhsp / 100,
    bit_len = bitscore / pmin(qlen, slen),
    bit_self = bitscore / self_bit,
    norm_bit_self = bit_self * qcovhsp / 100,
    rank = row_number()
  ) %>%
  arrange(qaccver, desc(bitscore)) %>% 
  left_join(tax)


outblast <- opt$outfile

orphans_file <- gsub("\\..*", "_orphans.tsv", outblast)

write_delim(out_blast_df, file = outblast, delim = "\t")
write_delim(orphans_df, file = orphans_file, delim = "\t")

# TODO:
# maybe find other way to summarize multiple query subject hits! BUT I don't think it is a big deal
