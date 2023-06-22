library("optparse")

option_list <- list(
  make_option(c("-i", "--input"),
    type = "character", default = NULL,
    help = "blast file", dest = "input"
  ),
  make_option(c("-k", "--kingdom"),
    type = "character", default = "Eukaryota",
    help = "kingdom of the seed species. Default: Eukaryota", dest = "king"
  ),
  make_option(c("-o", "--output"),
    type = "character", default = "",
    help = "output filename", metavar = "character", dest = "outfile"
  )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# opt <- NULL
# opt$input <- "output_hgt/Adeanei/blast/Adeanei_filtered.blast"
# opt$king <- "Eukaryota"

suppressMessages(library(tidyverse))

blast_results <- read_delim(opt$input)

selfblast <- blast_results %>% 
  group_by(qaccver) %>% 
  summarise(self_bit = min(self_bit))

AI_df <- blast_results %>%
  filter(!is.na(k)) %>%
  select(qaccver, k, n_hits, evalue, bitscore, bit_self) %>%
  mutate(ingroup = ifelse(k == opt$king, "G", "O")) %>%
  group_by(qaccver, ingroup, .drop = FALSE) %>%
  mutate(norm_bitscore = bitscore * exp(1)^(-10 * ((max(bitscore) - bitscore) / bitscore))) %>%
  summarise(bbh = min(evalue), 
            bs = max(bitscore), 
            sum_normbs = sum(norm_bitscore), 
            sum_norm_selfbs = sum(bit_self),
            out_pct = min(n() / n_hits * 100)) %>%
  ungroup() %>%
  # if no obs of other group set pct to 0 and best evalue to 1
  tidyr::complete(qaccver, ingroup, fill = list(bbh = 1, out_pct = 0, bs = 0, sum_normbs = 0, sum_norm_selfbs = 0)) %>%
  pivot_wider(names_from = ingroup, values_from = c(bbh, out_pct, bs, sum_normbs, sum_norm_selfbs)) %>%
  mutate(
    AI = log(bbh_G + 1e-200) - log(bbh_O + 1e-200),
    HGT_index = (bs_O - bs_G),
    AHS = sum_normbs_O - sum_normbs_G, 
    AHS2 = sum_norm_selfbs_O - sum_norm_selfbs_G) %>% 
  left_join(selfblast) %>% 
  mutate(AI2 = HGT_index / self_bit)

# AI index from https://doi.org/10.1126/science.1156407
# log(bbh_G+1e-200)-log(bbh_O+1e-200)
# AI from https://doi.org/10.1093/molbev/msw073
# (bbh_O/self_bit - bbh_G/self_bit)
# HGT_index from https://doi.org/10.1371/journal.pgen.1003035
# AHS from AvP 10.1371/journal.pcbi.1010686

# outAI <- paste0(opt$outfile, gsub("\\..*", "_AI.txt", basename(opt$input)))
outAI <- opt$outfile
write_delim(AI_df, file = outAI, delim = "\t")
