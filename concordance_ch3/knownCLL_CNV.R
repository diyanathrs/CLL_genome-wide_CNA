# test overlap with known CNVs in CLL
library(dplyr)
library(GenomicRanges)
library(tidyr)

source('functions_new.R')

cll_cnvs  <- tibble::tribble(
  ~cnv_name,   ~chr,   ~start,     ~end,       ~expected_type,
  "del13q",    "chr13", 48000000,   52000000,   "Loss",
  "tri12",     "chr12", 1,          133851895,  "Gain",
  "del11q",    "chr11", 105000000,  125000000,  "Loss",
  "del17p",    "chr17", 7500000,    7700000,   "Loss",
  "del6q",     "chr6",  90000000,   170000000,  "Loss",
  "gain2p",    "chr2",  1,          50000000,   "Gain",
  "gain8q",    "chr8",  117000000,  146364022,  "Gain",
  "del14q",    "chr14", 95000000,   107349540,  "Loss",
  "tri18",     "chr18", 1,          78077248,   "Gain"
)

target_gr <- GRanges(seqnames = cll_cnvs$chr,
  ranges = IRanges(start = cll_cnvs$start, end = cll_cnvs$end),
  cnv_name = cll_cnvs$cnv_name,
  expected_type = cll_cnvs$expected_type)

cnv_gr <- load.raw_gap(snp.ct = 5,len = 5, use.merged.cnv = T)

names(mcols(cnv_gr))[3] <- 'cnv_type'

cnv_df2 <- data.frame(cnv_gr)

#check overlap
hits <- findOverlaps(cnv_gr, target_gr)

overlap_df <- tibble(cnv_idx = queryHits(hits),
  target_idx = subjectHits(hits),
  sample_id = mcols(cnv_gr)$sample_id[queryHits(hits)],
  caller = mcols(cnv_gr)$caller[queryHits(hits)],
  cnv_type = mcols(cnv_gr)$cnv_type[queryHits(hits)],
  cnv_name = mcols(target_gr)$cnv_name[subjectHits(hits)],
  expected_type = mcols(target_gr)$expected_type[subjectHits(hits)],
  cnv_width = width(cnv_gr)[queryHits(hits)],
  target_width = width(target_gr)[subjectHits(hits)],
  overlap_bp = width(pintersect(cnv_gr[queryHits(hits)], target_gr[subjectHits(hits)]))) %>%  
  mutate(prop_target_covered = overlap_bp / target_width,
    prop_cnv_overlapped = overlap_bp / cnv_width,
    reciprocal_overlap = pmin(prop_target_covered, prop_cnv_overlapped))

# Filter to expected CNV type and overlap threshold
known_cnv_hits <- overlap_df %>%
  filter(cnv_type == expected_type,
    prop_target_covered >= 0.8 | prop_cnv_overlapped >= 0.8)

#summaries
sample_known_cnvs <- known_cnv_hits %>%
  group_by(sample_id, cnv_name, expected_type) %>%
  summarise(
    present = TRUE,
    n_callers = n_distinct(caller),
    callers = paste(sort(unique(caller)), collapse = ";"),
    max_prop_target_covered = max(prop_target_covered, na.rm = TRUE),
    max_prop_cnv_overlapped = max(prop_cnv_overlapped, na.rm = TRUE),
    max_reciprocal_overlap = max(reciprocal_overlap, na.rm = TRUE),
    .groups = "drop")

known_cnv_matrix <- sample_known_cnvs %>%
  dplyr::select(sample_id, cnv_name, present) %>%
  distinct() %>%
  pivot_wider(names_from = cnv_name,
    values_from = present,
    values_fill = FALSE)

# caller support
sample_known_cnvs_2caller <- sample_known_cnvs %>%
  filter(n_callers >= 2)

known_cnv_matrix_2caller <- sample_known_cnvs_2caller %>%
  dplyr::select(sample_id, cnv_name, present) %>%
  distinct() %>%
  pivot_wider(
    names_from = cnv_name,
    values_from = present,
    values_fill = FALSE
  )

# 10. freq table
known_cnv_freq <- sample_known_cnvs %>%
  group_by(cnv_name, expected_type) %>%
  summarise(n_samples = n_distinct(sample_id),
    .groups = "drop"
  ) %>%  mutate(percent = 100 * n_samples / n_distinct(cnv_df2$sample_id)
  ) %>% arrange(desc(n_samples))

known_cnv_freq

# 2 caller
known_cnv_freq_2caller <- sample_known_cnvs_2caller %>%
  group_by(cnv_name, expected_type) %>%
  summarise(
    n_samples = n_distinct(sample_id),
    .groups = "drop"
  ) %>%
  mutate(
    percent = 100 * n_samples / n_distinct(cnv_df2$sample_id)
  ) %>%
  arrange(desc(n_samples))

known_cnv_freq_2caller

# 11 plot
library(ggplot2)

ggplot(known_cnv_freq, aes(x = reorder(cnv_name, n_samples), y = percent)) +
  geom_col(aes(fill=cnv_name)) +
  coord_flip() +
  labs(
    x = "Known CLL CNV",
    y = "Samples with event (%)",
    title = "Frequency of known CLL CNVs in the cohort"
  ) + scale_fill_brewer(palette = 'Paired') + 
  theme_thesis()

