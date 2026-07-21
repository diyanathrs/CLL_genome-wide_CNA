# Use FISH annotations to check sensitivity of calls
# Check FISH results vs each caller
# FISH results vs RO results
library(dplyr)
library(tidyr)
library(stringr)

source('functions_new.R')

##########################
## method 1 - do not use #
##########################
fish_df <- read.csv('All CLL outcome_fish.csv')
names(fish_df)
fish_df <- fish_df %>% dplyr::select(Sample, FISH)
names(fish_df) <- c('sample_id', 'FISH_result')

fish_long <- fish_df %>%
  mutate(FISH_result = str_to_lower(FISH_result),
    FISH_result = str_squish(FISH_result)) %>%
  # Split combined annotations into separate rows
  separate_rows(FISH_result, sep = " ") %>%
  # Clean individual annotation labels
  mutate(FISH_result = str_squish(FISH_result),
    fish_test = case_when(str_detect(FISH_result, "13q|dleu|rb1") ~ "del13q",
      str_detect(FISH_result, "17p|tp53") ~ "del17p",
      str_detect(FISH_result, "11q|atm") ~ "del11q",
      str_detect(FISH_result, "del6q21") ~ "del6q",
      str_detect(FISH_result, "normal|negative|none|no abnormality") ~ "normal",
      TRUE ~ FISH_result)) %>%
  # Remove empty annotations
  filter(!is.na(fish_test),
         fish_test != "",
         fish_test != "na")  %>%
  mutate(expected_cnv = case_when(
    grepl("^del", fish_test) ~ "Loss",
    grepl("^tri", fish_test) ~ "Gain",
    fish_test == "normal" ~ "Normal",
    TRUE ~ NA_character_)) %>%
  distinct(sample_id, fish_test, expected_cnv)

head(fish_long)
table(fish_long$fish_test)
unique(fish_long$sample_id)

# ------------------------------------------------------------
# 1. Starting table
# ------------------------------------------------------------
# fish_long must already have:
# sample_id, fish_test, expected_cnv
#
# Example:
# fish_long
#   sample_id fish_test expected_cnv
#   896_AP    del11q    Loss
#   127_RAN   del11q    Loss


# ------------------------------------------------------------
# 2. Define FISH target regions
#    IMPORTANT: update coordinates to match your genome build
# ------------------------------------------------------------

fish_regions <- data.frame(
  fish_test = c("del11q", "del13q", "tri8", "del17p", "tri12",
                 "del6q", "del8p", "del12", "del14", "tri18"),
  chr = c("chr11", "chr13", "chr8", "chr17", "chr12",
          "chr6", "chr8", "chr12", "chr14", "chr18"),
  start = c(54644206, 19000001, 1, 1, 1,
            61830167, 1, 1, 1, 1),
  end = c(135006516, 115169878, 146364022, 22263006, 133851895,
          171115067, 43838887, 133851895, 107349540, 78077248),
  expected_cnv_region = c("loss", "loss", "gain", "loss", "gain",
           "loss", "loss", "loss", "loss", "gain"))

setdiff(fish_long$fish_test, fish_regions$fish_test)
# ------------------------------------------------------------
# 3. Attach genomic coordinates to FISH annotations
# ------------------------------------------------------------
fish_for_overlap <- fish_long %>%
  left_join(fish_regions, by = "fish_test") %>%
  mutate(expected_cnv = coalesce(expected_cnv, expected_cnv_region)) %>%
  dplyr::select(sample_id, fish_test, expected_cnv, chr, start, end) %>%
  filter(!is.na(chr), !is.na(start), !is.na(end))

table(fish_for_overlap$fish_test)
# ------------------------------------------------------------
# 4. Convert FISH and CNV calls to GRanges
# ------------------------------------------------------------
fish_gr <- GRanges(seqnames = fish_for_overlap$chr,
  ranges = IRanges(start = fish_for_overlap$start,
    end = fish_for_overlap$end),
  sample_id = fish_for_overlap$sample_id,
  fish_test = fish_for_overlap$fish_test,
  expected_cnv = fish_for_overlap$expected_cnv)

# load cll.cnv gr
cnv_gr <- load.raw_gap(snp.ct = 5, len = 5, use.merged.cnv = F)

# ------------------------------------------------------------
# 5. Find all CNV fragments overlapping FISH regions
# ------------------------------------------------------------

hits <- findOverlaps(fish_gr, cnv_gr)

fragment_hits <- tibble(fish_idx = queryHits(hits),
  cnv_idx = subjectHits(hits)) %>%
  mutate(sample_id = mcols(fish_gr)$sample_id[fish_idx],
    sample_cnv = mcols(cnv_gr)$sample_id[cnv_idx],
    
    fish_test = mcols(fish_gr)$fish_test[fish_idx],
    expected_cnv = mcols(fish_gr)$expected_cnv[fish_idx],
    
    caller = mcols(cnv_gr)$caller[cnv_idx],
    cnv_state = mcols(cnv_gr)$loss_gain[cnv_idx],
    
    fish_chr = as.character(seqnames(fish_gr))[fish_idx],
    fish_start = start(fish_gr)[fish_idx],
    fish_end = end(fish_gr)[fish_idx],
    fish_width = width(fish_gr)[fish_idx],
    
    cnv_chr = as.character(seqnames(cnv_gr))[cnv_idx],
    cnv_start = start(cnv_gr)[cnv_idx],
    cnv_end = end(cnv_gr)[cnv_idx],
    cnv_width = width(cnv_gr)[cnv_idx],
    
    overlap_bp = width(pintersect(fish_gr[fish_idx], cnv_gr[cnv_idx])),
    overlap_fraction_fish = overlap_bp / fish_width,
    overlap_fraction_cnv = overlap_bp / cnv_width
  ) %>%  filter(sample_id == sample_cnv,
    expected_cnv == cnv_state
  )

# ------------------------------------------------------------
# 6. Summarise fragmented overlaps per sample / FISH target / caller
# ------------------------------------------------------------

fragment_summary <- fragment_hits %>%
  group_by(sample_id, fish_test, expected_cnv, caller) %>%
  summarise(
    n_fragments = dplyr::n(),
    total_fragment_bp = sum(cnv_width),
    total_overlap_bp = sum(overlap_bp),
    max_single_overlap_fraction = max(overlap_fraction_fish),
    summed_overlap_fraction_fish = sum(overlap_bp) / fish_width[1],
    min_cnv_start = min(cnv_start),
    max_cnv_end = max(cnv_end),
    fragment_span_bp = max_cnv_end - min_cnv_start + 1,
    .groups = "drop"
  )


# ------------------------------------------------------------
# 7. Calculate unique FISH-region coverage by fragmented CNVs
#    This avoids double-counting overlapping CNV fragments
# ------------------------------------------------------------

unique_fragment_coverage <- lapply(seq_len(nrow(fish_for_overlap)), function(i) {
  h <- fragment_hits %>%
    filter(fish_idx == i)
  
  if (nrow(h) == 0) return(NULL)
  
  bind_rows(
    lapply(split(h, h$caller), function(hc) {
      
      cnv_sub <- cnv_gr[hc$cnv_idx]
      fish_sub <- fish_gr[i]
      
      intersections <- pintersect(cnv_sub, fish_sub)
      unique_bp <- sum(width(reduce(intersections)))
      
      tibble(
        sample_id = mcols(fish_sub)$sample_id,
        fish_test = mcols(fish_sub)$fish_test,
        expected_cnv = mcols(fish_sub)$expected_cnv,
        caller = unique(hc$caller),
        unique_overlap_bp = unique_bp,
        fish_width = width(fish_sub),
        unique_overlap_fraction_fish = unique_bp / width(fish_sub)
      )
    })
  )
}) %>% bind_rows()


# ------------------------------------------------------------
# 8. Build complete comparison table
#    Includes callers with no overlap as zero
# ------------------------------------------------------------

all_callers <- sort(unique(cnv_gr$caller))

comparison_frag <- fish_for_overlap %>%
  select(sample_id, fish_test, expected_cnv) %>%
  distinct() %>%
  crossing(caller = all_callers) %>%
  left_join(
    fragment_summary,
    by = c("sample_id", "fish_test", "expected_cnv", "caller")
  ) %>%
  left_join(
    unique_fragment_coverage,
    by = c("sample_id", "fish_test", "expected_cnv", "caller")
  ) %>%
  mutate(
    n_fragments = replace_na(n_fragments, 0),
    total_fragment_bp = replace_na(total_fragment_bp, 0),
    total_overlap_bp = replace_na(total_overlap_bp, 0),
    max_single_overlap_fraction = replace_na(max_single_overlap_fraction, 0),
    summed_overlap_fraction_fish = replace_na(summed_overlap_fraction_fish, 0),
    unique_overlap_bp = replace_na(unique_overlap_bp, 0),
    unique_overlap_fraction_fish = replace_na(unique_overlap_fraction_fish, 0),
    
    cnv_detected_any_fragment = ifelse(n_fragments > 0, 1, 0),
    fragmented_call = ifelse(n_fragments > 1, 1, 0),
    
    cnv_detected_50pct_unique = ifelse(unique_overlap_fraction_fish >= 0.50, 1, 0),
    cnv_detected_20pct_unique = ifelse(unique_overlap_fraction_fish >= 0.20, 1, 0)
  )

# ------------------------------------------------------------
# 9. Simple tables
# ------------------------------------------------------------

# Number of FISH regions with / without CNV detection by caller
comparison_frag %>%
  count(caller, cnv_detected_any_fragment, name = "n_fish_regions")

# Number of fragmented calls by caller
comparison_frag %>%
  count(caller, fragmented_call, name = "n_fish_regions")

# Combined detection and fragmentation table
comparison_frag %>%
  count(caller, cnv_detected_any_fragment, fragmented_call, name = "n_fish_regions")


# ------------------------------------------------------------
# 10. Summary by caller and FISH lesion
# ------------------------------------------------------------

caller_fish_summary <- comparison_frag %>%
  group_by(caller, fish_test, expected_cnv) %>%
  summarise(n_tested = dplyr::n(),
    n_detected_any_fragment = sum(cnv_detected_any_fragment),
    n_detected_20pct_unique = sum(cnv_detected_20pct_unique),
    n_detected_50pct_unique = sum(cnv_detected_50pct_unique),
    n_fragmented = sum(fragmented_call),
    
    detection_rate_any_fragment = n_detected_any_fragment / n_tested,
    detection_rate_20pct_unique = n_detected_20pct_unique / n_tested,
    detection_rate_50pct_unique = n_detected_50pct_unique / n_tested,
    
    fragmentation_rate_among_detected = n_fragmented / n_detected_any_fragment,
    
    mean_fragments_when_detected = mean(n_fragments[n_fragments > 0]),
    median_fragments_when_detected = median(n_fragments[n_fragments > 0]),
    mean_unique_overlap_fraction = mean(unique_overlap_fraction_fish),
    median_unique_overlap_fraction = median(unique_overlap_fraction_fish),
    
    .groups = "drop"
  )

caller_fish_summary


# ------------------------------------------------------------
# 11. Inspect fragmented FISH overlaps
# ------------------------------------------------------------

fragmented_examples <- comparison_frag %>%
  filter(fragmented_call == 1) %>%
  arrange(caller, fish_test, desc(n_fragments))

fragmented_examples


# ------------------------------------------------------------
# 12. Keep the raw fragment-level table for manual inspection
# ------------------------------------------------------------

fragment_hits %>%
  arrange(sample_id, fish_test, caller, cnv_start)

#####################################################
### FISH method 2 >> simple
##################################################

library(dplyr)
library(tidyr)
library(GenomicRanges)

source('functions_new.R')
# load cll.cnv gr
cnv_gr <- load.raw_gap(snp.ct = 5, len = 5, use.merged.cnv = T, only.common.sam = T)

fish_df <- read.csv('All CLL outcome_fish.csv')
names(fish_df)
fish_df <- fish_df %>% dplyr::select(Sample, FISH, Study)
names(fish_df) <- c('sample_id', 'FISH_result', 'Study')

fish_long <- fish_df %>%
  mutate(FISH_result = str_to_lower(FISH_result),
         FISH_result = str_squish(FISH_result)) %>%
  # Split combined annotations into separate rows
  separate_rows(FISH_result, sep = " ") %>%
  # Clean individual annotation labels
  mutate(FISH_result = str_squish(FISH_result),
         fish_test = case_when(str_detect(FISH_result, "13q|dleu|rb1") ~ "del13q",
                               str_detect(FISH_result, "17p|tp53") ~ "del17p",
                               str_detect(FISH_result, "11q|atm") ~ "del11q",
                               str_detect(FISH_result, "del6q21") ~ "del6q",
                               str_detect(FISH_result, "normal|negative|none|no abnormality") ~ "normal",
                               TRUE ~ FISH_result)) %>%
  # Remove empty annotations
  filter(!is.na(fish_test),
         fish_test != "",
         fish_test != "na")  %>%
  mutate(expected_cnv = case_when(
    grepl("^del", fish_test) ~ "Loss",
    grepl("^tri", fish_test) ~ "Gain",
    fish_test == "normal" ~ "Normal",
    TRUE ~ NA_character_)) %>%
  distinct(sample_id, fish_test,Study, expected_cnv)

head(fish_long)
table(fish_long$Study)
length(unique(fish_long$sample_id))

# ------------------------------------------------------------
# 1. Add FISH genomic coordinates
# ------------------------------------------------------------

fish_regions <- tibble::tribble(
  ~fish_test,   ~chr,   ~start,     ~end,       ~expected_cnv_region,
  "del13q",    "chr13", 48000000,   52000000,   "Loss",
  "tri12",     "chr12", 1,          133851895,  "Gain",
  "del11q",    "chr11", 107000000,  110000000,  "Loss",
  "del17p",    "chr17", 1,          7700000,   "Loss",
  "del6q",     "chr6",  110000000,   113000000,  "Loss",
  "gain2p",    "chr2",  1,          50000000,   "Gain",
  "gain8q",    "chr8",  117000000,  146364022,  "Gain",
  "del14q",    "chr14", 95000000,   107349540,  "Loss",
  "tri18", "chr18", 15000000, 18560000, "Gain")

fish_for_overlap <- fish_long %>%
  left_join(fish_regions, by = "fish_test") %>%
  mutate(expected_cnv = coalesce(expected_cnv, expected_cnv_region)) %>%
  dplyr::select(sample_id, fish_test, expected_cnv, chr, start, end) %>%
  filter(!is.na(chr), !is.na(start), !is.na(end))

fish_n <- fish_for_overlap %>%
  distinct(sample_id, fish_test, expected_cnv) %>%
  count(fish_test, expected_cnv, name = "n_fish_samples")

length(unique(fish_for_overlap$sample_id))
sum(fish_n$n_fish_samples)
# ------------------------------------------------------------
# 2. Convert to GRanges
# ------------------------------------------------------------

fish_gr <- GRanges(seqnames = fish_for_overlap$chr,
  ranges = IRanges(
    start = fish_for_overlap$start,
    end = fish_for_overlap$end
  ),
  sample_id = fish_for_overlap$sample_id,
  fish_test = fish_for_overlap$fish_test,
  expected_cnv = fish_for_overlap$expected_cnv
)


# ------------------------------------------------------------
# 3. Find overlapping fragments
# ------------------------------------------------------------

hits <- findOverlaps(fish_gr, cnv_gr)

fragment_hits <- tibble(
  fish_idx = queryHits(hits),
  cnv_idx = subjectHits(hits)) %>%
  mutate(sample_id = mcols(fish_gr)$sample_id[fish_idx],
    sample_cnv = mcols(cnv_gr)$sample_id[cnv_idx],
    
    fish_test = mcols(fish_gr)$fish_test[fish_idx],
    expected_cnv = mcols(fish_gr)$expected_cnv[fish_idx],
    
    caller = mcols(cnv_gr)$caller[cnv_idx],
    cnv_state = mcols(cnv_gr)$loss_gain[cnv_idx],
    
    fish_width = width(fish_gr)[fish_idx],
    cnv_width = width(cnv_gr)[cnv_idx],
    
    overlap_bp = width(pintersect(fish_gr[fish_idx], cnv_gr[cnv_idx]))
  ) %>%
  filter(
    sample_id == sample_cnv,
    expected_cnv == cnv_state
  )


# ------------------------------------------------------------
# 4. Simple fragmented-call summary
#    Sum all overlapping CNV fragment lengths
# ------------------------------------------------------------

fragment_length_summary <- fragment_hits %>%
  group_by(sample_id, fish_test, expected_cnv, caller) %>%
  summarise(n_fragments = dplyr::n(),
    fish_width = fish_width[1],
    
    summed_fragment_length_bp = sum(cnv_width),
    summed_overlap_bp = sum(overlap_bp),
    
    fraction_by_fragment_length = summed_fragment_length_bp / fish_width,
    fraction_by_overlap_length = summed_overlap_bp / fish_width,
    
    .groups = "drop")


# ------------------------------------------------------------
# 5. Add callers with no overlap as zero
# ------------------------------------------------------------

all_callers <- sort(unique(cnv_gr$caller))

comparison_simple_frag <- fish_for_overlap %>%
  dplyr::select(sample_id, fish_test, expected_cnv) %>%
  distinct() %>%
  tidyr::crossing(caller = all_callers) %>%
  left_join(
    fragment_length_summary,
    by = c("sample_id", "fish_test", "expected_cnv", "caller")
  ) %>% mutate(
    n_fragments = replace_na(n_fragments, 0),
    fish_width = replace_na(fish_width, 0),
    summed_fragment_length_bp = replace_na(summed_fragment_length_bp, 0),
    summed_overlap_bp = replace_na(summed_overlap_bp, 0),
    fraction_by_fragment_length = replace_na(fraction_by_fragment_length, 0),
    fraction_by_overlap_length = replace_na(fraction_by_overlap_length, 0),
    
    cnv_detected_any = ifelse(n_fragments > 0, 1, 0),
    fragmented_call = ifelse(n_fragments > 1, 1, 0),
    detected_5pct_overlap = ifelse(fraction_by_overlap_length >= 0.05, 1, 0),
    detected_20pct_overlap = ifelse(fraction_by_overlap_length >= 0.20, 1, 0),
    detected_50pct_overlap = ifelse(fraction_by_overlap_length >= 0.50, 1, 0)
  )

caller_fish_summary <- comparison_simple_frag %>%
  group_by(caller, fish_test, expected_cnv) %>%
  summarise(
    n_tested = dplyr::n(),
    n_detected_any = sum(cnv_detected_any),
    n_detected_5pct = sum(detected_5pct_overlap),
    n_detected_20pct = sum(detected_20pct_overlap),
    n_fragmented = sum(fragmented_call),
    detection_rate_any = n_detected_any / n_tested,
    detection_rate_5pct = n_detected_5pct / n_tested,
    detection_rate_20pct = n_detected_20pct / n_tested,
    mean_fraction_overlap = mean(fraction_by_overlap_length),
    median_fraction_overlap = median(fraction_by_overlap_length),
    .groups = "drop"
  ) %>%
  left_join(fish_n, by = c("fish_test", "expected_cnv")) %>%
  mutate(
    fish_test_label = paste0(fish_test, "\nFISH n=", n_fish_samples)
  )


table(comparison_simple_frag$fish_test)

comparison_simple_frag %>%
  count(caller, cnv_detected_any, fragmented_call, name = "n_fish_regions")

# summarize
comparison_simple_frag %>%
  group_by(caller, fish_test) %>%
  summarise(n_tested = dplyr::n(),
    n_detected_any = sum(cnv_detected_any),
    n_fragmented = sum(fragmented_call),
    mean_n_fragments = mean(n_fragments),
    median_n_fragments = median(n_fragments),
    mean_fraction_overlap = mean(fraction_by_overlap_length),
    median_fraction_overlap = median(fraction_by_overlap_length),
    n_detected_5pct = sum(detected_5pct_overlap),
    n_detected_20pct = sum(detected_20pct_overlap),
    n_detected_50pct = sum(detected_50pct_overlap),
    .groups = "drop" )

#check
test <- comparison_simple_frag %>% group_by(fish_test, sample_id) %>% summarise(dplyr::n())
range(test$`dplyr::n()`)

############################
### for thesis
############################
#Main table: FISH-region detection by caller
fish.table <- comparison_simple_frag %>%
  group_by(caller, fish_test) %>%
  summarise(n_tested = dplyr::n(),
    n_detected_any = sum(cnv_detected_any),
    n_fragmented = sum(fragmented_call),
    mean_n_fragments = mean(n_fragments),
    median_n_fragments = median(n_fragments),
    mean_fraction_overlap = mean(fraction_by_overlap_length),
    median_fraction_overlap = median(fraction_by_overlap_length),
    n_detected_5pct = sum(detected_5pct_overlap),
    n_detected_20pct = sum(detected_20pct_overlap),
    n_detected_50pct = sum(detected_50pct_overlap),
    .groups = "drop"
  )

# for thesis
fish.table <- fish.table %>% dplyr::select(-c(5,6,9)) %>% arrange(fish_test)

fish.table %>% dplyr::select(c(2,1,3,5,6,4,7,9)) %>% gt() %>%
  tab_header(title = "Caller concordance with FISH detected alterations"
  ) %>%
  fmt_number(
    columns = 5,
    decimals = 2
  ) %>%
  cols_align(
    align = "left",
    columns = everything())

comparison_simple_frag_plot <- comparison_simple_frag %>%
  left_join(fish_n, by = c("fish_test", "expected_cnv")) %>%
  mutate(
    fish_test_label = paste0(fish_test, "\nFISH n=", n_fish_samples)
  )

#2. Fragmentation table by caller
fragmentation_table <- comparison_simple_frag %>%
  group_by(caller, fish_test) %>%
  summarise(
    n_tested = dplyr::n(),
    n_detected = sum(cnv_detected_any),
    n_fragmented = sum(fragmented_call),
    fragmentation_rate_all = n_fragmented / n_tested,
    fragmentation_rate_detected = n_fragmented / n_detected,
    
    mean_fragments_when_detected = mean(n_fragments[n_fragments > 0], na.rm = TRUE),
    median_fragments_when_detected = median(n_fragments[n_fragments > 0], na.rm = TRUE),
    max_fragments = max(n_fragments),
    
    .groups = "drop" )

# per sample fish concordnace
sample_level_table <- comparison_simple_frag %>%
  arrange(sample_id, fish_test, caller) %>%
  dplyr::select(sample_id,
    fish_test,
    expected_cnv,
    caller,
    n_fragments,
    fraction_by_overlap_length,
    cnv_detected_any,
    fragmented_call )

#discordant table
missed_by_all_callers <- comparison_simple_frag %>%
  group_by(sample_id, fish_test, expected_cnv) %>%
  summarise(
    n_callers_detected = sum(detected_5pct_overlap),
    callers_detected = paste(caller[detected_5pct_overlap == 1], collapse = ";"),
    max_overlap_fraction = max(fraction_by_overlap_length),
    .groups = "drop"
  ) %>%
  filter(n_callers_detected == 0)

table(missed_by_all_callers$fish_test)

# only one caller 
single_caller_only <- comparison_simple_frag %>%
  group_by(sample_id, fish_test, expected_cnv) %>%
  summarise(
    n_callers_detected = sum(detected_5pct_overlap),
    callers_detected = paste(caller[detected_5pct_overlap == 1], collapse = ";"),
    max_overlap_fraction = max(fraction_by_overlap_length),
    .groups = "drop"
  ) %>%
  filter(n_callers_detected == 1)

table(single_caller_only$fish_test)

###############################
### fish figs
#############################
#library(ggplot2)
names(comparison_simple_frag)

comparison_simple_frag %>%
  ggplot(aes(x = caller, y = cnv_detected_any)) +
  geom_col() +
  facet_wrap(~ fish_test) +
#  ylim(0, 1) + 
  labs(x = "CNV caller",
    y = "Detection rate",
    title = "Detection of FISH-defined lesions by CNV caller",
    subtitle = "Detection defined as ≥50% overlap with FISH target region") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

names(caller_fish_summary)

# version2 - goes in thesis
fisn_detection.main <- caller_fish_summary %>%
  ggplot(aes(x = caller, y = detection_rate_5pct, fill=caller)) +
  geom_col() +
  facet_wrap(~ fish_test_label) +
  ylim(0, 1) +
  labs(
    x = "CNV caller",
    y = "Detection rate",
    title = "Detection of FISH-defined lesions by CNV caller",
    #subtitle = "Detection defined as ≥10% overlap with FISH target region"
  ) +
  theme_thesis() + theme(legend.position = 'none') +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) + scale_fill_manual(values = caller_cols)

ggsave(filename = "thesis_out/FISH_detection_main.png",
       plot = fisn_detection.main,
       width = 8, height = 6, dpi = 300)

#version3
comparison_simple_frag_plot %>%
  filter(cnv_detected_any == 1) %>%
  ggplot(aes(x = caller, y = n_fragments)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.4) +
  facet_wrap(~ fish_test_label) +
  labs(
    x = "CNV caller",
    y = "Number of overlapping CNV fragments",
    title = "Fragmentation of calls overlapping FISH-defined regions"
  ) +
  theme_bw()

#heatmap
comparison_simple_frag %>%
  mutate(fish_sample = paste(sample_id, fish_test, sep = "_")
  ) %>% 
  ggplot(aes(x = caller, y = fish_sample, fill = as.factor(detected_5pct_overlap))) +
  geom_tile() +
  facet_wrap(~ fish_test, scales = "free_y") +
  labs(
    x = "CNV caller",
    y = "Sample–FISH lesion",
    fill = "Detected",
    title = "Sample-level detection of FISH-defined lesions"
  ) +
  theme_bw() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

#5. Caller-support threshold plot

caller_support_fish <- comparison_simple_frag %>%
  group_by(sample_id, fish_test, expected_cnv) %>%
  summarise(
    n_callers_detected = sum(cnv_detected_any),
    n_callers_5pct = sum(detected_5pct_overlap),
    n_callers_50pct = sum(detected_50pct_overlap),
    .groups = "drop")

caller_support_summary <- caller_support_fish %>%
  tidyr::crossing(min_callers = 1:4) %>%
  mutate(detected_at_threshold = n_callers_detected >= min_callers
  ) %>%
  group_by(min_callers) %>%
  summarise(
    n_fish_regions = dplyr::n(),
    n_detected = sum(detected_at_threshold),
    detection_rate = n_detected / n_fish_regions,
    .groups = "drop")

caller_support_summary %>%
  ggplot(aes(x = min_callers, y = detection_rate)) +
  geom_line() +
  geom_point(size = 2) +
  scale_x_continuous(breaks = 1:4) +
  ylim(0, 1) +
  labs(
    x = "Minimum number of supporting callers",
    y = "Detection rate",
    title = "Detection of FISH-defined lesions by caller-support threshold"
  ) +
  theme_bw()

#6. Stacked bar plot: number of callers detecting each FISH lesion - for thesis

caller_support_fish %>%
  ggplot(aes(x = factor(n_callers_5pct), fill=factor(n_callers_5pct))) +
  geom_bar() +
  facet_wrap(~ fish_test) +
  labs(
    x = "Number of CNV callers detecting FISH lesion",
    y = "Number of FISH-positive sample-lesions",
    title = "Multi-caller support for FISH-defined lesions"
  ) +  theme_thesis() + theme(legend.position = 'none')#


#################
## Extra plots
################
#1. Dot plot: detection rate by caller and FISH lesion
caller_fish_summary %>%
  ggplot(aes(x = caller, y = detection_rate_5pct)) +
  geom_point(size = 3) +
  facet_wrap(~ fish_test_label) +
  ylim(0, 1) +
  labs(
    x = "CNV caller",
    y = "Detection rate",
    title = "Detection of FISH-defined lesions by CNV caller"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#2. Heatmap: caller detection rate per FISH lesion
caller_fish_summary %>%
  ggplot(aes(x = caller, y = fish_test_label, fill = detection_rate_5pct)) +
  geom_tile() +
  geom_text(aes(label = scales::percent(detection_rate_5pct, accuracy = 1))) +
  scale_fill_gradient(low = "white", high = "steelblue", limits = c(0, 1)) +
  labs(
    x = "CNV caller",
    y = "FISH lesion",
    fill = "Detection rate",
    title = "FISH-CNV concordance by caller"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#3. Heatmap: number of detected FISH-positive samples
caller_fish_summary %>%
  ggplot(aes(x = caller, y = fish_test_label, fill = n_detected_5pct)) +
  geom_tile() +
  geom_text(aes(label = paste0(n_detected_5pct, "/", n_fish_samples))) +
  labs(
    x = "CNV caller",
    y = "FISH lesion",
    fill = "Detected samples",
    title = "Number of FISH-positive samples detected by CNV callers"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#4. Caller-support plot: how many callers detect each FISH lesion
caller_support_fish <- comparison_simple_frag %>%
  group_by(sample_id, fish_test, expected_cnv) %>%
  summarise(
    n_callers_detected = sum(detected_5pct_overlap),
    .groups = "drop"
  ) %>%
  left_join(fish_n, by = c("fish_test", "expected_cnv")) %>%
  mutate(
    fish_test_label = paste0(fish_test, "\nFISH n=", n_fish_samples)
  )

# 2nd fig for thesis
caller.sup <- caller_support_fish %>%
  ggplot(aes(x = factor(n_callers_detected), fill = factor(n_callers_detected))) +
  geom_bar() +
  facet_wrap(~ fish_test_label) +
  labs(
    x = "Number of callers",
    y = "Number of FISH-positive samples",
    title = "Multi-caller support for FISH-defined lesions"
  ) +  theme_thesis() + theme(legend.position = 'none')+ scale_fill_brewer(palette = 'PuRd') #OrRd

caller.sup


first2 <- cowplot::plot_grid(fisn_detection.main, caller.sup, align = 'h', ncol = 2, labels = 'AUTO')
first2

ggsave(filename = "thesis_out/FISH_first2.png", 
       plot = first2,
       width = 12, height = 6.5, dpi = 300)

#5. Stacked bar: caller-support categories by lesion
caller_support_fish %>%
  count(fish_test_label, n_callers_detected) %>%
  group_by(fish_test_label) %>%
  mutate(prop = n / sum(n)) %>%
  ggplot(aes(x = fish_test_label, y = prop, fill = factor(n_callers_detected))) +
  geom_col() +
  scale_y_continuous(labels = scales::percent) +
  labs(
    x = "FISH lesion",
    y = "Proportion of FISH-positive samples",
    fill = "Number of callers",
    title = "Distribution of CNV caller support for FISH-positive lesions"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#6. Fragmentation rate heatmap
caller_fish_summary %>%
  mutate(
    fragmentation_rate = n_fragmented / n_detected_any
  ) %>%
  ggplot(aes(x = caller, y = fish_test_label, fill = fragmentation_rate)) +
  geom_tile() +
  geom_text(aes(label = scales::percent(fragmentation_rate, accuracy = 1))) +
  scale_fill_gradient(low = "white", high = "firebrick", na.value = "grey90") +
  labs(
    x = "CNV caller",
    y = "FISH lesion",
    fill = "Fragmentation rate",
    title = "Fragmentation of CNV calls overlapping FISH regions"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#7. Boxplot: overlap fraction by caller - need for thesis 3
fish3 <- comparison_simple_frag_plot %>%
  filter(cnv_detected_any == 1) %>%
  ggplot(aes(x = caller, y = fraction_by_overlap_length, fill=caller)) +
  geom_boxplot(outliers = F ) +
  geom_point(alpha=0.3, shape = 21, size = 2) +
 # geom_jitter(width = 0.15, alpha = 0.4)+
  facet_wrap(~ fish_test_label) +
  scale_y_continuous(labels = scales::percent) +
  labs(
    x = "CNV caller",
    y = "Fraction of FISH region overlapped",
    title = "CNV overlap coverage across FISH-defined regions"
  ) +
  theme_thesis() + scale_fill_manual(values = caller_cols) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = 'none')

#8. Scatter plot: number of fragments vs overlap fraction - need for thesis 4
fish4 <- comparison_simple_frag_plot %>%
  filter(cnv_detected_any == 1) %>%
  ggplot(aes(x = n_fragments, y = fraction_by_overlap_length, col=caller)) +
  geom_point(alpha = 0.5, size=2) +
  facet_grid(fish_test_label~caller) +
  scale_y_continuous(labels = scales::percent) +
  scale_x_continuous(n.breaks = 3)+
  labs(
    x = "Number of overlapping CNV fragments",
    y = "Fraction of FISH region overlapped",
    title = "CNV fragmentation and FISH-region coverage"
  ) + scale_color_manual(values = caller_cols)+
  theme_thesis()+ theme(legend.position = 'none', axis.text.x = element_text(angle = 45, hjust = 1))

fish4

#save 3 and 4
fish34 <- cowplot::plot_grid(fish3, fish4, ncol = 2, labels = 'AUTO', rel_widths = c(1,1), axis = 'b')
fish34

ggsave(filename = "thesis_out/FISH_second2.png", 
       plot = fish34,
       width = 12, height = 7, dpi = 300)

#9. Detection threshold curve
threshold_summary <- comparison_simple_frag %>%
  group_by(caller, fish_test) %>%
  summarise(
    any_overlap = mean(cnv_detected_any),
    overlap_5 = mean(detected_5pct_overlap),
    overlap_20 = mean(detected_20pct_overlap),
    overlap_50 = mean(detected_50pct_overlap),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = c(any_overlap, overlap_5, overlap_20, overlap_50),
    names_to = "threshold",
    values_to = "detection_rate"
  ) %>%
  mutate(
    threshold = factor(
      threshold,
      levels = c("any_overlap", "overlap_5","overlap_20", "overlap_50"),
      labels = c("Any overlap", "≥10% overlap","≥20% overlap", "≥50% overlap")
    )
  )

threshold_summary %>%
  ggplot(aes(x = threshold, y = detection_rate, group = caller, col = caller)) +
  geom_line() +
  geom_point(size = 2) +
  facet_wrap(~ fish_test) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(
    x = "CNV detection definition",
    y = "Detection rate",
    linetype = "Caller",
    title = "Effect of overlap threshold on FISH-CNV concordance"
  ) +  theme_thesis()+  scale_color_manual(values = caller_cols) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

##################################
### Plot BAF/LRR plots for all
##################################
# use missed_by=_all table
library(readr)
source('functions_new.R')

       
missed_by_all_callers$chr <- extract_numeric(missed_by_all_callers$fish_test)

for (i in seq_len(nrow(missed_by_all_callers))) {
plot_LRR.BAF(sample = missed_by_all_callers$sample_id[i], 
             chr = missed_by_all_callers$chr[i], 
             loss_gain = missed_by_all_callers$expected_cnv[i])
  }


lapply(seq_len(nrow(missed_by_all_callers)), function(i) {
  plot_LRR.BAF(sample = missed_by_all_callers$sample_id[i],
    chr = missed_by_all_callers$chr[i],
    loss_gain = missed_by_all_callers$expected_cnv[i]
  )
})

table(missed_by_all_callers$fish_test)
unique(missed_by_all_callers$sample_id)
