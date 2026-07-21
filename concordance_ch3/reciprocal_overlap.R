# Using reciprocal overlap because check_reduced.R showed
# reduced ranges approach isn't working
#########################################################
#CNVs are concordant only if:
#1. Same sample
#2. Same CNV type: loss with loss, gain with gain
#3. Different callers
#4. Reciprocal overlap >= 50%
#########################################################
library(GenomicRanges)
library(dplyr)
library(regioneR)
library(GenomeInfoDb)
library(parallel)
library(tidyverse)

source('functions_new.R')

##########################################################
### load raw CNVs or gap.merged
#######################################################
use.merged.cnv = T
gr.cll.cnv <- load.raw_gap(snp.ct = 50, len = 1000, use.merged.cnv = use.merged.cnv, only.common.sam = T)

################################################
### start RO workflow
#############################################
gr.cll.cnv$cnv_uid <- seq_along(gr.cll.cnv)
table(gr.cll.cnv$loss_gain)
# check uids
gr.cll.cnv[gr.cll.cnv$cnv_uid=='127627']

# split sample-wise
sam_split <- split(gr.cll.cnv, paste(gr.cll.cnv$sample_id, gr.cll.cnv$loss_gain, sep = "__"))

# run reciprocal overlap function
ro_pairs <- mclapply(sam_split, get_ro_pairs, mc.cores = 10) %>% bind_rows()

ro_threshold <- 0.7

concordant_pairs <- ro_pairs %>% filter(reciprocal_overlap >= ro_threshold)

#save RO_hit_df
saveRDS(ro_pairs, paste0("RO_pairs_", ifelse(use.merged.cnv,'gapMerged','raw'),".Rds"))

##########################
# start from here
# first combine RO pairs
##########################
ro_pairs <- readRDS(paste0("RO_pairs_", ifelse(use.merged.cnv,'gapMerged','raw'),".Rds"))

## add nested overlap
ro_pairs2 <- ro_pairs %>% mutate(smaller_width = pmin(width_query, width_subject),
    larger_width  = pmax(width_query, width_subject),
    prop_smaller_overlapped = overlap_bp / smaller_width,
    prop_larger_overlapped  = overlap_bp / larger_width,
    strict_ro_pass = reciprocal_overlap >= 0.5,
    nested_pass = prop_smaller_overlapped >= 0.8,
    relaxed_match = strict_ro_pass | nested_pass,
    size_ratio = larger_width / smaller_width)

ro_pairs2 %>% summarise(total_pairs = dplyr::n(),
    strict_ro_pairs = sum(strict_ro_pass, na.rm = TRUE),
    nested_only_pairs = sum(!strict_ro_pass & nested_pass, na.rm = TRUE),
    relaxed_pairs = sum(relaxed_match, na.rm = TRUE),
    median_size_ratio = median(size_ratio, na.rm = TRUE))

# by caller pair 
ro_pairs2 %>% group_by(caller_query, caller_subject) %>%
  summarise(total_pairs = dplyr::n(),
    strict_ro_pairs = sum(strict_ro_pass, na.rm = TRUE),
    nested_only_pairs = sum(!strict_ro_pass & nested_pass, na.rm = TRUE),
    relaxed_pairs = sum(relaxed_match, na.rm = TRUE),
    median_size_ratio = median(size_ratio, na.rm = TRUE),
    .groups = "drop")

range(ro_pairs2$size_ratio)
# plot - we already know this -no point
ggplot(ro_pairs2, aes(x = log10(size_ratio), y = reciprocal_overlap)) +
  geom_point(aes(col = loss_gain), alpha = 0.1) +
  geom_hline(yintercept = ro_threshold, linetype = 3) +
  facet_grid(~loss_gain) +
  labs(x = "log10(size ratio between larger and smaller CNA)",
    y = "Reciprocal overlap",
    title = "Effect of CNA size imbalance on reciprocal overlap") +
  theme_thesis() + scale_color_manual(values = cnv_cols)+ 
  theme(legend.position = 'none')

 ggplot(ro_pairs2, aes(x = log10(size_ratio), y = prop_smaller_overlapped)) +
  geom_point(alpha = 0.2) +
  geom_hline(yintercept = ro_threshold, linetype = "dashed") +
  labs(x = "log10(size ratio between larger and smaller CNA)",
    y = "Proportion of smaller CNA overlapped",
    title = "Nested overlap captures smaller CNAs within larger CNAs") +
  theme_bw()

#If many points remain high on this plot while reciprocal overlap is low, 
#that confirms a nested CNA problem rather than complete caller disagreement.
#Then summarise this numerically:

ro_pairs2 %>%
  mutate(
    size_ratio_bin = case_when(
      size_ratio < 2 ~ "<2x",
      size_ratio < 10 ~ "2–10x",
      size_ratio < 100 ~ "10–100x",
      TRUE ~ ">100x"
    )
  ) %>%
  group_by(size_ratio_bin) %>%
  summarise(
    n_pairs = dplyr::n(),
    median_ro = median(reciprocal_overlap, na.rm = TRUE),
    strict_ro_pass = mean(strict_ro_pass, na.rm = TRUE),
    nested_pass = mean(nested_pass, na.rm = TRUE),
    relaxed_match = mean(relaxed_match, na.rm = TRUE),
    .groups = "drop"
  )

#####################################################
## compare RO pairs
####################################################

ro_summary_overall <- ro_pairs2 %>%  dplyr::summarise(n_pairs =  dplyr::n(),
    median_ro = median(reciprocal_overlap, na.rm = TRUE),
    mean_ro = mean(reciprocal_overlap, na.rm = TRUE),
    q25_ro = quantile(reciprocal_overlap, 0.25, na.rm = TRUE),
    q75_ro = quantile(reciprocal_overlap, 0.75, na.rm = TRUE),
    min_ro = min(reciprocal_overlap, na.rm = TRUE),
    max_ro = max(reciprocal_overlap, na.rm = TRUE))

ro_summary_overall

ro_summary_by_type <- ro_pairs2 %>%
  group_by(loss_gain) %>% summarise(n_pairs =  dplyr::n(),
    median_ro = median(reciprocal_overlap, na.rm = TRUE),
    mean_ro = mean(reciprocal_overlap, na.rm = TRUE),
    q25_ro = quantile(reciprocal_overlap, 0.25, na.rm = TRUE),
    q75_ro = quantile(reciprocal_overlap, 0.75, na.rm = TRUE), .groups = "drop")

ro_summary_by_type

# summary by caller pair

ro_pairs2 <- ro_pairs %>% mutate(caller_1 = pmin(caller_query, caller_subject),
    caller_2 = pmax(caller_query, caller_subject),
    caller_pair = paste(caller_1, caller_2, sep = " vs "))

ro_summary_by_pair <- ro_pairs2 %>%  group_by(caller_pair) %>%  summarise(n_pairs =  dplyr::n(),
    median_ro = median(reciprocal_overlap, na.rm = TRUE),
    mean_ro = mean(reciprocal_overlap, na.rm = TRUE),
    q25_ro = quantile(reciprocal_overlap, 0.25, na.rm = TRUE),
    q75_ro = quantile(reciprocal_overlap, 0.75, na.rm = TRUE),
    median_overlap_kb = median(overlap_bp / 1000, na.rm = TRUE),
    median_width_query_kb = median(width_query / 1000, na.rm = TRUE),
    median_width_subject_kb = median(width_subject / 1000, na.rm = TRUE),
    .groups = "drop") %>%   arrange(desc(median_ro))

ro_summary_by_pair

# summary by caller pair and loss/gain
ro_summary_by_pair_type <- ro_pairs2 %>%
  group_by(caller_pair, loss_gain) %>% summarise(n_pairs =  dplyr::n(),
    median_ro = median(reciprocal_overlap, na.rm = TRUE),
    mean_ro = mean(reciprocal_overlap, na.rm = TRUE),
    q25_ro = quantile(reciprocal_overlap, 0.25, na.rm = TRUE),
    q75_ro = quantile(reciprocal_overlap, 0.75, na.rm = TRUE),
  #  ro_ge_0.5 = sum(reciprocal_overlap >= 0.5, na.rm = TRUE),
    ro_ge_thresh = sum(reciprocal_overlap >= ro_threshold, na.rm = TRUE),
    .groups = "drop") %>%  arrange(caller_pair, loss_gain)

ro_pass_pct <- round(ro_summary_by_pair_type$ro_ge_thresh/ro_summary_by_pair_type$n_pairs*100, 2)

ro_summary_by_pair_type$ro_ge_thresh_pcr <- paste0(ro_summary_by_pair_type$ro_ge_thresh,' (',ro_pass_pct,'%)')

ro_summary_by_pair_type

# use this plot in thesis
ro_summary_by_pair_type %>% dplyr::select(-c(8)) %>%
  gt() %>%
  tab_header(title = "Number of RO-supported pairs by caller combination and CNV type."
  ) %>%
  fmt_number(
    columns = where(is.numeric),
    decimals = 2
  ) %>%
  cols_align(
    align = "left",
    columns = everything()
  )

#box plot
ro_dist <- ro_pairs2 %>% ggplot(aes(reciprocal_overlap, caller_pair)) + 
  geom_boxplot(aes(fill=loss_gain),outliers = F) +
  theme_thesis() + scale_fill_manual(values=cnv_cols) +
  labs(title='Reciprocal overlap distribution of caller pairs',
       fill='CNV event',
       x='Reciprocal overlap', y='')

ro_dist

ggsave(filename = "thesis_out/RO_dist_callers_boxplt.png",
       plot = ro_dist,
       width = 7, height = 6, dpi = 300)

# Count pairs passing reciprocal-overlap thresholds
ro_threshold_summary <- ro_pairs2 %>%
  summarise(total_pairs =  dplyr::n(),
    ro_ge_0.1 = sum(reciprocal_overlap >= 0.1, na.rm = TRUE),
    ro_ge_0.3 = sum(reciprocal_overlap >= 0.3, na.rm = TRUE),
    ro_ge_0.5 = sum(reciprocal_overlap >= 0.5, na.rm = TRUE),
    ro_ge_0.7 = sum(reciprocal_overlap >= 0.7, na.rm = TRUE),
    ro_ge_0.9 = sum(reciprocal_overlap >= 0.9, na.rm = TRUE)) %>%
  pivot_longer(cols = starts_with("ro_ge"),
    names_to = "threshold", values_to = "n_pairs") %>%
  mutate(fraction = n_pairs / total_pairs)

ro_threshold_summary

# Threshold summary by caller pair
ro_threshold_by_pair <- ro_pairs2 %>%
  group_by(caller_pair) %>% summarise(total_pairs =  dplyr::n(),
    ro_ge_0.5 = sum(reciprocal_overlap >= 0.5, na.rm = TRUE),
    ro_ge_0.7 = sum(reciprocal_overlap >= 0.7, na.rm = TRUE),
    ro_ge_0.9 = sum(reciprocal_overlap >= 0.9, na.rm = TRUE),
    frac_ge_0.5 = ro_ge_0.5 / total_pairs,
    frac_ge_0.7 = ro_ge_0.7 / total_pairs,
    frac_ge_0.9 = ro_ge_0.9 / total_pairs,
    .groups = "drop") %>% arrange(desc(frac_ge_0.5))

ro_threshold_by_pair

# plot ro distribution
ggplot(ro_pairs2, aes(x = reciprocal_overlap)) +
  geom_histogram(bins = 50) + theme_bw() +
  labs(x = "Reciprocal overlap",
    y = "Number of caller-pair overlaps")

ggplot(ro_pairs2, aes(x = reciprocal_overlap)) +
  geom_histogram(bins = 50) +
  facet_wrap(~ loss_gain) +
  theme_bw() +  labs(x = "Reciprocal overlap",
    y = "Number of caller-pair overlaps")

ro_dist_caller <- ggplot(ro_pairs2, aes(x = reciprocal_overlap, fill = loss_gain)) +
  geom_histogram(bins = 30) +
  facet_wrap(~caller_pair) +
  scale_y_continuous(labels = comma) +
  geom_vline(xintercept = ro_threshold, linetype = 3, colour = 'orange') +
  labs(title = "Reciprocal overlap distribution by caller pair",
    x = "Reciprocal overlap", fill='CNV type',
    y = "Number of caller-pair overlaps") +
  theme_thesis() + scale_fill_manual(values = cnv_cols)
  #theme(legend.position = 'none')

ro_dist_caller

ggsave(filename = "thesis_out/RO_dist_callers.png",
       plot = ro_dist_caller,
       width = 9, height = 6, dpi = 300)


##########################################################
#### Get concordance based on RO - get event support
##########################################################
library(dplyr)
library(igraph)

ro_threshold <- 0.7

ro_pairs_filt <- ro_pairs2 %>%
  filter(reciprocal_overlap >= ro_threshold)

nrow(ro_pairs_filt) / nrow(ro_pairs)*100

# get RO clusters
ro_event_support <- ro_pairs_filt %>%
  group_by(sample_id, loss_gain) %>%
  group_modify(~ get_ro_clusters(.x)) %>%
  ungroup()

nrow(ro_event_support)
# summarize 
support_summary <- ro_event_support %>%
  count(n_callers, name = "n_events") %>% arrange(n_callers)

support_summary
sum(support_summary$n_events)
nrow(ro_event_support)

support_summary_by_type <- ro_event_support %>%
  count(loss_gain, n_callers, name = "n_events") %>%
  group_by(loss_gain) %>%
  mutate(fraction = n_events / sum(n_events)) %>%
  ungroup()

support_summary_by_type %>% gt()
sum(support_summary_by_type$n_events)

# check concordance groups
table(ro_event_support$callers)

data.frame(table(ro_event_support$callers)) %>% ggplot(aes(Var1, Freq)) + 
  geom_bar(stat = 'identity') + theme(axis.text.x = element_text(angle=90))

# plot
library(ggplot2)

ro.clusters <- ggplot(support_summary_by_type,
       aes(x = factor(n_callers), y = n_events, fill = loss_gain)) +
  geom_col(position = "dodge") + theme_bw() + 
  scale_y_continuous(labels = comma)+
  labs(x = "Number of supporting callers",
    y = "Number of reciprocal-overlap CNV events",
    title= "Number of clusters with multicaller support",
    fill = "CNV type") + theme_thesis() +
  scale_fill_manual(values = cnv_cols)

ro.clusters

ggsave(filename = "thesis_out/RO_clusters.png",
       plot = ro.clusters,
       width = 7, height = 5, dpi = 300)


##################################
### use probe counts
##################################
probe_summary_overall <- ro_pairs %>%  summarise(n_pairs = dplyr::n(),
    median_min_probes = median(min_probes, na.rm = TRUE),
    median_max_probes = median(max_probes, na.rm = TRUE),
    median_mean_probes = median(mean_probes, na.rm = TRUE),
    q25_min_probes = quantile(min_probes, 0.25, na.rm = TRUE),
    q75_min_probes = quantile(min_probes, 0.75, na.rm = TRUE),
    min_min_probes = min(min_probes, na.rm = TRUE),
    max_max_probes = max(max_probes, na.rm = TRUE))

probe_summary_overall

# by caller pair
probe_summary_by_pair <- ro_pairs2 %>% group_by(caller_pair) %>%
  summarise(n_pairs = dplyr::n(),
    median_min_probes = median(min_probes, na.rm = TRUE),
    median_max_probes = median(max_probes, na.rm = TRUE),
    median_mean_probes = median(mean_probes, na.rm = TRUE),
    q25_min_probes = quantile(min_probes, 0.25, na.rm = TRUE),
    q75_min_probes = quantile(min_probes, 0.75, na.rm = TRUE),
    .groups = "drop") %>% arrange(desc(median_min_probes))

probe_summary_by_pair

# by loss /gain
probe_summary_by_type <- ro_pairs %>%
  group_by(loss_gain) %>%
  summarise(n_pairs = dplyr::n(),
    median_min_probes = median(min_probes, na.rm = TRUE),
    median_mean_probes = median(mean_probes, na.rm = TRUE),
    median_max_probes = median(max_probes, na.rm = TRUE),
    .groups = "drop")

probe_summary_by_type

# Count RO pairs by probe-support bins
ro_pairs2 <- ro_pairs2 %>%  mutate(min_probe_bin = case_when(
      min_probes < 10 ~ "<10",
      min_probes >= 10 & min_probes < 50 ~ "10-50",
      min_probes >= 50 & min_probes < 200 ~ "50-200",
      min_probes >= 200 & min_probes < 1000 ~ "200-1K",
      min_probes >= 1000 ~ ">=1K",
      TRUE ~ NA_character_),  min_probe_bin = factor(min_probe_bin,
      levels = c("<10", "10-50", "50-200", "200-1K", ">=1K")))

ro_pairs2 <- ro_pairs2 %>%  mutate(max_probe_bin = case_when(
  max_probes < 10 ~ "<10",
  max_probes >= 10 & max_probes < 50 ~ "10-50",
  max_probes >= 50 & max_probes < 200 ~ "50-200",
  max_probes >= 200 & max_probes < 1000 ~ "200-1K",
  max_probes >= 1000 ~ ">=1K",
  TRUE ~ NA_character_),  max_probe_bin = factor(max_probe_bin,
                                                 levels = c("<10", "10-50", "50-200", "200-1K", ">=1K")))

probe_bin_summary <- ro_pairs2 %>% count(max_probe_bin, name = "n_pairs") %>%
  mutate(fraction = n_pairs / sum(n_pairs) )

probe_bin_summary

#Plot number of RO pairs by minimum probe bin
ggplot(probe_bin_summary,  aes(x = max_probe_bin, y = n_pairs)) +
  geom_col() +  theme_bw() +
  facet_grid()+
  labs( x = "Minimum probe count in RO pair",
        y = "Number of reciprocal-overlap pairs",
        title = "Distribution of probe support among RO pairs")


# split by loss/gain
probe_bin_by_type <- ro_pairs2 %>%
  count(loss_gain, min_probe_bin, name = "n_pairs") %>%
  group_by(loss_gain) %>%
  mutate(fraction = n_pairs / sum(n_pairs)) %>%
  ungroup()

probe_bin_by_type

# plot min probs by caller pair
ggplot(ro_pairs2, aes(x = caller_pair, y = min_probes)) +
  geom_boxplot(outlier.alpha = 0.2) +
  scale_y_log10() +  coord_flip() + theme_bw() +
  labs(x = "Caller pair",
    y = "Minimum probe count in RO pair",
    title = "Probe support of reciprocal-overlap CNV pairs" )

# RO vs min probes
ggplot(ro_pairs2, aes(x = min_probes, y = reciprocal_overlap)) +
  geom_point(alpha = 0.3) +
  scale_x_log10() +  theme_bw() +
  labs(x = "Minimum probe count in RO pair",
    y = "Reciprocal overlap",
    title = "Relationship between probe support and reciprocal overlap")

# Concordant pairs passing probe thresholds

probe_threshold_summary <- ro_pairs %>% summarise(total_pairs = dplyr::n(),
    ro_ge_0.5 = sum(reciprocal_overlap >= 0.5, na.rm = TRUE),
    ro_ge_0.5_min_probe_ge_5 =
      sum(reciprocal_overlap >= 0.5 & min_probes >= 5, na.rm = TRUE),
    ro_ge_0.5_min_probe_ge_10 =
      sum(reciprocal_overlap >= 0.5 & min_probes >= 10, na.rm = TRUE),
    ro_ge_0.5_min_probe_ge_20 =
      sum(reciprocal_overlap >= 0.5 & min_probes >= 20, na.rm = TRUE)) %>%
  mutate(frac_ro_ge_0.5 = ro_ge_0.5 / total_pairs,
    frac_ro_ge_0.5_min_probe_ge_5 =
      ro_ge_0.5_min_probe_ge_5 / total_pairs,
    frac_ro_ge_0.5_min_probe_ge_10 =
      ro_ge_0.5_min_probe_ge_10 / total_pairs,
    frac_ro_ge_0.5_min_probe_ge_20 =
      ro_ge_0.5_min_probe_ge_20 / total_pairs)

probe_threshold_summary

###################
### prob diff vs RO
###################
head(ro_pairs2)
ro_pairs2 <- ro_pairs2 %>% mutate(probe_diff = max_probes - min_probes,
    probe_ratio = max_probes / min_probes)

probe_diff_summary <- ro_pairs2 %>%
  summarise(n_pairs = dplyr::n(),
    median_probe_diff = median(probe_diff, na.rm = TRUE),
    mean_probe_diff = mean(probe_diff, na.rm = TRUE),
    q25_probe_diff = quantile(probe_diff, 0.25, na.rm = TRUE),
    q75_probe_diff = quantile(probe_diff, 0.75, na.rm = TRUE),
    max_probe_diff = max(probe_diff, na.rm = TRUE),
    median_probe_ratio = median(probe_ratio, na.rm = TRUE))

probe_diff_summary

# by caller pair
probe_diff_by_pair <- ro_pairs2 %>%
  group_by(caller_pair) %>% summarise(n_pairs = dplyr::n(),
    median_probe_diff = median(probe_diff, na.rm = TRUE),
    median_probe_ratio = median(probe_ratio, na.rm = TRUE),
    median_ro = median(reciprocal_overlap, na.rm = TRUE),
    .groups = "drop") %>%
  arrange(desc(median_probe_diff))

probe_diff_by_pair

# summary by prob diff bin
ro_pairs2 <- ro_pairs2 %>%
  mutate(probe_diff_bin = case_when(
      probe_diff == 0 ~ "0",
      probe_diff <= 5 ~ "1-5",
      probe_diff <= 10 ~ "6-10",
      probe_diff <= 20 ~ "11-20",
      probe_diff <= 50 ~ "21-50",
      probe_diff > 50 ~ ">50",
      TRUE ~ NA_character_),
    probe_diff_bin = factor(
      probe_diff_bin, levels = c("0", "1-5", "6-10", "11-20", "21-50", ">50"))
  )

probe_diff_bin_summary <- ro_pairs2 %>%  group_by(probe_diff_bin) %>%
  summarise(n_pairs = dplyr::n(),
    median_ro = median(reciprocal_overlap, na.rm = TRUE),
    mean_ro = mean(reciprocal_overlap, na.rm = TRUE),
    .groups = "drop")

probe_diff_bin_summary

# plot probe-count difference vs reciprocal overlap
ggplot(ro_pairs2, aes(x = probe_diff + 1,
           y = reciprocal_overlap)) +
  geom_point(alpha = 0.3) +
  scale_x_log10() +
  theme_bw() +
  labs(x = "Probe count difference between matched CNVs (max - min + 1, log10)",
    y = "Reciprocal overlap",
    title = "Probe-count difference versus reciprocal overlap")

# Boxplot of reciprocal overlap by probe-difference bin
ggplot(ro_pairs2, aes(x = probe_diff_bin,
           y = reciprocal_overlap)) +
  geom_boxplot(outlier.alpha = 0.2) +
  theme_bw() +
  labs(x = "Probe count difference between matched CNVs",
    y = "Reciprocal overlap",
    title = "Reciprocal overlap by probe-count difference")


############################
### Caller pair heatmap
#########################
ro_threshold <- 0.7

ro_pairs2 <- ro_pairs %>%
  mutate(
    caller_1 = pmin(caller_query, caller_subject),
    caller_2 = pmax(caller_query, caller_subject),
    caller_pair = paste(caller_1, caller_2, sep = " vs "),
    ro_pass = reciprocal_overlap >= ro_threshold
  )

pairwise_concordance <- ro_pairs2 %>%
  group_by(caller_1, caller_2) %>%
  summarise(
    n_overlapping_pairs = dplyr::n(),
    n_concordant_pairs = sum(ro_pass, na.rm = TRUE),
    fraction_concordant = n_concordant_pairs / n_overlapping_pairs,
    median_ro = median(reciprocal_overlap, na.rm = TRUE),
    .groups = "drop"
  )

pairwise_concordance

ggplot(pairwise_concordance,
       aes(x = caller_1, y = caller_2, fill = n_concordant_pairs)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = n_concordant_pairs), size = 4) +
  theme_thesis() +
  labs(
    x = "Caller",
    y = "Caller",
    fill = "Concordant\nRO pairs",
    title = paste0("Pairwise CNV caller concordance, RO ≥ ", ro_threshold)
  ) 

ggplot(pairwise_concordance,
       aes(x = caller_1, y = caller_2, fill = fraction_concordant)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = round(fraction_concordant, 2)), size = 4) +
  scale_fill_gradient(low = "white", high = "steelblue", limits = c(0, 1)) +
  theme_bw() +
  labs(
    x = "Caller",
    y = "Caller",
    fill = "Fraction",
    title = paste0("Fraction of overlapping CNV pairs passing RO ≥ ", ro_threshold)
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

pairwise_concordance_type <- ro_pairs2 %>%
  group_by(loss_gain, caller_1, caller_2) %>%
  summarise(
    n_overlapping_pairs = dplyr::n(),
    n_concordant_pairs = sum(ro_pass, na.rm = TRUE),
    fraction_concordant = n_concordant_pairs / n_overlapping_pairs,
    median_ro = median(reciprocal_overlap, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(pairwise_concordance_type,
       aes(x = caller_1, y = caller_2, fill = n_concordant_pairs)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = n_concordant_pairs), size = 3.5) +
  facet_wrap(~ loss_gain) +
  theme_bw() +
  labs(
    x = "Caller",
    y = "Caller",
    fill = "Concordant\nRO pairs",
    title = paste0("Pairwise CNV caller concordance by CNV type, RO ≥ ", ro_threshold)
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggplot(pairwise_concordance,
       aes(x = caller_1, y = caller_2, fill = median_ro)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = round(median_ro, 2)), size = 4) +
  scale_fill_gradient(low = "white", high = "steelblue", limits = c(0, 1)) +
  theme_bw() +
  labs(
    x = "Caller",
    y = "Caller",
    fill = "Median RO",
    title = "Median reciprocal overlap by caller pair"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

##########################
### length vs RO
###########################

library(dplyr)
library(ggplot2)

ro_pairs_len <- ro_pairs %>%
  mutate(
    min_width_kb  = pmin(width_query, width_subject, na.rm = TRUE) / 1000,
    max_width_kb  = pmax(width_query, width_subject, na.rm = TRUE) / 1000,
    mean_width_kb = rowMeans(cbind(width_query, width_subject), na.rm = TRUE) / 1000,
    size_ratio = pmin(width_query, width_subject, na.rm = TRUE) /
      pmax(width_query, width_subject, na.rm = TRUE),
    caller_1 = pmin(caller_query, caller_subject),
    caller_2 = pmax(caller_query, caller_subject),
    caller_pair = paste(caller_1, caller_2, sep = " vs ")
  )


  #B. Same plot by loss/gain

ro_pairs_len[ro_pairs_len$mean_width_kb  >= 10000,]

#2. Bin CNVs by length and summarise RO
ro_pairs_len <- ro_pairs_len %>%
  mutate(
    size_bin = case_when(
      mean_width_kb < 2 ~ "<2 kb",
      mean_width_kb >= 2 & mean_width_kb < 10 ~ "2-10 kb",
      mean_width_kb >= 10 & mean_width_kb < 100 ~ "10-100 kb",
      mean_width_kb >= 100 & mean_width_kb < 500 ~ "100-500 kb",
      mean_width_kb >= 500 & mean_width_kb < 1000 ~ "500 kb-1 Mb",
      mean_width_kb >= 1000 & mean_width_kb < 10000 ~ "1 Mb - 10 Mb",
      mean_width_kb >= 10000 ~ ">10 Mb",
      TRUE ~ NA_character_
    ),
    size_bin = factor(
      size_bin, 
      levels = c("<2 kb", "2-10 kb", "10-100 kb",
                 "100-500 kb", "500 kb-1 Mb", "1 Mb - 10 Mb", ">10 Mb")
    )
  )


ro_pairs_len %>% filter(reciprocal_overlap > 0.7) %>% count(caller_pair)
  #summarise(mean(reciprocal_overlap, na.rm=T))
  
ro_by_size <- ro_pairs_len %>% 
  group_by(size_bin, loss_gain) %>%
  summarise(
    n_pairs = dplyr::n(),
    median_ro = median(reciprocal_overlap, na.rm = TRUE),
    mean_ro = mean(reciprocal_overlap, na.rm = TRUE),
    q25_ro = quantile(reciprocal_overlap, 0.25, na.rm = TRUE),
    q75_ro = quantile(reciprocal_overlap, 0.75, na.rm = TRUE),
    frac_ro_ge_0.5 = mean(reciprocal_overlap >= 0.5, na.rm = TRUE),
    frac_ro_ge_0.7 = mean(reciprocal_overlap >= 0.7, na.rm = TRUE),
    .groups = "drop"
  )

ro_by_size %>% gt() %>% 
  fmt_number(
    columns = where(is.numeric),
    decimals = 2
  ) %>%
  cols_align(
    align = "left",
    columns = everything()
  )

# plot
ggplot(ro_by_size,
       aes(x = size_bin, y = mean_ro)) +
  geom_col() +
  theme_bw() +
  labs(
    x = "CNV length bin",
    y = "Median reciprocal overlap",
    title = "Median reciprocal overlap by CNV length"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggplot(ro_by_size,
       aes(x = size_bin, y = frac_ro_ge_0.7)) +
  geom_col() +
  theme_bw() +
  labs(
    x = "CNV length bin",
    y = "Fraction of pairs with RO ≥ 0.7",
    title = "Caller concordance by CNV length"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ro_by_size_type <- ro_pairs_len %>%  filter(reciprocal_overlap >= 0.7) %>%
  group_by(loss_gain, size_bin, caller_pair) %>%
  summarise(
    n_pairs = dplyr::n(),
    mean_ro = mean(reciprocal_overlap, na.rm = TRUE),
    frac_ro_ge_0.5 = mean(reciprocal_overlap >= 0.5, na.rm = TRUE),
    frac_ro_ge_0.7 = mean(reciprocal_overlap >= 0.7, na.rm = TRUE),
    .groups = "drop"
  )

view(ro_by_size_type)


ro_length <- ggplot(ro_by_size_type,
       aes(y = n_pairs , x = size_bin, fill = caller_pair)) +
  geom_col(position = "fill") +  
  #scale_y_log10() +
  facet_wrap(~loss_gain) +scale_fill_brewer(palette = 'Dark2')+
  theme_thesis() +
  labs(
    x = "CNV length bin",
    y = "Fraction of concordant pairs (RO ≥ 0.7)",
    fill = "Caller pair",
    title = "Concordant pairs by segment size"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ro_length
# save
ggsave(filename = "thesis_out/RO_vs_length.png",
       plot = ro_length,
       width = 9, height = 6, dpi = 300)
