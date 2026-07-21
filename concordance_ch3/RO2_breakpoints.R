# Use this part to investigate breakpoints
# of CNA clusters selected by RO
library(GenomicRanges)
library(dplyr)
library(regioneR)
library(GenomeInfoDb)
library(parallel)
library(tidyverse)
library(plyranges)

source('functions_new.R')
##################################
# load raw CNV or gap.merged
##################################
use.merged.cnv = T
gr.cll.cnv <- load.raw_gap(snp.ct = 50, len = 1000, use.merged.cnv = use.merged.cnv, only.common.sam = T)

gr.cll.cnv$cnv_uid <- seq_along(gr.cll.cnv)

##################
# load RO data and
# combine RO pairs
##########################
ro_pairs <- readRDS(paste0("RO_pairs_",ifelse(use.merged.cnv,'gapMerged','raw'),".Rds"))

# filter 
ro_threshold <- 0.7

ro_pairs_filt <- ro_pairs %>% filter(reciprocal_overlap >= ro_threshold)

# make RO clusters
ro_clusters <- make_ro_clusters(ro_pairs_filt)
nrow(ro_clusters)
length(unique(ro_clusters$ro_cluster_id))


# add pos 
cnv_pos <- data.frame(cnv_uid = gr.cll.cnv$cnv_uid,
  chr = as.character(seqnames(gr.cll.cnv)),
  start = start(gr.cll.cnv), end = end(gr.cll.cnv),
  width = width(gr.cll.cnv),  caller = gr.cll.cnv$caller,
  sample_id = gr.cll.cnv$sample_id,
  loss_gain = gr.cll.cnv$loss_gain)

cluster_members <- ro_clusters %>%
  left_join(cnv_pos, by = c("cnv_uid", "sample_id", "loss_gain"))

# now calculate breakpoint var per cluster
breakpoint_variance <- cluster_members %>%
  group_by(sample_id, loss_gain, ro_cluster_id, chr) %>%
  summarise(n_calls = n(), n_callers = n_distinct(caller),
    callers = paste(sort(unique(caller)), collapse = ";"),
    min_start = min(start),
    max_start = max(start),
    median_start = median(start),
    mean_start = mean(start),
    start_range_bp = max(start) - min(start),
    start_sd_bp = sd(start),
    start_iqr_bp = IQR(start),
    min_end = min(end),
    max_end = max(end),
    median_end = median(end),
    mean_end = mean(end),
    end_range_bp = max(end) - min(end),
    end_sd_bp = sd(end),
    end_iqr_bp = IQR(end),
    union_start = min(start),
    union_end = max(end),
    union_width_bp = max(end) - min(start) + 1,
    intersection_start = max(start),
    intersection_end = min(end),
    intersection_width_bp = pmax(0, min(end) - max(start) + 1),
    median_call_width_bp = median(width), .groups = "drop")

# Add normalized breakpoint variance
breakpoint_variance <- breakpoint_variance %>%
  mutate(start_range_fraction = start_range_bp / union_width_bp,
    end_range_fraction = end_range_bp / union_width_bp,
    total_breakpoint_range_bp = start_range_bp + end_range_bp,
    total_breakpoint_range_fraction = total_breakpoint_range_bp / union_width_bp)

# save this for future 
saveRDS(breakpoint_variance, 'breakpoint_variance_out.Rds')


## other plots ##

# summary
breakpoint_summary_by_support <- breakpoint_variance %>%
  group_by(loss_gain, n_callers) %>% summarise(n_clusters = dplyr::n(),
    median_start_range_bp = median(start_range_bp, na.rm = TRUE),
    median_end_range_bp = median(end_range_bp, na.rm = TRUE),
    median_total_breakpoint_range_bp = median(total_breakpoint_range_bp, na.rm = TRUE),
    
    iqr_start_range_bp = IQR(start_range_bp, na.rm = TRUE),
    iqr_end_range_bp = IQR(end_range_bp, na.rm = TRUE),
    
    median_union_width_bp = median(union_width_bp, na.rm = TRUE),
    median_intersection_width_bp = median(intersection_width_bp, na.rm = TRUE),
    .groups = "drop")

breakpoint_summary_by_support %>% gt()

###### initial plots#################
#how many had exact breakpoints
##############################
names(breakpoint_variance)

#1. RO cluster size distribution
ggplot(breakpoint_variance,
       aes(x = median_call_width_bp)) +
  geom_histogram(bins = 50, fill = "steelblue", color = "black") +
  scale_x_log10(labels = comma) +
 # facet_grid(loss_gain ~ n_callers) +
  labs(
    title = "Distribution of RO cluster sizes",
    x = "Median CNV size (bp, log10)",
    y = "Number of RO clusters"
  ) +
  theme_bw()

agreement_df <- breakpoint_variance %>%
  mutate(
    perfect = start_range_bp == 0 & end_range_bp == 0
  ) %>%
  count(loss_gain, n_callers,callers, perfect) %>%
  group_by(loss_gain, n_callers) %>%
  mutate(percent = 100 * n / sum(n))

#2. Percentage of perfect breakpoint agreement
ggplot(agreement_df,
       aes(x = factor(n_callers),
           y = percent,
           fill = perfect)) +
  geom_col(position = "stack") +
  labs(
    title = "Perfect breakpoint agreement across RO clusters",
    x = "Number of callers",
    y = "Percentage of RO clusters",
    fill = "Perfect agreement"
  ) +
  facet_wrap(~loss_gain) +
  theme_bw()

names(breakpoint_variance)

#3. Distribution of breakpoint variance
ggplot(breakpoint_variance,
       aes(x = total_breakpoint_range_bp + 1)) +
  geom_histogram(bins = 50, 
                 fill = "tomato",
                 color = "black") +
  scale_x_log10(labels = comma) +
  scale_y_log10(labels = comma)+
  facet_grid(~n_callers) +
  labs(
    title = "Distribution of breakpoint disagreement",
    x = "Total breakpoint range (bp, log10)",
    y = "Number of RO clusters"
  ) +
  theme_bw()

#4. Breakpoint variance vs CNV size <<< use
break.size <- ggplot(breakpoint_variance,
       aes(x = median_call_width_bp,
           y = total_breakpoint_range_bp + 1,
           colour = factor(n_callers))) +
  geom_point(alpha = 0.6, size = 1.2) +
  scale_x_log10(
    breaks = c(1, 100, 1e4, 1e6, 100e6),
    labels = c("0","0.1 kb", "10 kb", "1 Mb", "100 Mb")
  ) +
  scale_y_log10(
  breaks = c(1, 100, 1e4, 1e6),
  labels = c("0","0.1 kb", "10 kb", "1 Mb")
) +
  facet_wrap(~loss_gain) +
  labs(
    title = "Breakpoint variance versus RO cluster size",
    x = "Median CNV size (bp)",
    y = "Breakpoint variability (bp)",
    colour = "Callers"
  ) + scale_color_brewer(palette = 'Set2')+
  theme_thesis()

break.size

ggsave(filename = "thesis_out/breakpoint_vs_size.png",
       plot = break.size,
       width = 9, height = 6, dpi = 300)

# union width vs callers
ggplot(breakpoint_variance,
       aes(x = factor(n_callers),
           y = intersection_width_bp + 1,
           fill = loss_gain)) +
  geom_boxplot(outlier.alpha = 0.2) +
  scale_y_log10() +
  labs(
    x = "Number of supporting callers",
    y = "Union width (bp, log10 scale)",
    fill = "CNV type",
    title = "Union width of RO-supported CNV clusters"
  ) +
  scale_fill_manual(values = cnv_cols) +
  theme_thesis()

ggplot(breakpoint_variance,
       aes(x = union_width_bp + 1,
           colour = loss_gain)) +
  stat_ecdf(linewidth = 1) +
  scale_x_log10() +
  labs(
    x = "Union width (bp)",
    y = "Cumulative proportion",
    colour = "CNV type",
    title = "Cumulative distribution of RO cluster union widths"
  ) +
  scale_color_manual(values = cnv_cols) +
  theme_thesis()

#5. Relative breakpoint variance
ggplot(breakpoint_variance,
       aes(x = factor(n_callers),
           y = total_breakpoint_range_fraction,
           fill = loss_gain)) +
  geom_boxplot(outlier.alpha = 0.3) +
  scale_y_continuous(labels = percent_format()) +
  labs(
    title = "Relative breakpoint disagreement",
    x = "Number of callers",
    y = "Breakpoint range / CNV size"
  ) +
  theme_bw()

#6. Start vs end breakpoint variability
plot_df <- breakpoint_variance %>%
  dplyr::select(start_range_fraction,
         end_range_fraction) %>%
  tidyr::pivot_longer(
    everything(),
    names_to = "Breakpoint",
    values_to = "Fraction"
  )

ggplot(plot_df,
       aes(Breakpoint, Fraction)) +
  geom_violin(fill = "grey85") +
  geom_boxplot(width = 0.15, outlier.alpha = 0.3) +
  scale_y_continuous(labels = percent_format()) +
  labs(
    title = "Start and end breakpoint variability",
    y = "Breakpoint range / CNV size"
  ) +
  theme_bw()

#7. Breakpoint agreement by caller support
ggplot(breakpoint_variance,
       aes(x = factor(n_callers),
           y = total_breakpoint_range_bp + 1,
           fill = factor(n_callers))) +
  geom_boxplot(outlier.alpha = 0.3) +
  scale_y_log10(labels = comma) +
  labs(
    title = "Breakpoint disagreement by caller support",
    x = "Number of callers",
    y = "Breakpoint range (bp)"
  ) +
  theme_bw()

# or perfect agreement
agreement_summary <- breakpoint_variance %>%
  mutate(perfect = total_breakpoint_range_bp == 0) %>%
  group_by(n_callers) %>%
  summarise(
    proportion = mean(perfect),
    .groups = "drop"
  )

ggplot(agreement_summary,
       aes(factor(n_callers), proportion)) +
  geom_col(fill = "steelblue") +
  scale_y_continuous(labels = percent_format()) +
  labs(
    title = "Perfect breakpoint agreement by caller support",
    x = "Number of callers",
    y = "Proportion of RO clusters"
  ) +
  theme_bw()


#additional <<<< add
ecdf_breakpoint <- ggplot(breakpoint_variance,
       aes(total_breakpoint_range_bp + 1,
           colour = loss_gain)) +
  stat_ecdf(linewidth = 1) +
  facet_wrap(~n_callers)+
  scale_x_log10(
    breaks = c(1, 100, 1e4, 1e6),
    labels = c("0","0.1 kb", "10 kb", "1 Mb")) +
  labs(
    title = "Cumulative distribution of breakpoint disagreement",
    x = "Breakpoint variability (bp)",
    y = "Cumulative proportion",
    color='CNV type'
  ) + scale_color_manual(values = cnv_cols) +
  theme_thesis()

ecdf_breakpoint

ggsave(filename = "thesis_out/breakpoint_ECDF.png",
       plot = ecdf_breakpoint,
       width = 8, height = 5, dpi = 300)

######################################
## plots
#######################################

library(ggplot2)

range(breakpoint_variance$total_breakpoint_range_bp)

ggplot(breakpoint_variance,
  aes(x = factor(n_callers),
    y = total_breakpoint_range_bp+1)) +
  geom_boxplot(outlier.alpha = 0.2) +
  scale_y_log10() +
  facet_wrap(~ loss_gain) +
  labs(x = "Number of supporting callers",
    y = "Total breakpoint range, bp, log10 scale",
    title = "Breakpoint variability of reciprocal-overlap CNV clusters") +
  theme_bw()

# normalized
ggplot(breakpoint_variance,
  aes(x = factor(n_callers),
    y = total_breakpoint_range_fraction)) +
  geom_boxplot(outlier.alpha = 0.2) +
  facet_wrap(~ loss_gain) +
  labs(x = "Number of supporting callers",
    y = "Breakpoint range / cluster union width",
    title = "Relative breakpoint variability of RO-supported CNV clusters") +
  theme_bw()

##################################
## Check breakpoints for RO pairs not clusters
##################################
cnv_gr <- gr.cll.cnv

cnv_pos <- data.frame(
  cnv_uid = cnv_gr$cnv_uid,
  chr = as.character(seqnames(cnv_gr)),
  start = start(cnv_gr),
  end = end(cnv_gr),
  width_bp = width(cnv_gr),
  caller = cnv_gr$caller,
  sample_id = cnv_gr$sample_id,
  loss_gain = cnv_gr$loss_gain
)

ro_pairs_pos <- ro_pairs %>%
  left_join(
    cnv_pos %>%
      rename(
        query_uid = cnv_uid,
        chr_query = chr,
        start_query = start,
        end_query = end,
        width_query_pos = width_bp,
        caller_query_pos = caller
      ),
    by = c("query_uid", "sample_id", "loss_gain")
  )

ro_pairs_pos <- ro_pairs_pos %>%
  left_join(
    cnv_pos %>%
      rename(
        subject_uid = cnv_uid,
        chr_subject = chr,
        start_subject = start,
        end_subject = end,
        width_subject_pos = width_bp,
        caller_subject_pos = caller
      ),
    by = c("subject_uid", "sample_id", "loss_gain")
  )

#2. Calculate pairwise breakpoint variability
ro_pairs_breakpoints <- ro_pairs_pos %>%
  rowwise() %>%
  mutate(
    caller_pair = paste(sort(c(caller_query, caller_subject)), collapse = "_")
  ) %>%
  ungroup() %>%
  mutate(
    start_diff_bp = abs(start_query - start_subject),
    end_diff_bp = abs(end_query - end_subject),
    total_breakpoint_diff_bp = start_diff_bp + end_diff_bp,
    
    pair_union_start = pmin(start_query, start_subject),
    pair_union_end = pmax(end_query, end_subject),
    pair_union_width_bp = pair_union_end - pair_union_start + 1,
    
    pair_intersection_start = pmax(start_query, start_subject),
    pair_intersection_end = pmin(end_query, end_subject),
    pair_intersection_width_bp = pmax(0, pair_intersection_end - pair_intersection_start + 1),
    
    relative_breakpoint_diff =
      total_breakpoint_diff_bp / pair_union_width_bp
  )


#3. Compare across different RO thresholds
ro_thresholds <- c(0.7, 0.9)

bp_by_threshold <- lapply(ro_thresholds, function(thresh) {
  
  ro_pairs_breakpoints %>%
    filter(reciprocal_overlap >= thresh) %>%
    mutate(ro_threshold = thresh)
  
}) %>%
  bind_rows()

#4. Summary table by caller pair and RO threshold 
bp_summary_caller_pair <- bp_by_threshold %>%
  group_by(loss_gain, caller_pair, ro_threshold) %>%
  summarise(
    n_pairs = dplyr::n(),
    
    median_start_diff_bp = median(start_diff_bp, na.rm = TRUE),
    median_end_diff_bp = median(end_diff_bp, na.rm = TRUE),
    median_total_bp_diff = median(total_breakpoint_diff_bp, na.rm = TRUE),
    
    iqr_total_bp_diff = IQR(total_breakpoint_diff_bp, na.rm = TRUE),
    median_relative_bp_diff = median(relative_breakpoint_diff, na.rm = TRUE),
    iqr_relative_bp_diff = IQR(relative_breakpoint_diff, na.rm = TRUE),
    
    median_ro = median(reciprocal_overlap, na.rm = TRUE),
    median_pair_union_width_bp = median(pair_union_width_bp, na.rm = TRUE),
    
    .groups = "drop"
  )

#5. Plot breakpoint variability by caller pair <<<use this

ggplot(
  bp_by_threshold,
  aes(
    x = caller_pair,
    y = total_breakpoint_diff_bp
  )
) +
  geom_boxplot(outlier.alpha = 0.15) +
 # scale_y_log10() +
  facet_grid(~ ro_threshold) +
  labs(
    x = "Caller pair",
    y = "Total breakpoint difference, bp, log10 scale",
    title = "Pairwise breakpoint variability by caller pair and RO threshold"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


# normalized
ggplot(
  bp_by_threshold,
  aes(
    x = caller_pair,
    y = relative_breakpoint_diff
  )
) +
  geom_boxplot(outlier.alpha = 0.15) +
  facet_grid(loss_gain ~ ro_threshold) +
  labs(
    x = "Caller pair",
    y = "(Start difference + end difference) / pair union width",
    title = "Relative pairwise breakpoint variability by caller pair and RO threshold"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

#6. Cleaner line plot of median variability across RO thresholds
ggplot(bp_summary_caller_pair,
  aes(
    x = ro_threshold,
    y = median_relative_bp_diff,
    group = caller_pair,
    colour = caller_pair
  )
) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_wrap(~ loss_gain) +
  labs(
    x = "Reciprocal overlap threshold",
    y = "Median relative breakpoint difference",
    colour = "Caller pair",
    title = "Effect of RO threshold on pairwise breakpoint variability"
  ) +
  theme_thesis() + theme(legend.position = 'right')

#2. Main figure: relative breakpoint variability by caller pair

ggplot(bp_by_threshold,
       aes(x = caller_pair,
           y = relative_breakpoint_diff)) +
  geom_boxplot(outlier.alpha = 0.1) +
  facet_grid(loss_gain ~ ro_threshold) +
  labs(
    x = "Caller pair",
    y = "(Start difference + end difference) / pair union width",
    title = "Relative pairwise breakpoint variability by caller pair and RO threshold"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#3. Threshold sensitivity line plot
ggplot(bp_summary_caller_pair,
       aes(x = ro_threshold,
           y = median_relative_bp_diff,
           group = caller_pair,
           colour = caller_pair)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_wrap(~ loss_gain) +
  labs(
    x = "Reciprocal overlap threshold",
    y = "Median relative breakpoint difference",
    colour = "Caller pair",
    title = "Breakpoint agreement improves with stricter RO thresholds"
  ) +
  theme_bw()

# 4. Raw breakpoint difference plot
ggplot(bp_by_threshold,
       aes(x = caller_pair,
           y = total_breakpoint_diff_bp + 1)) +
  geom_boxplot(outlier.alpha = 0.1) +
  scale_y_log10() +
  facet_grid(loss_gain ~ ro_threshold) +
  labs(
    x = "Caller pair",
    y = "Total breakpoint difference + 1 bp, log10 scale",
    title = "Absolute pairwise breakpoint differences"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#6. Relationship between CNV size and breakpoint difference <<<<< use this

size.vs.breakpoint <- ggplot(bp_by_threshold,
       aes(x = pair_union_width_bp+1,
           y = total_breakpoint_diff_bp+1, col=loss_gain)) +
  geom_point(alpha = 0.15, size = 0.5) +
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap(~ loss_gain) +
  labs(x = "Pair union width, log10 scale",
    y = "Total breakpoint difference, log10 scale",
    title = "Relationship between CNV size and breakpoint variability"
  ) +  theme_thesis() + scale_colour_manual(values = cnv_cols) +
  theme(legend.position = 'none') 

size.vs.breakpoint

ggsave(filename = "thesis_out/CNV size vs breakpoint var.png",
       plot = size.vs.breakpoint,
       width = 8, height = 5, dpi = 300)


#7. Proportion of tightly matched breakpoint pairs
bp_tight_summary <- bp_by_threshold %>%
  group_by(loss_gain, caller_pair, ro_threshold) %>%
  summarise(
    n_pairs = dplyr::n(),
    pct_within_10kb = mean(total_breakpoint_diff_bp <= 10000, na.rm = TRUE) * 100,
    pct_within_50kb = mean(total_breakpoint_diff_bp <= 50000, na.rm = TRUE) * 100,
    pct_within_100kb = mean(total_breakpoint_diff_bp <= 100000, na.rm = TRUE) * 100,
    .groups = "drop"
  )

bp_tight_summary

####################################
# breakpoint agreement analysis << RO is correlated with breakpoint agreement
###################################
library(dplyr)
library(tidyr)
library(ggplot2)

bp_agreement_summary <- bp_by_threshold %>%
  mutate(
    both_exact   = start_diff_bp == 0 & end_diff_bp == 0,
    start_exact  = start_diff_bp == 0,
    end_exact    = end_diff_bp == 0,
    within_10kb  = total_breakpoint_diff_bp <= 10000,
    within_50kb  = total_breakpoint_diff_bp <= 50000,
    within_100kb = total_breakpoint_diff_bp <= 100000
  ) %>%
  group_by(loss_gain, caller_pair, ro_threshold) %>%
  summarise(
    n_pairs = dplyr::n(),
    pct_both_exact = mean(both_exact, na.rm = TRUE) * 100,
    pct_start_exact = mean(start_exact, na.rm = TRUE) * 100,
    pct_end_exact = mean(end_exact, na.rm = TRUE) * 100,
    pct_within_10kb = mean(within_10kb, na.rm = TRUE) * 100,
    pct_within_50kb = mean(within_50kb, na.rm = TRUE) * 100,
    pct_within_100kb = mean(within_100kb, na.rm = TRUE) * 100,
    .groups = "drop"
  )


bp_agreement_counts <- bp_by_threshold %>%
  mutate(
    both_exact   = start_diff_bp == 0 & end_diff_bp == 0,
    start_exact  = start_diff_bp == 0,
    end_exact    = end_diff_bp == 0,
    within_10kb  = total_breakpoint_diff_bp <= 10000,
    within_50kb  = total_breakpoint_diff_bp <= 50000,
    within_100kb = total_breakpoint_diff_bp <= 100000
  ) %>%
  group_by(loss_gain, caller_pair, ro_threshold) %>%
  summarise(
    n_pairs = dplyr::n(),
    n_both_exact = sum(both_exact, na.rm = TRUE),
    n_start_exact = sum(start_exact, na.rm = TRUE),
    n_end_exact = sum(end_exact, na.rm = TRUE),
    n_within_10kb = sum(within_10kb, na.rm = TRUE),
    n_within_50kb = sum(within_50kb, na.rm = TRUE),
    n_within_100kb = sum(within_100kb, na.rm = TRUE),
    .groups = "drop"
  )


# plot
ggplot(bp_agreement_summary,
       aes(x = caller_pair,
           y = pct_both_exact)) +
  geom_col() +
#  facet_grid(loss_gain ~ ro_threshold) +
  labs(
    x = "Caller pair",
    y = "Pairs within 50 kb breakpoint difference (%)",
    title = "Breakpoint agreement by caller pair and RO threshold"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#Alternative: plot several agreement cut-offs together.
bp_agreement_long <- bp_agreement_summary %>%
  dplyr::select(loss_gain, caller_pair, ro_threshold,
         pct_both_exact, pct_within_10kb,
         pct_within_50kb, pct_within_100kb) %>%
  pivot_longer(
    cols = starts_with("pct_"),
    names_to = "agreement_category",
    values_to = "percent"
  ) %>%
  mutate(
    agreement_category = recode(
      agreement_category,
      pct_both_exact = "Exact",
      pct_within_10kb = "Within 10 kb",
      pct_within_50kb = "Within 50 kb",
      pct_within_100kb = "Within 100 kb"
    )
  )

ggplot(bp_agreement_long,
       aes(x = caller_pair,
           y = percent,
           fill = agreement_category)) +
  geom_col(position = "dodge") +
  scale_y_log10(labels = comma) +
  facet_grid(~loss_gain) +
  labs(
    x = "Caller pair",
    y = "RO pairs (%)",
    fill = "Agreement",
    title = "Breakpoint agreement categories by caller pair"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# abs numbers ##############
breakpoints.RO <- ggplot(bp_agreement_counts,
       aes(x = caller_pair,
           y = n_both_exact,
           fill=caller_pair)) +
  geom_col() +
  scale_y_continuous(labels = comma) +
  facet_grid(~ ro_threshold) +
  labs(
    x = "Caller pair",
    y = "Number of RO pairs within 50 kb",
    title = "Number of breakpoint-concordant RO pairs by caller pair"
  ) +
  theme_thesis() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8), 
        legend.position = 'none')+
  scale_fill_brewer(palette = 'Set3')

breakpoints.RO

# save this
ggsave(filename = "thesis_out/Number of breakpoint-concordant RO pairs by caller pair.png",
       plot = breakpoints.RO,
       width = 8, height = 6, dpi = 300)

# type 2
bp_agreement_counts_long <- bp_agreement_counts %>%
  dplyr::select(loss_gain, caller_pair, ro_threshold,
         n_both_exact, n_within_10kb,
         n_within_50kb, n_within_100kb) %>%
  pivot_longer(
    cols = starts_with("n_"),
    names_to = "agreement_category",
    values_to = "n_pairs"
  ) %>%
  mutate(
    agreement_category = recode(
      agreement_category,
      n_both_exact = "Exact",
      n_within_10kb = "Within 10 kb",
      n_within_50kb = "Within 50 kb",
      n_within_100kb = "Within 100 kb"
    )
  )

ggplot(bp_agreement_counts_long,
       aes(x = caller_pair,
           y = n_pairs,
           fill = agreement_category)) +
  geom_col(position = 'dodge') +
  facet_grid(~loss_gain ) +
  scale_y_continuous(labels = comma)+
  labs(
    x = "Caller pair",
    y = "Number of RO pairs",
    fill = "Agreement",
    title = "Breakpoint agreement counts by caller pair"
  ) +
  theme_thesis() + scale_fill_brewer(palette = 'Set2') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
 