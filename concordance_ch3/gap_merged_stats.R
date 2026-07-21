# get gap merged stats and plots
cnv_cols <- c("Gain" = "#2C7BB6",
              "Loss" = "#D7191C")


caller_cols <- c("PennCNV"   = "#1B9E77",
                 "QuantiSNP" = "#D95F02",
                 "iPattern"  = "#7570B3",
                 "Nexus"     = "#E7298A")

stage_cols <- c("Raw" = "#D95F02",
                "Gap-merged" = "#1B9E77")


### load raw CNVs or gap.merged
source('functions.R')
source('functions_new.R')

# load all raw CNVs regardless of length or prob counts
cnv_gr <- load.raw_gap(snp.ct = 5, len = 5, use.merged.cnv = F) 
cnv_gr.gap <- load.raw_gap(snp.ct = 5, len = 5, use.merged.cnv = T)


library(dplyr)
library(ggplot2)
library(scales)

cnv_all <- bind_rows(data.frame(cnv_gr) %>% mutate(stage = "Raw"),
  data.frame(cnv_gr.gap) %>% mutate(stage = "Gap-merged"))

cnv_all <- cnv_all %>% 
  mutate(width_kb = width / 1000, width_mb = width / 1e6)

cnv_all$stage <- factor(cnv_all$stage, levels = c('Raw','Gap-merged') )

#check
cnv_all %>%
  group_by(stage, caller, loss_gain) %>%
  summarise(
    n = n(),
    median_width_kb = median(width_kb, na.rm = TRUE),
    mean_width_kb = mean(width_kb, na.rm = TRUE),
    max_width_kb = max(width_kb, na.rm = TRUE),
    .groups = "drop"
  )

cnv_all %>%
  filter(stage == "Gap-merged", n_merged_segments > 1) %>%
  summarise(
    median_width_kb = median(width_kb, na.rm = TRUE),
    median_n_segments = median(n_merged_segments, na.rm = TRUE)
  )

cnv_counts <- cnv_all %>%
  group_by(stage, caller, loss_gain) %>%
  summarise(n_cnvs = n(),
    .groups = "drop")

raw.gap.calls <- ggplot(cnv_counts, aes(x = caller, y = n_cnvs, fill = stage)) +
  geom_col(position = position_dodge(width = 0.8), alpha = 1) +
  facet_wrap(~ loss_gain) +
  scale_y_continuous(labels = comma) +
  labs(x = "CNV caller",
    y = "Number of CNV calls",
    fill = "Stage", 
    title = "CNV call counts before and after gap merging"
  ) + theme_thesis() + scale_fill_manual(values = stage_cols)

raw.gap.calls

ggsave(filename = "thesis_out/raw_gap_counts.png",
       plot = raw.gap.calls,
       width = 9,  height = 5,
       dpi = 300)

# size
ggplot(cnv_all, aes(x = stage, y = width_kb , col = stage)) +
  geom_boxplot(outlier.alpha = 0.15) +
  facet_grid(loss_gain ~ caller) +
  scale_y_log10(labels = comma) +
  labs(x = NULL,
    y = "CNV size (kb, log10 scale)",
    fill = "Processing stage",
    title = "CNV size distributions before and after gap merging") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

##################################################
# Plot C: CNVs per sample before vs after merging
#################################################
cnvs_per_sample <- cnv_all %>%
  group_by(stage, sample_id, caller, loss_gain) %>%
  summarise(n_cnvs = dplyr::n(),
    total_width_mb = sum(width_mb, na.rm = TRUE),
    .groups = "drop")

ggplot(cnvs_per_sample, aes(x = stage, y = n_cnvs, col = stage)) +
  geom_boxplot(outlier.alpha = 0.4) + facet_grid(loss_gain ~ caller) +
  scale_y_continuous(labels = comma) +  labs(x = NULL,  y = "CNVs per sample",
    fill = "Processing stage",
    title = "Per-sample CNV burden before and after gap merging") +
  theme_bw() + scale_fill_manual(values=stage_cols)

# summary table
summary_before_after <- cnv_all %>%
  group_by(stage, caller, loss_gain) %>%
  summarise(n_cnvs = dplyr::n(),
    n_samples = n_distinct(sample_id),
    
    median_cnvs_per_sample = {
      tmp <- cur_data() %>%
        count(sample_id)
      median(tmp$n, na.rm = TRUE)
    },
    
    mean_cnvs_per_sample = {
      tmp <- cur_data() %>%
        count(sample_id)
      mean(tmp$n, na.rm = TRUE)
    },
    
    median_width_kb = median(width_kb, na.rm = TRUE),
    q25_width_kb = quantile(width_kb, 0.25, na.rm = TRUE),
    q75_width_kb = quantile(width_kb, 0.75, na.rm = TRUE),
    
    median_probes = median(no_of_probes, na.rm = TRUE),
    q25_probes = quantile(no_of_probes, 0.25, na.rm = TRUE),
    q75_probes = quantile(no_of_probes, 0.75, na.rm = TRUE),
    
    total_cnv_burden_mb = sum(width_mb, na.rm = TRUE),
    .groups = "drop"
  )

summary_before_after

# clearner
summary_before_after_thesis <- summary_before_after %>%
  mutate(
    width_kb_IQR = paste0(
      round(median_width_kb, 1), " (",
      round(q25_width_kb, 1), "-",
      round(q75_width_kb, 1), ")"
    ),
    probes_IQR = paste0(
      round(median_probes, 1), " (",
      round(q25_probes, 1), "-",
      round(q75_probes, 1), ")"
    ),
    total_cnv_burden_mb = round(total_cnv_burden_mb, 1),
    mean_cnvs_per_sample = round(mean_cnvs_per_sample, 2),
    median_cnvs_per_sample = round(median_cnvs_per_sample, 2)
  ) %>%
  dplyr::select(
    stage,
    caller,
    loss_gain,
    n_cnvs,
    n_samples,
    median_cnvs_per_sample,
    mean_cnvs_per_sample,
    width_kb_IQR,
    probes_IQR,
    total_cnv_burden_mb
  )

summary_before_after_thesis

library(dplyr)
library(tidyr)

change_summary <- summary_before_after %>%
  dplyr::select(
    stage, caller, loss_gain,
    n_cnvs,
    median_width_kb,
    total_cnv_burden_mb,
    median_probes
  ) %>%
  pivot_wider(
    names_from = stage,
    values_from = c(
      n_cnvs,
      median_width_kb,
      total_cnv_burden_mb,
      median_probes
    )
  ) %>%
  mutate(
    cnv_count_change_pct =
      100 * (`n_cnvs_Gap-merged` - n_cnvs_Raw) / n_cnvs_Raw,
    
    median_size_change_pct =
      100 * (`median_width_kb_Gap-merged` - median_width_kb_Raw) / median_width_kb_Raw,
    
    total_burden_change_pct =
      100 * (`total_cnv_burden_mb_Gap-merged` - total_cnv_burden_mb_Raw) / total_cnv_burden_mb_Raw,
    
    median_probe_change_pct =
      100 * (`median_probes_Gap-merged` - median_probes_Raw) / median_probes_Raw
  ) %>%
  mutate(
    across(where(is.numeric), ~ round(.x, 1))
  )

# add percentage
change_summary_compact <- change_summary %>%
  transmute(Caller = caller,
    Type = loss_gain,
    `CNVs raw` = n_cnvs_Raw,
    `CNVs merged` = `n_cnvs_Gap-merged`,
    `Δ CNVs (%)` = cnv_count_change_pct,
    `Median size raw (kb)` = median_width_kb_Raw,
    `Median size merged (kb)` = `median_width_kb_Gap-merged`,
    `Δ size (%)` = median_size_change_pct,
    `Δ burden (%)` = total_burden_change_pct
  )

# use this table
library(gt)

change_summary_compact %>%
  gt() %>%
  tab_header(
    title = "Change in CNV attributes after gap merging"
  ) %>%
  fmt_number(
    columns = where(is.numeric),
    decimals = 1
  ) %>%
  cols_align(
    align = "center",
    columns = everything()
  )

## plot 
percent_raw <- ggplot(change_summary, aes(x = caller, y = cnv_count_change_pct, fill = loss_gain)) +
  geom_col(position = position_dodge(width = 0.8)) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  labs(x = "CNV caller",
    y = "Change in CNV count after merging (%)",
    fill = "CNV type",
    title = "Percentage change in CNV call count after gap merging"
  ) +  theme_bw()

percent_raw

ggsave(filename = "thesis_out/raw_gap_percent_count.png",
       plot = percent_raw,
       width = 6,  height = 4,
       dpi = 300)

# distribution 2

cnv_all2 <- cnv_all %>%
  filter(width_kb > 0) %>%
  mutate(log10_width_kb = log10(width_kb))


# combined
cnv_all <- cnv_all %>%
  filter(width_kb > 0) %>%
  mutate(log10_width_kb = log10(width_kb))

ggplot(cnv_all, aes(x = log10_width_kb, fill = stage)) +
  geom_histogram(
    position = "identity",
    alpha = 0.45,
    bins = 50
  ) +
  facet_wrap(~ caller, ncol = 2) +
  scale_x_continuous(
    breaks = log10(c(1, 10, 100, 1000, 10000, 100000)),
    labels = c("1", "10", "100", "1,000", "10,000", "100,000")
  ) +
  labs(
    x = "CNV size (kb, log10 scale)",
    y = "Number of CNVs",
    fill = "Stage",
    title = "CNV size distribution before and after gap merging"
  ) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme_bw()


ggplot(cnv_all, aes(x = log10_width_kb, fill = stage)) +
  geom_histogram(
    position = "identity",
    alpha = 0.45,
    bins = 40
  ) +
  facet_grid(loss_gain ~ caller) +
  scale_x_continuous(
    breaks = log10(c(1, 10, 100, 1000, 10000, 100000)),
    labels = c("1", "10", "100", "1,000", "10,000", "100,000")
  ) +
  labs(
    x = "CNV size (kb, log10 scale)",
    y = "Number of CNVs",
    fill = "Stage",
    title = "CNV size distributions before and after gap merging by caller and CNV type"
  ) +  theme_bw()

# best plt

gap_merg_dist <- ggplot(cnv_all, aes(x = log10_width_kb, colour = stage, fill = stage)) +
  geom_histogram(
    position = "identity",
    alpha = 0.25,
    bins = 50) +
  facet_wrap(~ caller, ncol = 2) +
  scale_x_continuous(
    breaks = log10(c(1, 10, 100, 1000, 10000, 100000)),
    labels = c("1kb", "10kb", "100kb", "1Mb", "10Mb", "100Mb")
  ) +
  labs(x = "CNV size (log10 scale)",
    y = "Number of CNVs",
    colour = "Stage", fill = "Stage",
    title = "CNV size histograms before and after gap merging"
  ) + theme_thesis() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_fill_manual(values = stage_cols)

gap_merg_dist

ggsave(filename = "thesis_out/raw_gap_histogram.png",
       plot = gap_merg_dist, width = 8, 
       height = 6, dpi = 300)

########################
### Final thesis ready
########################

theme_thesis <- function() {
  theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 11),
      axis.title = element_text(face = "bold"),
      strip.background = element_rect(fill = "grey95", colour = "black"),
      strip.text = element_text(face = "bold"),
      legend.position = "top",
      legend.title = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.2, colour = "grey90")
    )
}

# hist
p_hist_count <- ggplot(cnv_all, aes(x = log10_width_kb, fill = stage, colour = stage)) +
  geom_histogram(
    position = "identity",
    alpha = 0.35,
    bins = 45,
    linewidth = 0.25
  ) +
  facet_wrap(~ caller, ncol = 2) +
  scale_fill_manual(values = stage_cols) +
  scale_colour_manual(values = stage_cols) +
  scale_x_continuous(
    breaks = log10(c(1, 10, 100, 1000, 10000, 100000)),
    labels = c("1kb", "10kb", "100kb", "1Mb", "10Mb", "100Mb")
  ) +
  labs(
    title = "Distribution of CNV sizes before and after gap merging",
  #  subtitle = "Overlayed histograms shown separately for each caller",
    x = "CNV size (log10 scale)",
    y = "Number of CNVs",
    fill = "Processing stage",
    colour = "Processing stage"
  ) +
  theme_thesis()+ theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

p_hist_count


ggsave(filename = "thesis_out/cnv_size_hist_overlay_by_caller.png",
  plot = p_hist_count,
  width = 10,
  height = 7,
  dpi = 300
)

## more merged insights
library(dplyr)
library(ggplot2)
library(scales)

# manual
cnv_all %>%
  filter(n_merged_segments  > 2) %>%
  group_by(caller) %>%
  summarise(n_merged_cnvs = n())

merged_insight_summary <- cnv_all %>%
  filter(stage == "Gap-merged") %>%
  group_by(caller, loss_gain) %>%
  summarise(
    n_merged_cnvs = n(),
    
    n_single_segment_cnvs = sum(n_merged_segments == 1, na.rm = TRUE),
    n_multi_segment_cnvs = sum(n_merged_segments > 1, na.rm = TRUE),
    
    pct_multi_segment_cnvs = 100 * mean(n_merged_segments > 1, na.rm = TRUE),
    
    median_n_segments = median(n_merged_segments, na.rm = TRUE),
    q25_n_segments = quantile(n_merged_segments, 0.25, na.rm = TRUE),
    q75_n_segments = quantile(n_merged_segments, 0.75, na.rm = TRUE),
    max_n_segments = max(n_merged_segments, na.rm = TRUE),
    
    median_merged_width_kb = median(merged_width_bp / 1000, na.rm = TRUE),
    q25_merged_width_kb = quantile(merged_width_bp / 1000, 0.25, na.rm = TRUE),
    q75_merged_width_kb = quantile(merged_width_bp / 1000, 0.75, na.rm = TRUE),
    
    .groups = "drop"
  )

merged_insight_summary

# out table
merged_insight_table <- merged_insight_summary %>%
  mutate(
    `Merged CNVs` = n_merged_cnvs,
    `Multi-segment CNVs (%)` = round(pct_multi_segment_cnvs, 1),
    `Segments per merged CNV, median (IQR)` = paste0(
      round(median_n_segments, 1), " (",
      round(q25_n_segments, 1), "-",
      round(q75_n_segments, 1), ")"
    ),
    `Max segments` = max_n_segments,
    `Merged size kb, median (IQR)` = paste0(
      round(median_merged_width_kb, 1), " (",
      round(q25_merged_width_kb, 1), "-",
      round(q75_merged_width_kb, 1), ")"
    )
  ) %>%
  dplyr::select(caller,
    loss_gain,
    `Merged CNVs`,
    `Multi-segment CNVs (%)`,
    `Segments per merged CNV, median (IQR)`,
    `Max segments`,
    `Merged size kb, median (IQR)`
  )

merged_insight_table

#3. Plot the distribution of n_merged_segments
p_n_segments <- cnv_all %>%
  filter(stage == "Gap-merged") %>%
  mutate(
    n_merged_segments_capped = pmin(n_merged_segments, 10),
    segment_group = ifelse(
      n_merged_segments_capped == 10,
      "10+",
      as.character(n_merged_segments_capped)
    ),
    segment_group = factor(
      segment_group,
      levels = c(as.character(1:9), "10+")
    )
  ) %>%
  ggplot(aes(x = segment_group, fill = loss_gain)) +
  geom_bar(position = "dodge") +
  facet_wrap(~ caller, ncol = 2) +
  scale_y_continuous(labels = comma) +
  labs(
    x = "Number of raw segments contributing to merged CNV",
    y = "Number of merged CNVs",
    fill = "CNV type",
    title = "Number of raw segments collapsed into each gap-merged CNV"
  ) + scale_fill_manual(values = cnv_cols) + theme_thesis()

p_n_segments

# 4. Plot only CNVs that were actually merged
p_true_merges <- cnv_all %>%
  filter(stage == "Gap-merged", n_merged_segments > 1) %>%
  mutate(
    n_merged_segments_capped = pmin(n_merged_segments, 10),
    segment_group = ifelse(
      n_merged_segments_capped == 10,
      "10+",
      as.character(n_merged_segments_capped)
    ),
    segment_group = factor(
      segment_group,
      levels = c(as.character(2:9), "10+")
    )
  ) %>%
  ggplot(aes(x = segment_group, fill = loss_gain)) +
  geom_bar(position = "dodge") +
  facet_wrap(~ caller, ncol = 2) +
  scale_y_continuous(labels = comma) +
  labs(
    x = "Number of raw segments contributing to merged CNV",
    y = "Number of truly merged CNVs",
    fill = "CNV type",
    title = "Complexity of true gap-merging events"
  ) + scale_fill_manual(values = cnv_cols) + theme_thesis()

p_true_merges

# 5. Compare merged width against number of merged segments

p_segments_vs_width <- cnv_all %>%
  filter(stage == "Gap-merged", n_merged_segments > 1) %>%
  mutate(
    merged_width_kb = merged_width_bp / 1000
  ) %>%
  ggplot(aes(x = n_merged_segments, y = merged_width_kb, col=loss_gain)) +
  geom_point(alpha = 0.3, size = 0.8) +
  facet_grid(loss_gain ~ caller) +
  scale_y_log10(labels = comma) +
  scale_x_continuous(breaks = pretty_breaks()) +
  labs(
    x = "Number of raw segments contributing to merged CNV",
    y = "Merged CNV size (kb, log10 scale)",
    title = "Relationship between merging complexity and merged CNV size"
  ) + scale_color_manual(values = cnv_cols) + theme_thesis() + theme(legend.position = 'none')
# use this
p_segments_vs_width


ggsave(filename = "thesis_out/merged_complex_cnv_size.png",
       plot = p_segments_vs_width,
       width = 8,
       height = 5,
       dpi = 300)

#6. 2. Summarise fragmentation using n_merged_segments

merge_fragmentation_summary <- cnv_all %>%
  filter(stage == "Gap-merged") %>%
  group_by(caller, loss_gain) %>%
  summarise(
    n_final_cnvs = dplyr::n(),
    
    n_unchanged_cnvs = sum(n_merged_segments == 1, na.rm = TRUE),
    n_true_merged_cnvs = sum(n_merged_segments > 1, na.rm = TRUE),
    
    pct_true_merged_cnvs = 100 * mean(n_merged_segments > 1, na.rm = TRUE),
    
    median_segments_per_cnv = median(n_merged_segments, na.rm = TRUE),
    q25_segments_per_cnv = quantile(n_merged_segments, 0.25, na.rm = TRUE),
    q75_segments_per_cnv = quantile(n_merged_segments, 0.75, na.rm = TRUE),
    max_segments_per_cnv = max(n_merged_segments, na.rm = TRUE),
    
    median_width_kb = median(width / 1000, na.rm = TRUE),
    q25_width_kb = quantile(width / 1000, 0.25, na.rm = TRUE),
    q75_width_kb = quantile(width / 1000, 0.75, na.rm = TRUE),
    
    .groups = "drop"
  )

merge_fragmentation_summary

# thesis ready table
merge_fragmentation_table <- merge_fragmentation_summary %>%
  transmute(
    Caller = caller,
    Type = loss_gain,
    `Final CNVs` = n_final_cnvs,
    `True merged CNVs` = n_true_merged_cnvs,
    `True merged (%)` = round(pct_true_merged_cnvs, 1),
    `Segments per CNV, median (IQR)` = paste0(
      round(median_segments_per_cnv, 1), " (",
      round(q25_segments_per_cnv, 1), "-",
      round(q75_segments_per_cnv, 1), ")"
    ),
    `Max segments` = max_segments_per_cnv,
    `Final size kb, median (IQR)` = paste0(
      round(median_width_kb, 1), " (",
      round(q25_width_kb, 1), "-",
      round(q75_width_kb, 1), ")"
    )
  )

merge_fragmentation_table

#4. Compare size of unchanged vs truly merged CNVs
cnv_all2 <- cnv_all %>%
  filter(stage == "Gap-merged") %>%
  mutate(
    merge_status = if_else(
      n_merged_segments > 1,
      "Multi-segment merged",
      "Single-segment unchanged"
    ),
    width_kb = width / 1000
  )

merge_status_size_summary <- cnv_all2 %>%
  group_by(caller, loss_gain, merge_status) %>%
  summarise(
    n_cnvs = n(),
    median_width_kb = median(width_kb, na.rm = TRUE),
    q25_width_kb = quantile(width_kb, 0.25, na.rm = TRUE),
    q75_width_kb = quantile(width_kb, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

merge_status_size_summary

# 5. Plot size distributions by merge status
ggplot(cnv_all2, aes(x = width_kb, fill = merge_status)) +
  geom_histogram(
    bins = 50,
    position = "identity",
    alpha = 0.4
  ) +
  facet_grid(loss_gain ~ caller) +
  scale_x_log10(
    labels = scales::comma,
    breaks = c(1, 10, 100, 1000, 10000, 100000)
  ) +
  labs(
    x = "Final CNV size (kb, log10 scale)",
    y = "Number of CNVs",
    fill = "Merge status",
    title = "Size distribution of unchanged and multi-segment gap-merged CNVs"
  ) +
  theme_bw()

#6. Plot number of merged segments
ggplot(
  cnv_all %>%
    filter(stage == "Gap-merged") %>%
    mutate(
      n_segments_capped = pmin(n_merged_segments, 10),
      segment_group = ifelse(n_segments_capped == 10, "10+", as.character(n_segments_capped)),
      segment_group = factor(segment_group, levels = c(as.character(1:9), "10+"))
    ),
  aes(x = segment_group, fill = loss_gain)
) +
  geom_bar(position = "dodge") +
  facet_wrap(~ caller, ncol = 2) +
  scale_y_continuous(labels = scales::comma) +
  labs(
    x = "Number of raw segments contributing to final CNV",
    y = "Number of final CNVs",
    fill = "CNV type",
    title = "Fragmentation collapsed by gap merging"
  ) +
  theme_bw()

#############
# final table
#############
library(dplyr)
library(tidyr)

gap_merge_summary <- cnv_all %>%
  group_by(stage, caller, loss_gain) %>%
  summarise(
    n_cnvs = n(),
    mean_width_kb = mean(width_kb, na.rm = TRUE),
    median_width_kb = median(width_kb, na.rm = TRUE),
    max_width_kb = max(width_kb, na.rm = TRUE),
    total_width_mb = sum(width_kb, na.rm = TRUE) / 1000,
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = stage,
    values_from = c(
      n_cnvs,
      mean_width_kb,
      median_width_kb,
      max_width_kb,
      total_width_mb
    )
  ) %>%
  mutate(
    cnv_reduction_pct =
      100 * (n_cnvs_Raw - `n_cnvs_Gap-merged`) / n_cnvs_Raw,
    
    mean_size_change_pct =
      100 * (`mean_width_kb_Gap-merged` - mean_width_kb_Raw) / mean_width_kb_Raw,
    
    median_size_change_pct =
      100 * (`median_width_kb_Gap-merged` - median_width_kb_Raw) / median_width_kb_Raw,
    
    total_width_change_pct =
      100 * (`total_width_mb_Gap-merged` - total_width_mb_Raw) / total_width_mb_Raw
  )

gap_merge_table <- gap_merge_summary %>%
  transmute(
    Caller = caller,
    Type = loss_gain,
    
    `Raw CNVs` = n_cnvs_Raw,
    `Merged CNVs` = `n_cnvs_Gap-merged`,
    `CNV reduction (%)` = round(cnv_reduction_pct, 1),
    
    `Mean size raw (kb)` = round(mean_width_kb_Raw, 1),
    `Mean size merged (kb)` = round(`mean_width_kb_Gap-merged`, 1),
    `Mean size change (%)` = round(mean_size_change_pct, 1),
    
    `Median size raw (kb)` = round(median_width_kb_Raw, 1),
    `Median size merged (kb)` = round(`median_width_kb_Gap-merged`, 1),
    
    `Total burden change (%)` = round(total_width_change_pct, 1)
  )

gap_merge_table

library(flextable)

ft <- flextable(gap_merge_table) %>%
  autofit()

ft

# bulk reduction
cnv_all %>%
  count(stage) %>%
  tidyr::pivot_wider(names_from = stage, values_from = n) %>%
  mutate(reduction_pct = round(100 * (Raw - `Gap-merged`) / Raw, 1)
  )

