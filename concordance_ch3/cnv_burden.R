# generate data and figs for cnv burden
# sample cnv burden by caller (by calls, by length)
# sample cnv burden by chr before and after gap merging

library(dplyr)
library(ggplot2)
library(stringr)
library(forcats)
library(patchwork)
library(scales)


source('functions_new.R')

raw_cnv_gr      <- load.raw_gap(snp.ct = 5, len = 5, use.merged.cnv = F, only.common.sam = F, qc.pass.only = T)
merged_cnv_gr   <- load.raw_gap(snp.ct = 5, len = 5, use.merged.cnv = T, only.common.sam = F, qc.pass.only = T)

# Consistent colours for thesis
caller_cols <- c("PennCNV"   = "#1B9E77",
  "QuantiSNP" = "#D95F02",
  "iPattern"  = "#7570B3",
  "Nexus"     = "#E7298A")

gr_to_cnv_df <- function(gr, dataset_label = "raw") {
  
  df <- as.data.frame(gr)
  
  # Standardise CNV type column
  if ("cnv_type" %in% names(df)) {
    df <- df %>%
      mutate(cnv_state = cnv_type)
  } else if ("loss_gain" %in% names(df)) {
    df <- df %>%
      mutate(cnv_state = loss_gain)
  } else {
    stop("Could not find either 'cnv_type' or 'loss_gain' in GRanges metadata.")
  }
  
  # Basic required metadata checks
  required_cols <- c("seqnames", "start", "end", "width", "sample_id", "caller", "cnv_state")
  missing_cols <- setdiff(required_cols, names(df))
  
  if (length(missing_cols) > 0) {
    stop(
      "Missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  df %>%
    mutate(
      dataset = dataset_label,
      chr = as.character(seqnames),
      width_bp = width,
      width_mb = width_bp / 1e6,
      cnv_state = case_when(
        cnv_state %in% c("Loss", "loss", "DEL", "Deletion", "del") ~ "Loss",
        cnv_state %in% c("Gain", "gain", "DUP", "Amplification", "dup") ~ "Gain",
        TRUE ~ as.character(cnv_state)
      )
    ) %>%
    dplyr::select(
      dataset,
      sample_id,
      caller,
      chr,
      start,
      end,
      width_bp,
      width_mb,
      cnv_state,
      everything()
    )
}

raw_cnv_df <- gr_to_cnv_df(gr = raw_cnv_gr, dataset_label = "Raw CNV calls")
#raw cnv stats
raw_cnv_df %>% group_by(loss_gain, caller) %>% summarise(count=n())

gapmerged_cnv_df <- gr_to_cnv_df(gr = merged_cnv_gr,  dataset_label = "Gap-merged CNV calls")
# merged stats
gapmerged_cnv_df %>% group_by(loss_gain, caller) %>% summarise(count=n())
range(gapmerged_cnv_df$merged_width_bp)

###########################
### for gen coverage
## do reduce
##########################
raw_cnv_gr_cov <- GenomicRanges::reduce(raw_cnv_gr)

cnv_all <- gapmerged_cnv_df

# use only common samples?
# Identify samples present in all 4 callers
common_samples <- cnv_all %>% distinct(sample_id, caller) %>%
  dplyr::count(sample_id, name = "n_callers") %>%
  filter(n_callers == 4) %>% pull(sample_id)

# Filter CNV table to only those samples
cnv_all <- cnv_all %>% filter(sample_id %in% common_samples)

# check
glimpse(cnv_all)

range(cnv_all$merged_width_bp, na.rm = T)
table(cnv_all$caller)
table(cnv_all$cnv_state)

cnv_all <- cnv_all %>%
  mutate(chr = case_when(
      chr %in% as.character(1:22) ~ paste0("chr", chr),
      chr == "X" ~ "chrX",
      chr == "Y" ~ "chrY",
      TRUE ~ chr),
    chr = factor(chr, levels = c(paste0("chr", 1:22), "chrX", "chrY")))

make_diverging_burden <- function(df, value_col) {
  df %>%    mutate(
      cnv_state = case_when(
        cnv_state %in% c("Loss", "loss", "DEL", "Deletion", "del") ~ "Loss",
        cnv_state %in% c("Gain", "gain", "DUP", "Amplification", "dup") ~ "Gain",
        TRUE ~ as.character(cnv_state)
      ),
      burden_signed = case_when(
        cnv_state == "Loss" ~ -abs(.data[[value_col]]),
        cnv_state == "Gain" ~  abs(.data[[value_col]]),
        TRUE ~ .data[[value_col]]
      )
    )
}


chr_burden_n <- cnv_all %>%
  group_by(dataset, chr, cnv_state) %>%
  summarise(
    n_cnvs = n(),
    n_samples = n_distinct(sample_id),
    n_callers = n_distinct(caller),
    .groups = "drop")

# plot chr-wise
p_chr_n <- ggplot(chr_burden_n, aes(x = chr, y = n_cnvs, fill = cnv_state)) +
  geom_col(position = "stack") +
  #facet_wrap(~ dataset, ncol = 1, scales = "free_y") +
  scale_y_continuous(labels = comma) +
  labs(
    title = "Chromosome-wise CNV burden by number of CNV calls",
    x = "Chromosome",
    y = "Number of CNV calls",
    fill = "CNV type"
  ) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  ) + scale_fill_manual(values = cnv_cols)

p_chr_n

chr_burden_n_div <- chr_burden_n %>%
  make_diverging_burden(value_col = "n_cnvs")

p_chr_n_div <- ggplot(chr_burden_n_div, aes(x = chr, y = burden_signed, fill = cnv_state)) +
  geom_col() +
  geom_hline(yintercept = 0, linewidth = 0.4) +
#  facet_wrap(~ dataset, ncol = 1, scales = "free_y") +
  scale_y_continuous(
    labels = function(x) comma(abs(x))
  ) +
  labs(
    title = "A. Chromosome-wise CNV burden by number of calls",
    # subtitle = "Gains are shown above zero; losses are shown below zero",
    x = "Chromosome",
    y = "Number of CNV calls",
    fill = "CNV type"
  ) + theme_thesis() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank())+ scale_fill_manual(values = cnv_cols)


p_chr_n_div

###################
# by CNV coverage
##################
chr_burden_coverage <- cnv_all %>%
  group_by(dataset, chr, cnv_state) %>%
  summarise(
    total_cnv_bp = sum(width_bp, na.rm = TRUE),
    total_cnv_mb = total_cnv_bp / 1e6,
    mean_cnv_mb = mean(width_mb, na.rm = TRUE),
    median_cnv_mb = median(width_mb, na.rm = TRUE),
    n_cnvs = n(),
    .groups = "drop"
  )

p_chr_cov <- ggplot(chr_burden_coverage, aes(x = chr, y = total_cnv_mb, fill = cnv_state)) +
  geom_col(position = "stack") +
  facet_wrap(~ dataset, ncol = 1, scales = "free_y") +
  scale_y_continuous(labels = comma) +
  labs(
    title = "Chromosome-wise CNV burden by total CNV coverage",
    x = "Chromosome",
    y = "Total CNV coverage (Mb)",
    fill = "CNV type"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position='none',
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  )+ scale_fill_manual(values = cnv_cols)

p_chr_cov

chr_burden_coverage_div <- chr_burden_coverage %>%
  make_diverging_burden(value_col = "total_cnv_mb")

p_chr_cov_div <- ggplot(chr_burden_coverage_div, aes(x = chr, y = burden_signed, fill = cnv_state)) +
  geom_col() +
  geom_hline(yintercept = 0, linewidth = 0.4) +
 # facet_wrap(~ dataset, ncol = 1, scales = "free_y") +
  scale_y_continuous(
    labels = function(x) comma(abs(x))
  ) +
  labs(
    title = "B. Chromosome-wise CNV burden by total coverage",
   # subtitle = "Gains are shown above zero; losses are shown below zero",
    x = "Chromosome",
    y = "Total CNV coverage (Mb)",
    fill = "CNV type"
  ) + theme_thesis() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), 
    panel.grid.minor = element_blank(), legend.position = 'none',
  ) + scale_fill_manual(values = cnv_cols)

p_chr_cov_div

# together
p_chr_n_div / p_chr_cov_div 

ggsave(filename = "thesis_out/chromosome_wise_cnv_burden2.png",
  plot = p_chr_n_div / p_chr_cov_div,
  width = 9,
  height = 8,
  dpi = 300)

##################################
##Part B: Caller-wise CNV burden
##################################

caller_burden_n <- cnv_all %>%
  group_by(dataset, caller, cnv_state) %>%
  summarise(
    n_cnvs = n(),
    n_samples = n_distinct(sample_id),
    n_chromosomes = n_distinct(chr),
    .groups = "drop"
  ) %>% mutate(burden_signed = case_when(
    cnv_state == "Loss" ~ -abs(n_cnvs ),
    cnv_state == "Gain" ~  abs(n_cnvs )
  ))

# new div plots
p_caller_n_div <- ggplot(caller_burden_n, aes(x = caller, y = burden_signed, fill = cnv_state)) +
  geom_col() +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  #facet_wrap(~ dataset, ncol = 1, scales = "free_y") +
  scale_y_continuous(
    labels = function(x) comma(abs(x))
  ) +
  labs(
    title = "Chromosome-wise CNV burden by total coverage",
    # subtitle = "Gains are shown above zero; losses are shown below zero",
    x = "CNV caller",
    y = "Total CNV coverage (Mb)",
    fill = "CNV type"
  ) +
  theme_thesis() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), 
    panel.grid.minor = element_blank()
  ) + scale_fill_manual(values = cnv_cols)

p_caller_n_div


# by cov
caller_burden_coverage <- cnv_all %>%
  group_by(dataset, caller, cnv_state) %>%
  summarise(
    total_cnv_bp = sum(width_bp, na.rm = TRUE),
    total_cnv_mb = total_cnv_bp / 1e6,
    mean_cnv_mb = mean(width_mb, na.rm = TRUE),
    median_cnv_mb = median(width_mb, na.rm = TRUE),
    n_cnvs = n(),
    .groups = "drop"
  ) %>% mutate(burden_signed = case_when(
    cnv_state == "Loss" ~ -abs(total_cnv_mb),
    cnv_state == "Gain" ~  abs(total_cnv_mb)
  ))

# new div plots
p_caller_cov_div <- ggplot(caller_burden_coverage, aes(x = caller, y = burden_signed, fill = cnv_state)) +
  geom_col() +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  #facet_wrap(~ dataset, ncol = 1, scales = "free_y") +
  scale_y_continuous(
    labels = function(x) comma(abs(x))
  ) +
  labs(
    title = "Chromosome-wise CNV burden by total coverage",
    # subtitle = "Gains are shown above zero; losses are shown below zero",
    x = "CNV caller",
    y = "Total CNV coverage (Mb)",
    fill = "CNV type"
  ) +
  theme_thesis() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), 
    panel.grid.minor = element_blank()
  ) + scale_fill_manual(values = cnv_cols)

p_caller_cov_div

p_caller_n_div / p_caller_cov_div

ggsave(filename = "thesis_out/caller_wise_cnv_burden_raw_vs_gapmerged.png",
  plot = p_caller_n / p_caller_cov,
  width = 10, height = 8, dpi = 300)


p_caller_n_div

########################################
##Part C: Stratify by caller and chromosome together
########################################
chr_caller_burden_n <- cnv_all %>%
  group_by(dataset, caller, chr, cnv_state) %>%
  summarise(n_cnvs = n(),
    n_samples = n_distinct(sample_id),
    .groups = "drop")

# by cov
chr_caller_burden_n_div <- cnv_all %>%
  group_by(dataset, caller,chr, cnv_state) %>%
  summarise(n_cnvs = n(),
            n_samples = n_distinct(sample_id),
            .groups = "drop")

  summarise(total_cnv_bp = sum(width_bp, na.rm = TRUE),
    total_cnv_mb = total_cnv_bp / 1e6,
    mean_cnv_mb = mean(width_mb, na.rm = TRUE),
    median_cnv_mb = median(width_mb, na.rm = TRUE),
    n_cnvs = n(),
    .groups = "drop"
  )


p_chr_caller_n <- ggplot(chr_caller_burden_n, aes(x = chr, y = n_cnvs, fill = cnv_state)) +
  geom_col(position = "stack") +
  facet_grid(dataset ~ caller, scales = "free_y") +
  scale_y_continuous(labels = comma) +
  labs(title = "Chromosome-wise CNV count stratified by caller",
    x = "Chromosome",
    y = "Number of CNV calls",
    fill = "CNV type") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank())

p_chr_caller_n

# new div plots
p_chr_caller_n_div <- ggplot(chr_caller_burden_n_div, aes(x = chr, y = n_cnvs, fill = cnv_state)) +
  geom_col() +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  #facet_wrap(~ dataset, ncol = 1, scales = "free_y") +
  scale_y_continuous(
    labels = function(x) comma(abs(x))
  ) +
  labs(
    title = "Chromosoomal coverage",
    # subtitle = "Gains are shown above zero; losses are shown below zero",
    x = "CNV caller",
    y = "Total CNV coverage (Mb)",
    fill = "CNV type"
  ) +
  theme_thesis() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), 
    panel.grid.minor = element_blank()
  ) + scale_fill_manual(values = cnv_cols)

p_chr_caller_n_div

chr_caller_burden_cov <- cnv_all %>%
  group_by(dataset, caller, chr, cnv_state) %>%
  summarise(total_cnv_bp = sum(width_bp, na.rm = TRUE),
    total_cnv_mb = total_cnv_bp / 1e6,
    n_cnvs = n(), .groups = "drop")

p_chr_caller_cov <- ggplot(chr_caller_burden_cov, aes(x = chr, y = total_cnv_mb, fill = cnv_state)) +
  geom_col(position = "stack") +
  facet_grid(dataset ~ caller, scales = "free_y") +
  scale_y_continuous(labels = comma) +
  labs(title = "Chromosome-wise CNV coverage stratified by caller",
    x = "Chromosome",
    y = "Total CNV coverage (Mb)",
    fill = "CNV type") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank())

p_chr_caller_cov

ggsave(filename = "thesis_out/chromosome_by_caller_cnv_count.png",
  plot = p_chr_caller_n,
  width = 14,  height = 8,
  dpi = 300)

ggsave(
  filename = "thesis_out/chromosome_by_caller_cnv_coverage.png",
  plot = p_chr_caller_cov,
  width = 14,
  height = 8,
  dpi = 300
)

##################################
##Part E: Per-sample CNV burden
#################################
cnv_all <- bind_rows(raw_cnv_df, gapmerged_cnv_df)

sample_caller_burden <- cnv_all %>%
  group_by(dataset, sample_id, caller, cnv_state) %>%
  summarise(n_cnvs = n(),
    total_cnv_mb = sum(width_mb, na.rm = TRUE),
    mean_cnv_mb = mean(width_mb, na.rm = TRUE),
    median_cnv_mb = median(width_mb, na.rm = TRUE),
    .groups = "drop"
  )

# boxplot
p_sample_n <- ggplot(sample_caller_burden, aes(x = caller, y = n_cnvs, color = cnv_state)) +
  geom_boxplot(outlier.alpha = 0.25) +
  facet_wrap(~ dataset, ncol = 2) +
  scale_y_continuous(labels = comma) +
  labs(
    title = "Per-sample CNV burden by caller",
    x = "CNV caller",
    y = "Number of CNVs per sample",
    col = "CNV type"
  ) +
  theme_thesis() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  ) + scale_color_manual(values = cnv_cols)

p_sample_n

# coverage
p_sample_cov <- ggplot(sample_caller_burden, aes(x = caller, y = total_cnv_mb, col = cnv_state)) +
  geom_boxplot(outlier.alpha = 0.25) +
  facet_wrap(~ dataset, ncol = 1, scales = "free_y") +
  scale_y_continuous(labels = comma) +
  labs(
    title = "Per-sample CNV coverage burden by caller",
    x = "CNV caller",
    y = "Total CNV coverage per sample (Mb)",
    color = "CNV type"
  ) +
  theme_bw(base_size = 12) 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank() + ylim(0 , 199)
  )

p_sample_cov

p_sample_n | p_sample_cov

##############
# for focal vs broad
##################

cnv_df <- cnv_all %>%
  mutate(width_bp = end - start,
         width_mb = width_bp / 1e6,
         cnv_class = ifelse(width_mb <= 1, "Focal", "Broad")
  )
cnv_df$dataset <- factor(cnv_df$dataset, levels = c('Raw CNV calls', 'Gap-merged CNV calls')) 

sample_burden <- cnv_df %>%
  group_by(sample_id, caller, loss_gain, cnv_class, dataset ) %>%
  summarise(
    burden_mb = sum(width_mb),
    .groups = "drop"
  )

sample_burden <- cnv_df %>%
  group_by(sample_id, caller, loss_gain, cnv_class, dataset ) %>%
  summarise(
  burden = n(),
  .groups = "drop"
)

ggplot(sample_burden,
       aes(x = caller,
           y = burden,
           col = loss_gain)) +
  geom_boxplot(outlier.alpha = 0.2) +
  facet_grid(cnv_class~dataset, scales = "free_y") +
    labs(
    title = "Per-sample CNV Burden by CNV Class",
    x = "",
    y = "CNV burden per sample", color='CNV type',
  ) +
  theme_thesis() + scale_color_manual(values = cnv_cols)
