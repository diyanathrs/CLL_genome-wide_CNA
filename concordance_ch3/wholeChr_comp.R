# Investigate whole CHR alterations detected by algos
library(GenomicRanges)
library(GenomeInfoDb)
library(rtracklayer)
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)

#-----------------------------
# Get hg19 chromosome sizes
#-----------------------------
chrom_info <- getChromInfoFromUCSC("hg19") |>
  as.data.frame() |>
  filter(grepl("^chr([1-9]|1[0-9]|2[0-2]|X|Y)$", chrom)) |>
  transmute(
    chr = chrom,
    chr_start = 1,
    chr_end = size
  )

#-----------------------------
# Get hg19 cytobands
#-----------------------------
cyto <- getTable(ucscTableQuery("hg19", table="cytoBand")) |>
  as.data.frame()

# UCSC tables are 0-based, convert to 1-based GRanges-style coordinates
centromeres <- cyto |>
  filter(gieStain == "acen") |>
  group_by(chrom) |>
  summarise(
    cent_start = min(chromStart) + 1,
    cent_end   = max(chromEnd),
    .groups = "drop"
  ) |>
  rename(chr = chrom)

#-----------------------------
# Build chromosome-arm targets
#-----------------------------
arms_df <- chrom_info |>
  inner_join(centromeres, by = "chr") |>
  transmute(
    chr,
    p_start = chr_start,
    p_end   = cent_start - 1,
    q_start = cent_end + 1,
    q_end   = chr_end
  ) |>
  pivot_longer(
    cols = c(p_start, p_end, q_start, q_end),
    names_to = c("arm", ".value"),
    names_pattern = "([pq])_(start|end)"
  ) |>
  mutate(
    target_id = paste0(chr, arm),
    target_type = "arm"
  ) |>
  dplyr::select(chr, start, end, target_id, target_type)

chrom_df <- chrom_info |>
  transmute(
    chr,
    start = chr_start,
    end = chr_end,
    target_id = chr,
    target_type = "chromosome"
  )

min(chrom_df$end/1e6)

targets_df <- bind_rows(arms_df, chrom_df) |>
  mutate(target_width = end - start + 1)

targets_gr <- GRanges(
  seqnames = targets_df$chr,
  ranges = IRanges(targets_df$start, targets_df$end),
  target_id = targets_df$target_id,
  target_type = targets_df$target_type,
  target_width = targets_df$target_width
)

### load CNV data
source('functions_new.R')
cnv_gr <- load.raw_gap(snp.ct = 5, len = 5, use.merged.cnv = T, only.common.sam = T)

cnv_gr <- keepSeqlevels(cnv_gr, value = paste0("chr", c(1:22)),
  pruning.mode = "coarse")

cnv_gr <- trim(cnv_gr)

# Keep required metadata only
mcols(cnv_gr)$sample_id <- as.character(cnv_gr$sample_id)
mcols(cnv_gr)$caller <- as.character(cnv_gr$caller)
mcols(cnv_gr)$loss_gain <- as.character(cnv_gr$loss_gain)

#3. Reduce calls within each sample/caller/CNV-type
group_key <- paste(
  cnv_gr$sample_id,
  cnv_gr$caller,
  cnv_gr$loss_gain,
  sep = "___"
)

cnv_split <- split(cnv_gr, group_key)
cnv_reduced_list <- GenomicRanges::reduce(cnv_split)

cnv_reduced <- unlist(cnv_reduced_list, use.names = FALSE)

group_info <- data.frame(group_key = rep(names(cnv_reduced_list), 
                                         lengths(cnv_reduced_list))) |>
  separate(group_key,
    into = c("sample_id", "caller", "loss_gain"),
    sep = "___",
    remove = FALSE
  )

mcols(cnv_reduced)$sample_id <- group_info$sample_id
mcols(cnv_reduced)$caller <- group_info$caller
mcols(cnv_reduced)$loss_gain <- group_info$loss_gain

#4. Calculate chromosome/arm coverage per caller
hits <- findOverlaps(cnv_reduced, targets_gr)

overlap_df <- data.frame(cnv_idx = queryHits(hits),
  target_idx = subjectHits(hits),
  overlap_bp = width(pintersect(cnv_reduced[queryHits(hits)], targets_gr[subjectHits(hits)]))
)

event_coverage <- overlap_df |>
  mutate(sample_id = cnv_reduced$sample_id[cnv_idx],
    caller = cnv_reduced$caller[cnv_idx],
    loss_gain = cnv_reduced$loss_gain[cnv_idx],
    target_id = targets_gr$target_id[target_idx],
    target_type = targets_gr$target_type[target_idx],
    target_width = targets_gr$target_width[target_idx]
  ) |>
  group_by(sample_id, caller, loss_gain, target_id, target_type, target_width) |>
  summarise(
    covered_bp = sum(overlap_bp),
    covered_fraction = covered_bp / first(target_width),
    .groups = "drop"
  )

#5. Define large arm/chromosome-level CNVs
arm_threshold <- 0.70
chrom_threshold <- 0.90

large_events <- event_coverage |>
  mutate(
    is_large_event = case_when(
      target_type == "arm" &
        covered_fraction >= arm_threshold ~ TRUE,
      
      target_type == "chromosome" &
        covered_fraction >= chrom_threshold ~ TRUE,
      
      TRUE ~ FALSE
    )
  ) |>
  filter(is_large_event)

#6. Count large chromosome/arm CNVs per caller
large_event_counts <- large_events |>
  group_by(caller, loss_gain, target_type) |>
  summarise(
    n_events = dplyr::n(),
    n_samples = n_distinct(sample_id),
    n_targets = n_distinct(target_id),
    .groups = "drop"
  ) |>
  arrange(target_type, loss_gain, desc(n_events))

large_event_counts

# target
large_event_counts_by_target <- large_events |>
  group_by(caller, loss_gain, target_type, target_id) |>
  summarise(
    n_events = dplyr::n(),
    n_samples = n_distinct(sample_id),
    .groups = "drop"
  ) |>
  arrange(target_type, target_id, loss_gain, caller)

large_event_counts_by_target

# plot by caller
ggplot(large_event_counts,
       aes(x = caller, y = n_events, fill = loss_gain)) +
  geom_col(position = "dodge") +
  facet_wrap(~ target_type, scales = "free_y") +
  theme_bw() +
  labs(
    x = "Caller",
    y = "Number of large chromosome/arm-level CNVs",
    fill = "CNV type",
    title = "Large chromosome and chromosome-arm level CNVs by caller"
  )

# target
ggplot(large_event_counts_by_target,
       aes(x = caller, y = target_id, fill = n_events)) +
  geom_tile() +
  facet_grid(loss_gain ~ target_type, scales = "free_y", space = "free_y") +
  theme_bw() +
  labs(
    x = "Caller",
    y = "Chromosome / arm",
    fill = "N events",
    title = "Distribution of large CNVs across chromosome arms and chromosomes"
  )

#8. Test whether Nexus large events are detected by other callers
nexus_large <- large_events |>
  filter(caller == "Nexus") |>
  dplyr::select(
    sample_id,
    loss_gain,
    target_id,
    target_type,
    nexus_covered_fraction = covered_fraction
  )

other_callers <- setdiff(unique(event_coverage$caller), "Nexus")

nexus_support <- nexus_large |>
  tidyr::crossing(caller = other_callers) |>
  left_join(
    event_coverage |>
      dplyr::select(sample_id, caller, loss_gain, target_id, target_type, covered_fraction),
    by = c("sample_id", "caller", "loss_gain", "target_id", "target_type")
  ) |>
  mutate(
    covered_fraction = replace_na(covered_fraction, 0),
    support_threshold = ifelse(target_type == "arm", arm_threshold, chrom_threshold),
    supported = covered_fraction >= support_threshold
  )

# summarize
nexus_support_summary <- nexus_support |>
  group_by(target_type, loss_gain, caller) |>
  summarise(
    n_nexus_events = dplyr::n(),
    n_supported = sum(supported),
    pct_supported = 100 * n_supported / n_nexus_events,
    median_other_caller_coverage = median(covered_fraction),
    .groups = "drop"
  ) |>
  arrange(target_type, loss_gain, caller)

nexus_support_summary

#9. Count how many of the 3 other callers support each Nexus event
nexus_event_support_count <- nexus_support |>
  group_by(sample_id, loss_gain, target_id, target_type, nexus_covered_fraction) |>
  summarise(
    n_other_callers_supporting = sum(supported),
    supporting_callers = paste(caller[supported], collapse = "; "),
    max_other_caller_coverage = max(covered_fraction),
    mean_other_caller_coverage = mean(covered_fraction),
    .groups = "drop"
  ) |>
  mutate(
    support_class = case_when(
      n_other_callers_supporting == 0 ~ "Nexus only",
      n_other_callers_supporting == 1 ~ "Supported by 1 other caller",
      n_other_callers_supporting == 2 ~ "Supported by 2 other callers",
      n_other_callers_supporting == 3 ~ "Supported by all 3 other callers"
    )
  )

nexus_event_support_count

#summary
nexus_support_class_summary <- nexus_event_support_count |>
  group_by(target_type, loss_gain, support_class) |>
  summarise(
    n_events = dplyr::n(),
    pct_events = 100 * n_events / sum(n_events),
    .groups = "drop"
  )

nexus_support_class_summary

# plot
ggplot(nexus_support_class_summary,
       aes(x = support_class, y = n_events, fill = loss_gain)) +
  geom_col(position = "dodge") +
  facet_wrap(~ target_type, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    x = "Support among non-Nexus callers",
    y = "Number of Nexus large events",
    fill = "CNV type",
    title = "Support for Nexus chromosome/arm-level CNVs by other callers"
  )

#10.identify Nexus only large events
nexus_only_events <- nexus_event_support_count |>
  filter(n_other_callers_supporting == 0) |>
  arrange(target_type, loss_gain, target_id, sample_id)

nexus_only_events

nexus_only_partial_signal <- nexus_support |>
  semi_join(
    nexus_only_events,
    by = c("sample_id", "loss_gain", "target_id", "target_type")
  ) |>
  arrange(sample_id, target_id, loss_gain, caller)

nexus_only_partial_signal

#11. Useful diagnostic plot: Nexus vs other caller coverage
ggplot(nexus_support,
       aes(x = nexus_covered_fraction,
           y = covered_fraction,
           col=loss_gain)) +
  geom_point(alpha = 0.4) +
  geom_hline(yintercept = arm_threshold, linetype = 2, linewidth = 0.3) +
  facet_grid(loss_gain ~ caller) +
  theme_thesis() +
  labs(
    x = "Nexus covered fraction of chromosome/arm",
    y = "Other caller covered fraction of same chromosome/arm",
    title = "Coverage support for Nexus large CNVs by other callers"
  ) + scale_color_manual(values = cnv_cols)

#For chromosome-only events, use the chromosome threshold:
ggplot(nexus_support |> filter(target_type == "chromosome"),
       aes(x = nexus_covered_fraction,
           y = covered_fraction,
           col=loss_gain)) +
  geom_point(alpha = 0.4) +
  geom_hline(yintercept = chrom_threshold, linetype = 2, linewidth = 0.3) +
  facet_grid(loss_gain ~ caller) +
  theme_thesis() +
  labs(
    x = "Nexus chromosome coverage",
    y = "Other caller chromosome coverage",
    title = "Support for Nexus whole-chromosome CNVs"
  ) +  scale_color_manual(values = cnv_cols)

#A strong evidence table would be:
nexus_absence_summary <- nexus_event_support_count |>
  mutate(
    other_signal_class = case_when(
      max_other_caller_coverage == 0 ~ "No overlap in other callers",
      max_other_caller_coverage < 0.20 ~ "<20% covered by other callers",
      max_other_caller_coverage < 0.50 ~ "20-50% covered by other callers",
      max_other_caller_coverage < 0.70 ~ "50-70% covered by other callers",
      TRUE ~ "Subthreshold but substantial coverage"
    )
  ) |>
  group_by(target_type, loss_gain, other_signal_class) |>
  summarise(
    n_events = dplyr::n(),
    .groups = "drop"
  )

nexus_absence_summary

