# load RO output and do nested as well
library(dplyr)
library(tidyr)
library(stringr)

use.merged.cnv <- T

ro_pairs <- readRDS(paste0("RO_pairs_", ifelse(use.merged.cnv,'gapMerged','raw'),".Rds"))

# -------------------------------
# 1. Set thresholds
# -------------------------------

ro_cutoff <- 0.50
nested_cutoff <- 0.80

# -------------------------------
# 2. Add nested-overlap statistics
# -------------------------------
all_edges <- ro_pairs %>%
  mutate( query_uid = as.character(query_uid),
    subject_uid = as.character(subject_uid),
    
    smaller_width = pmin(width_query, width_subject),
    larger_width  = pmax(width_query, width_subject),
    
    prop_smaller_overlapped = overlap_bp / smaller_width,
    prop_larger_overlapped  = overlap_bp / larger_width,
    
    size_ratio = larger_width / smaller_width,
    
    strict_ro_pass = reciprocal_overlap >= ro_cutoff,
    
    nested_pass = prop_smaller_overlapped >= nested_cutoff,
    
    nested_only = nested_pass & !strict_ro_pass,
    
    strict_or_nested = strict_ro_pass | nested_pass,
    
    caller_pair = paste(
      pmin(caller_query, caller_subject),
      pmax(caller_query, caller_subject),
      sep = "_vs_"))

# Quick check
all_edges %>%
  summarise(n_pairs = dplyr::n(),
    n_strict_ro = sum(strict_ro_pass, na.rm = TRUE),
    n_nested = sum(nested_pass, na.rm = TRUE),
    n_nested_only = sum(nested_only, na.rm = TRUE),
    prop_strict_ro = mean(strict_ro_pass, na.rm = TRUE),
    prop_nested = mean(nested_pass, na.rm = TRUE),
    prop_nested_only = mean(nested_only, na.rm = TRUE))

summary(all_edges$reciprocal_overlap)
table(all_edges$strict_ro_pass, useNA = "ifany")

# build OR clusters using RO threshold 
# -------------------------------
# 3. Build RO clusters from strict RO edges
# -------------------------------

strict_edges <- all_edges %>% filter(strict_ro_pass) %>%
  dplyr::select(sample_id, loss_gain, query_uid,
    subject_uid) %>%  distinct()

make_ro_clusters <- function(df) {
  
  graph_df <- df %>%
    dplyr::select(from = query_uid, to = subject_uid)
  
  g <- igraph::graph_from_data_frame(graph_df, directed = FALSE)
  comps <- igraph::components(g)$membership
  
  tibble(event_uid = names(comps),
    ro_component = as.integer(comps)
  )
}

ro_cluster_map <- strict_edges %>%
  group_by(sample_id, loss_gain) %>%
  group_modify(~ make_ro_clusters(.x)) %>%
  ungroup() %>%
  mutate(ro_cluster_id = paste(sample_id, loss_gain, ro_component, sep = "__")
  ) %>%
  dplyr::select(event_uid,
    sample_id,
    loss_gain,
    ro_cluster_id)

# -------------------------------
# 4. Create event table
# -------------------------------
query_events <- all_edges %>%
  transmute(event_uid = query_uid,
    sample_id = sample_id,
    loss_gain = loss_gain,
    caller = caller_query,
    width = width_query)

subject_events <- all_edges %>%
  transmute(
    event_uid = subject_uid,
    sample_id = sample_id,
    loss_gain = loss_gain,
    caller = caller_subject,
    width = width_subject )

event_table <- bind_rows(query_events, subject_events) %>%
  distinct()

ro_event_support <- event_table %>%
  inner_join(ro_cluster_map,
    by = c("event_uid", "sample_id", "loss_gain"))

head(ro_event_support)

# -------------------------------
# 5. Add RO cluster IDs to both sides of all pairwise edges
# -------------------------------

edges_with_clusters <- all_edges %>%
  left_join(ro_cluster_map %>%
      rename(query_uid = event_uid,
        query_cluster = ro_cluster_id),
    by = c("query_uid", "sample_id", "loss_gain")) %>%
  left_join(ro_cluster_map %>%
      rename(subject_uid = event_uid,
        subject_cluster = ro_cluster_id),
    by = c("subject_uid", "sample_id", "loss_gain"))

within_ro_cluster_edges <- edges_with_clusters %>%
  filter(!is.na(query_cluster),
    !is.na(subject_cluster),
    query_cluster == subject_cluster) %>%
  mutate(ro_cluster_id = query_cluster)

head(within_ro_cluster_edges)

# -------------------------------
# 6. Cluster-level nested-overlap summary
# -------------------------------

nested_stats_by_cluster <- within_ro_cluster_edges %>%
  group_by(ro_cluster_id, sample_id, loss_gain) %>%
  summarise(
    n_internal_pairs = dplyr::n(),
    
    n_strict_ro_pairs = sum(strict_ro_pass, na.rm = TRUE),
    n_nested_pairs = sum(nested_pass, na.rm = TRUE),
    n_nested_only_pairs = sum(nested_only, na.rm = TRUE),
    
    has_nested_pair = any(nested_pass, na.rm = TRUE),
    has_nested_only_pair = any(nested_only, na.rm = TRUE),
    
    median_prop_smaller_overlapped =
      median(prop_smaller_overlapped, na.rm = TRUE),
    
    max_prop_smaller_overlapped =
      max(prop_smaller_overlapped, na.rm = TRUE),
    
    median_size_ratio =
      median(size_ratio, na.rm = TRUE),
    
    max_size_ratio =
      max(size_ratio, na.rm = TRUE),
    
    .groups = "drop"
  )

head(nested_stats_by_cluster)

# -------------------------------
# 7. Event and caller support per RO cluster
# -------------------------------

ro_cluster_event_stats <- ro_event_support %>%
  group_by(ro_cluster_id, sample_id, loss_gain) %>%
  summarise(n_segments = dplyr::n(),
    n_callers = n_distinct(caller),
    callers = paste(sort(unique(caller)), collapse = ";"),
    median_segment_width = median(width, na.rm = TRUE),
    max_segment_width = max(width, na.rm = TRUE),
    min_segment_width = min(width, na.rm = TRUE),
    cluster_size_ratio = max_segment_width / min_segment_width,
    .groups = "drop")

head(ro_cluster_event_stats)
nrow(ro_cluster_event_stats)

# -------------------------------
# 8. Final annotated RO cluster table
# -------------------------------

ro_cluster_summary_nested <- ro_cluster_event_stats %>%
  left_join(nested_stats_by_cluster,
    by = c("ro_cluster_id", "sample_id", "loss_gain") ) %>%
  mutate(across(c(n_internal_pairs,
        n_strict_ro_pairs,
        n_nested_pairs, n_nested_only_pairs),
      ~ replace_na(.x, 0)),
    has_nested_pair = replace_na(has_nested_pair, FALSE),
    has_nested_only_pair = replace_na(has_nested_only_pair, FALSE))

head(ro_cluster_summary_nested)

# -------------------------------
# 9. Overall nested-overlap summary across RO clusters
# -------------------------------

overall_nested_cluster_summary <- ro_cluster_summary_nested %>%
  summarise(n_ro_clusters = dplyr::n(),
    
    clusters_with_nested_overlap =
      sum(has_nested_pair, na.rm = TRUE),
    
    clusters_with_nested_only_overlap =
      sum(has_nested_only_pair, na.rm = TRUE),
    
    prop_clusters_with_nested_overlap =
      clusters_with_nested_overlap / n_ro_clusters,
    
    prop_clusters_with_nested_only_overlap =
      clusters_with_nested_only_overlap / n_ro_clusters,
    
    median_nested_pairs_per_cluster =
      median(n_nested_pairs, na.rm = TRUE),
    
    median_nested_only_pairs_per_cluster =
      median(n_nested_only_pairs, na.rm = TRUE),
    
    median_cluster_size_ratio =
      median(cluster_size_ratio, na.rm = TRUE),
    
    max_cluster_size_ratio =
      max(cluster_size_ratio, na.rm = TRUE))

overall_nested_cluster_summary

# -------------------------------
# 10. Caller-pair nested-overlap summary within RO clusters
# -------------------------------

nested_by_caller_pair <- within_ro_cluster_edges %>%
  group_by(caller_pair) %>%
  summarise(n_pairs = dplyr::n(),
    
    n_strict_ro_pairs = sum(strict_ro_pass, na.rm = TRUE),
    n_nested_pairs = sum(nested_pass, na.rm = TRUE),
    n_nested_only_pairs = sum(nested_only, na.rm = TRUE),
    
    prop_strict_ro = n_strict_ro_pairs / n_pairs,
    prop_nested = n_nested_pairs / n_pairs,
    prop_nested_only = n_nested_only_pairs / n_pairs,
    
    median_size_ratio = median(size_ratio, na.rm = TRUE),
    max_size_ratio = max(size_ratio, na.rm = TRUE),
    
    .groups = "drop") %>%
  arrange(desc(prop_nested_only), desc(n_nested_only_pairs))

nested_by_caller_pair

# sum up by gain_loss
nested_by_loss_gain <- ro_cluster_summary_nested %>%
  group_by(loss_gain) %>%
  summarise(n_ro_clusters = dplyr::n(),
    
    clusters_with_nested_overlap =
      sum(has_nested_pair, na.rm = TRUE),
    
    clusters_with_nested_only_overlap =
      sum(has_nested_only_pair, na.rm = TRUE),
    
    prop_clusters_with_nested_overlap =
      clusters_with_nested_overlap / n_ro_clusters,
    
    prop_clusters_with_nested_only_overlap =
      clusters_with_nested_only_overlap / n_ro_clusters,
    
    median_cluster_size_ratio =
      median(cluster_size_ratio, na.rm = TRUE),
    
    .groups = "drop")

nested_by_loss_gain

# plot
library(ggplot2)

ggplot(ro_cluster_summary_nested,
  aes(x = log10(cluster_size_ratio),
    fill = has_nested_only_pair)) +
  geom_histogram(bins = 50, alpha = 0.8, position = "identity") +
  facet_wrap() +
  labs(
    x = "log10(cluster-level size ratio)",
    y = "Number of RO clusters",
    fill = "Contains nested-only pair",
    title = "Nested-only overlap within 0.5 reciprocal-overlap clusters") +
  theme_bw()

# -------------------------------
# 11. Percentage of smaller calls by caller within RO clusters
# -------------------------------

smaller_call_by_cluster <- within_ro_cluster_edges %>%
  mutate(
    smaller_caller = case_when(
      width_query < width_subject ~ caller_query,
      width_subject < width_query ~ caller_subject,
      TRUE ~ "equal_width"
    ),
    
    larger_caller = case_when(
      width_query > width_subject ~ caller_query,
      width_subject > width_query ~ caller_subject,
      TRUE ~ "equal_width"
    ),
    
    smaller_uid = case_when(
      width_query < width_subject ~ query_uid,
      width_subject < width_query ~ subject_uid,
      TRUE ~ NA_character_
    ),
    
    larger_uid = case_when(
      width_query > width_subject ~ query_uid,
      width_subject > width_query ~ subject_uid,
      TRUE ~ NA_character_
    )
  ) %>% filter(smaller_caller != "equal_width")

# summarize 
smaller_call_percent_by_cluster <- smaller_call_by_cluster %>%
  group_by(ro_cluster_id, sample_id, loss_gain, smaller_caller) %>%
  summarise(
    n_smaller_pair_instances = dplyr::n(),
    n_unique_smaller_calls = n_distinct(smaller_uid),
    n_nested_smaller_instances = sum(nested_pass, na.rm = TRUE),
    n_nested_only_smaller_instances = sum(nested_only, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(ro_cluster_id, sample_id, loss_gain) %>%
  mutate(
    total_smaller_pair_instances = sum(n_smaller_pair_instances),
    total_unique_smaller_calls = sum(n_unique_smaller_calls),
    
    percent_smaller_pair_instances =
      100 * n_smaller_pair_instances / total_smaller_pair_instances,
    
    percent_unique_smaller_calls =
      100 * n_unique_smaller_calls / total_unique_smaller_calls,
    
    percent_nested_smaller_instances =
      100 * n_nested_smaller_instances / sum(n_nested_smaller_instances),
    
    percent_nested_only_smaller_instances =
      100 * n_nested_only_smaller_instances / sum(n_nested_only_smaller_instances)
  ) %>%
  ungroup()

smaller_call_percent_by_cluster

# check for each caller
nested_smaller_call_percent_by_cluster <- smaller_call_by_cluster %>%
  filter(nested_pass) %>%
  group_by(ro_cluster_id, sample_id, loss_gain, smaller_caller) %>%
  summarise(
    n_nested_smaller_pair_instances = dplyr::n(),
    n_unique_nested_smaller_calls = n_distinct(smaller_uid),
    .groups = "drop"
  ) %>% group_by(ro_cluster_id, sample_id, loss_gain) %>%
  mutate(
    total_nested_smaller_pair_instances =
      sum(n_nested_smaller_pair_instances),
    
    total_unique_nested_smaller_calls =
      sum(n_unique_nested_smaller_calls),
    
    percent_nested_smaller_pair_instances =
      100 * n_nested_smaller_pair_instances /
      total_nested_smaller_pair_instances,
    
    percent_unique_nested_smaller_calls =
      100 * n_unique_nested_smaller_calls /
      total_unique_nested_smaller_calls
  ) %>%
  ungroup()

## Overall summary across all RO clusters
#To report this in your thesis, you probably also want a global summary:
  
overall_smaller_call_percent_by_caller <- smaller_call_by_cluster %>%
  group_by(smaller_caller) %>%
  summarise(
    n_smaller_pair_instances = dplyr::n(),
    n_unique_smaller_calls = n_distinct(smaller_uid),
    n_nested_smaller_instances = sum(nested_pass, na.rm = TRUE),
    n_nested_only_smaller_instances = sum(nested_only, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    percent_smaller_pair_instances =
      100 * n_smaller_pair_instances / sum(n_smaller_pair_instances),
    
    percent_unique_smaller_calls =
      100 * n_unique_smaller_calls / sum(n_unique_smaller_calls),
    
    percent_nested_smaller_instances =
      100 * n_nested_smaller_instances / sum(n_nested_smaller_instances),
    
    percent_nested_only_smaller_instances =
      100 * n_nested_only_smaller_instances / sum(n_nested_only_smaller_instances)
  ) %>%
  arrange(desc(percent_unique_smaller_calls))

overall_smaller_call_percent_by_caller

#And the nested-only version:
overall_nested_only_smaller_percent_by_caller <- smaller_call_by_cluster %>%
  filter(nested_only) %>%
  group_by(smaller_caller) %>%
  summarise(
    n_nested_only_pair_instances = dplyr::n(),
    n_unique_nested_only_smaller_calls = n_distinct(smaller_uid),
    .groups = "drop"
  ) %>%
  mutate(
    percent_nested_only_pair_instances =
      100 * n_nested_only_pair_instances /
      sum(n_nested_only_pair_instances),
    
    percent_unique_nested_only_smaller_calls =
      100 * n_unique_nested_only_smaller_calls /
      sum(n_unique_nested_only_smaller_calls)
  ) %>%
  arrange(desc(percent_unique_nested_only_smaller_calls))

overall_nested_only_smaller_percent_by_caller

#Add this back to your RO cluster summary
smaller_caller_composition <- nested_smaller_call_percent_by_cluster %>%
  mutate(smaller_caller_percent = paste0(smaller_caller,
      ": ",
      round(percent_unique_nested_smaller_calls, 1),
      "%"
    )) %>%
  group_by(ro_cluster_id, sample_id, loss_gain) %>%
  summarise(
    nested_smaller_caller_composition =
      paste(smaller_caller_percent, collapse = "; "),
    .groups = "drop")

ro_cluster_summary_nested2 <- ro_cluster_summary_nested %>%
  left_join(
    smaller_caller_composition,
    by = c("ro_cluster_id", "sample_id", "loss_gain"))
