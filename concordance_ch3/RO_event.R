# 4. Event-level: CNV length vs caller support
library(igraph)

get_ro_clusters_with_coords <- function(df) {
  
  if (nrow(df) == 0) return(NULL)
  
  edges <- df %>%
    transmute(
      from = paste0(caller_query, "_", query),
      to   = paste0(caller_subject, "_", subject)
    )
  
  g <- graph_from_data_frame(edges, directed = FALSE)
  comp <- components(g)$membership
  
  cluster_df <- data.frame(
    node = names(comp),
    cluster_id = as.integer(comp)
  ) %>%
    mutate(
      caller = sub("_[^_]+$", "", node),
      raw_index = as.integer(sub("^.*_", "", node))
    )
  
  cluster_df %>%
    left_join(
      df %>%
        dplyr::select(query, subject, width_query, width_subject) %>%
        pivot_longer(
          cols = c(query, subject),
          names_to = "role",
          values_to = "raw_index"
        ) %>%
        mutate(width_bp = ifelse(role == "query", width_query, width_subject)) %>%
        distinct(raw_index, width_bp),
      by = "raw_index"
    ) %>%
    group_by(cluster_id) %>%
    summarise(
      n_callers = n_distinct(caller),
      callers = paste(sort(unique(caller)), collapse = ";"),
      n_cnvs = dplyr::n(),
      median_width_kb = median(width_bp / 1000, na.rm = TRUE),
      mean_width_kb = mean(width_bp / 1000, na.rm = TRUE),
      min_width_kb = min(width_bp / 1000, na.rm = TRUE),
      max_width_kb = max(width_bp / 1000, na.rm = TRUE),
      .groups = "drop"
    )
}

ro_threshold <- 0.7

ro_events_len <- ro_pairs_len %>%
  filter(reciprocal_overlap >= ro_threshold) %>%
  group_by(sample_id, loss_gain) %>%
  group_modify(~ get_ro_clusters_with_coords(.x)) %>%
  ungroup()


# summarise
ro_events_len <- ro_events_len %>%
  mutate(
    event_size_bin = case_when(
      median_width_kb < 10 ~ "<10 kb",
      median_width_kb >= 10 & median_width_kb < 50 ~ "10-50 kb",
      median_width_kb >= 50 & median_width_kb < 100 ~ "50-100 kb",
      median_width_kb >= 100 & median_width_kb < 500 ~ "100-500 kb",
      median_width_kb >= 500 & median_width_kb < 1000 ~ "500 kb-1 Mb",
      median_width_kb >= 1000 & median_width_kb < 10000 ~ "1 Mb - 10 Mb",
      median_width_kb >= 10000 & median_width_kb < 100000 ~ "10 Mb - 100 Mb",
      median_width_kb >= 100000 ~ ">100 Mb",
      TRUE ~ NA_character_
    ),
    event_size_bin = factor(
      event_size_bin,
      levels = c("<10 kb", "10-50 kb", "50-100 kb",
                 "100-500 kb", "500 kb-1 Mb", "1 Mb - 10 Mb", "10 Mb - 100 Mb", ">100 Mb")
    )
  )

ro_events_len <- ro_events_len %>%
  mutate(
    event_size_bin = case_when(
      max_width_kb < 10 ~ "<10 kb",
      max_width_kb >= 10 & max_width_kb < 50 ~ "10-50 kb",
      max_width_kb >= 50 & max_width_kb < 100 ~ "50-100 kb",
      max_width_kb >= 100 & max_width_kb < 500 ~ "100-500 kb",
      max_width_kb >= 500 & max_width_kb < 1000 ~ "500 kb-1 Mb",
      max_width_kb >= 1000 & max_width_kb < 10000 ~ "1 Mb - 10 Mb",
      max_width_kb >= 10000 & max_width_kb < 100000 ~ "10 Mb - 100 Mb",
      max_width_kb >= 100000 ~ ">100 Mb",
      TRUE ~ NA_character_
    ),
    event_size_bin = factor(
      event_size_bin,
      levels = c("<10 kb", "10-50 kb", "50-100 kb",
                 "100-500 kb", "500 kb-1 Mb", "1 Mb - 10 Mb", "10 Mb - 100 Mb", ">100 Mb")
    )
  )



support_by_size <- ro_events_len %>% group_by(loss_gain) %>%
  count(event_size_bin, n_callers, name = "n_events") %>%
  group_by(event_size_bin, loss_gain) %>%
  mutate(
    fraction = n_events / sum(n_events)
  ) %>%
  ungroup()

support_by_size

# plot
ro.event.len <- ggplot(support_by_size,
       aes(y = event_size_bin, x = n_events, fill = factor(n_callers))) +
  facet_wrap(~loss_gain) +
  #geom_bar(stat='identity')+
  #scale_y_log10(labels = comma)
  geom_col(position = "fill") +
  theme_thesis() + scale_fill_brewer(palette = 'Set2')+
  labs(
    y = "Event size bin",
    x = "Fraction of concordant events",
    fill = "Supporting callers",
    title = "Caller support by CNV event size"
  ) 
 # theme(axis.text.x = element_text(angle = 45, hjust = 1))

ro.event.len
# save 
ggsave(filename = "thesis_out/RO_event_length.png",
       plot = ro.event.len,
       width = 8, height = 5, dpi = 300)
