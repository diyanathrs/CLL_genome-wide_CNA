############################################################
## 0. Packages
############################################################

library(tidyverse)
library(GenomicRanges)
library(GenomeInfoDb)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(patchwork)


source('functions_new.R')
cnv_type_cols <- c("Loss" = "#D55E00","Gain" = "#0072B2")
caller_cols <- c("PennCNV"   = "#1B9E77", "QuantiSNP" = "#D95F02",
                 "iPattern"  = "#7570B3", "Nexus"     = "#E7298A")
# load CNV data

cnv_gr <- load.raw_gap(snp.ct = 5, len = 5, use.merged.cnv = T, only.common.sam = T)

############################################################
## 3. Get hg19 gene coordinates
############################################################
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
gene_gr <- genes(txdb)

# Keep standard chromosomes only
gene_gr <- keepSeqlevels(gene_gr, paste0("chr", c(1:22, "X", "Y")), pruning.mode = "coarse")

# Add Entrez IDs
gene_gr$entrez_id <- names(gene_gr)

# Map Entrez IDs to gene symbols
gene_symbols <- mapIds(org.Hs.eg.db,
  keys = gene_gr$entrez_id,
  column = "SYMBOL",
  keytype = "ENTREZID",
  multiVals = "first")

gene_gr$gene_symbol <- unname(gene_symbols)

# If symbol is missing, use Entrez ID as fallback
gene_gr$gene_name <- ifelse(
  is.na(gene_gr$gene_symbol),
  paste0("ENTREZ_", gene_gr$entrez_id),
  gene_gr$gene_symbol
)

gene_gr$gene_width_bp <- width(gene_gr)

gene_gr

############################################################
## 4. Find CNV-gene overlaps
############################################################

hits <- findOverlaps(cnv_gr, gene_gr, ignore.strand = TRUE)

cnv_gene_overlap <- tibble(cnv_index  = queryHits(hits),
  gene_index = subjectHits(hits)) %>%
  mutate(cnv_uid   = mcols(cnv_gr)$cnv_uid[cnv_index],
    sample_id = mcols(cnv_gr)$sample_id[cnv_index],
    caller    = mcols(cnv_gr)$caller[cnv_index],
    cnv_type  = mcols(cnv_gr)$loss_gain[cnv_index],
    cnv_chr   = as.character(seqnames(cnv_gr))[cnv_index],
    cnv_start = start(cnv_gr)[cnv_index],
    cnv_end   = end(cnv_gr)[cnv_index],
    cnv_width_bp = width(cnv_gr)[cnv_index],
    
    gene_chr  = as.character(seqnames(gene_gr))[gene_index],
    gene_start = start(gene_gr)[gene_index],
    gene_end   = end(gene_gr)[gene_index],
    gene_width_bp = width(gene_gr)[gene_index],
    entrez_id = mcols(gene_gr)$entrez_id[gene_index],
    gene_symbol = mcols(gene_gr)$gene_symbol[gene_index],
    gene_name = mcols(gene_gr)$gene_name[gene_index])

# Add exact overlap width
overlap_ranges <- pintersect(cnv_gr[queryHits(hits)], gene_gr[subjectHits(hits)])

cnv_gene_overlap <- cnv_gene_overlap %>%
  mutate(overlap_bp = width(overlap_ranges),
    cnv_fraction_overlapping_gene = overlap_bp / cnv_width_bp,
    gene_fraction_covered_by_cnv = overlap_bp / gene_width_bp)

head(cnv_gene_overlap)

saveRDS(cnv_gene_overlap, "cnv_gene_overlap_hg19.rds")

############################################################
## 5. Gene-level burden summary
############################################################
cnv_gene_overlap <- readRDS('cnv_gene_overlap_hg19.rds')

n_total_samples <- n_distinct(cnv_gr$sample_id)

gene_sample_caller <- cnv_gene_overlap %>%
  group_by(gene_name,
    gene_symbol,
    entrez_id,
    gene_chr,
    gene_start,
    gene_end,
    gene_width_bp,
    sample_id,
    caller,
    cnv_type
  ) %>% summarise(n_cnv_calls = dplyr::n(),
    total_overlap_bp = sum(overlap_bp, na.rm = TRUE),
    max_gene_fraction_covered = max(gene_fraction_covered_by_cnv, na.rm = TRUE),
    .groups = "drop") %>%
  mutate(# Prevent fragmented or overlapping calls from producing impossible fractions >1
    total_gene_fraction_covered = pmin(total_overlap_bp / gene_width_bp, 1)
  )


# i want only >50% of the genes covered
gene_sample_caller <- gene_sample_caller %>% filter(max_gene_fraction_covered > 0.5)

gene_sample_summary <- gene_sample_caller %>%
  group_by(gene_name,
    gene_symbol,
    entrez_id,
    gene_chr,
    gene_start,
    gene_end,
    gene_width_bp,
    sample_id,
    cnv_type
  ) %>%
  summarise(n_callers_overlapping_gene = n_distinct(caller),
    callers_overlapping_gene = paste(sort(unique(caller)), collapse = ";"),
    n_total_cnv_calls = sum(n_cnv_calls),
    max_gene_fraction_covered = max(max_gene_fraction_covered, na.rm = TRUE),
    total_gene_fraction_covered = max(total_gene_fraction_covered, na.rm = TRUE),
    .groups = "drop"
  )

saveRDS(gene_sample_summary, 'gene_sample_summary.rds')
gene_sample_summary <- readRDS('gene_sample_summary.rds')

gene_burden_summary <- gene_sample_summary %>%
  group_by(gene_name,
    gene_symbol,
    entrez_id,
    gene_chr,
    gene_start,
    gene_end,
    gene_width_bp,
    cnv_type
  ) %>%
  summarise(n_samples_affected = n_distinct(sample_id),
    prop_samples_affected = n_samples_affected / n_total_samples,
    
    mean_callers_per_affected_sample = mean(n_callers_overlapping_gene, na.rm = TRUE),
    median_callers_per_affected_sample = median(n_callers_overlapping_gene, na.rm = TRUE),
    max_callers_in_any_sample = max(n_callers_overlapping_gene, na.rm = TRUE),
    
    n_total_cnv_calls = sum(n_total_cnv_calls),
    mean_gene_fraction_covered = mean(total_gene_fraction_covered, na.rm = TRUE),
    median_gene_fraction_covered = median(total_gene_fraction_covered, na.rm = TRUE),
    
    .groups = "drop") %>%
  arrange(desc(n_samples_affected), desc(mean_callers_per_affected_sample))

head(gene_burden_summary, 20)
dim(gene_burden_summary)

############################################################
## 6. Caller-specific support per gene
############################################################

gene_caller_summary <- gene_sample_caller %>%
  group_by(
    gene_name,
    gene_symbol,
    entrez_id,
    gene_chr,
    gene_start,
    gene_end,
    cnv_type,
    caller
  ) %>%
  summarise(
    n_samples_by_caller = n_distinct(sample_id),
    n_cnv_calls_by_caller = sum(n_cnv_calls),
    mean_gene_fraction_covered_by_caller = mean(total_gene_fraction_covered, na.rm = TRUE),
    .groups = "drop"
  )

gene_caller_wide <- gene_caller_summary %>%
  dplyr::select(
    gene_name,
    gene_symbol,
    entrez_id,
    gene_chr,
    gene_start,
    gene_end,
    cnv_type,
    caller,
    n_samples_by_caller
  ) %>%
  pivot_wider(
    names_from = caller,
    values_from = n_samples_by_caller,
    values_fill = 0,
    names_prefix = "n_samples_"
  )

gene_burden_with_callers <- gene_burden_summary %>%
  left_join(
    gene_caller_wide,
    by = c(
      "gene_name",
      "gene_symbol",
      "entrez_id",
      "gene_chr",
      "gene_start",
      "gene_end",
      "cnv_type"
    )
  )

head(gene_burden_with_callers, 20)

write.csv(gene_burden_with_callers,
  "gene_cnv_burden_with_caller_support_hg19.csv",
  row.names = FALSE)

gene_burden_with_callers <- read_csv('gene_cnv_burden_with_caller_support_hg19.csv')

############################################################
## 7. High-confidence gene-level CNV calls
############################################################

gene_sample_highconf <- gene_sample_summary %>%
  mutate(
    highconf_2caller = n_callers_overlapping_gene >= 2,
    highconf_3caller = n_callers_overlapping_gene >= 3,
    highconf_4caller = n_callers_overlapping_gene >= 4
  )

highconf_gene_summary <- gene_sample_highconf %>%
  group_by(
    gene_name,
    gene_symbol,
    entrez_id,
    gene_chr,
    gene_start,
    gene_end,
    gene_width_bp,
    cnv_type
  ) %>%
  summarise(
    n_samples_any_caller = n_distinct(sample_id),
    n_samples_2plus_callers = n_distinct(sample_id[highconf_2caller]),
    n_samples_3plus_callers = n_distinct(sample_id[highconf_3caller]),
    n_samples_4_callers = n_distinct(sample_id[highconf_4caller]),
    
    prop_samples_any_caller = n_samples_any_caller / n_total_samples,
    prop_samples_2plus_callers = n_samples_2plus_callers / n_total_samples,
    prop_samples_3plus_callers = n_samples_3plus_callers / n_total_samples,
    prop_samples_4_callers = n_samples_4_callers / n_total_samples,
    
    .groups = "drop"
  ) %>%
  arrange(desc(n_samples_2plus_callers), desc(n_samples_any_caller))

head(highconf_gene_summary, 20)

write.csv(highconf_gene_summary,
  "highconfidence_gene_cnv_summary_hg19.csv",
  row.names = FALSE)

highconf_gene_summary <- read_csv('highconfidence_gene_cnv_summary_hg19.csv')

############################################################
## 8. Plot top genes by affected samples
############################################################
top_n_genes <- 50

plot_top_genes <- gene_burden_summary %>%
  group_by(gene_name, cnv_type) %>%
  summarise(
    n_samples_affected = max(n_samples_affected),
    prop_samples_affected = max(prop_samples_affected),
    mean_callers_per_affected_sample = max(mean_callers_per_affected_sample),
    .groups = "drop"
  ) %>%
  arrange(desc(n_samples_affected)) %>%
  slice_head(n = top_n_genes) %>%
  mutate(
    gene_name = fct_reorder(gene_name, n_samples_affected)
  )

p1 <- ggplot(plot_top_genes,
  aes(x = gene_name,
    y = n_samples_affected,
    fill = cnv_type
  )
) +
  geom_col(position = "dodge") +
  coord_flip() +
  labs(x = "Gene",
    y = "Number of affected samples",
    fill = "CNV type",
    title = "Top genes affected by CNVs",
    subtitle = "Ranked by number of samples with overlapping CNVs"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  ) + scale_fill_manual(values = cnv_cols)

p1

############################################################
## Plot 1B. Top N genes separately for losses and gains,
## stacked by caller
############################################################

top_n_genes <- 50

# Select top N genes separately within each CNV type
# For example: top 30 losses and top 30 gains
top_genes_by_unique_samples <- gene_burden_summary %>%
  group_by(cnv_type, gene_name) %>%
  summarise(
    unique_affected_samples = max(n_samples_affected),
    .groups = "drop"
  ) %>%
  group_by(cnv_type) %>%
  arrange(desc(unique_affected_samples), .by_group = TRUE) %>%
  slice_head(n = top_n_genes) %>%
  ungroup()

# Get caller-specific sample counts for those top genes
plot_top_genes_stacked_by_caller <- gene_caller_summary %>%
  semi_join(
    top_genes_by_unique_samples,
    by = c("gene_name", "cnv_type")
  ) %>%
  left_join(
    top_genes_by_unique_samples,
    by = c("gene_name", "cnv_type")
  ) %>%
  group_by(cnv_type) %>%
  mutate(
    gene_name = forcats::fct_reorder(gene_name, unique_affected_samples)
  ) %>%
  ungroup()

p1_stacked_by_caller <- ggplot(
  plot_top_genes_stacked_by_caller,
  aes(
    x = gene_name,
    y = n_samples_by_caller,
    fill = caller
  )
) +
  geom_col(position = "stack") +
  coord_flip() +
  facet_wrap(~ cnv_type, scales = "free_y") +
  labs(
    x = "Gene",
    y = "Affected samples counted per caller",
    fill = "CNV caller",
    title = paste0("Top ", top_n_genes, " genes affected by CNVs"),
    subtitle = "Top genes selected separately for losses and gains; stacked bars show caller-specific affected sample counts"
  ) +
  theme_thesis() 

p1_stacked_by_caller

############################################################
## Plot 1C. Top N genes separately for losses and gains,
## stacked by caller + unique affected sample marker
############################################################

unique_sample_marker <- top_genes_by_unique_samples %>%
  group_by(cnv_type) %>%
  mutate(
    gene_name = forcats::fct_reorder(gene_name, unique_affected_samples)
  ) %>%
  ungroup()

p1_stacked_by_caller_with_marker <- ggplot(
  plot_top_genes_stacked_by_caller,
  aes(
    x = gene_name,
    y = n_samples_by_caller,
    fill = caller
  )
) +
  geom_col(position = "stack") +
  geom_point(
    data = unique_sample_marker,
    aes(
      x = gene_name,
      y = unique_affected_samples
    ),
    inherit.aes = FALSE,
    colour = "black",
    size = 2.2
  ) +
  coord_flip() +
  facet_wrap(~ cnv_type, scales = "free_y") +
  labs(
    x = "Gene",
    y = "Affected samples counted per caller",
    fill = "CNV caller",
    title = paste0("Top ", top_n_genes, " genes affected by CNVs"),
    subtitle = "Top genes selected separately for losses and gains; black points show unique affected samples"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold")
  )

p1_stacked_by_caller_with_marker

# correct gene order

library(tidytext)

top_n_genes <- 50

top_genes_by_unique_samples <- gene_burden_summary %>%
  group_by(cnv_type, gene_name) %>%
  summarise(
    unique_affected_samples = max(n_samples_affected),
    .groups = "drop"
  ) %>%
  group_by(cnv_type) %>%
  arrange(desc(unique_affected_samples), .by_group = TRUE) %>%
  slice_head(n = top_n_genes) %>%
  ungroup()

plot_top_genes_stacked_by_caller <- gene_caller_summary %>%
  semi_join(
    top_genes_by_unique_samples,
    by = c("gene_name", "cnv_type")
  ) %>%
  left_join(
    top_genes_by_unique_samples,
    by = c("gene_name", "cnv_type")
  ) %>%
  mutate(
    gene_name_ordered = reorder_within(
      gene_name,
      unique_affected_samples,
      cnv_type
    )
  )

unique_sample_marker <- top_genes_by_unique_samples %>%
  mutate(
    gene_name_ordered = reorder_within(
      gene_name,
      unique_affected_samples,
      cnv_type
    )
  )

p1_stacked_by_caller_with_marker <- ggplot(
  plot_top_genes_stacked_by_caller,
  aes(
    x = gene_name_ordered,
    y = n_samples_by_caller,
    fill = caller
  )
) +
  geom_col(position = "stack") +
  geom_point(
    data = unique_sample_marker,
    aes(
      x = gene_name_ordered,
      y = unique_affected_samples
    ),
    inherit.aes = FALSE,
    colour = "white",
    size = 1.5
  ) +
  coord_flip() +
  facet_wrap(~ cnv_type, scales = "free_y") +
  scale_x_reordered() +
  labs(
    x = "Gene",
    y = "Affected samples counted per caller",
    fill = "CNV caller",
    title = paste0("Top ", top_n_genes, " genes affected by CNVs"),
    #subtitle = "Top genes selected separately for losses and gains; black points show unique affected samples"
  ) +
  theme_thesis() + scale_fill_manual(values = caller_cols)

p1_stacked_by_caller_with_marker

ggsave("thesis_out/plot1_top_genes_by_affected_samples_stacked_by_caller.png",
  p1_stacked_by_caller_with_marker,
  width = 9,
  height = 10,
  dpi = 300)

############################################################
## 10. Scatter plot: recurrence vs caller support
############################################################
library(ggrepel)

# Define recurrence threshold
recurrence_cutoff <- 250
# or use:
# recurrence_cutoff <- quantile(gene_burden_summary$n_samples_affected, 0.90)

label_df <- gene_burden_summary %>%
  filter(
    n_samples_affected >= recurrence_cutoff,
    median_callers_per_affected_sample >= 3
  )

p3 <- gene_burden_summary %>%
  filter(n_samples_affected > 0) %>%
  ggplot(
    aes(
      x = n_samples_affected,
      y = mean_callers_per_affected_sample,
      colour = cnv_type
    )
  ) +
  geom_point(alpha = 0.6, size = 2) +
  geom_text_repel(
    data = label_df,
    aes(label = gene_name),
    size = 3.5,
    box.padding = 0.4,
    point.padding = 0.5,
    max.overlaps = Inf, nudge_y = 0.1,
    show.legend = FALSE
  ) +
#  scale_x_log10() +
  labs(
    x = "Number of affected samples (log10 scale)",
    y = "Mean number of callers per affected sample",
    colour = "CNV type",
    title = "Gene-level CNV recurrence versus caller support",
    subtitle = paste0(
      "Labels indicate genes detected by ≥3.5 callers on average and affecting ≥",
      recurrence_cutoff,
      " samples"
    )
  ) +
  theme_thesis()+ scale_color_manual(values = cnv_cols)


p3


############################################################
## 11. Heatmap-style caller support plot
############################################################

top_gene_names <- gene_burden_summary %>%
  arrange(desc(n_samples_affected)) %>%
  distinct(gene_name) %>%
  slice_head(n = 30) %>%
  pull(gene_name)

plot_caller_heatmap <- gene_caller_summary %>%
  filter(gene_name %in% top_gene_names) %>%
  mutate(
    gene_name = factor(gene_name, levels = rev(top_gene_names))
  )

p4 <- ggplot(
  plot_caller_heatmap,
  aes(
    x = caller,
    y = gene_name,
    fill = n_samples_by_caller
  )
) +
  geom_tile(colour = "white") +
  facet_wrap(~ cnv_type) +
  labs(
    x = "Caller",
    y = "Gene",
    fill = "Affected samples",
    title = "Caller-specific gene CNV burden",
    subtitle = "Top genes ranked by number of affected samples"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold")
  )

p4



####################
# frequency of callers
######################

caller_combo_df <- gene_sample_summary %>%
  mutate(
    # split, sort, recombine → ensures "A;B" == "B;A"
    caller_set = callers_overlapping_gene %>%
      str_split(";") %>%
      map(~ sort(.x)) %>%
      map_chr(~ paste(.x, collapse = ";"))
  )

caller_combo_freq <- caller_combo_df %>%
  count(cnv_type, caller_set, sort = TRUE) %>%
  rename(n_gene_sample_events = n)


caller_combo_freq <- caller_combo_freq %>%
  group_by(cnv_type) %>%
  mutate(
    prop = n_gene_sample_events / sum(n_gene_sample_events),
    perc = prop*100
  ) %>% ungroup()

top_combos <- caller_combo_freq %>%
  group_by(cnv_type) %>%
  slice_max(n_gene_sample_events, n = 10) %>%
  ungroup()

caller_combo_freq <- caller_combo_freq %>%
  mutate(
    n_callers = str_count(caller_set, ";") + 1,
    combo_class = case_when(
      n_callers == 1 ~ "Single caller",
      n_callers == 2 ~ "Pairwise agreement",
      n_callers == 3 ~ "Three-way agreement",
      n_callers >= 4 ~ "Full consensus"
    )
  )

caller_combo_freq %>% dplyr::select(-4) %>% arrange(cnv_type) %>% mutate(
  across(where(is.numeric), ~ round(.x, 2))
) %>%  gt() 

dim(caller_combo_df)
sum(caller_combo_freq$n_gene_sample_events)

# plot
library(ggplot2)

caller_combo_freq$caller_set <- gsub('Nexus','N',caller_combo_freq$caller_set)
caller_combo_freq$caller_set <- gsub('PennCNV','P',caller_combo_freq$caller_set)
caller_combo_freq$caller_set <- gsub('QuantiSNP','Q',caller_combo_freq$caller_set)
caller_combo_freq$caller_set <- gsub('iPattern','I',caller_combo_freq$caller_set)

p_combo <- caller_combo_freq %>%
  group_by(cnv_type) %>%
#  slice_max(n_gene_sample_events, n = 10) %>%
  ungroup() %>%
  mutate(
    caller_set = fct_reorder(caller_set, n_gene_sample_events)
  ) %>%
  ggplot(aes(x = caller_set, y = n_gene_sample_events, fill = as.factor(n_callers))) +
  geom_col(show.legend = T) +
  coord_flip() +
  facet_wrap(~ cnv_type, scales = "free_y") +
  labs(
    x = "Caller combination",
    y = "Number of gene-sample CNV events",
    fill= "Supporting callers",
#    subtitle = "Based on gene-level CNV overlap events"
  ) +
  theme_thesis() + scale_fill_brewer(palette = 'Reds')

p_combo

caller_combo_freq$combo_class <- factor(caller_combo_freq$combo_class, levels 
                                        = c("Single caller", "Pairwise agreement" ,"Three-way agreement" ,"Full consensus"))

p_combo_class <- caller_combo_freq %>%
  group_by(cnv_type, combo_class) %>%
  summarise(total = sum(n_gene_sample_events), .groups = "drop") %>%
  ggplot(aes(x = combo_class, y = total, fill = cnv_type)) +
  geom_col(position = "dodge") +
  labs(
    x = " ",
    y = "Gene-sample CNV events",
   title = "CNV caller agreement structure",
    fill='CNV type'
  ) +
  theme_thesis() + scale_fill_manual(values = cnv_cols)

p_combo_class/p_combo


ggsave("thesis_out/gene_concordance.png",
       plot = p_combo_class/p_combo,
       width = 8,  height = 9,
       dpi = 300)


######################
## new summary
##################

library(dplyr)

gene_summary <- gene_burden_with_callers %>%
  mutate(
    recurrence_class = case_when(
      n_samples_affected == 1 ~ "Singleton",
      n_samples_affected <= 20 ~ "Rare (2-20)",
      n_samples_affected <= 100 ~ "Moderately recurrent (20-100)",
      TRUE ~ "Highly recurrent (>100)"
    )
  )

recurrence_summary <- gene_summary %>%
  count(cnv_type, recurrence_class) %>%
  group_by(cnv_type) %>%
  mutate(percent = round(100 * n / sum(n), 1))

recurrence_summary



gene_summary2 <- gene_summary %>%
  mutate(
    caller_support = case_when(
      median_callers_per_affected_sample < 1.5 ~ "Mean ~1 caller",
      median_callers_per_affected_sample < 2.5 ~ "Mean ~2 callers",
      median_callers_per_affected_sample < 3.5 ~ "Mean ~3 callers",
      TRUE                                  ~ "Mean ~4 callers"
    )
  )

gene_summary2 %>%
  count(cnv_type,
        caller_support
  ) %>%
  group_by(cnv_type) %>%
  mutate(percent = round(100*n/sum(n),1)) %>% ungroup()

caller_support_summary <-
  gene_summary2 %>%
  count(cnv_type,
    recurrence_class,
    caller_support
  ) %>%
  group_by(cnv_type) %>%
  mutate(percent = round(100*n/sum(n),1)) %>% ungroup()



sum(caller_support_summary$n)

caller_support_summary %>% gt()

########## get genes

gene_summary2 %>% filter(recurrence_class=='Highly recurrent (>100)' & caller_support=="Mean ~4 callers") %>%
  group_by(cnv_type) %>% summarise(dplyr::n())


gene_summary2 %>% filter(recurrence_class=='Highly recurrent (>100)') %>%
  group_by(cnv_type, caller_support) %>% summarise(dplyr::n())

###

#4. Highly recurrent genes detected by ALL FOUR callers
highly_recurrent_all4 <-
  gene_summary2 %>%
  filter(
    recurrence_class=="Moderately recurrent (20-100)",
    caller_support=="Mean ~4 callers"
  ) %>%
  arrange(desc(n_samples_affected))

highly_recurrent_all4

highly_recurrent_all4$gene_name

#Highly recurrent singleton-caller genes
highly_recurrent_singletons <-
  gene_summary2 %>%
  filter(
    recurrence_class=="Highly recurrent (>200)",
    max_callers_in_any_sample==1
  ) %>%
  arrange(desc(n_samples_affected))

#6. Genes always supported by all four callers
genes_all4 <-  highconf_gene_summary %>%
  filter(n_samples_any_caller ==
      n_samples_4_callers
  ) %>%
  arrange(desc(n_samples_any_caller))

#7. Genes always singleton
genes_singleton_only <-
  highconf_gene_summary %>%
  filter(
    n_samples_2plus_callers==0
  ) %>%
  arrange(desc(n_samples_any_caller))

#8. Percentage of samples with consensus
gene_consensus <-
  highconf_gene_summary %>%
  mutate(
    pct_consensus =
      100 *
      n_samples_4_callers /
      n_samples_any_caller
  )

gene_consensus <-
  gene_consensus %>%
  mutate(
    consensus_class =
      case_when(
        pct_consensus==100 ~ "Always 4 callers",
        pct_consensus>=75 ~ "Mostly 4 callers",
        pct_consensus>=50 ~ "Half 4 callers",
        TRUE ~ "Rare 4 caller support"
      )
  )

#9. Which caller contributes unique genes?
caller_singletons <-
  gene_summary %>%
  filter(max_callers_in_any_sample==1) %>%
  mutate(
    singleton_caller =
      case_when(
        n_samples_Nexus>0 ~ "Nexus",
        n_samples_PennCNV>0 ~ "PennCNV",
        n_samples_QuantiSNP>0 ~ "QuantiSNP",
        n_samples_iPattern>0 ~ "iPattern"
      )
  ) %>%
  count(singleton_caller)

#10. Exclusive detection patterns
gene_patterns <-
  gene_summary %>%
  mutate(
    
    Nexus =
      n_samples_Nexus>0,
    
    PennCNV =
      n_samples_PennCNV>0,
    
    QuantiSNP =
      n_samples_QuantiSNP>0,
    
    iPattern =
      n_samples_iPattern>0
  ) %>%
  mutate(
    pattern =
      paste0(
        ifelse(Nexus,"N","-"),
        ifelse(PennCNV,"P","-"),
        ifelse(QuantiSNP,"Q","-"),
        ifelse(iPattern,"I","-")
      )
  )

pattern_summary <-
  gene_patterns %>%
  count(pattern, sort=TRUE)

highconf_gene_summary %>%
  filter(
    n_samples_any_caller >= recurrence_threshold,
    n_samples_any_caller == n_samples_4_callers
  )
