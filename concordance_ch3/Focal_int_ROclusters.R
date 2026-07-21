# interrogate focal alterations using RO clusters

library(GenomicRanges)
library(dplyr)
source('functions_new.R')
library(karyoploteR)
library(gt)

###############
### plot RO pairs
##################
breakpoint_variance <- readRDS('breakpoint_variance_out.Rds')

# stats
caller_combo_by_type <- breakpoint_variance %>%
  count(n_callers, loss_gain, callers, name = "n") %>%
  group_by(loss_gain) %>%
  mutate(
    fraction = n / sum(n),
    percent = 100 * fraction
  ) %>%
  arrange(loss_gain, desc(n)) %>% ungroup

# get table similar to gene level concordance
caller_combo_by_type %>% dplyr::select(2,3,4,6,1) %>% arrange(loss_gain) %>% mutate(
  across(where(is.numeric), ~ round(.x, 2))) %>%  gt() 

sum(caller_combo_by_type$n)

caller_combo_support <- breakpoint_variance %>%
  count(n_callers, loss_gain, callers, name = "n") %>%
  group_by(n_callers) %>%
  mutate(
    fraction = n / sum(n),
    percent = 100 * fraction
  ) %>%
  arrange(n_callers, desc(n))

caller_combo_support

# plot
caller_combo_support$callers <- gsub('Nexus','N',caller_combo_support$callers)
caller_combo_support$callers <- gsub('PennCNV','P',caller_combo_support$callers)
caller_combo_support$callers <- gsub('QuantiSNP','Q',caller_combo_support$callers)
caller_combo_support$callers <- gsub('iPattern','I',caller_combo_support$callers)

p_combo <- caller_combo_support %>%
  group_by(loss_gain) %>%
  #  slice_max(n_gene_sample_events, n = 10) %>%
  ungroup() %>%
  mutate(
    caller_set = fct_reorder(callers, n)
  ) %>%
  ggplot(aes(x = caller_set, y = n, fill = as.factor(n_callers))) +
  geom_col(show.legend = T) +
  coord_flip() +
  scale_y_continuous(labels = comma)+
  facet_wrap(~ loss_gain, scales = "free_y") +
  labs(
    x = "Caller combination",
    y = "Number of RO-based CNV clusters",
    #title = "Most frequent CNV caller combinations",
    fill= "Supporting callers",
    #    subtitle = "Based on gene-level CNV overlap events"
  ) +
  theme_thesis() + scale_fill_brewer(palette = 'Reds')

p_combo

# plot 2
caller_combo_support <- caller_combo_support %>%
  mutate(caller_set = case_when(
      n_callers == 1 ~ "Single caller",
      n_callers == 2 ~ "Pairwise agreement",
      n_callers == 3 ~ "Three-way agreement",
      n_callers >= 4 ~ "Full consensus"
    )
  )

caller_combo_support$caller_set <- factor(caller_combo_support$caller_set, levels 
                                        = c("Single caller", "Pairwise agreement" ,"Three-way agreement" ,"Full consensus"))

p_combo_class <- caller_combo_support %>%
  group_by(loss_gain, caller_set) %>%
  summarise(total = sum(n), .groups = "drop") %>%
  ggplot(aes(x = caller_set, y = total, fill = loss_gain)) +
  geom_col(position = "dodge") +
  scale_y_continuous(labels = comma)+
  labs(
   x = "",
    y = "RO-based CNV clusters",
    title = "CNV caller agreement structure",
    fill='CNV type'
  ) +
  theme_thesis() + scale_fill_manual(values = cnv_cols)

p_combo_class
#library(patchwork)
p_combo_class / p_combo

ggsave("thesis_out/RO_cluster_concordance.png",
       plot = p_combo_class/p_combo,
       width = 7,  height = 8,
       dpi = 300)


#######################
# RO pair manual
######################
breakpoint_variance <- readRDS('breakpoint_variance_out.Rds')
nrow(breakpoint_variance)

breakpoint_variance$callers <- gsub('Nexus','N',breakpoint_variance$callers)
breakpoint_variance$callers <- gsub('PennCNV','P',breakpoint_variance$callers)
breakpoint_variance$callers <- gsub('QuantiSNP','Q',breakpoint_variance$callers)
breakpoint_variance$callers <- gsub('iPattern','I',breakpoint_variance$callers)


breakpoint_variance <- breakpoint_variance %>%
  mutate(width_mb = union_width_bp / 1e6,
         cnv_class = ifelse(width_mb <= 1, "Focal", "Broad")
  )

breakpoint_variance %>% group_by(loss_gain, cnv_class, callers ) %>% 
  count() %>% pivot_wider(
    names_from = cnv_class,
    values_from = n,
    values_fill = 0,
    names_prefix = "count_"
  ) %>%
  arrange(loss_gain, desc(count_Focal)) %>% ungroup %>% gt()

#RO_breakpont.var.fil <- RO_breakpont.var %>% filter(callers == 'iPattern;PennCNV;QuantiSNP')

# plot RO clusters using lrr.baf

#dir.create('thesis_out/LrrBaf_plts_focal')

clusters.4 <- RO_breakpont.var.fil %>% filter(n_callers==4)
nrow(clusters.4)
hist(clusters.4$median_call_width_bp)

plt.dat <- clusters.4[500,]

plot_LRR.BAF.RO(clusters.4[500,], span = 2, p.size = 0.5)

grp <- 'Focal'
RO_breakpont.var.fil <- breakpoint_variance %>% filter(loss_gain =='Gain' & cnv_class == grp) %>% filter(n_callers == 2)
nrow(RO_breakpont.var.fil)

mclapply(seq_len(nrow(RO_breakpont.var.fil)), function(i) {
  plot_LRR.BAF.RO(RO_breakpont.var.fil[i,], span = 1, p.size = 0.7, group = grp)
}, mc.cores = 4)


###############################
### PLOT gene concordant events
##############################
gene_sample_summary <- readRDS('gene_sample_summary.rds')
gene_cnv_burden <- read.csv('gene_cnv_burden_with_caller_support_hg19.csv')
dim(gene_sample_summary)

gene <- 'OR4K2'
#search

to.plt <- gene_cnv_burden[grep(gene, gene_cnv_burden$gene_name),]
sam.plt <- gene_sample_summary[grep(gene, gene_sample_summary$gene_name),]

unique(gene_sample_summary$callers_overlapping_gene)
gene.sam.fil <- gene_sample_summary %>% filter(callers_overlapping_gene=='Nexus')
plot(table(gene.sam.fil$gene_chr))
nrow(gene.sam.fil)


plot_LRR.BAF.gene(sam.plt[8,], start.add = 1, end.add = 1, p.size = 0.5, tick = 15)
