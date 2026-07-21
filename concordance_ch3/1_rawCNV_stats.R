# This should be the first step
# Remove problematic samples before working on raw CNV calls
# Plot LRR SD and BAF drift, Waviness factor and remove
# Then use no of CNVs/sample - plot first - maybe merging adjacent would reduce this

library(GenomicRanges)
library(dplyr)
library(tidyr)
library(ggplot2)
library(gt)
library(scales)
library(flextable)
library(officer)

###############
### initial QC
###############
cnv_cols <- c("Gain" = "#2C7BB6",
              "Loss" = "#D7191C")

caller_cols <- c("PennCNV"   = "#1B9E77",
                 "QuantiSNP" = "#D95F02",
                 "iPattern"  = "#7570B3",
                 "Nexus"     = "#E7298A")


### load raw CNVs or gap.merged
source('functions.R')
source('functions_new.R')

# load all raw CNVs regardless of length or prob counts
cnv_gr <- load.raw_gap(snp.ct = 5, len = 50, use.merged.cnv = F, only.common.sam = F, qc.pass.only = T)
length(cnv_gr)
data.frame(mcols(cnv_gr)) %>%  group_by(study_centre) %>%
  summarise(n_samples = n_distinct(sample_id), .groups = "drop")

cnv_gr.gap <- load.raw_gap(snp.ct = 50, len = 1000, use.merged.cnv = T)

###################################
## KM plots
##################################
library(karyoploteR)

# change plot params based on sample length
pp <- getDefaultPlotParams(plot.type=2) # 60 and 120
pp$ideogramheight <- 100

pdf(file = 'thesis_out/KP_gap_CNV.pdf', width = 10, height = 10)

kp <- plotKaryotype(plot.type = 2, chromosomes = 'autosomal', cex=1, plot.params = pp) #chromosomes = 'chr6
#kpPlotRegions(kp, data = raw_cnv_gr_cov, data.panel = 1, col="#2C7BB6", border="#2C7BB6", r0=0.05, r1=1)

#kpDataBackground(kp, color = "#FFFFFFAA")
kpPlotDensity(kp, data = subset(cnv_gr.gap, cnv_gr.gap$loss_gain=='Gain'), data.panel = 1, ymax = 1800,
              col="#2C7BB6", border="#2C7BB6", r0=0.05, r1=1)
kpPlotDensity(kp, data = subset(cnv_gr.gap, cnv_gr.gap$loss_gain=='Loss'), data.panel = 2, ymax = 3000,
              col="#D7191C", border="#D7191C", r0=0.05, r1=1)
dev.off()

kpPlotDensity(kp, data = subset(cnv_gr, cnv_gr$loss_gain=='Gain'), 
              window.size = 0.5e6, data.panel = 2, col="#2C7BB6", border="#2C7BB6", r0=0.5, r1=1.2)
kpPlotDensity(kp, data = subset(cnv_gr, cnv_gr$loss_gain=='Loss'), 
              window.size = 0.5e6, data.panel = 2, col="#D7191C", border="#D7191C", r0=0.5, r1=-0.2)

kpPlotCoverage(kp, data = subset(cnv_gr, cnv_gr$loss_gain=='Loss'), col = cnv_type_cols[1], data.panel = 1)
kpPlotCoverage(kp, data = subset(cnv_gr, cnv_gr$loss_gain=='Gain'), col = cnv_type_cols[2], data.panel = 2)
kpPlotRegions(kp, data=gr.cll.cnv.sam[[caller.idx]], col='blue', r1 = 0.5)
kpPlotRegions(kp, data=gr.cll.gap.sam[[caller.idx]], col='red', r0 = 0.5)


########################
## bulk CNV tables - stats
########################

# present these data in early chapter 3 ###################################
# table with below columns
# Caller	Total CNVs	CNVs/sample median	Median CNV size	Median probes	Losses	Gains
# Purpose: shows whether one caller is systematically liberal or conservative.
###############################################

cnv_gr <- cnv_gr %>% filter(!sample_id %in% qc.fail) # use pre merge
length(cnv_gr)

cnv_df <- data.frame(sample_id    = cnv_gr$sample_id,
                     chr = seqnames(cnv_gr),
                     caller       = cnv_gr$caller,
                     loss_gain    = cnv_gr$loss_gain,
                     copy_number  = cnv_gr$CN,
                     no_of_probes = cnv_gr$no_of_probes,
                     width_bp     = width(cnv_gr)) %>%
  mutate(width_kb = width_bp / 1000, width_mb = width_bp / 1e6)

nrow(cnv_df)
length(unique(cnv_df$sample_id))

# only keep overlapping samples
# Identify samples present in all 4 callers
common_samples <- cnv_df %>% distinct(sample_id, caller) %>%
  count(sample_id, name = "n_callers") %>%
  filter(n_callers == 4) %>% pull(sample_id)

# Filter CNV table to only those samples
cnv_df <- cnv_df %>% filter(sample_id %in% common_samples)

cnvs_per_sample <- cnv_df %>% group_by(caller, sample_id) %>%
  summarise(n_cnvs = dplyr::n(), total_cnv_mb = sum(width_mb, na.rm = TRUE),
            .groups = "drop")

overview_table <- cnv_df %>%  group_by(caller) %>%
  summarise(`Total CNVs` = dplyr::n(),
            `Samples with CNVs` = n_distinct(sample_id),
            `Median size (kb)` = median(width_kb, na.rm = TRUE),
            `Mean size (kb)` = mean(width_kb, na.rm = TRUE),
            `Median probes` = median(no_of_probes, na.rm = TRUE),
            `Losses` = sum(loss_gain == "Loss", na.rm = TRUE),
            `Gains` = sum(loss_gain == "Gain", na.rm = TRUE),
            .groups = "drop") %>% 
  left_join(cnvs_per_sample %>% group_by(caller) %>%
              summarise(`Median CNVs/sample` = median(n_cnvs, na.rm = TRUE),
                        `Mean CNVs/sample` = mean(n_cnvs, na.rm = TRUE),
                        `Median burden/sample (Mb)` = median(total_cnv_mb, na.rm = TRUE),
                        `Mean burden/sample (Mb)` = mean(total_cnv_mb, na.rm = TRUE),
                        .groups = "drop"),  by = "caller")


overview_table_clean <- overview_table %>%
  mutate(across(where(is.numeric), ~ round(.x, 2)))

overview_gt <- overview_table_clean %>%
  gt() %>%  tab_header(title = "Overview of raw CNV calls by caller") %>% 
  cols_label(caller = "CNV caller")

overview_gt

#save
gtsave(overview_gt, "raw_CNV_caller_overview.html")


raw_cnv_by_caller <- cnv_df %>%
  group_by(caller) %>%
  summarise(n_cnvs = dplyr::n(),
    #n_samples = n_distinct(sample_id),
    median_cnvs_per_sample = median(table(sample_id)),
    mean_width_kb = mean(width_kb, na.rm = TRUE),
    median_width_kb = median(width_kb, na.rm = TRUE),
    min_width_kb = min(width_kb, na.rm = TRUE),
    max_width_mb = max(width_mb, na.rm = TRUE),
    median_probes = median(no_of_probes , na.rm = TRUE),
    mean_probes = mean(no_of_probes , na.rm = TRUE),
    .groups = "drop")

raw_cnv_by_caller


# cnv type
cnv_type_by_caller <- cnv_df %>%
  group_by(caller, loss_gain) %>%
  summarise(n_cnvs = dplyr::n(),
    n_samples = n_distinct(sample_id),
    median_width_kb = median(width_kb, na.rm = TRUE),
    median_probes = median(no_of_probes, na.rm = TRUE),
    .groups = "drop") %>%
  arrange(caller, loss_gain)

cnv_type_by_caller

# print table
ft <- flextable::flextable(raw_cnv_by_caller) %>%  flextable::autofit()

#########################
### PLOTS
#########################
library(ggplot2)

# CNV per sample
cnvs_per_sample <- cnv_df %>% group_by(sample_id, caller, loss_gain) %>%
  summarise(n_cnvs = length(sample_id), .groups = "drop")

# summary
cnvs_per_sample %>% group_by(caller) %>% summarise(mean_cnvs = mean(n_cnvs),
                                                   median_cnvs = median(n_cnvs), sd_cnvs = sd(n_cnvs),
                                                   min_cnvs = min(n_cnvs), max_cnvs = max(n_cnvs))

ggplot(cnvs_per_sample, aes(x = caller, y = n_cnvs, col = loss_gain)) +
  geom_boxplot(outlier.alpha = 0.3, position = position_dodge(width = 0.8)) +
  labs(x = "CNV caller", y = "Number of CNVs per sample",
    title = "CNV burden per sample by caller",
    #subtitle = "Separated by copy-number losses and gains"
    ) +
  theme_thesis() + scale_color_manual(values = cnv_cols)+
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))


###########################
## gains vs losses
###########################
gain_loss_summary <- cnv_df %>%  group_by(caller, loss_gain) %>% 
  summarise(n = length(sample_id), .groups = "drop")

# plot
ggplot(gain_loss_summary, aes(x = caller, y = n, fill = loss_gain)) +
  geom_col(position = "dodge") + theme_bw() +
 scale_fill_manual(values = cnv_cols) 

##############################
### CNV size distribution
#############################
cnv_df$caller <- factor(cnv_df$caller, levels = c('PennCNV','QuantiSNP','iPattern','Nexus' ))

size_hist <- ggplot(cnv_df, aes(x = width_bp, fill = caller)) +
  geom_histogram(bins = 75, alpha = 0.6, position = "identity") +
  scale_x_log10() +   theme_bw(base_size = 12) +
  labs(x = "CNV size (kb, log scale)", y = "Count", fill='CNV Caller')+
  scale_fill_manual(values = caller_cols) 

size_hist

ggsave(filename = "thesis_out/cnv_size_hist.png",
  plot = size_hist, width = 8,
  height = 5, dpi = 300)


ggplot(cnv_df, aes(x = width_bp, fill = caller)) +
  geom_density(alpha = 0.8, position = "identity") +
  scale_x_log10() +   theme_bw() +
  labs(x = "CNV size (kb, log scale)", y = "Count")+
  scale_fill_brewer(palette = 'Set2')

##############################
# CNV size versus probe count
##############################
ggplot(cnv_df, aes(x = no_of_probes,
           y = width_kb, colour = caller)) +
  geom_point(alpha = 0.3) +
  scale_x_log10() +  scale_y_log10() +
  theme_bw()

ggplot(cnv_df, aes(x = width_mb , y = probes)) +
  geom_point(alpha = 0.25, size = 0.8) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.8) +
  scale_x_log10(labels = label_number()) +
  scale_y_log10(labels = label_number()) +
  scale_colour_manual(values = cnv_type_cols) +
  facet_wrap(~ caller) +
  labs(x = "CNA length (Mb, log10 scale)",
    y = "Number of probes (log10 scale)",
    colour = "CNA type",
    title = "Relationship between CNA length and probe count",
    subtitle = "Each point represents one CNA call") +
  theme_bw()

#### into size bins
cnv_df <- cnv_df %>% mutate(size_bin = cut(width_bp, breaks = c(1, 1e3, 1e4, 1e5, 1e6, 1e7, Inf),
                        labels = c("<1kb","1-10kb", "10-100kb",
                                   "100kb-1Mb", "1-10Mb", ">10 Mb"),  right = FALSE))

# plots for raw CNV calls
cnv_df %>% ggplot(aes(size_bin, fill=caller)) + geom_bar(position='dodge')

###################################
### Study-centre effects
####################################
centre_summary <- cnv_df %>%
  group_by(study_centre, sample_id, caller) %>%
  summarise(n_cnvs = length(sample_id),
    .groups = "drop")

ggplot(centre_summary,
       aes(x = study_centre, y = n_cnvs, fill = caller)) +
  geom_boxplot() +  theme_bw() + coord_flip()




########################################
## New plots
#########################################
library(dplyr)
library(ggplot2)
library(scales)
library(tidyr)

## -------------------------------------------------------
## Manual colours
## -------------------------------------------------------
cnv_gr <- load.raw_gap(snp.ct = 5, len = 5, use.merged.cnv = F, only.common.sam = F, qc.pass.only = T)

cnv_df <- data.frame(sample_id    = cnv_gr$sample_id,
                     caller       = cnv_gr$caller,
                     loss_gain    = cnv_gr$loss_gain,
                     copy_number  = cnv_gr$CN,
                     no_of_probes = cnv_gr$no_of_probes,
                     width_bp     = width(cnv_gr)) %>%
  mutate(width_kb = width_bp / 1000, width_mb = width_bp / 1e6)

# no of CNVs
plot_df <- cnv_df %>%
  dplyr::count(caller, loss_gain) %>%
  mutate(n = ifelse(loss_gain == "Loss", -n, n))

no.cnv <- ggplot(plot_df,
       aes(caller, n, fill = loss_gain)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = cnv_cols) +
  scale_y_continuous(labels = abs) +
  labs(x = "",
       y = "Number of CNVs",
       fill = "") + theme_thesis() +
  theme(legend.position = 'none')
no.cnv
# Total genomic coverage

plot_df <- cnv_df %>%
  group_by(caller, loss_gain) %>%
  summarise(total_mb = sum(width_mb),
            .groups = "drop") %>%
  mutate(total_mb = ifelse(loss_gain == "Loss",
                           -total_mb,
                           total_mb))

coverage <- ggplot(plot_df,
                   aes(caller,
                       total_mb,
                       fill = loss_gain)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = cnv_cols) +
  scale_y_continuous(labels = abs) +
  labs(y = "Total CNV coverage (Mb)",
       x = "") + theme_thesis() +
  theme(legend.position = 'none')

#8. Chromosome-wise burden (counts)

chr_df <- data.frame(
  chr = as.character(seqnames(cnv_gr)),
  caller = cnv_gr$caller,
  loss_gain = cnv_gr$loss_gain)

plot_df <- chr_df %>%
  dplyr::count(chr, caller, loss_gain) %>%
  mutate(n = ifelse(loss_gain == "Loss",
                    -n,
                    n)) %>%
  mutate(chr = case_when(
    chr %in% as.character(1:22) ~ paste0("chr", chr),
    chr == "X" ~ "chrX",
    chr == "Y" ~ "chrY",
    TRUE ~ chr),
    chr = factor(chr, levels = c(paste0("chr", 1:22), "chrX", "chrY")))

chr_cnv <- ggplot(plot_df,
       aes(chr,
           n,
           fill = loss_gain)) +
  geom_col() +
  #facet_wrap(~caller, ncol = 1) +
  scale_fill_manual(values = cnv_cols) +
  scale_y_continuous(labels = abs) +
  labs(x = "Chromosome",
       y = "Number of CNVs") + theme_thesis()+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  , legend.position = 'none')
  
# chr wise coverage
chr_cov <- cnv_df %>%
  mutate(chr = as.character(seqnames(cnv_gr))) %>%   # if not already in df
  group_by(chr, caller, loss_gain) %>%
  summarise(total_mb = sum(width_mb), .groups = "drop") %>%
  mutate(total_mb = ifelse(loss_gain == "Loss",
                           -total_mb,
                           total_mb)) %>%
  mutate(chr = case_when(
    chr %in% as.character(1:22) ~ paste0("chr", chr),
    chr == "X" ~ "chrX",
    chr == "Y" ~ "chrY",
    TRUE ~ chr),
    chr = factor(chr, levels = c(paste0("chr", 1:22), "chrX", "chrY")))

chr_cov.plt <- ggplot(chr_cov,
       aes(x = chr,
           y = total_mb,
           fill = loss_gain)) +
  geom_col() +
 # facet_wrap(~caller, ncol = 1) +
  scale_fill_manual(values = cnv_cols) +
  scale_y_continuous(labels = abs) +
  labs(
   # title = "Chromosome-wise CNV Genomic Coverage by Caller",
    x = "Chromosome",
    y = "Total CNV coverage (Mb)",
    fill = ""
  ) + theme_thesis() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)  ,legend.position = 'none')



first.plt <- (no.cnv | coverage) / (chr_cnv | chr_cov.plt)
first.plt

ggsave(filename = "thesis_out/1_cnvs_caller_chr.png",
       plot = first.plt, width = 10,
       height = 7, dpi = 300)


#3. CNV size distribution
size.box <- ggplot(cnv_df,
       aes(caller,
           width_mb,
           fill = caller)) +
  geom_boxplot(outlier.alpha = 0.15) +
  geom_jitter()+
  scale_fill_manual(values = caller_cols) +
  scale_y_log10() +
  labs(y = "CNV size (Mb, log10)",
       x = "")

size.box

#6. Probe count distribution
ggplot(cnv_df,
       aes(caller,
           no_of_probes,
           fill = caller)) +
  geom_boxplot(outlier.alpha = 0.2) +
  scale_fill_manual(values = caller_cols) +
  scale_y_log10() +
  labs(y = "Number of probes (log10)",
       x = "")

#7. Gain/Loss proportions
plot_df <- cnv_df %>%
  count(caller, loss_gain) %>%
  group_by(caller) %>%
  mutate(prop = n / sum(n))

ggplot(plot_df,
       aes(caller,
           prop,
           fill = loss_gain)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = cnv_cols) +
  scale_y_continuous(labels = percent) +
  labs(y = "Percentage of CNVs",
       x = "")


#9. CNV size density
ggplot(cnv_df,
       aes(width_mb,
           colour = caller)) +
  geom_density() +
  scale_colour_manual(values = caller_cols) +
  scale_x_log10() +
  labs(x = "CNV size (Mb)",
       y = "Density")


############
# size
#############
library(forcats)
cnv_gr <- load.raw_gap(snp.ct = 5, len = 5, use.merged.cnv = F, only.common.sam = T)

cnv_df <- data.frame(sample_id    = cnv_gr$sample_id,
                     caller       = cnv_gr$caller,
                     loss_gain    = cnv_gr$loss_gain,
                     copy_number  = cnv_gr$CN,
                     no_of_probes = cnv_gr$no_of_probes,
                     width_bp     = width(cnv_gr)) %>%
  mutate(width_kb = width_bp / 1000, width_mb = width_bp / 1e6)

size_breaks <- c(0, 0.05, 0.1, 0.5, 1, 5, Inf)
size_labels <- c("<50 kb",
                 "50–100 kb",
                 "100–500 kb",
                 "0.5–1 Mb",
                 "1–5 Mb",
                 ">5 Mb")

plot_df <- cnv_df %>%
  mutate(size_bin = cut(width_mb,
                        breaks = size_breaks,
                        labels = size_labels,
                        right = FALSE)) %>%
  dplyr::count(caller, loss_gain, size_bin) %>%
  group_by(caller, loss_gain) %>%
  mutate(percent = 100 * n / sum(n))

ggplot(plot_df,
       aes(size_bin, percent, fill = caller)) +
  geom_col(position = position_dodge(width = 0.8)) +
  facet_wrap(~loss_gain) +
  scale_fill_manual(values = caller_cols) +
  labs( title = "Distribution of CNV Sizes by Caller",
        x = "CNV size",
       y = "CNVs (%)",
       fill = "") +
  theme_thesis() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

### new plot CNV burden stacked ##################

cnv_df <- data.frame(sample_id    = cnv_gr$sample_id,
                     chr = gsub('chr','',seqnames(cnv_gr)),
                     start= start(cnv_gr),
                     end  = end(cnv_gr),
                     caller       = cnv_gr$caller,
                     loss_gain    = cnv_gr$loss_gain,
                     copy_number  = cnv_gr$CN,
                     no_of_probes = cnv_gr$no_of_probes,
                     width_bp     = width(cnv_gr)) %>%
  mutate(width_kb = width_bp / 1000, width_mb = width_bp / 1e6)


chr_burden <- cnv_df %>%
  group_by(chr, caller, loss_gain) %>%
  summarise(cnv_count = dplyr::n(),
            .groups = "drop"
  )


#library(ggplot2)

cnv.chr <- ggplot(chr_burden, aes(x = chr, y = cnv_count, fill = caller)) +
  geom_bar(stat = "identity") +
  labs(title = "Chromosome-wise CNV burden by caller",
    x = "Chromosome",
    y = "Number of CNVs",
    fill = "Caller"
  ) + facet_grid(~loss_gain)+ scale_y_continuous(labels = comma)+
  theme_thesis() + scale_fill_manual(values = caller_cols) 
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

chr_burden <- cnv_df %>%
  mutate(width = end - start) %>%
  group_by(chr, caller, loss_gain) %>%
  summarise(
    total_mb = sum(width) / 1e6,
    .groups = "drop"
  )

# coverage
cov_chr <- ggplot(chr_burden, aes(x = chr, y = total_mb, fill = caller)) +
  geom_bar(stat = "identity") +
  labs(
   # title = "Chromosome-wise CNV burden by caller",
    x = "Chromosome",
    y = "Total CNV coverage (Mb)",
    fill = "Caller"
  ) + facet_grid(~loss_gain)+
  scale_y_continuous(labels = comma)+
  theme_thesis() + scale_fill_manual(values = caller_cols) +
  theme(legend.position = 'none')

cnv.chr / cov_chr

ggsave(filename = "thesis_out/2_stacked_caller_chr.png",
       plot = cnv.chr / cov_chr, width = 14,
       height = 9, dpi = 300)

####################
### focal vs broad CNVs
########################
library(dplyr)

cnv_df <- cnv_df %>%
  mutate(width_bp = end - start,
    width_mb = width_bp / 1e6,
    cnv_class = ifelse(width_mb <= 1, "Focal", "Broad")
  )

overall_table <- cnv_df %>%
  group_by(cnv_class) %>%
  summarise(
    n_cnvs = n(),
    total_mb = sum(width_mb),
    mean_size_mb = mean(width_mb),
    median_size_mb = median(width_mb),
    .groups = "drop"
  )

overall_table

# caller wise
caller_table <- cnv_df %>%
  group_by(caller, cnv_class) %>%
  summarise(
    n_cnvs = n(),
    total_mb = sum(width_mb),
    mean_size_mb = mean(width_mb),
    .groups = "drop"
  ) %>%
  arrange(caller, cnv_class)

caller_table

caller_prop <- cnv_df %>%
  group_by(caller, cnv_class, chr) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(caller) %>%
  mutate(prop = n / sum(n))

# plot
#library(ggplot2)

ggplot(caller_prop, aes(x = chr, y = n, fill = cnv_class)) +
  geom_bar(stat = "identity") +
  facet_wrap(~caller)+
  labs(
    title = "Focal vs Broad CNV burden by caller",
    x = "Caller",
    y = "CNV count",
    fill = "CNV class"
  ) + theme_thesis()

ggplot(caller_prop, aes(x = caller, y = prop, fill = cnv_class)) +
  geom_bar(stat = "identity") +
  labs(
    title = "Proportion of focal vs broad CNVs per caller",
    x = "Caller",
    y = "Proportion",
    fill = "CNV class"
  ) +
  theme_minimal()

# dist
overall_plot_df <- cnv_df %>%
  group_by(cnv_class) %>%
  summarise(n = n(), .groups = "drop")

ggplot(overall_plot_df, aes(x = cnv_class, y = n, fill = cnv_class)) +
  geom_bar(stat = "identity") +
  labs(
    title = "Overall focal vs broad CNV burden",
    x = "CNV class",
    y = "Count"
  ) +
  theme_minimal()

sample_focal <- cnv_df %>%
  group_by(sample_id, cnv_class) %>%
  summarise(n_cnvs = n(), .groups = "drop")

ggplot(sample_focal, aes(x = sample_id, y = n_cnvs, fill = cnv_class)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  theme(axis.text.x = element_blank()) +
  labs(
    title = "Sample-level focal vs broad CNV burden",
    x = "Sample",
    y = "CNV count"
  )
