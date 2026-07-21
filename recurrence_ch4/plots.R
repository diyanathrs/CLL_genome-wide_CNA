# generate plots and tables for thesis

library(ggplot2)
library(patchwork)
library(gt)

###########################
# disjoint stats and figs
###########################
# read all.djs
df.dj.loss <- readRDS('CLL_filtered_disj_cox_Broad_WTscore_5Loss_OS.Rds') %>% mutate(type='Loss') 
df.dj.gain <- readRDS('CLL_filtered_disj_cox_Broad_WTscore_3Gain_OS.Rds') %>% mutate(type='Gain')
df.cll.djs <- rbind(df.dj.gain, df.dj.loss)
#saveRDS(df.dj.all, 'cll_disjoints_all.Rds')

df.cll.djs <- readRDS('cll_disjoints_all.Rds')
df.cll.djs <- df.cll.djs %>% mutate(single_pct=n_samples_ge1/sample_support*100) 
df.cll.djs <- df.cll.djs %>% mutate(three_pct=n_samples_ge3/sample_support*100) 
df.cll.djs <- df.cll.djs %>% mutate(above3_pct=(n_samples_ge3+n_samples_ge4)/sample_support*100) 
#for text
range(df.cll.djs$single_pct[df.cll.djs$single_pct > 0], na.rm = TRUE)
range(df.cll.djs$above3_pct[df.cll.djs$above3_pct > 0], na.rm = TRUE)
table(df.cll.djs$type)

## caller sup

summary_support <- df.cll.djs %>%
  group_by(type) %>%
  summarise(single_caller = sum(n_samples_ge1)/sum(sample_support),
    three_callers = sum(n_samples_ge3)/sum(sample_support),
    four_callers  = sum(n_samples_ge4)/sum(sample_support),
    above3_callers = sum(n_samples_ge3 + n_samples_ge4)/sum(sample_support) )

summary_support <- df.cll.djs %>%
  group_by(type) %>% summarise(single_caller = sum(n_samples_ge1 > 0),
            three_callers = sum(n_samples_ge3 > 0),
            four_callers  = sum(n_samples_ge4 > 0),
            above3_callers = sum(n_samples_ge3 > 0 | n_samples_ge4 > 0))


summary_support

# filtered
df.cll.djs.fil <- df.cll.djs %>% filter(weighted_score > 15)

df.cll.djs.fil %>%
  group_by(type) %>% 
  summarise(single_caller = sum(n_samples_ge1 > 0),
                               three_callers = sum(n_samples_ge3 > 0),
                               four_callers  = sum(n_samples_ge4 > 0),
                               above3_callers = sum(n_samples_ge3 > 0 | n_samples_ge4 > 0))
  
min(df.cll.djs.fil$width)
# summary
type_summary <- df.cll.djs %>%  group_by(type) %>%
  summarise(Total_disjoints = n(),
    Median_width_kb = median(width)/1000,
    Mean_width_kb = mean(width)/1000,
    Median_sample_support = mean(sample_support),
    Mean_sample_support = mean(sample_support),
    Median_weighted_score = median(weighted_score),
    Mean_weighted_score = mean(weighted_score),
    Median_mean_callers = median(mean_callers),
    Mean_mean_callers = mean(mean_callers),
    Median_confidence = median(mean_confidence),
    Mean_confidence = mean(mean_confidence) )

data.frame(t(type_summary)) %>% gt(rownames_to_stub = T)

# Width distribution << go in thesis
size <- ggplot(df.cll.djs, aes(width/1000, fill=type)) +
  geom_histogram(bins=100, position = "identity",
                 alpha = 0.7) + scale_x_log10() +
  labs(x="Disjoint width (kb)",
       y="Number of disjoints", fill='CNV type') + scale_fill_manual(values = cnv_cols) + 
  theme_thesis() + theme(legend.position = 'blank')

# Sample support distribution
ggplot(df.cll.djs,
       aes(sample_support,
           fill=type))+ 
  geom_histogram(binwidth=4, position = "identity")+
  labs(x="Sample support",
       y="Number of disjoints")+ scale_fill_manual(values = cnv_cols) + 
  theme_thesis() 

#Weighted score distribution << take
weighted_score <- ggplot(df.cll.djs, aes(weighted_score,
           fill=type))+
  geom_histogram(bins=75, position = "identity",
                 alpha = 0.7)+ geom_vline(xintercept = 15, linetype = 3, linewidth = 0.6, color='white') +
  labs(x="Weighted score", y="No of disjoints", fill='CNV type')+ scale_fill_manual(values = cnv_cols) + 
  theme_thesis() 

weighted_score

#mean callers
mean.cllrs <- ggplot(df.cll.djs, aes(mean_callers,
           fill=type))+
  geom_histogram(binwidth=0.05, position = "identity",
                 alpha = 0.7)+
  labs(x="Mean callers supporting disjoints", 
       y= "Number of disjoints", fill='CNV type')+ scale_fill_manual(values = cnv_cols) + 
  theme_thesis() + theme(legend.position = 'blank')

#conf
ggplot(df.cll.djs, aes(mean_confidence,
           fill=type))+
  geom_histogram(bins=50)+
  labs(x="Mean confidence") + scale_fill_manual(values = cnv_cols) + 
  theme_thesis()

#Sample support vs weighted score << goes
sup.vs.conf <- ggplot(df.cll.djs,  aes(sample_support,
                         mean_confidence,
           colour=type))+
  geom_point(alpha=0.4)+
  labs(color='CNV type', x='Sample support', y='Mean confidence per disjoint',) + 
  geom_smooth(method="lm", se=FALSE, linetype = 1, linewidth = 0.5) + scale_color_manual(values = cnv_cols) + 
  theme_thesis()

fst.plt <- (size | mean.cllrs) / (sup.vs.conf | weighted_score) + plot_annotation(tag_levels = 'A')

# save to one
ggsave(filename = "plots/disj_stats.png",
       plot = fst.plt,
       width = 10, height = 8,
       dpi = 300)


# chr-wise
chr_counts <- df.cll.djs %>%  count(seqnames, type)

ggplot(chr_counts,
       aes(seqnames,
           n,
           fill=type))+
  geom_col(position="stack")+
  theme_bw()+
  theme(axis.text.x=element_text(angle=90,vjust=.5))

# chr-wise score distribution
df.cll.djs %>%   group_by(seqnames, type) %>%
  summarise(mean_confidence = mean(above3_pct, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(seqnames = factor(seqnames,
                           levels = paste0("chr", c(1:22, "X", "Y")))) %>%
  ggplot(aes(seqnames, mean_confidence, fill = type)) +
  geom_col(position = "dodge") +
  theme_bw() +
  labs(x = "Chromosome",
       y = "Mean confidence score") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

#Top recurrent disjoints
top_support <- df.cll.djs %>%
  arrange(desc(sample_support)) %>%
  dplyr::select(seqnames,start,end,width,type,
         sample_support,
         weighted_score,
         mean_callers,
         mean_confidence) %>%
  head(20)

top_support

#Caller agreement
ggplot(df.cll.djs,    aes(max_callers,
           fill=type))+
  geom_bar()+
  labs(x="Maximum callers supporting interval") +
  scale_fill_manual(values = cnv_cols)+ theme_thesis()


#Correlation matrix
library(GGally)

df.cll.djs %>%
  dplyr::select(sample_support,
         weighted_score,
         mean_callers,
         mean_confidence,
         max_callers) %>%
  ggpairs()


###############
#plot whole chr
###############
kp <- plotKaryotype(plot.type = 2, chromosomes = 'autosomal') #chromosomes = 'chr6'
kpPlotDensity(kp, data=gr.whole_chr[gr.whole_chr$loss_gain=='Loss'], col='red', data.panel = 1, r0 = 0.2, r1 = 1)
kpPlotDensity(kp, data=gr.whole_chr[gr.whole_chr$loss_gain=='Gain'], col='blue', data.panel = 2, r0 = 0.2, r1 = 1)

whole.df <- as.data.frame(gr.whole_chr) %>%
  count(caller, loss_gain, seqnames) %>%
  mutate(n = ifelse(loss_gain == "Loss", -n, n)
  )

whole <- ggplot(whole.df, aes(x = seqnames, y = n, fill = caller)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = 0, colour = "black", linewidth = 0.3) +
  scale_fill_manual(values = caller_cols) +
  scale_y_continuous(labels = abs) + ylim(c(-100, 600)) +
  labs(x = "Chromosome",
    y = "Number of broad calls",
    fill = "CNV Type") + theme_thesis() + 
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

whole
#plot non whole chr
kp <- plotKaryotype(plot.type = 2, chromosomes = 'autosomal') #chromosomes = 'chr6'
kpPlotDensity(kp, data=gr.non.whole_chr[gr.non.whole_chr$loss_gain=='Loss'],  col='red', data.panel = 1, r0 = 0.4, r1 = 1)
kpPlotDensity(kp, data=gr.non.whole_chr[gr.non.whole_chr$loss_gain=='Gain'], col='blue', data.panel = 2, r0 = 0.4, r1 = 1)

focal.df <- as.data.frame(gr.non.whole_chr) %>%
  count(caller, loss_gain, seqnames) %>%
  mutate(n = ifelse(loss_gain == "Loss", -n, n)
  )

focal <- ggplot(focal.df, aes(x = seqnames, y = n, fill = caller)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = 0, colour = "black", linewidth = 0.3) +
  scale_fill_manual(values = caller_cols) +
  scale_y_continuous(labels = abs) +
  labs(
    x = "Chromosome",
    y = "Number of focal calls",
    fill = "CNV Type"
  ) + theme_thesis() + theme(axis.text.x = element_text(angle = 45, hjust = 1))

focal

scnd.plt <- whole / focal + plot_annotation(tag_levels = 'A')

# save 
ggsave(filename = "plots/separated_cnv_chr.png",
       plot = scnd.plt,
       width = 8,
       height = 8,
       dpi = 300)




#djs per chr
df.cll_dj %>% group_by(seqnames) %>% count %>% ggplot(aes(seqnames, n)) + geom_bar(stat='identity')

df.cll_dj %>% group_by(seqnames) %>% summarise(sum=sum(sample_support))
df.cll_dj %>% group_by(seqnames) %>% summarise(mean_score=mean(mean_confidence), mean_callers=mean(mean_callers)) %>% 
  ggplot(aes(seqnames, mean_callers)) + geom_bar(stat = 'identity')

