## get cytobands
library(GenomicRanges)
library(rtracklayer)
library(ggrepel)

source('../concordance_ch3/functions.R')
source('../concordance_ch3/functions_new.R')

#######################################
## pre plots
##################################
# load genome-wide cox results and model
#load cox results for OS
df.cox.gain <- readRDS('CLL_filtered_disj_cox_Broad_WTscore_3Gain_OS.Rds') %>% mutate(cnv_type='Gain')
df.cox.loss <- readRDS('CLL_filtered_disj_cox_Broad_WTscore_3Loss_OS.Rds')  %>% mutate(cnv_type='Loss')
#comb
df.cox.OS.sig <- rbind(df.cox.gain, df.cox.loss) #%>% filter(FDR <= 0.05 & n_samples > 10)

dim(df.cox.OS.sig)
names(df.cox.OS.sig)
table(df.cox.OS.sig$type)
range(df.cox.OS.sig$n_samples)

df.cox.OS.sig %>% ggplot(aes(HR, -log10(FDR), fill=cnv_type, size = n_samples)) + 
  geom_point(alpha = 0.4, shape = 21, col='black') + 
  geom_hline(yintercept = -log10(0.05), ,colour = 'black', linetype = 2, linewidth = 0.4)+ 
  geom_vline(xintercept = 1, colour = 'black', linetype = 2, linewidth = 0.4)+
  scale_fill_manual(values = cnv_cols)+ theme_thesis()

#  convert back to GR
gr.cox.OS.sig <- df.cox.OS.sig[df.cox.OS.sig$type == "disjoint",] %>%
  tidyr::separate(feature, into = c("chr", "coords"), sep = ":") %>%
  tidyr::separate(coords, into = c("start", "end"), sep = "-") %>%
  mutate(start = as.numeric(start),
         end = as.numeric(end)
  ) %>% makeGRangesFromDataFrame(
    keep.extra.columns = TRUE,
    seqnames.field = "chr",
    start.field = "start",
    end.field = "end"
  )

# kp plot
kp <- plotKaryotype(plot.type = 2, chromosomes = 'autosomal') #chromosomes = 'chr6'
kpPlotCoverage(kp, data=gr.cox.OS.sig[gr.cox.OS.sig$cnv_type=='Loss'], col='red', data.panel = 1, r1 = 0.5)
kpPlotRegions(kp, data=gr.cox.OS.sig[gr.cox.OS.sig$cnv_type=='Gain'], col='blue', data.panel = 2, r1 = 0.5)

########################################
## remove adjacent results within loci
########################################
source('functions_feature_selection.R')

surv.lst <- c('OS','TTFT','PT')

df.cox.out <- list()
gr.cox.out <- GRangesList()
gr.cox.features <- GRangesList()
whole.cox.combined  <- list()

for (sv in surv.lst){
  print(sv)
  #cox_all_features.Rds
  #load cox results
  df.cox.gain <- readRDS(paste0('CLL_filtered_disj_cox_Broad_WTscore_3Gain_',sv,'.Rds')) %>% filter(FDR < 0.05) %>% mutate(cnv_type='Gain')
  df.cox.loss <- readRDS(paste0('CLL_filtered_disj_cox_Broad_WTscore_3Loss_',sv,'.Rds'))  %>% filter(FDR < 0.05) %>% mutate(cnv_type='Loss')
  # comb
  df.cox.out[[sv]] <- rbind(df.cox.gain, df.cox.loss) 
  #names(df.cox.out[[sv]])
  print(table( df.cox.out[[sv]]$type))
  # get whole chr
  whole.cox.combined[[sv]] <- df.cox.out[[sv]] %>% filter(type == "whole_chr")
  
  # convert to granges
  gr.cox.out[[sv]] <- df.cox.out[[sv]][df.cox.out[[sv]]$type == "disjoint",] %>%
    tidyr::separate(feature, into = c("chr", "coords"), sep = ":") %>%
    tidyr::separate(coords, into = c("start", "end"), sep = "-") %>%
    mutate(start = as.numeric(start),
           end = as.numeric(end)) %>%
    makeGRangesFromDataFrame(
      keep.extra.columns = TRUE,
      seqnames.field = "chr",
      start.field = "start",
      end.field = "end")
  print(elementNROWS(gr.cox.out[sv]))
  # reduce and pick top
  gr.cox.features[[sv]] <- reduce_cox_hits(gr = gr.cox.out[[sv]], priority = 'FDR', max_gap = 0)
  gr.cox.features[[sv]]$feature <- names(gr.cox.features[[sv]])
}

# combine features into one granges
# flatten gr list
df.cox.features <- GenomicRanges::as.data.frame(gr.cox.features)
# create all df.list and save
df.whole.cox.combined <- bind_rows(whole.cox.combined, .id = 'group_name')
df.whole.cox.combined$FDR <- df.whole.cox.combined$p # replace FDR

common.cols <- intersect(names(df.cox.features), names(df.whole.cox.combined))

features.cox.all.df <- rbind(df.cox.features[common.cols], df.whole.cox.combined[common.cols])
range(features.cox.all.df$FDR)
table(features.cox.all.df$type)

write.csv(features.cox.all.df, paste0('cox_all_features.csv'), quote = F, row.names = F)

########################
### All disjoint stats
########################
for (sv in surv.lst){
  print(sv)
  df.cox.gain <- readRDS(paste0('CLL_filtered_disj_cox_Broad_WTscore_3Gain_',sv,'.Rds')) %>% mutate(cnv_type='Gain')
  df.cox.loss <- readRDS(paste0('CLL_filtered_disj_cox_Broad_WTscore_3Loss_',sv,'.Rds'))  %>% mutate(cnv_type='Loss')
  # comb
  df.cox.out[[sv]] <- rbind(df.cox.gain, df.cox.loss) 
 }
df.cox.out.comb <- rbindlist(df.cox.out, idcol = 'surv')
names(df.cox.out.comb)
#stats
df.cox.out.comb %>% filter(FDR < 0.05) %>% group_by(surv, cnv_type) %>% count()
test <- df.cox.out.comb %>%  filter(FDR < 0.05) %>% group_by(cnv_type, feature, cytoband, n_samples) %>% count() 

common_features <- df.cox.out.comb %>%
  filter(FDR < 0.05) %>%
  count(feature) %>%
  filter(n == 3) %>%          # significant in all 3 survival types
  pull(feature)

# use manually selected features

# plotting
plot.df <- df.cox.out.comb %>%
  filter(feature %in% common_features[1:100]) %>%
  mutate(
    lower = exp(log(HR) - 1.96*sqrt(1/cox_n_sam)),
    upper = exp(log(HR) + 1.96*sqrt(1/cox_n_sam)))

ggplot(plot.df, aes(HR,
           reorder(feature, HR),
           colour = surv)) +
  
  geom_vline(xintercept = 1,
             linetype = 2,
             colour = "grey50") +
  
  geom_point(position = position_dodge(width = 0.6),
             size = 3) +
  
  geom_errorbarh(aes(xmin = lower,
                     xmax = upper),
                 position = position_dodge(width = 0.6),
                 height = 0.2) +
  
  scale_x_log10() +
  
  labs(x = "Hazard ratio",
       y = "",
       colour = "Outcome") +
  
  theme_bw()

################################
## correlations of cna disjoints
################################
library(Matrix)
library(pheatmap)

#get features
features <- read.csv('cox_all_features.csv') %>% filter(n_samples > 15)
table(features$group_name)
# load CNA mat
CNA.mat.gain <- readRDS(paste0('CNA_mat_Broad_WTscore_3Gain.Rds'))
CNA.mat.loss <- readRDS(paste0('CNA_mat_Broad_WTscore_3Loss.Rds'))  
# keep sig gain
sig.feat.gain <- intersect(features$feature, colnames(CNA.mat.gain))
CNA.mat.gain <- as.matrix(CNA.mat.gain[, sig.feat.gain])
#colnames(CNA.mat.gain) <- paste0(colnames(CNA.mat.gain))
# losses
sig.feat.loss <- intersect(features$feature, colnames(CNA.mat.loss))
CNA.mat.loss <- as.matrix(CNA.mat.loss[, sig.feat.loss])
#colnames(CNA.mat.loss) <- paste0(colnames(CNA.mat.loss))
#comb
all.samples <- union(rownames(CNA.mat.loss), rownames(CNA.mat.gain))

# Empty matrices
gain.full <- Matrix(0, nrow = length(all.samples),
                    ncol = ncol(CNA.mat.gain),
                    sparse = TRUE,
                    dimnames = list(all.samples, colnames(CNA.mat.gain)))

loss.full <- Matrix(0,  nrow = length(all.samples),
                    ncol = ncol(CNA.mat.loss),
                    sparse = TRUE,
                    dimnames = list(all.samples, colnames(CNA.mat.loss)))

# Fill in existing values
gain.full[rownames(CNA.mat.gain), ] <- CNA.mat.gain
loss.full[rownames(CNA.mat.loss), ] <- CNA.mat.loss
sum(duplicated(c(colnames(gain.full), colnames(loss.full))))

# Combine
cna.mat <- cbind(gain.full, loss.full)
corr.mat <- cor(as.matrix(cna.mat))

surv <- 'OS'
# Plot heatmap
#features <- features %>% filter(group_name == 'surv')
features$feature2 <- features$feature
features$feature2 <- gsub('_',':',features$feature2)
chr.spt <- features %>% tidyr::separate(feature2, into = c("chr", "coords"), sep = ":") %>% 
  tidyr::separate(coords, into = c("start", "end"), sep = "-") 
features$chr <- chr.spt$chr
features$start <- chr.spt$start
features$start <- gsub('Loss|Gain','1',features$start)
# order chr
chr.levels <- paste0("chr", c(1:22, "X", "Y"))
features$chr <- factor(features$chr, levels = chr.levels)
features <- features[order(features$chr, features$start),]
ord <- features$feature
idx <- match(ord, rownames(corr.mat))
corr.mat <- corr.mat[idx, idx]

ann <- data.frame(CNV = features$cnv_type, row.names = features$feature)
ann <- data.frame(CNV = features$cnv_type, row.names = features$feature)
ann <- ann[rownames(corr.mat), , drop = FALSE]
ann.colors <- list(CNV = c(Gain = cnv_cols[[1]], Loss = cnv_cols[[2]]))

#library(RColorBrewer)
my_colors <- colorRampPalette((brewer.pal(n = 7, name = "OrRd")))(50)

pheatmap(corr.mat, clustering_method = "complete",
         border_color = NA,  fontsize = 7,  cluster_cols = F, cluster_rows = F,
        # annotation_row = ann, annotation_col = ann,
         annotation_colors = ann.colors, 
         main = paste("Correlation of Significant CNA Features-",surv))


corr.df <- reshape2::melt(corr.mat)

corr.df <- corr.df |>
  dplyr::rename(feature1 = Var1,
                feature2 = Var2,
                cor = value) |>
  dplyr::filter(feature1 != feature2) |>      # remove self-correlations
  dplyr::filter(abs(cor) >= 0.9)              # choose your threshold

# summarize corr
hist(corr.df$cor)
corr.summary <- bind_rows(corr.df %>% 
                            transmute(feature = feature1,
                partner = feature2,
                cor), corr.df %>%
      transmute(feature = feature2,
                partner = feature1,
                cor)) %>% group_by(feature) %>% summarise(
    max_cor = max(abs(cor)), n_correlated = dplyr::n(),
    correlated_features = paste(partner, collapse = "; "),
    .groups = "drop")

# merge to feature table
cox.sig <- left_join(df.cox.OS.sig, corr.summary, by = c("feature"))

cox.sig$keep <- TRUE

for(i in seq_len(nrow(corr.df))){
  f1 <- corr.df$feature1[i]
  f2 <- corr.df$feature2[i]
  
  i1 <- match(f1, cox.sig$feature)
  i2 <- match(f2, cox.sig$feature)
  
  if(cox.sig$FDR[i1] > cox.sig$FDR[i2])
    cox.sig$keep[i1] <- FALSE
  else
    cox.sig$keep[i2] <- FALSE
}

View(cox.sig)
table(cox.sig$keep)
cox.sig$feature[cox.sig$keep]
#save
write.csv(cox.sig, paste0('cox_all_features_with_corr_',surv,'.csv'), quote = F, row.names = F)

#############################
## volcano plot
############################
#load cox results for OS
df.cox.loss <- readRDS(paste0('CLL_filtered_disj_cox_Broad_WTscore_3Loss_',surv,'.Rds')) %>% mutate(cnv_type='Loss')
df.cox.gain <- readRDS(paste0('CLL_filtered_disj_cox_Broad_WTscore_3Gain_',surv,'.Rds')) %>% mutate(cnv_type='Gain')
#comb
df.cox.OS.sig <- rbind(df.cox.gain, df.cox.loss) %>% filter(cox_n_sam > 10) #filter(FDR <= 0.05 & n_samples > 10)
annot <- read.csv(paste0('cox_all_features_',surv,'.csv'))

# volcano plot
df.cox.OS.sig %>% ggplot(aes(HR, -log10(FDR))) + 
  geom_point(aes(size = n_samples, fill=cnv_type),alpha = 0.4, shape = 21, col='black') + 
  geom_label_repel(data = annot[1:80,], aes(label = cytoband) , size = 2) +
  geom_hline(yintercept = -log10(0.05), ,colour = 'black', linetype = 2, linewidth = 0.4)+ 
  geom_vline(xintercept = 1, colour = 'black', linetype = 2, linewidth = 0.4)+
  scale_fill_manual(values = cnv_cols)+ theme_thesis()

###################################
## multivariable cox and KP
###################################
surv <- 'OS'
source('../concordance_ch3/functions.R')
source('../concordance_ch3/functions_new.R')

features <- read.csv(paste0('cox_all_features_with_corr_',surv,'.csv'))
#features <- read.csv('cox_all_features_for_plt.csv')
features <- features[features$keep,]
range(features$n_samples)
# add more hits
in.man <- data.frame(feature=c('chr6:155647047-157404463','chr3:47448767-47488255'),
                     n_samples=c(31,27),
                     cytoband=c('chr6q25.3 (ARID1B)', 'chr3p21.31 (SCAP)'))

features <- rbind(features[c(1,6,9)], in.man)
#View(features)

clinical <- read.csv('../survival_mdr/All CLL outcome_spss.csv')
clinical$Study <- gsub('^Hull.*','Hull',clinical$Study)
clinical$Study <- gsub('^Newcastle.*','Newcastle',clinical$Study)
clinical$Study <- gsub('^Oxford.*','Oxford',clinical$Study)

table(clinical$Study)

# combine CNA cox mat#../survival_mdr/All CLL outcome_spss.csv combine CNA cox mat
CNA.mat.gain <- readRDS(paste0('CNA_mat_Broad_WTscore_3Gain_',surv,'.Rds'))
# keep sig gain
sig.feat.gain <- intersect(features$feature, colnames(CNA.mat.gain))
CNA.mat.gain <- as.matrix(CNA.mat.gain[, sig.feat.gain])
colnames(CNA.mat.gain) <- paste0(colnames(CNA.mat.gain))
# losses
CNA.mat.loss <- readRDS(paste0('CNA_mat_Broad_WTscore_3Loss_',surv,'.Rds'))
sig.feat.loss <- intersect(features$feature, colnames(CNA.mat.loss))
CNA.mat.loss <- as.matrix(CNA.mat.loss[, sig.feat.loss])
colnames(CNA.mat.loss) <- paste0(colnames(CNA.mat.loss))
#comb
all.samples <- union(rownames(CNA.mat.loss), rownames(CNA.mat.gain))
# Empty matrices
gain.full <- Matrix(0, nrow = length(all.samples),
                    ncol = ncol(CNA.mat.gain),
                    sparse = TRUE,
                    dimnames = list(all.samples, colnames(CNA.mat.gain)))

loss.full <- Matrix(0,  nrow = length(all.samples),
                    ncol = ncol(CNA.mat.loss),
                    sparse = TRUE,
                    dimnames = list(all.samples, colnames(CNA.mat.loss)))

# Fill in existing values
gain.full[rownames(CNA.mat.gain), ] <- CNA.mat.gain
loss.full[rownames(CNA.mat.loss), ] <- CNA.mat.loss

# Combine
cna.mat <- cbind(gain.full, loss.full)
cnv.df.feat <- as.data.frame(as.matrix(cna.mat))

# add to clinical feat
setdiff(clinical$Sample, rownames(cnv.df.feat))
cnv_clin <- cnv.df.feat %>%  tibble::rownames_to_column("Sample") %>%
  left_join(clinical, by = "Sample")

length(cnv_clin$Sample)
names(cnv_clin)

source('functions_feature_selection.R')
rm(dat)
# plot in a loop
pdf(paste0("surv_plots_auto/",surv,"_KM.pdf"), width = 6.5, height = 5.5, onefile = T)
for (i in 1:nrow(features)) {
  feat <- features$feature[i]
  #tit <- bquote( "Deletion chr10q24.32" ~ italic(.("(SFR1, CFAP43)")))
  tit <-  paste(features$cnv_type[i],'at',features$cytoband[i])# SFR1;CFAP43 chr10q24.32
  tit <- sub("Gain_"," whole", tit)
  # use KM function
  plt <- plot_km_feature(df = cnv_clin, title = tit, time = "OSdays", #OS_PTstatus
                         status = "OSstatus",
                         feature = feat, show_hr = T, conf.int = F)
  # save plot
  print(plt)
  }
dev.off()

#########################
## multivariable cox
#########################
time_var <- 'OSdays/365.25' #OS_PTdays
event_var <- 'OSstatus' #OS_PTstatus
feat <- "`chr3:128080357-128741641`"

rhs <- paste(c(feat, "Age","Sex", "Binet", "VH"), collapse = " + ")
rhs <- paste(c(feat, "Study"), collapse = " + ")

# cox adjusted
cox.formula <- as.formula(paste0("Surv(",time_var,",",event_var,") ~ ",rhs))
cox.fit <- coxph(cox.formula, data = cnv_clin)
summary(cox.fit)




#### removed
for (i in names(df.cox.in)) {
  df.cox <- df.cox.in[[i]]
  gr.cox.out <- df.cox[df.cox$type == "disjoint",] %>% filter(FDR < 0.05) %>%
    mutate(feature2=feature) %>% tidyr::separate(feature2, into = c("chr", "coords"), sep = ":") %>%
    tidyr::separate(coords, into = c("start", "end"), sep = "-") %>%
    mutate(start = as.numeric(start),
           end = as.numeric(end)) %>% 
    makeGRangesFromDataFrame(keep.extra.columns = TRUE,
                             seqnames.field = "chr", start.field = "start", end.field = "end")
  length(gr.cox.out)
  range(gr.cox.out$n_samples)
  
  gr.cox.features[[i]] <- reduce_cox_hits(gr = gr.cox.out, priority = 'FDR', max_gap = 0)
}


# combine features into one granges
# flatten gr list
df.cox.OS.features <- GenomicRanges::as.data.frame(gr.cox.features)

df.cox.OS.feature <- unlist(GRangesList(gr.cox.features), use.names = F)
length(df.cox.OS.feature)

View(data.frame(gr.cox.features))

# create all df.list and save
whole.cox.combined <- bind_rows(lapply(df.cox.in, filter, type == "whole_chr"))
whole.cox.combined$FDR <- whole.cox.combined$p # replace FDR
cox.all.df <- rbind(data.frame(df.cox.OS.feature)[6:15], whole.cox.combined)

write.csv(cox.all.df, paste0('cox_all_features_',surv,'.csv'), quote = F, row.names = F)





