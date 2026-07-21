## New 2024 Oct ##
library(survival)
library(survminer)
library(dplyr)
library(karyoploteR)
library(GenomeInfoDb)

setwd('~/cnv_cll/survival_mdr/')
#source('functions.R')

type <- 'Loss'
algo <- 'Nexus'
snp_ct <- 20
len <- 0.1e6

# file name
save.file <- paste0('cna_mat_',algo,'_',type, snp_ct,'snps_',len/1e6,'Mb.rds')
if (!file.exists(paste0(getwd(),'/cox_mat/',save.file))) {
  
# use length insted of probes - > 500kb
gr.cll.cnv <- read.csv('All_cnvs_with_ox.csv') %>% 
  filter(length > len & numSNP > snp_ct) %>%  #& !study %in% c('Oxford','Oxford-ARC','Oxford-ADM')
  filter(!chr %in% c('X','Y','23')) %>% filter(CNV_type==type) %>% toGRanges() 
if (!algo =='all') {
  gr.cll.cnv <- gr.cll.cnv %>% filter(method==algo)
}
if (length(gr.cll.cnv) < 1) {stop("No CNAs left after algorithm filtering, check algo argument..")}

seqlevelsStyle(gr.cll.cnv) <- 'UCSC'
head(gr.cll.cnv)
range(gr.cll.cnv$length)
range(gr.cll.cnv$numSNP)
table(gr.cll.cnv$method)

# some stats of raw cnv calls
length(unique(gr.cll.cnv$new_sam_id))
plot(density(gr.cll.cnv$numSNP))

#load outcome data
cll.outcome <- read.csv('All CLL outcome_spss.csv') %>%  filter(!Study %in% c('Oxford-ARC', 'Oxford-ADM'))  #to remove oxford samples 
length(unique(cll.outcome$Sample))
#from the surv remove them from the outcome table
table(cll.outcome$Study)
setdiff(gr.cll.cnv$new_sam_id, cll.outcome$Sample)

# start the pipeline
count_samples.algo <- function(algo, cutoff){
  if (algo=='all') {
    gr.cnv.disjoints <- disjoin(gr.cll.cnv)
    count <- count_overlaps(gr.cnv.disjoints, gr.cll.cnv)
  } else {
  # return disjs as granges list
    gr.cnv.disjoints <- disjoin(gr.cll.cnv[gr.cll.cnv$method==algo])
  # create out columns from the disjoints
  #out_table$seqnames <- as.numeric(out_table$seqnames)
    gr_cnv_algo <- gr.cll.cnv %>% filter(method==algo) 
  #calculate overlaps
    count <- count_overlaps(gr.cnv.disjoints, gr_cnv_algo)
  }
  # add ovlps to disjoints object
  # add disjoint id
  gr.cnv.disjoints$count <- count
  # filter cutoff
  gr.cnv.disjoints <- filter(gr.cnv.disjoints, count > cutoff)

  return(gr.cnv.disjoints)
}

count_samples <- function(cutoff){
  gr.cnv.disjoints <- disjoin(gr.cll.cnv)
  count <- count_overlaps(gr.cnv.disjoints, gr.cll.cnv)
  # add ovlps to disjoints object
  # add disjoint id
  gr.cnv.disjoints$count <- count
  # filter cutoff
  gr.cnv.disjoints <- filter(gr.cnv.disjoints, count > cutoff)
  return(gr.cnv.disjoints)
}

ranges.out <- count_samples(10)
length(ranges.out)
# kp plot
kp <- plotKaryotype(plot.type = 1, chromosomes = 'autosomal') #chromosomes = 'chr6'
kpPlotCoverage(kp, data=ranges.out, show.0.cov = T, col=ifelse(type=='Loss','red','blue'))


ranges.out$disj_id <- seq_along(ranges.out)
head(ranges.out)

##################################
# create CNA matrix for algo/all #
##################################
# create CNA matrix - goal is to add all ranges with 1/0 for alt
cna.mat <- mclapply(seq_along(ranges.out), function(i) {
  #apply per range
  outcome <- cll.outcome[1]
  ovlps.nex <- find_overlaps(gr.cll.cnv, range(ranges.out[i])) #filter_by_overlaps
  dat.ovlps <- data.frame(Sample=unique(ovlps.nex$new_sam_id), alteration=1) #for nexus
  outcome <- left_join(outcome, dat.ovlps,by='Sample')
  outcome$alteration[is.na(outcome$alteration)] <- 0  
  setnames(outcome, 'alteration', paste0(as.character(ranges.out[i])))
  }, mc.cores = 10)

#set kay for fast join
cna.mat <- lapply(cna.mat, as.data.table)
lapply(cna.mat, setkey, Sample)

cna.mat <- Reduce(function(x, y) cbind(x, y[, -1, with = FALSE]), cna.mat)
#merge surv data
names(cll.outcome)
cna.mat <- merge(cll.outcome[c(1,3:13)], cna.mat, by='Sample', all=T)
head(cna.mat)

# save cna.mat as rds
saveRDS(cna.mat, paste0('cox_mat/', save.file))
} else {
  cna.mat <- readRDS(paste0('cox_mat/', save.file))
}
head(cna.mat)

# below goes for chapter 2- main methods
##################################
#### cox function from mat ###
###################################
plot_surv <- function(surv_type, cutoff) {
# Define survival object
cna.mat$TreatmentCenter <- factor(cna.mat$TreatmentCenter)
cna.mat$Array <- factor(cna.mat$Array)
cna.mat$Study <- factor(cna.mat$Study)
cna.mat$Sex <- factor(cna.mat$Sex)

cna.mat.out <- cna.mat[, grep("chr", names(cna.mat))]
head(names(cna.mat))
feature_names <- colnames(cna.mat.out)
counts.vec <- sapply(feature_names, function(f) sum(cna.mat[[f]] == 1, na.rm=TRUE))
valid_features <- feature_names[
  sapply(feature_names, function(f) sum(cna.mat[[f]] == 1, na.rm=TRUE) >= cutoff)]
feature_names <- valid_features
# Apply univariate Cox model to each CNA
# univ_results <- lapply(feature_names, function(feature) {
#   formula <- as.formula(paste("Surv(OSdays, OSstatus) ~", feature, "+ TreatmentCenter"))
#   #model <- coxph(surv_obj ~ feature + TreatmentCenter + Array + Sex, data = cna.mat) #add centre as covariate? what else?
#   model <- coxph(formula, data = cna.mat)
#   coef_info <- summary(model)$coefficients[1, c("coef", "Pr(>|z|)")]
#   return(coef_info)
# })
if (surv_type=='OS'){
univ_results <- mclapply(feature_names, function(feat) {
  x <- cna.mat[[feat]]
  model <- coxph(Surv(OSdays, OSstatus) ~ x + TreatmentCenter + Study + Array, data = cna.mat)
  summary(model)$coefficients[1, c("exp(coef)", "Pr(>|z|)")]
}, mc.cores = 8)
} else if (surv_type=='TTFT') {
  univ_results <- mclapply(feature_names, function(feat) {
    x <- cna.mat[[feat]]
    model <- coxph(Surv(TTFTdays, TTFTstatus) ~ x + TreatmentCenter + Array, data = cna.mat)
    summary(model)$coefficients[1, c("exp(coef)", "Pr(>|z|)")]
  }, mc.cores = 8)
} else {
  univ_results <- mclapply(feature_names, function(feat) {
    x <- cna.mat[[feat]]
    model <- coxph(Surv(OS_PTdays, OS_PTstatus) ~ x + TreatmentCenter + Array, data = cna.mat)
    summary(model)$coefficients[1, c("exp(coef)", "Pr(>|z|)")]
  }, mc.cores = 8)
}

# Convert to data frame
univ_results <- do.call(rbind, univ_results)
rownames(univ_results) <- feature_names
univ_results <- as.data.frame(univ_results)
univ_results$FDR <- p.adjust(univ_results$`Pr(>|z|)`, method = "BH")

##################
## plot mahnatton
##################
cox_pval <- GRanges(rownames(univ_results))
seqlevelsStyle(cox_pval) <- 'UCSC'
cox_pval$FDR <- univ_results$FDR

#add chr col
chr.col <- as.numeric(gsub('chr','', seqnames(cox_pval)))
cox_pval$chr.col <- chr.col%%2
table(cox_pval$chr.col)
cox_pval$chr.col <- gsub('^1','#FFBD07AA', cox_pval$chr.col)
cox_pval$chr.col <- gsub('^0','#00A6EDAA', cox_pval$chr.col)
cox_pval$pos <- (end(cox_pval) + start(cox_pval))/2

range(-log10(cox_pval$FDR))

chrs <- paste0('chr',seq(1:22))
#png(paste(title,'.png'),width = 12,height = 7,units = 'in',res = 180)
kp <- plotKaryotype(plot.type = 4, chromosomes = chrs, labels.plotter = NULL,main = surv_type)
kpAddChromosomeNames(kp,srt=45,cex=0.6)
#kpAddBaseNumbers(kp)
#kpPlotDensity(kp, data=gr.nex, r0=0, r1=0.3, window.size = 1000, col="orchid")
kpPoints(kp, data=cox_pval, y = -log10(cox_pval$FDR), pch=21, col='black', cex = 1,
         lwd = 0.2, ymax = max(-log10(cox_pval$FDR),na.rm = T), bg=cox_pval$chr.col)

}
plot_surv('OS', 10)

kpText(kp, data = GRanges(test), labels = paste0(test$cyto,'\n',test$gene,'','\n','N=',test$nex), y = test$logp, 
       ymax = max(del_pvals$logp),pos=4, cex=0.55, col='black')

range(del_pvals$nex)/500
#plot(density(del_pvals.noox$os.nex.log))
kpAxis(kp,ymin = 0,ymax = max(del_pvals$logp), cex=0.8,tick.pos = c(0,1.30103,max(del_pvals$logp)))
kpAddLabels(kp,'-log10(FDR)',srt = 90,cex = 0.8, r0 = 0,side = 'left')

kpAbline(kp, h=-log10(0.05), lty=3, ymax=max(del_pvals$logp), ymin =min(del_pvals$logp),cex=1.5,
         col='red')
#dev.off()

###################
### KP plots ######
###################
cox_reg <- function(surv_type, cutoff) {
  # Define survival object
  cna.mat$TreatmentCenter <- factor(cna.mat$TreatmentCenter)
  cna.mat$Array <- factor(cna.mat$Array)
  cna.mat$Study <- factor(cna.mat$Study)
  cna.mat$Sex <- factor(cna.mat$Sex)
  
  cna.mat.out <- cna.mat[, grep("chr", names(cna.mat))]
  head(names(cna.mat))
  feature_names <- colnames(cna.mat.out)
  counts.vec <- sapply(feature_names, function(f) sum(cna.mat[[f]] == 1, na.rm=TRUE))
  valid_features <- feature_names[
    sapply(feature_names, function(f) sum(cna.mat[[f]] == 1, na.rm=TRUE) >= cutoff)]
  feature_names <- valid_features

  if (surv_type=='OS'){
    univ_results <- mclapply(feature_names, function(feat) {
      x <- cna.mat[[feat]]
      model <- coxph(Surv(OSdays, OSstatus) ~ x + TreatmentCenter + Study + Array, data = cna.mat)
      summary(model)$coefficients[1, c("exp(coef)", "Pr(>|z|)")]
    }, mc.cores = 8)
  } else if (surv_type=='TTFT') {
    univ_results <- mclapply(feature_names, function(feat) {
      x <- cna.mat[[feat]]
      model <- coxph(Surv(TTFTdays, TTFTstatus) ~ x + TreatmentCenter + Array, data = cna.mat)
      summary(model)$coefficients[1, c("exp(coef)", "Pr(>|z|)")]
    }, mc.cores = 8)
  } else {
    univ_results <- mclapply(feature_names, function(feat) {
      x <- cna.mat[[feat]]
      model <- coxph(Surv(OS_PTdays, OS_PTstatus) ~ x + TreatmentCenter + Array, data = cna.mat)
      summary(model)$coefficients[1, c("exp(coef)", "Pr(>|z|)")]
    }, mc.cores = 8)
  }
  
  # Convert to data frame
  univ_results <- do.call(rbind, univ_results)
  rownames(univ_results) <- feature_names
  univ_results <- as.data.frame(univ_results)
  univ_results$FDR <- p.adjust(univ_results$`Pr(>|z|)`, method = "BH") 
  return(univ_results)
  }

cox_out <- cox_reg('PT', 10)

idx <- grep('chr3:46336781-46697991', names(cna.mat))
x <- cna.mat[[idx]]
survfit(Surv(OSdays/365.25, OSstatus) ~x, data = cna.mat) %>% 
  ggsurvplot(conf.int=F, pval=TRUE, risk.table=TRUE, 
             xlab='Time in years',  legend=c(0.8,0.8),
            palette=c("#66CD00", "#B22222", "#4682B4"), 
             risk.table.y.text = FALSE)
##############################################################
################### Feature selection #######################
# first Univariate Cox
# then LASSO-Cox
# later Random Survival Forest
# Compare prediction performance using
# c=index, time-dependent AUC
##############################################################
#test with surv data and typed meatadata
#cna.mat <- read.csv('All CLL outcome_spss.csv')
#cna.mat <- as.matrix(cna.mat)
# start with lasso cox
library(glmnet)

#remove NA samples from cna.mat.nex
cna.mat.lasso <- cna.mat[!is.na(cna.mat$OSdays),]
#remove low/high recurrence
prevalence <- rowMeans(cna.mat.lasso[13:ncol(cna.mat.lasso)] != 0, na.rm = TRUE)
plot(prevalence)
cna.mat.lasso <- cna.mat.lasso[, prevalence >= 0.03 & prevalence <= 0.97]
# Prepare input
x <- as.matrix(cna.mat.lasso[, grep("chr", names(cna.mat.lasso))])
#x <- as.matrix(cna.mat.lasso[22:length(cna.mat)])
y <- Surv(cna.mat.lasso$OSdays, cna.mat.lasso$OSstatus)

# Fit LASSO
fit <- cv.glmnet(x, y, family = "cox", alpha = 1)

# Extract selected features
selected_features <- rownames(coef(fit, s = "lambda.min"))[which(coef(fit, s = "lambda.min") != 0)]
selected_features

#count total
fet.recurrance <- data.frame(colSums(cna.mat[22:ncol(cna.mat)]))
range(fet.recurrance)
#get selected features
fet.recurrance[rownames(fet.recurrance) %in% selected_features,]

# extract lasso selected
# Coefficients at the chosen lambda
coef_lasso <- coef(fit, s = "lambda.min")
# Alternatively:
# coef_lasso <- coef(cvfit, s = "lambda.min")

coef_lasso <- as.matrix(coef_lasso)

selected_features <- rownames(coef_lasso)[coef_lasso[, 1] != 0]
selected_beta <- coef_lasso[selected_features, 1]

selected_features
selected_beta
length(selected_features)

# Calculate each patient’s LASSO prognostic score
# Ensure columns are in exactly the same order as the coefficients
x_selected <- as.matrix(cna.mat[, selected_features, drop = FALSE])

# Weighted prognostic score
cna.mat$lasso_score <- as.numeric(x_selected %*% selected_beta)

summary(cna.mat$lasso_score)

#####################################################################
# Random survival forest
#####################################################################
# rows = patients
# columns = clinical + CNA features
# OSdays, OSstatus, age, sex, treatment, chr11q_del, chr13q_del, chr17p_del, ...
######################################################################
library(randomForestSRC)
library(survival)

set.seed(1)

names(cna.mat)

rsf_fit <- rfsrc(Surv(OSdays, OSstatus) ~ Age + Sex + TreatmentCenter  +
    chr11q_del + chr13q_del + chr17p_del + chr12_gain,
  data = cna.mat, ntree = 1000,
  importance = TRUE)

print(rsf_fit)

######################################
## check correlation
######################################
cna_correlation <- cor(selected_features,
  use = "pairwise.complete.obs",
  method = "spearman")

round(cna_correlation, 2)


##############################
###comb cox model
##############################
selected_features <- selected_features[-c(7,8,9)]
combined_formula <- reformulate(selected_features,
  response = "Surv(OSdays, OSstatus)")

combined_cna_model <- coxph(combined_formula,
  data = cna.mat,  x = TRUE)

summary(combined_cna_model)

model_summary <- summary(combined_cna_model)

combined_results <- data.frame(
  feature = rownames(model_summary$coefficients),
  beta = model_summary$coefficients[, "coef"],
  HR = model_summary$coefficients[, "exp(coef)"],
  lower_95_CI = model_summary$conf.int[, "lower .95"],
  upper_95_CI = model_summary$conf.int[, "upper .95"],
  p_value = model_summary$coefficients[, "Pr(>|z|)"],
  row.names = NULL)

combined_results

# forest
ggforest( combined_cna_model,
  data = cna.mat,
  main = "Multivariable Cox model of selected recurrent CNAs",
  cpositions = c(0.02, 0.25, 0.45),
  fontsize = 0.8,
  refLabel = "Reference",
  noDigits = 2)

###################################
## plot KM for selected features ##
###################################
# Create the survival object
surv_obj <- Surv(time = cna.mat$OSdays/365.25, event = cna.mat$OSstatus)
fit <- survfit(surv_obj ~ cna.mat$Del_10q24.32, data = cna.mat)

# Plot the KM curve
ggsurvplot(fit,
           data = cna.mat,
           pval = TRUE,
           conf.int = TRUE,
           risk.table = TRUE,
           legend.labs = c("No del", "del"),
           title = "CLL Survival",
           xlab = "Time (years)",
           ylab = "Survival Probability")

### check cox for age sex etc.
head(cll.outcome)
# Create the survival object
surv_obj <- Surv(time = cll.outcome$OSdays/365.25, event = cll.outcome$OSstatus)
fit <- survfit(surv_obj ~ cll.outcome$TreatmentCenter, data = cll.outcome)

# Plot the KM curve
ggsurvplot(fit,
           data = cll.outcome,
           pval = TRUE,
           conf.int = F,
           risk.table = TRUE,
           title = "CLL Survival",
           xlab = "Time (years)",
           ylab = "Survival Probability")

# cox for cll.outcome 
names(cll.outcome)
fml <- as.formula(paste("Surv(OSdays,OSstatus) ~",paste(names(cll.outcome)[23:30], collapse = ' + ')))

cox.old <- coxph(fml, data = cll.outcome)
summary(cox.old)
ggforest(cox.old)

### cox model for new cnas
cur_fet <- selected_features
fml <- as.formula(paste("Surv(OSdays,OSstatus) ~", paste0('`',cur_fet,'`', collapse = ' + ')))

cox_model <- coxph(fml, data = cna.mat)
summary(cox_model)
ggforest(cox_model)

#####################
## Idea 2- dividing groups and KM
#####################
# divide patients into low risk and high risk based on features 
# if any of 24 detected its high risk?


