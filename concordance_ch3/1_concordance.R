# This is the immediate analysis after calling CNVs
# This is the major part of the chapter 3
# For each sample - take CNVs and look for concordance across algos
# resulting concordant CNV disjoints (> 3) can be used in cox_model.R script for recurrence
# workflow. rawCNVs > split gains/del > split algo > merge gaps > cnv QC based on size and probes > concordance
library(dplyr)
library(GenomicRanges)
library(parallel)
library(karyoploteR)
library(data.table)

source('functions.R')
source('functions_new.R')

# load all raw CNVs regardless of length or prob counts
gr.cll.cnv <- load.raw_gap(snp.ct = 5, len = 5, use.merged.cnv = F)

names(mcols(gr.cll.cnv))
range(width(gr.cll.cnv))
range(gr.cll.cnv$no_of_probes)

# crude kp plot
kp <- plotKaryotype(plot.type = 1, chromosomes = 'autosomal') #chromosomes = 'chr6'
kpPlotCoverage(kp, data=gr.cll.cnv)

######################################################
### Merge adjacent CNVs - don't use this
# first merge adjacent CNVs in same sample if 
# they have gaps less than 20% of total comb length
######################################################
gap.frac <- 0.2
gap.file <- paste0('Gap_',gap.frac,'_merged.rds')

#### new DT version
gap.merged.cnv <- merge_by_gap_fraction_dt(gr.cll.cnv, max_gap_fraction = gap.frac)

range(width(gap.merged.cnv))
range(gap.merged.cnv$no_of_probes)

# save gap merged
saveRDS(gap.merged.cnv, gap.file) 

###### old ver ###
merge_adjacent=T

gap.frac <- 0.2

if (isTRUE(merge_adjacent)) {
  
gap.file <- paste0('Gap_',gap.frac,'_merged.rds')
harmonized.file <- paste0('harmonized_disj.rds')
if (file.exists(gap.file)) {
  gap.merged.cnv <- readRDS(gap.file)
} else {
  #merge adjacent for each algo
  algo.split <- split(gr.cll.cnv, paste(gr.cll.cnv$caller, gr.cll.cnv$loss_gain, sep = "__"))
  names(algo.split)

  # do a loop for each algo as mcapply used in the function
  gap.merged.cnv <- GRangesList()
  for (i in seq_along(algo.split)) {
    print(names(algo.split[i]))
    gap.merged.cnv[[i]] <- merge_by_gap_fraction(algo.split[[i]], max_gap_fraction = gap.frac, ncores = 8)
    names(gap.merged.cnv[i]) <- names(algo.split[i])
  }
  # save gap merged
  saveRDS(gap.merged.cnv, gap.file) 
  # remove harmonized file if this step is redone
  file.remove(harmonized.file)
} 

# combine all merged
gap.merged.cnv <- unlist(GRangesList(gap.merged.cnv), use.names = FALSE)
range(gr.cll.cnv$numSNP)
range(gap.merged.cnv$numSNP)
hist(width(gr.cll.cnv))
hist(width(gap.merged.cnv))
} else {
  gap.merged.cnv <- gr.cll.cnv
  harmonized.file <- paste0('harmonized_disj.rds')
}
# plot
cowplot::plot_grid(length.plt(gr.cll.cnv, bin.size = 100), length.plt(gap.merged.cnv, bin.size = 100), ncol =1)
#view(data.frame(gap.merged.cnv))

# crude kp plot
kp <- plotKaryotype(plot.type = 1, chromosomes = 'chr6') #chromosomes = 'chr6'
kpPlotCoverage(kp, data=gr.cll.cnv, col='blue')
kpPlotCoverage(kp, data=gap.merged.cnv, col='red')

#####################
### Second QC on CNVs
#####################
# QC on numSNPs and length
# remove problematic regions such as centromeres, telemores - use UCSC track
# remove CNVs overlapping segdups
# use confidence score?

# algo wise stats

############################################################
# WITHIN-SAMPLE HARMONIZATION
############################################################
# Goal:
# Use disjoints for concordance as well
############################################################
source('functions.R')
# Split by sample + state
#gr_split <- split(gap.merged.cnv, paste0(gap.merged.cnv$new_sam_id, "_", gap.merged.cnv$CNV_type))
if (!file.exists(harmonized.file)) {
  ## harmonize sutdy names and sample ids
  gap.merged.cnv$study <- gsub('^Bourn.*','Bournemouth', gap.merged.cnv$study)
  gap.merged.cnv$study <- gsub('^Hull.*','Hull', gap.merged.cnv$study)
  gap.merged.cnv$study <- gsub('^Newcastle.*','Newcastle', gap.merged.cnv$study)  
  gap.merged.cnv$study <- gsub('^Oxford.*','Oxford', gap.merged.cnv$study) 
  gap.merged.cnv$study <- gsub('^South.*','Southampton', gap.merged.cnv$study) 
  gap.merged.cnv$new_sam_id <- gsub('1536_HRH-HOE', '1536_HRH_HOE', gap.merged.cnv$new_sam_id)
  
gr_split <- split(gap.merged.cnv, gap.merged.cnv$new_sam_id)
table(gap.merged.cnv$study)

harmonized_list <- mclapply(gr_split, harmonize_sample, mc.cores = 10)

harmonized_gr <- unlist(GRangesList(harmonized_list))
mcols(harmonized_gr)
table(harmonized_gr$n_callers)
# save gap merged
saveRDS(harmonized_gr, harmonized.file)
} else {
  harmonized_gr <- readRDS(harmonized.file)
}

# plot with karyoplotR
kp <- plotKaryotype(plot.type = 1, chromosomes = 'chr1') #chromosomes = 'chr6'
kpPlotCoverage(kp, data=gap.merged.cnv, show.0.cov = F, col = 'red')
kpPlotCoverage(kp, data=harmonized_gr, show.0.cov = F)

#harmonized_dat <- data.frame(harmonized_gr)

#########################################################
# goto 2_recurrence.R after this point
#########################################################
## Concordance only analysis for chap 3
# create reduced ranges and check with kp
# for each reduced range check concordance and report?
# maybe remove large alterations (arm or whole chr)
#########################################################
####################################################
## Load for raw CNVs
####################################################
# load all raw CNVs regardless of length or prob counts
gr.cll.cnv <- read.csv('../survival_mdr/All_cnvs_with_ox.csv') %>%  #& !study %in% c('Oxford','Oxford-ARC','Oxford-ADM')
  filter(!chr %in% c('X','Y','23')) %>% filter(CNV_type==type & numSNP > 10) %>% toGRanges()  
seqlevelsStyle(gr.cll.cnv) <- 'UCSC'
if (length(gr.cll.cnv) < 1) {stop("No CNAs left after algorithm filtering, check algo argument..")}
# filter out larger alt
#gr.cll.cnv <- gr.cll.cnv %>% filter(width(gr.cll.cnv) < 5e6)
range(width(gr.cll.cnv))

# check reduced size before merging adjacent
table(gr.cll.cnv$method)
gr_split.raw <- split(gr.cll.cnv, gr.cll.cnv$new_sam_id)

#check stats
sam.len <- data.frame(max(width(gr_split.raw)))
idx <- which(names(gr_split.raw)=='B9039')
kp <- plotKaryotype(plot.type = 2, chromosomes = 'chr9') #chromosomes = 'chr6'
kpPlotRegions(kp, data=gr_split.raw[[idx]], col='red', data.panel = 1)

# now reduce within sam and plot
red_raw.sam <- reduce(gr_split.raw)
red_all.sam <- unlist(red_raw.sam, use.names = F)

kp <- plotKaryotype(plot.type = 2, chromosomes = 'chr4') #chromosomes = 'chr6'
kpPlotCoverage(kp, data=red_all.sam, col='red', data.panel = 1, show.0.cov = F) # sam reduced
kpPlotCoverage(kp, data=harmonized_gr[harmonized_gr$n_callers>=3], col='blue', data.panel = 1, show.0.cov = F) # sam disj

############################
## check gap merged cnvs with KP plt
############################
type <- 'Loss'
gap.frac <- 0.2

# load gap merged cnv data
gap.file <- paste0('Gap_', gap.frac,'_merged_', type,'.rds')
gap.merged.cnv <- readRDS(gap.file)
gap.merged.unlist <- unlist(GRangesList(gap.merged.cnv), use.names = FALSE)

#check stats
range(width(gap.merged.cnv))

# split samples
gap_split.sam <- split(gap.merged.unlist, gap.merged.unlist$new_sam_id)
barplot(min(width(gap_split.sam)))

#check stats
gap.sam.len <- data.frame(max(width(gap_split.sam)))
gap.idx <- which(names(gap_split.sam)=='B9039')
#kp <- plotKaryotype(plot.type = 1, chromosomes = 'autosomal') #chromosomes = 'chr6'
kpPlotRegions(kp, data=gap_split.sam[[gap.idx]], col='blue', data.panel = 2)


#check reduced size
red.sam <- reduce(gr_split.raw)
sum(elementNROWS(red.sam))
red.sam.unlist <- unlist(GRangesList(red.sam), use.names = FALSE)
red.all <- reduce(gr.cll.cnv)
length(red.all)

# KP plots
# for all
kp <- plotKaryotype(plot.type = 2, chromosomes = 'chr3') #chromosomes = 'chr6'
kpPlotCoverage(kp, data=red.all, col='pink', data.panel = 1, show.0.cov = T)
kpPlotCoverage(kp, data=red.sam.unlist, col='red', data.panel = 1, show.0.cov = T)

#for sample wise
kp <- plotKaryotype(plot.type = 2, chromosomes = 'chr17') #chromosomes = 'chr6'
kpPlotRegions(kp, data=gr_split.raw[[1]], col='pink', data.panel = 1)
kpPlotRegions(kp, data=red.sam[[1]], col='blue', data.panel = 2)
kpPlotRegions(kp, data=red.sam.unlist, show.0.cov = T, col = 'red')


