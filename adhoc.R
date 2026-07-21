library(survival)
library(survminer)
library(stringr)
library(dplyr)
library(plyranges)
library(GenomicRanges)
library(regioneR)
library(parallel)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(karyoploteR)
library(parallel)

setwd('/home/dean91/cnv_cll/survival_mdr/')
# Loading data and params ####
type='Loss'
# for losses numSNP = 15, for gains numSNP was 50
# new snp cutoff is 20
snp_ct <- 5
gr.cll.cnv <- read.csv('All_cnvs_with_ox.csv') %>% filter(numSNP > snp_ct) %>%  
  filter(!chr %in% c('X','Y','23')) %>% filter(CNV_type==type) %>% toGRanges() 

# some stats of raw cnv calls
length(unique(gr.cll.cnv$new_sam_id))
plot(density(gr.cll.cnv$length))

gr.cll.cnv$sm_algo <- paste0(gr.cll.cnv$new_sam_id,'_',gr.cll.cnv$method)
test <- as.data.frame(gr.cll.cnv[,c(4,11)])
snp_avf <- function()
max.snp <- test %>% aggregate(numSNP~sm_algo,FUN=max)
max.snp %>% filter(numSNP<100)

