library(survival)
library(survminer)
library(stringr)
library(dplyr)
library(plyranges)
library(GenomicRanges)
library(regioneR)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(karyoploteR)
library(parallel)
library(GenomeInfoDb)

# load functions
source('functions.R')

##################
###Process CNVRs #
##################
setwd('/home/dean91/cnv_cll/survival_mdr/')
# Loading data and params ###
type='Loss'
# for losses numSNP = 15, for gains numSNP was 50
# new snp cutoff is 20
snp_ct <- 200
gr.cll.cnv <- read.csv('All_cnvs_with_ox.csv') %>% filter(numSNP > snp_ct & !study %in% c('Oxford','Oxford-ARC','Oxford-ADM')) %>%  
  filter(!chr %in% c('X','Y','23')) %>% filter(CNV_type==type) %>% toGRanges() 
seqlevelsStyle(gr.cll.cnv) <- 'UCSC'

# some stats of raw cnv calls
length(unique(gr.cll.cnv$new_sam_id))
plot(density(gr.cll.cnv$length))

# convert rpt to normal from ONLY Leicester sample ids
#dat.cll.cnv <- data.frame(gr.cll.cnv)
#dat.cll.cnv$new_sam_id <- dat.cll.cnv$new_sam_id %>% gsub('_rpt','',.)

#load outcome data
cll.outcome <- read.csv('All CLL outcome_spss.csv') %>%  filter(!Study %in% c('Oxford-ARC', 'Oxford-ADM'))  #to remove oxford samples 
length(unique(cll.outcome$Sample))
#from the surv remove them from the outcome table

# start the pipeline
ranges.out<- count_samples()
# remove low recurrent disjs
ranges.out <- ranges.out %>% filter(pen+qsnp+ipn+nex > 10)
ranges.out <- ranges.out %>% filter(nex > 15) # use cutoff of nexus numbers instead
ranges.out$disj_id <- seq_along(ranges.out)

# load genes
genes.data <- makeGenesDataFromTxDb(txdb = TxDb.Hsapiens.UCSC.hg19.knownGene,
                                    plot.transcripts = F,plot.transcripts.structure = F)
genes.data <- addGeneNames(genes.data)

# Run lapply with 20 cores 
loss.out.all <- mclapply(seq_along(ranges.out), cal_survival, mc.cores = 10)
#loss.out <- mclapply(1:500, cal_survival,mc.cores = 20)
loss.out.all <- as.data.frame(do.call(rbind,loss.out.all))
loss.out.all <- loss.out.all[,c(1,2,3,4,5,6,7,8,10,9,11,14,12,15,13,16,17)]
# Saving the data
saveRDS(loss.out.all, paste0('cll_disj_',type,'_snp',snp_ct,'noOx.Rds'))

# Post-Hoc testing on selected disjs ##
# load dels and gains p-vals from data 
load_pvals <- readRDS('cll_disj_Loss_snp200noOx.Rds')
# get sample names for particular disj and plot surv
disj_id <- 5000
cal_survival(disj_id)
plot_survival(disj_id,'nex')
# what are the sample ids that are detected by others but not nex?
setdiff(outcomes$Sample[outcomes$alteration_all==1],outcomes$Sample[outcomes$alterati-on_nex==1])

###################
# Post analysis ###
###################

# get reduced table ###
reduced.os <- reduce_pval.list(load_pvals,'os')
reduced.pt <- reduce_pval.list(load_pvals,'pt')
reduced.ttfs <- reduce_pval.list(load_pvals,'ttfs')

# filter sig p-val
pvals.os <- load_pvals %>% filter(total > 15 & os.p.nex < 0.05) 
# remove chr12 due to trisomy 12

# get sig genes
# separate genes into rows
base.dat_filtered_curated <- sig_p.ranges %>% separate_rows(genes,sep = ",")
base.dat_filtered_curated$gene_2 <- paste0(base.dat_filtered_curated$seqnames,'_',
                                           base.dat_filtered_curated$genes)

sig_genes <- data.frame(table(base.dat_filtered_curated$genes)) # get a median p val
# add gene location to sig genes
sig_genes <- left_join(x = sig_genes,y = data.frame(genes.data$genes),join_by(Var1==name))
# collapse sig_ranges and get 
gr.collapsed_sig <- reduce(toGRanges(base.dat_filtered))

# post analysis losses - TTFT
del_pvals.ttfs <- del_pvals %>% filter(total > 12 & ttfs.p.nex < 0.05) 


# filter 
loss.reduced.os2 <- loss.reduced.os %>% filter(os.p.nex < 0.05 & total > 15)
loss.reduced.ttfs2 <- loss.reduced.ttfs %>% filter(ttfs.p.nex < 0.05 & total > 15)
loss.reduced.pt2 <- loss.reduced.pt %>% filter(pt.p.nex < 0.05 & total > 15)

loss.out.all$os.p.nex_minus <- -(loss.out.all$os.p.nex)

loss.out.all %>% ggplot(aes(disj_id,os.p.nex_minus)) + geom_point(size=0.5) + scale_y_log10()
  scale_y_continuous(trans='log10')

range(loss.reduced.ttfs2$ttfs.p.nex)

####################
## P-val Plotting###
####################
setwd('~/cnv_cll/survival_mdr/')
loss.out.all <- read.csv('cll_loss20noOX_pvals.csv')
loss.out.all$adj.p <- p.adjust(loss.out.all$os.p.nex,method = 'fdr')

range(loss.out.all$ttfs.p.nex,na.rm = T)
range(loss.out.all$adj.p,na.rm = T)

# p val plot
del_pvals <- loss.out.all %>% filter(nex > 15 & adj.p > 1e-15 )
# convert p val to -log10
del_pvals$adj.p[is.na(del_pvals$adj.p)] <- 1
del_pvals$log.adj.p <- -log10(del_pvals$adj.p)
range(del_pvals$adj.p)
range(del_pvals$log.adj.p)
#del_pvals %>% filter(os.nex.log > 177)
#add chr col
chr.col <- as.numeric(gsub('chr','',del_pvals$seqnames))
del_pvals$chr.col <- chr.col%%2

table(del_pvals$chr.col)
# convert to Granges
del_pvals$chr.col <- gsub('^1','#FFBD07AA',del_pvals$chr.col)
del_pvals$chr.col <- gsub('^0','#00A6EDAA',del_pvals$chr.col)
del_pvals$pos <- (del_pvals$end+del_pvals$start)/2
del_pvals <- GRanges(seqnames = del_pvals$seqnames,ranges = IRanges(del_pvals$pos,width = 1),nex=del_pvals$nex, 
                     logp=del_pvals$log.adj.p, colour=del_pvals$chr.col, gene=del_pvals$genes, disj=del_pvals$disj_id)
seqlevelsStyle(del_pvals) <- 'UCSC'
range(del_pvals$logp)

# get top hits by reducing
#top.hits <- reduce_pval.list(loss.out.all,'os')
top.hits <- as.data.frame(del_pvals) %>% filter(logp > 1.3) %>% arrange(-logp)
top.hits <- find_overlaps(GRanges(top.hits),kp$cytobands)
top.hits$cyto <- paste0(seqnames(top.hits),'.',top.hits$name)
top.hits <- data.frame(top.hits)
cyto.lst <- unique(top.hits$cyto)
gene.lst <- unique(top.hits$gene)

top.genes <- function(x){
  index <- match(x,top.hits$gene)
  hit <- top.hits[index,]
  return(data.frame(hit))
}
top.cyto <- function(x){
  index <- match(x,top.hits$cyto)
  hit <- top.hits[index,]
  return(data.frame(hit))
}

#get genes list
test.genes <- lapply(gene.lst[seq_along(gene.lst)],top.genes)
test.genes <- do.call(rbind,test.genes)
test.genes <- test.genes %>% group_by(seqnames) %>% slice(c(1))
#get cyto bands
test <- lapply(cyto.lst[seq_along(cyto.lst)],top.cyto)
test <- do.call(rbind,test)
#test <- test[-which(test$gene==''),]
test <- test %>% group_by(seqnames) %>% slice(1)
#remove artifact names 
#test <- test[-4,]
#test <- find_overlaps(GRanges(test),kp$cytobands)
test <- test.genes
test <- test[-c(6,11,16),]
test <- rbind(test,top.hits[c(71),])
# test$gene[11] <- 'ATM'
# test$nex[11] <- 232
# test <- test[-7,]
# test$cyto[9] <- 'chr10.q24'

# test <- rbind(test,top.hits[20117,])
#saveRDS(test,'genes_TTFT.rds')

# plot
#title <- 'Overall Survival'
#title <- 'Time to First Treatment'
title <- 'Post treatment survival'
chrs <- paste0('chr',seq(1:22))
#png(paste(title,'.png'),width = 12,height = 7,units = 'in',res = 180)
kp <- plotKaryotype(plot.type = 4, chromosomes = chrs,labels.plotter = NULL,main = title)
kpAddChromosomeNames(kp,srt=45,cex=0.6)
#kpAddBaseNumbers(kp)
#kpPlotDensity(kp, data=gr.nex, r0=0, r1=0.3, window.size = 1000, col="orchid")
kpPoints(kp, data=del_pvals,y = del_pvals$logp, pch=21, 
         cex=del_pvals$nex/500,bg=del_pvals$colour,col='black',lwd = 0.1,ymax = max(del_pvals$logp))

#test <- readRDS('genes_TTFT.rds')
kpText(kp, data = GRanges(test), labels = paste0(test$cyto,'\n',test$gene,'','\n','N=',test$nex), y = test$logp, 
       ymax = max(del_pvals$logp),pos=4, cex=0.55, col='black')

range(del_pvals$nex)/500
#plot(density(del_pvals.noox$os.nex.log))
kpAxis(kp,ymin = 0,ymax = max(del_pvals$logp), cex=0.8,tick.pos = c(0,1.30103,max(del_pvals$logp)))
kpAddLabels(kp,'-log10(FDR)',srt=90,cex=0.8, r0 = 0,side = 'left')

kpAbline(kp, h=-log10(0.05), lty=3, ymax=max(del_pvals$logp), ymin =min(del_pvals$logp),cex=1.5,
         col='red')
#dev.off()
###################
## plot survival ##
###################
library(ggsurvfit)

gr.cll.cnv.noox <- read.csv('All_cnvs_with_ox.csv') %>% filter(numSNP > snp_ct & !study %in% 
                                                                  c('Oxford','Oxford-ARC','Oxford-ADM')) %>%  
  filter(!chr %in% c('X','Y','23')) %>% filter(CNV_type==type) %>% toGRanges() 


seqlevelsStyle(gr.cll.cnv) <- 'UCSC'
cll.outcome <- read.csv('All CLL outcome_spss.csv') %>%  filter(!Study %in% c('Oxford-ARC', 'Oxford-ADM'))

unique(cll.outcome$Study)

samples_from_range <- function(i) {
  outcomes <- cll.outcome
  ovlps <- find_overlaps(gr.cll.cnv, range(GRanges(top.hits)[i]))#filter_by_overlaps
  dat.ovlps <- data.frame(Sample=unique(ovlps$new_sam_id),alteration_all=1)
  outcomes <- left_join(outcomes,dat.ovlps,by='Sample')
  outcomes$alteration_all[is.na(outcomes$alteration_all)] <- 0
  # adding nexus
  ovlps_nex <- ovlps %>% filter(method=='Nexus')
  if(length(ovlps_nex) > 1) {
    dat.ovlps.nex <- data.frame(Sample=unique(ovlps_nex$new_sam_id),deletion=1)
    outcomes <- left_join(outcomes,dat.ovlps.nex,by='Sample')
    outcomes$deletion[is.na(outcomes$deletion)] <- 0  } 
  return(outcomes)
}

i <- 4154
outcomes <- samples_from_range(i)
fitos <- survfit(Surv(OSdays/365.25, OSstatus) ~alteration_all, data = outcomes)
fitttfs <- survfit(Surv(TTFTdays/365.25, TTFTstatus) ~alteration_all, data = outcomes)
fitpt <- survfit(Surv(OS_PTdays/365.25, OS_PTstatus) ~alteration_all, data = outcomes)
fit <- list(All_os=fitos,All_ttfs=fitttfs,All_pt=fitpt)
# for nex
fitos.n <- survfit(Surv(OSdays/365.25, OSstatus) ~deletion, data = outcomes)
fitttfs.n <- survfit(Surv(TTFTdays/365.25, TTFTstatus) ~deletion, data = outcomes)
fitpt.n <- survfit(Surv(OS_PTdays/365.25, OS_PTstatus) ~deletion, data = outcomes)
fit.n <- list(Nexus_os=fitos.n,Nexus_ttfs=fitttfs.n,Nexus_pt=fitpt.n)

outcomes <- cll.outcome
tri.os <- survfit(Surv(OSdays/365.25, OSstatus) ~TreatmentCenter, data = outcomes)


table(outcomes$Triploid)
res <- ggsurvplot(tri.os, pval = TRUE, conf.int = F,pval.coord = c(20, 1),add.all = F,
                   risk.table = TRUE # Add risk table
              ) 

res$plot <- res$plot + labs(title = paste0(top.hits$cyto[i]),subtitle = top.hits$gene[i])
print(res)

#####
boxplot(outcomes$Age + outcomes$OSdays/365.25)
out.shrt <- outcomes
out.shrt$OSyears <- out.shrt$OSdays/365.25
out.shrt <- reshape2::melt(out.shrt[,c(1,7,51)])
ggplot(out.shrt) + geom_bar(aes(x=Sample,y=value,fill=variable),stat = 'identity',position = 'stack')

#######################################
## Interrogate p-vals from diff.algos##
######################################

loss.out.all$os.p.log <- -log10(loss.out.all$os.p)
loss.out.all$os.p.nex.log <- -log10(loss.out.all$os.p.nex)

loss.out.nex.sig <- loss.out.all %>% filter(os.p.nex < 0.05 & os.p > 0.05)

range(loss.out.nex.sig$os.p.nex)

loss.out.nex.sig %>% ggplot(aes(disj_id))  + geom_bar(aes(y=os.p.nex.log),stat = 'identity', col='blue') + 
  geom_bar(aes(y=os.p.log),stat = 'identity',col='red',alpha=0.5) +geom_hline(yintercept=-log10(0.05), linetype="dashed", 
                                                                           color = "green", size=0.5)
                                                                       