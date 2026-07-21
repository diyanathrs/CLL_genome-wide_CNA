# Calculate survival p-vals for all the MDRs?
# Need to include oxford data as well so use 3 algos that we used oxford data with
# quickly create disjoints from ranges and do surv
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
snp_ct <- 20
gr.cll.cnv <- read.csv('All_cnvs_with_ox.csv') %>% filter(numSNP > snp_ct) %>%  
  filter(!chr %in% c('X','Y','23')) %>% filter(CNV_type==type) %>% toGRanges() 

# some stats of raw cnv calls
length(unique(gr.cll.cnv$new_sam_id))
plot(density(gr.cll.cnv$length))

seqlevelsStyle(gr.cll.cnv) <- 'UCSC'
# convert rpt to normal from ONLY Leicester sample ids
#dat.cll.cnv <- data.frame(gr.cll.cnv)
#dat.cll.cnv$new_sam_id <- dat.cll.cnv$new_sam_id %>% gsub('_rpt','',.)

#load outcome data
cll.outcome <- read.csv('All CLL outcome_spss.csv') #%>%  filter(!Study %in% c('Oxford-ARC', 'Oxford-ADM'))  to remove oxford samples 
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

# Run lapply with 20 cores ####
loss.out.all <- mclapply(seq_along(ranges.out), cal_survival,mc.cores = 20)
#loss.out <- mclapply(1:500, cal_survival,mc.cores = 20)
loss.out.all <- as.data.frame(do.call(rbind,loss.out.all))
loss.out.all <- loss.out.all[,c(1,2,3,4,5,6,7,8,10,9,11,14,12,15,13,16,17)]
# Saving the data
write.csv(loss.out.all, paste0('cll_loss',snp_ct,'_pvals.csv'),quote = F,row.names = F)

# Post-Hoc testing on selected disjs ####
# load dels and gains p-vals from data 
del_pvals <- read.csv('cll_del_pvals.csv')
del_pvals.noox <- read.csv('cll_loss_pvals_noox.csv')
gain_pvals <- read.csv('cll_gains_pvals.csv')
# get sample names for particular disj and plot surv
disj_id <- 81641
cal_survival(disj_id)
plot_survival(disj_id,'nex')
# what are the sample ids that are detected by others but not nex?
setdiff(outcomes$Sample[outcomes$alteration_all==1],outcomes$Sample[outcomes$alteration_nex==1])

####################
# Post analysis ####
####################

# get reduced table ###
loss.reduced.os <- reduce_pval.list(loss.out.all,'os')
loss.reduced.pt <- reduce_pval.list(del_pvals.noox,'pt')
loss.reduced.ttfs <- reduce_pval.list(del_pvals.noox,'ttfs')

# filter sig p-val
pvals.os <- loss.out.all %>% filter(total > 15 & os.p.nex < 0.05) 
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

#############
##Adjust P###
#############
loss.out.all <- read.csv('cll_loss20_pvals.csv')
loss.out.all$os.adj.p <- p.adjust(loss.out.all$os.p.nex,method = 'fdr')

range(loss.out.all$os.p.nex,na.rm = T)
range(loss.out.all$os.adj.p,na.rm = T)

##################
## Plotting ######
##################
library(karyoploteR)
gr.nex <- read.csv('All_cnvs_with_ox.csv') %>% filter(numSNP > 5 & method=='Nexus') %>%  
  filter(!chr %in% c('X','Y','23')) %>% filter(CNV_type==type) %>% toGRanges()

# p val plot
del_pvals <- loss.out.all %>% filter(nex > 15 )#& os.p.nex <0.05)

# convert p val to -log10
del_pvals$os.p.nex[is.na(del_pvals$os.adj.p)] <- 1
del_pvals$os.nex.log <- -log10(del_pvals$os.adj.p)
range(del_pvals$os.nex.log)
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
                     logp=del_pvals$os.nex.log, colour=del_pvals$chr.col, gene=del_pvals$genes)

# get top hits by reducing
#top.hits <- reduce_pval.list(loss.out.all,'os')
top.hits <- as.data.frame(del_pvals) %>% filter(logp > 1.3) %>% arrange(-logp)
gene.lst <- unique(top.hits$gene)

top.genes <- function(x){
  index <- match(x,top.hits$gene)
  hit <- top.hits[index,]
  return((hit))
}

test <- lapply(gene.lst[seq_along(gene.lst)],top.genes)
test <- do.call(rbind,test)
test <- test[-2,]
test <- test %>% group_by(seqnames) %>% slice(1)
test <- data.frame(test)

# plot
title <- 'Overall Survival'
chrs <- paste0('chr',seq(1:22))
kp <- plotKaryotype(plot.type = 4, chromosomes = chrs,labels.plotter = NULL,main = title)
kpAddChromosomeNames(kp,srt=45,cex=0.6)
#kpAddBaseNumbers(kp)

#kpPlotDensity(kp, data=gr.nex, r0=0, r1=0.3, window.size = 1000, col="orchid")

kpPoints(kp, data=del_pvals,y = del_pvals$logp, pch=21, 
         cex=del_pvals$nex/400,bg=del_pvals$colour,col='black',lwd = 0.1,ymax = max(del_pvals$logp))

kpText(kp, data = toGRanges(test), labels = paste(test$gene,'\n','N=',test$nex), y = test$logp, 
       ymax = max(del_pvals$logp),pos=4, cex=0.5, col='black')

range(del_pvals$nex)/500
#plot(density(del_pvals.noox$os.nex.log))
kpAxis(kp,ymin = 0,ymax = max(del_pvals$logp), cex=0.8,tick.pos = c(0,1.30103,max(del_pvals$logp)))
kpAddLabels(kp,'-log10(FDR)',srt=90,cex=0.8, r0 = 0,side = 'left')

kpAbline(kp, h=-log10(0.05), lty=3, ymax=max(del_pvals$logp), ymin =min(del_pvals$logp),cex=1.5,
         col='red')

#plot(density(del_pvals$os.p.nex,na.rm = T))


### Plot survival##




# manually checked list - chr1_LINC00184** ,chr1_SPOCD1*, chr3_CAV3, chr3_SUMF1*,
# chr4_EVC(1-10mb),chr6_COL19A1,chr6_RWDD2A(78-100mb detected by sarah), chr8p,
# chr9_MLLT3(20-22mb), chr9_GABBR2(100-114), chr10_BTRC**,chr10_STK32C*,chr11_MYBPC3*, chr11_IGSF9B(81-134mb)
# chr12_WASHC4*, chr13_delu1/2, chr14_LGMN(68-105mb), chr15_EHD4*, chr17_DOC2B(can't find but low.p)
# chr17_VPS53(p.arm is very low.p), chr18_RNMT(p-arm), chr19_ZNF826P**, chr19_PSG6*-,
# chr20_DEFB125(p-arm is high sig), chr22_ZNF280B* (22-25mb - high sig)
#  Do the same for TTFS
#
# interesting disjs - loss - disj_6361
# Delete ####
dat.cll.cnv <-  read.csv('All_cnvs_with_ox.csv') 
sam_list <- data.frame(table(dat.cll.cnv$new_sam_id))

# first arc with 1_s
arc2_list <- data.frame(dat.cll.cnv$Sample_ID[grep(pattern ='arc',dat.cll.cnv$Sample_ID,ignore.case = T)])
arc_list <- data.frame(table(arc2_list))
# remove after first _ delim
str_extract(arc_list$dat.cll.cnv.Sample_ID.grep.pattern....arc...dat.cll.cnv.Sample_ID..,"[^_]*_[^_]*") # Works

# add arc_list_ed back to dat.cll.cnv
dat.cll.cnv$Sample_ID[grep(pattern ='arc',dat.cll.cnv$Sample_ID,ignore.case = T)] <- arc_list_ed

ox <- dat.cll.cnv %>% filter(study %in% c('Oxford','Oxford-ADM','Oxford-ARC'))

list.replacing <- dat.cll.cnv$Sample_ID
# use a loop to replace ox names
rep.nex.pen <- read.csv('data/nex_to_pen.csv')

for (i in seq_along(rep.nex.pen$from)) {
  
  list.replacing <- gsub(rep.nex.pen$from[i],rep.nex.pen$to[i],list.replacing)
  
}

dat.cll.cnv$new_sam_id <- list.replacing

#### - ####


# curves for each study