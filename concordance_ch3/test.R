library(GenomicRanges)
library(dplyr)



#############
## test
red <- reduce(gr_split[[1]])
overlaps <- findOverlaps(red, gr_split[[1]])

caller_support <- tapply(subjectHits(overlaps), queryHits(overlaps),
  function(x) {length(unique(gr_split[[1]]$method[x])) })

# Add support metadata
mcols(red)$n_callers <- caller_support
mcols(red)$sample <- unique(gr_split[[1]]$new_sam_id)
mcols(red)$cohort <- unique(gr_split[[1]]$study)
mcols(red)$type <- unique(gr_split[[1]]$type)


# check if count_overlaps same as findoverlaps - 
for (i in seq_along(gr_split)) {
  overlaps2 <- count_overlaps(reduce(gr_split[[i]]), gr_split[[i]])
  print(overlaps2)
}

####################################
### test gap merge
####################################
source('functions.R')

gr.cll.cnv <- read.csv('../survival_mdr/All_cnvs_with_ox.csv') %>%  #& !study %in% c('Oxford','Oxford-ARC','Oxford-ADM')
  filter(!chr %in% c('X','Y','23')) %>% filter(numSNP > snp.ct) %>% toGRanges()

gr.cll.cnv <- gr.cll.cnv %>% filter(new_sam_id=='1536_HRH_HOE')
table(gr.cll.cnv$method)

gap.frac <- 0.2

#merge adjacent for each algo
algo.split <- split(gr.cll.cnv, paste(gr.cll.cnv$method, gr.cll.cnv$CNV_type, sep = "__"))
names(algo.split)
# do a loop for each algo as mcapply used in the function
gap.merged.cnv <- GRangesList()

for (i in seq_along(algo.split)) {
  print(names(algo.split[i]))
  gap.merged.cnv[[i]] <- merge_by_gap_fraction(algo.split[[i]],max_gap_fraction = gap.frac, ncores = 5)
  names(gap.merged.cnv[i]) <- names(algo.split[i])
  }


############################################
### plot Kp for sample
############################################
library(plyranges)
library(GenomicRanges)
library(karyoploteR)

source('functions_new.R')
# load all raw CNVs regardless of length or prob counts
use.merged.cnv <- F

### load raw CNVs or gap.merged
gr.cll.cnv <- load.raw_gap(snp.ct = 5, len = 5, use.merged.cnv = F)
gr.cll.gap <- load.raw_gap(snp.ct = 5, len = 5, use.merged.cnv = T)

#check
range(gr.cll.cnv$no_of_probes)
median(gr.cll.cnv$no_of_probes)
mean(gr.cll.cnv$no_of_probes)

#gr.cll.cnv.sam <- gr.cll.cnv
# filter sam
gr.cll.cnv.sam <- gr.cll.cnv %>% filter(sample_id=='B8013')
gr.cll.cnv.sam <-  split(gr.cll.cnv.sam, gr.cll.cnv.sam$caller)
elementNROWS(gr.cll.cnv.sam)
gr.cll.gap.sam <- gr.cll.gap %>% filter(sample_id=='B8013')
gr.cll.gap.sam <-  split(gr.cll.gap.sam, gr.cll.gap.sam$caller)

caller.idx <- 3
# crude kp plot
kp <- plotKaryotype(plot.type = 1, chromosomes = 'chr2') #chromosomes = 'chr6'
kpPlotRegions(kp, data=gr.cll.cnv.sam[[caller.idx]], col='blue', r1 = 0.5)
kpPlotRegions(kp, data=gr.cll.gap.sam[[caller.idx]], col='red', r0 = 0.5)


# check each caller for all cnv
gr.cll.cnv.caller <-  split(gr.cll.cnv, gr.cll.cnv$caller)
elementNROWS(gr.cll.cnv.caller)

# crude kp plot
kp <- plotKaryotype(plot.type = 1, chromosomes = 'chr3') #chromosomes = 'chr6'
kpPlotCoverage(kp, data=gr.cll.cnv.caller[[1]], col='blue', show.0.cov = T)
kpPlotCoverage(kp, data=gr.cll.cnv.caller[[2]], col='red', show.0.cov = T)
kpPlotCoverage(kp, data=gr.cll.cnv.caller[[3]], col='green', show.0.cov = T)
kpPlotCoverage(kp, data=gr.cll.cnv.caller[[4]], col='orange', show.0.cov = T)


################################
##kp plot coverage all sam
################################
pp <- getDefaultPlotParams(plot.type= 1)
pp$ideogramheight <- 30
### change params for kp plot ####
gr.cll.cnv <- load.raw_gap(snp.ct = 5, len = 5, use.merged.cnv = T)
# other params
point.size <- 0.4
# margins
top_lrr_r1 <- 1.56
top_lrr_r <- 1.1
top_baf_r1 <- 0.94
top_baf_r <- 0.48

sam <- '2279_GMJ'
chr.plt <- 'chr17'
type <- 'Loss'
### load raw CNVs or gap.merged
gr.cll.cnv.sam <- gr.cll.cnv %>% filter(sample_id==sam & loss_gain==type)
#gr.cll.cnv.sam <- gr.cll.cnv %>% filter(sample_id==sam)
gr.cll.cnv.caller <-  split(gr.cll.cnv.sam, gr.cll.cnv.sam$caller)
# load  sig intensity
prb.file <- read.delim(paste0('../sig_intensities/',sam,'.txt'), comment.char = '#')
names(prb.file) <- c('snp', 'LRR','BAF', 'chr', 'pos')
#remove loh- Homozygous Copy Loss
prb.file <- GRanges(ranges = IRanges(start= prb.file$pos, width = 1), 
                    seqnames = prb.file$chr, lrr= prb.file$LRR, baf=prb.file$BAF)
seqlevelsStyle(prb.file) <- 'UCSC'

#target_region <- 'chr18:1-19560000'
# crude kp plot
png(filename = paste0('thesis_out/lrr_baf',sam,'.png'), width = 8, height = 5, res = 300, units = 'in')
kp <- plotKaryotype(plot.type = 1, chromosomes = chr.plt, plot.params = pp)#, zoom = target_region) #chromosomes = 'chr6'
kpAddBaseNumbers(kp, tick.dist = 30e6, tick.col="black", cex=0.7, add.units = T)
try(kpPlotRegions(kp, data=gr.cll.cnv.caller[[1]], col='#E7298A', r0 =0, r1 = 0.1, avoid.overlapping = F))
kpAddLabels(kp,labels = 'Nexus',r0 =0, r1 = 0.1, cex=0.7, pos = 1, label.margin = 0.08)
try(kpPlotRegions(kp, data=gr.cll.cnv.caller[[2]], col='#1B9E77', r0 =0.12, r1 = 0.22, avoid.overlapping = F))
kpAddLabels(kp,labels = 'Pcnv', r0 =0.12, r1 = 0.22, cex=0.7, pos = 1,label.margin = 0.08)
try(kpPlotRegions(kp, data=gr.cll.cnv.caller[[3]], col='#D95F02', r0 =0.24, r1 = 0.34, avoid.overlapping = F))
kpAddLabels(kp,labels = 'Qsnp', r0 =0.24, r1 = 0.34, cex=0.7, pos = 1,label.margin = 0.08)
try(kpPlotRegions(kp, data=gr.cll.cnv.caller[[4]], col='#7570B3', r0 =0.36, r1 = 0.46, avoid.overlapping = F))
kpAddLabels(kp,labels = 'Ipn', r0 =0.36, r1 = 0.46, cex=0.7, pos = 1,label.margin = 0.08)

## plot probe
#kp <- plotKaryotype(genome = 'hg19',plot.type = 1, chromosomes = 'chr6')
#kpAddCytobandLabels(karyoplot = kp,clipping = T,cex = 0.8) #srt=90
#kpAddMainTitle(kp,main =paste0(sam),cex=1)
# BAF
kpPoints(kp, data=prb.file, y=prb.file$baf, r0=top_baf_r, r1=top_baf_r1, cex=point.size, data.panel = 1,col='gray')#slategray3
kpAxis(kp, ymax = 1, ymin = 0, r0 = top_baf_r, r1 = top_baf_r1, cex=0.8, data.panel = 1)
kpAddLabels(kp,labels = 'BAF', r0 = top_baf_r, r1 = top_baf_r1, cex=0.8, srt=90, pos = 1, label.margin = 0.08)
# LRR
range(prb.file$lrr, na.rm = T)
y.max <- 2
y.min <- -2
kpPoints(kp, data=prb.file, y=prb.file$lrr, r0=top_lrr_r, r1=top_lrr_r1,cex=point.size , col='gray', ymin = y.min, ymax = y.max)
#kpPoints(kp, data=cnv, y=cnv$lrr, r0=top_lrr_r, r1=top_lrr_r1,cex=point.size , col='blue', ymin = y.min, ymax = y.max)
kpAxis(kp,ymin = y.min, ymax = y.max, r0=top_lrr_r, r1=top_lrr_r1, cex=0.8)
kpAddLabels(kp,labels = 'LRR',r0=top_lrr_r, r1=top_lrr_r1, cex=0.8,srt=90,pos = 1,label.margin = 0.08)
dev.off()
