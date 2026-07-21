# step 2- chapter 4 work
# carry forward from the concordance solved disjs
############################################################
# STEP 2
# CREATE COHORT-WIDE DISJOINT SEGMENTS - Recurrence
# Goal:
# Get the concordance for each disj but don't filter
# Filter after looking at the recurrence
# preserve focal structure
# # SCORE EACH DISJOINT SEGMENT
# create disjoints for recurrence from >3 concordant djs
# then count overlaps of all the concordant djs
# this way recurrent ranges are made using highly concordant djs
# but they capture all the djs regardless of concordance
############################################################
library(Matrix)
library(survival)

source('functions.R')
source('functions_new.R')

type <- 'Loss'


gr.cll.cnv <- load.raw_gap(snp.ct = 50, len = 500, use.merged.cnv = use.merged.cnv, only.common.sam = T)
gr.cll.cnv <- gr.cll.cnv %>% filter(loss_gain == type)

# create disjoints for recurrence using highly concordant djs
dj <- disjoin(gr.cll.cnv)
length(dj)

hits <- findOverlaps(dj, gr.cll.cnv)

# Number of samples overlapping each segment
sample_support <- tapply(subjectHits(hits), queryHits(hits),
                         function(x) {length(unique(gr.cll.cnv$sample_id[x])) } )

# Mean caller support
mean_callers <- tapply(subjectHits(hits), queryHits(hits), function(x) {
  sample_callers <- split(gr.cll.cnv$caller[x], gr.cll.cnv$sample_id[x])
  mean(sapply(sample_callers, function(z) length(unique(z))))
})

n_callers <- tapply(subjectHits(hits), queryHits(hits), function(x) {
  length(unique(gr.cll.cnv$caller[x]))
})


# Attach metadata
mcols(dj)$n_samples <- 0
mcols(dj)$mean_callers <- 0

mcols(dj)$n_samples[as.numeric(names(sample_support))] <- sample_support
mcols(dj)$mean_callers[as.numeric(names(mean_callers))] <- mean_callers


hist(dj$mean_callers)

# filter recurrence disjoints
dj_filtered <- dj[dj$n_samples >= 10]
# plot
kp <- plotKaryotype(plot.type = 1, chromosomes = 'autosomal') #chromosomes = 'chr6'
kpPlotCoverage(kp, data=dj_filtered, show.0.cov = T)
kpPlotCoverage(kp, data=gr.cll.cnv, show.0.cov = T, col='red')

# create disj matrix for cox
gr.cox.out <- cna.mat.cox(dj_filtered, gr.cll.cnv, surv_type = 'OS', 10)
range(gr.cox.out$HR)
range(gr.cox.out$FDR)
range(gr.cox.out$pval)
range(gr.cox.out$n_samples)

# plot
kp <- plotKaryotype(plot.type = 1, chromosomes = 'autosomal') #chromosomes = 'chr6'
kpPlotCoverage(kp, data=dj_filtered, show.0.cov = T, col = 'pink')
kpPlotCoverage(kp, data=gr.cox.out, show.0.cov = T)

unique(seqnames(gr.cox.out))
##############################
### plotting
##############################
## plot mahnatton
#add chr col
plot_manhatton <- function(gr.cox_results) {
chr.col <- as.numeric(gsub('chr','', seqnames(gr.cox_results)))
gr.cox_results$chr.col <- chr.col%%2
table(gr.cox_results$chr.col)
gr.cox_results$chr.col <- gsub('^1','#FFBD07AA', gr.cox_results$chr.col)
gr.cox_results$chr.col <- gsub('^0','#00A6EDAA', gr.cox_results$chr.col)
gr.cox_results$pos <- (end(gr.cox_results) + start(gr.cox_results))/2

range(-log10(gr.cox_results$FDR), na.rm = T)

chrs <- paste0('chr',seq(1:22))
#png(paste(title,'.png'),width = 12,height = 7,units = 'in',res = 180)
kp <- plotKaryotype(plot.type = 4, chromosomes = chrs, labels.plotter = NULL)
kpAddChromosomeNames(kp,srt=45,cex=0.6)
#kpAddBaseNumbers(kp)
#kpPlotDensity(kp, data=gr.nex, r0=0, r1=0.3, window.size = 1000, col="orchid")
kpPoints(kp, data=gr.cox_results, y = -log10(gr.cox_results$FDR), pch=21, col='black', cex = 0.6,
         lwd = 0.2, ymax = max(-log10(gr.cox_results$FDR),na.rm = T), bg=gr.cox_results$chr.col)
}

plot_manhatton(gr.cox.out)


kpText(kp, data = GRanges(test), labels = paste0(test$cyto,'\n',test$gene,'','\n','N=',test$nex), y = test$logp, 
       ymax = max(del_pvals$logp),pos=4, cex=0.55, col='black')

range(del_pvals$nex)/500
#plot(density(del_pvals.noox$os.nex.log))
kpAxis(kp,ymin = 0,ymax = max(del_pvals$logp), cex=0.8,tick.pos = c(0,1.30103,max(del_pvals$logp)))
kpAddLabels(kp,'-log10(FDR)',srt = 90,cex = 0.8, r0 = 0,side = 'left')

kpAbline(kp, h=-log10(0.05), lty=3, ymax=max(del_pvals$logp), ymin =min(del_pvals$logp),cex=1.5,
         col='red')
