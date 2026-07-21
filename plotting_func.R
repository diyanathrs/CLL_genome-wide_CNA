# Processing the sum_loss and sum_gain outputs created in earlier stages
# load packages
# Quickly plot the cnvrs and cnvs in them
library(ggplot2)
library(dplyr)
library(cowplot)
library(data.table)
library(karyoploteR)
library(regioneR)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(plyranges)
library(GenomicRanges)


# plotting function for reduced ranges ####
gr.cll.cnv %>% reduce()
 # create reduced ranges from loss out
gr.ranges.reduced  <- reduce(makeGRangesFromDataFrame(ranges.out))
start(gr.ranges.reduced) <- start(gr.ranges.reduced)-100000
end(gr.ranges.reduced) <- end(gr.ranges.reduced)+100000

# plot by a disjoint
disj_id
point_size=0.5
i <- 13
# plot the lowest p val for first red.range
t <- filter_by_overlaps(toGRanges(del_pvals.noox),ranges.out[disj_id])
t <- t %>% filter(total > 15)
t <- t[order(t$os.p.nex),]
length(t)

# now plot the t[1] ?
j=5
# need to go back to cnv table and get sample ids
cnvr_sam <- unique(filter_by_overlaps(gr.cll.cnv,ranges.out[disj_id])$Sample_ID)
#load snp data
sig <- read.delim(paste0('/home/dean91/cnv_cll/sig_intensities/',cnvr_sam[j],'.txt'))
names(sig) <- c('name','LRR','BAF','Chr','Position')
snp.data <- toGRanges(sig[,c("Chr", "Position", "Position", "BAF", "LRR")])
seqlevelsStyle(snp.data) <- "UCSC"

# separate individual cnvs
gr_cnvs_sam_n <- gr.cll.cnv %>% filter(method=='Nexus' & Sample_ID==cnvr_sam[j])
gr_cnvs_sam_i <- gr.cll.cnv %>% filter(method=='iPattern' & Sample_ID==cnvr_sam[j])
gr_cnvs_sam_p <- gr.cll.cnv %>% filter(method=='PennCNV' & Sample_ID==cnvr_sam[j]) 
gr_cnvs_sam_q <- gr.cll.cnv %>% filter(method=='QuantiSNP' & Sample_ID==cnvr_sam[j]) 

#png(filename = paste0(folder,'/',cnvr_sam[j],'.png'),width =1000 ,height = 600,res = 120)

kp <- plotKaryotype(genome = 'hg19',plot.type = 1,zoom=ranges.out[disj_id],plot.params = pp)
kpAddCytobandLabels(kp, cex=0.7,force.all = T) #srt=90\
kpAddBaseNumbers(kp,add.units = T,tick.dist = 0.5e5,units = 'mb')

# plot dgv bars too
kpPlotCoverage(kp,data = gr_dgv_cnv,r0 = 0.86,r1=1,col='#c542f5',clipping = T,border = F,show.0.cov = T)
kpAddLabels(kp,'DGV',r0 = 0.86,r1 = 1,srt=90,cex=0.6,label.margin = 0.04,pos = 1)

# add cnv rect
#kpPlotRegions(kp,data = gr_cnvs_sam_n,r0 = 0.57,r1=0.59,col='orange',clipping = T,border = NULL)
kpPlotRegions(kp,data = gr_cnvs_sam_i,r0 = 0.47,r1=0.48,col='blue',clipping = T,border = NULL)
kpPlotRegions(kp,data = gr_cnvs_sam_p,r0 = 0.49,r1=0.50,col='red',clipping = T,border = NULL)
kpPlotRegions(kp,data = gr_cnvs_sam_q,r0 = 0.51,r1=0.52,col='green',clipping = T,border = NULL)
kpPlotRegions(kp,data = gr_cnvs_sam_n,r0 = 0.53,r1=0.54,col='orange',clipping = T,border = NULL)

#kpAddCytobandsAsLine(kp)
kpAddMainTitle(kp,paste0('Sample: ',cnvr_sam[j],': CNVR_',i,':',granges(gr_red_loss)[i]),cex=0.8)
kpAxis(kp, r0=0.55, r1=0.85,cex=0.6,ymin = -2,ymax=1,tick.pos = c(1,0,-2))
kpPoints(kp, data=snp.data, y=snp.data$LRR, r0=0.75,col='black',pch = 21,bg='#f5b942',lwd=0.5,cex = point_size)
kpAddLabels(kp,'LRR',r0 = 0.55,r1 = 0.85,srt=90,cex=0.6,label.margin = 0.04,pos = 1)

kpAxis(kp, r0=0.15, r1=0.45,cex=0.6,ymax = 1,ymin = 0,tick.pos = c(0,1))
kpPoints(kp, data=snp.data, y=snp.data$BAF, r0=0.15, r1=0.45,col='black',ymax = 1,pch=21,bg='#42f5e3',lwd=0.5,cex = point_size)
kpAddLabels(kp,'BAF',r0 = 0.15,r1 = 0.45,srt=90,cex=0.6,label.margin = 0.04,pos = 1)

#kp <- plotKaryotype(genome = 'hg19',plot.type = 4,chromosomes = 'chr3')

kpPlotGenes(kp, data=genes.data, add.transcript.names = F,r1=0.18, plot.transcripts = F,
            gene.name.position = "right",data.panel = 1,gene.name.cex = 0.4,avoid.overlapping = T,
            clipping = F,gene.col = 'lightblue',gene.border.col = NA)
    




#end ####
#plot only dat_red_loss[x]
# need a way to find what samples are in each range of no_nex? 
# load all cnvs
# Load all cnvs (-ox)
dir.create('custom')

cnvs_loss$Sample_ID <- gsub('1536_HRH-HOE','1536_HRH_HOE',cnvs_loss$Sample_ID)

#add chr string to chr
cnvs_loss$chr <- paste0('chr',cnvs_loss$chr)
# convert all cnvs to gr object without subsampling to losses or gains
gr_cnvs_loss <- cnvs_loss  %>% toGRanges()
#see overlaping samples for particular cnvr
## new addition - reduce ranges for nonex

# plot probs for each nonex cnvr on a loop
# plot params
pp <- getDefaultPlotParams(plot.type=1)
pp$ideogramheight <- 5
pp$data1inmargin <- 2
pp$bottommargin <- 15
pp$topmargin <- 6

i= 2 # custom range from the list

folder <- paste0('custom/range_1algo_',seqnames(gr_red_loss)[i],'_',start(gr_red_loss)[i],'_',end(gr_red_loss)[i])
dir.create(path = folder)
# pulling genes
print('Pulling genes for custom CNVR')

# convert mrd p-val tabale to granges
gr.loss.reduced.noox <- reduce(toGRanges(loss.out.noox))
# use reduced ranges of loss.out
kp <- plotKaryotype(genome = 'hg19',plot.type = 1,zoom=gr.loss.reduced.noox[2],plot.params = pp)
genes.data <- makeGenesDataFromTxDb(txdb = TxDb.Hsapiens.UCSC.hg19.knownGene,karyoplot = kp,
                                    plot.transcripts = T,plot.transcripts.structure = F)
genes.data <- addGeneNames(genes.data)

# now get samples for each cnvr and cnvr coords
cnvr_sam <- unique(filter_by_overlaps(gr.cll.cnv,gr.loss.reduced.noox[1])$Sample_ID)


j <- 6
#for(j in 1:length(cnvr_sam)) {
  print(paste('Plotting sample:',cnvr_sam[j],'probes of CNVR_',i))
  
  #load snp data
  sig <- read.delim(paste0('/home/dean91/cnv_cll/sig_intensities/',cnvr_sam[j],'.txt'))   
  names(sig) <- c('name','LRR','BAF','Chr','Position')
  snp.data <- toGRanges(sig[,c("Chr", "Position", "Position", "BAF", "LRR")])
  seqlevelsStyle(snp.data) <- "UCSC"
  
  # separate individual cnvs
  gr_cnvs_sam_n <- gr.cll.cnv %>% filter(method=='Nexus' & Sample_ID==cnvr_sam[j])
  gr_cnvs_sam_i <- gr.cll.cnv %>% filter(method=='iPattern' & Sample_ID==cnvr_sam[j])
  gr_cnvs_sam_p <- gr.cll.cnv %>% filter(method=='PennCNV' & Sample_ID==cnvr_sam[j]) 
  gr_cnvs_sam_q <- gr.cll.cnv %>% filter(method=='QuantiSNP' & Sample_ID==cnvr_sam[j]) 
  
  #png(filename = paste0(folder,'/',cnvr_sam[j],'.png'),width =1000 ,height = 600,res = 120)
  
  kp <- plotKaryotype(genome = 'hg19',plot.type = 1,zoom=gr.loss.reduced.noox[2],plot.params = pp)
  kpAddCytobandLabels(kp, cex=0.7,force.all = T) #srt=90\
  kpAddBaseNumbers(kp,add.units = T,tick.dist = 2e6,units = 'mb')
  
  # plot dgv bars too
  kpPlotCoverage(kp,data = gr_dgv_cnv,r0 = 0.86,r1=1,col='#c542f5',clipping = T,border = F,show.0.cov = T)
  kpAddLabels(kp,'DGV',r0 = 0.86,r1 = 1,srt=90,cex=0.6,label.margin = 0.04,pos = 1)
  
  # add cnv rect
  #kpPlotRegions(kp,data = gr_cnvs_sam_n,r0 = 0.57,r1=0.59,col='orange',clipping = T,border = NULL)
  kpPlotRegions(kp,data = gr_cnvs_sam_i,r0 = 0.47,r1=0.48,col='blue',clipping = T,border = NULL)
  kpPlotRegions(kp,data = gr_cnvs_sam_p,r0 = 0.49,r1=0.50,col='red',clipping = T,border = NULL)
  kpPlotRegions(kp,data = gr_cnvs_sam_q,r0 = 0.51,r1=0.52,col='green',clipping = T,border = NULL)
  kpPlotRegions(kp,data = gr_cnvs_sam_n,r0 = 0.53,r1=0.54,col='orange',clipping = T,border = NULL)
  
  #kpAddCytobandsAsLine(kp)
  kpAddMainTitle(kp,paste0('Sample: ',cnvr_sam[j],': CNVR_',i,':',granges(gr_red_loss)[i]),cex=0.8)
  kpAxis(kp, r0=0.55, r1=0.85,cex=0.6,ymin = -2,ymax=1,tick.pos = c(1,0,-2))
  kpPoints(kp, data=snp.data, y=snp.data$LRR, r0=0.75,col='black',pch = 21,bg='#f5b942',lwd=0.5,cex = point_size)
  kpAddLabels(kp,'LRR',r0 = 0.55,r1 = 0.85,srt=90,cex=0.6,label.margin = 0.04,pos = 1)
  
  kpAxis(kp, r0=0.15, r1=0.45,cex=0.6,ymax = 1,ymin = 0,tick.pos = c(0,1))
  kpPoints(kp, data=snp.data, y=snp.data$BAF, r0=0.15, r1=0.45,col='black',ymax = 1,pch=21,bg='#42f5e3',lwd=0.5,cex = point_size)
  kpAddLabels(kp,'BAF',r0 = 0.15,r1 = 0.45,srt=90,cex=0.6,label.margin = 0.04,pos = 1)
  
  #kp <- plotKaryotype(genome = 'hg19',plot.type = 4,chromosomes = 'chr3')
  
  kpPlotGenes(kp, data=genes.data, add.transcript.names = F,r1=0.18, plot.transcripts = F,
              gene.name.position = "right",data.panel = 1,gene.name.cex = 0.4,avoid.overlapping = T,
              clipping = F,gene.col = 'lightblue',gene.border.col = NA)
  
  #dev.off()


####

# need a way to find what samples are in each range of no_nex? ####
# load all cnvs
# Load all cnvs (-ox)
dir.create('custom')
cnvs_loss <- read.csv('All_cnvs_noox_latest.csv') %>% filter(numSNP>=20) %>% 
  dplyr::select(1:5,9:12) %>% filter(!chr %in% c('X','Y','23') & CNV_type=='Loss')
#add chr string to chr
cnvs_loss$chr <- paste0('chr',cnvs_loss$chr)
# convert all cnvs to gr object without subsampling to losses or gains
gr_cnvs_loss <- cnvs_loss  %>% toGRanges()
#see overlaping samples for particular cnvr
## new addition - reduce ranges for nonex
gr_nonex <- toGRanges(no_nex)
gr_nonex <- reduce(gr_nonex)
# plot probs for each nonex cnvr on a loop
# plot params
pp <- getDefaultPlotParams(plot.type=1)
pp$ideogramheight <- 5
pp$data1inmargin <- 2
pp$bottommargin <- 15
pp$topmargin <- 6

i= 1
j=1
for (i in 1:length(gr_nonex)) {
  print(paste('Creating folder for reduced CNVR_',i,granges(gr_nonex)[i]))
  folder <- paste0('custom/cnvr_',i,'_',seqnames(gr_nonex)[i],'_',start(gr_nonex)[i],'_',end(gr_nonex)[i])
  dir.create(path = folder)
  print(paste('Plotting probes for reduced CNVR_',i,granges(gr_nonex)[i]))
  # pulling genes
  print('Pulling genes for CNVR')
  kp <- plotKaryotype(genome = 'hg19',plot.type = 1,zoom=granges(gr_nonex)[i]+0.5e6,plot.params = pp)
  genes.data <- makeGenesDataFromTxDb(txdb = TxDb.Hsapiens.UCSC.hg19.knownGene,karyoplot = kp,
                                      plot.transcripts = T,plot.transcripts.structure = F)
  genes.data <- addGeneNames(genes.data)
  
  # now get samples for each cnvr and cnvr coords
  cnvr_sam <- unique(filter_by_overlaps(gr_cnvs_loss,gr_nonex[i])$Sample_ID)
  
  for(j in 1:length(cnvr_sam)) {
    print(paste('Plotting sample:',cnvr_sam[j],'probes of CNVR_',i))
    
    #load snp data
    sig <- read.delim(paste0('../sig_intensities/',cnvr_sam[j],'.txt'))   
    names(sig) <- c('name','LRR','BAF','Chr','Position')
    snp.data <- toGRanges(sig[,c("Chr", "Position", "Position", "BAF", "LRR")])
    seqlevelsStyle(snp.data) <- "UCSC"
    
    # separate individual cnvs
    gr_cnvs_sam_n <- gr_cnvs_loss %>% filter(method=='Nexus' & Sample_ID==cnvr_sam[j])
    gr_cnvs_sam_i <- gr_cnvs_loss %>% filter(method=='iPattern' & Sample_ID==cnvr_sam[j])
    gr_cnvs_sam_p <- gr_cnvs_loss %>% filter(method=='PennCNV' & Sample_ID==cnvr_sam[j]) 
    gr_cnvs_sam_q <- gr_cnvs_loss %>% filter(method=='QuantiSNP' & Sample_ID==cnvr_sam[j]) 
    
    png(filename = paste0(folder,'/',cnvr_sam[j],'.png'),width =1000 ,height = 600,res = 120)
    
    kp <- plotKaryotype(genome = 'hg19',plot.type = 1,zoom=granges(gr_nonex)[i]+0.5e6,plot.params = pp)
    kpAddCytobandLabels(kp, cex=0.7,force.all = T) #srt=90\
    kpAddBaseNumbers(kp,add.units = T,tick.dist = 1e5,units = 'mb')
    
    # plot dgv bars too
    kpPlotCoverage(kp,data = gr_dgv_cnv,r0 = 0.86,r1=1,col='#c542f5',clipping = T,border = F,show.0.cov = T)
    kpAddLabels(kp,'DGV',r0 = 0.86,r1 = 1,srt=90,cex=0.6,label.margin = 0.04,pos = 1)
    
    # add cnv rect
    #kpPlotRegions(kp,data = gr_cnvs_sam_n,r0 = 0.57,r1=0.59,col='orange',clipping = T,border = NULL)
    kpPlotRegions(kp,data = gr_cnvs_sam_i,r0 = 0.46,r1=0.48,col='blue',clipping = T,border = NULL)
    kpPlotRegions(kp,data = gr_cnvs_sam_p,r0 = 0.49,r1=0.51,col='red',clipping = T,border = NULL)
    kpPlotRegions(kp,data = gr_cnvs_sam_q,r0 = 0.52,r1=0.54,col='green',clipping = T,border = NULL)
    kpPlotRegions(kp,data = gr_cnvs_sam_n,r0 = 0.53,r1=0.55,col='green',clipping = T,border = NULL)
    
    #kpAddCytobandsAsLine(kp)
    kpAddMainTitle(kp,paste0('Sample: ',cnvr_sam[j],': CNVR_',i,':',granges(gr_nonex)[i]),cex=0.8)
    kpAxis(kp, r0=0.55, r1=0.85,cex=0.6,ymin = -2,ymax=1,tick.pos = c(1,0,-2))
    kpPoints(kp, data=snp.data, y=snp.data$LRR, r0=0.75,col='black',pch = 21,bg='#f5b942',lwd=0.5)
    kpAddLabels(kp,'LRR',r0 = 0.55,r1 = 0.85,srt=90,cex=0.6,label.margin = 0.04,pos = 1)
    
    kpAxis(kp, r0=0.15, r1=0.45,cex=0.6,ymax = 1,ymin = 0,tick.pos = c(0,1))
    kpPoints(kp, data=snp.data, y=snp.data$BAF, r0=0.15, r1=0.45,col='black',ymax = 1,pch=21,bg='#42f5e3',lwd=0.5)
    kpAddLabels(kp,'BAF',r0 = 0.15,r1 = 0.45,srt=90,cex=0.6,label.margin = 0.04,pos = 1)
    
    #kp <- plotKaryotype(genome = 'hg19',plot.type = 4,chromosomes = 'chr3')
    
    kpPlotGenes(kp, data=genes.data, add.transcript.names = F,r1=0.15, plot.transcripts = F,
                gene.name.position = "right",data.panel = 1,gene.name.cex = 0.4,avoid.overlapping = T,
                clipping = F,gene.col = 'lightblue',gene.border.col = NA)
    
    dev.off()
  }
}


#genes.data <- mergeTranscripts(genes.data)
#####

## Input2
ch=3
Start= sum_loss$start[105958]
End = sum_loss$end[105958]
Start= 162.1e6
End = 162.2e6
Start <-Start-100000
End <- End+100000

# subset cnvs
cnvs_chr <- cnvs %>% filter(chr==ch & CNV_type=='Loss')
#first filter based on target region and then number
#cnvs_chr <- cnvs_chr%>% filter(posStart >= start)
cnvs_chr <- cnvs_chr%>% filter(posStart >= Start  & posEnd <= End)
cnvs_chr <- arrange(cnvs_chr,length)
cnvs_chr <- arrange(cnvs_chr,method)

# add number for each row in cnv_chr table
## Important!! generate numbers with intervals- seq_len don't do this
len <- length(cnvs_chr$chr)
cnv_ids <- seq(1,(len*3),by=3)
cnvs_chr <- cbind(cnvs_chr,cnv_ids)

plt2 <- cnvs_chr %>% ggplot(aes(xmin=posStart/1000,xmax=posEnd/1000,ymin=cnv_ids,ymax=cnv_ids+2))+
  geom_rect(aes(fill=method))+ xlim(Start/1000,End/1000)+ theme_test()+
  theme(axis.text.y = element_blank(),
        axis.ticks.y=element_blank(),
        legend.position = 'top')+labs(fill='Algorithm')+ theme(plot.margin = unit(c(0,0,0,0),"cm"))+
  xlab(paste0('Chr',ch,' (kb)\n'))+
  scale_fill_manual(values = c('blue','orange','red','green'))

#facet_grid(CNV_type~.)

cnvrs_chr <- sum_loss%>% filter(seqnames==paste0('chr',ch)) 

plt2+ geom_vline(xintercept=cnvrs_chr$start/1000,color='gray',linetype='dashed',alpha=0.7,size=0.01)+
  geom_vline(xintercept=cnvrs_chr$end/1000,color='gray',linetype='dashed',alpha=0.7,size=0.01)

# How to process sum_loss table? ####
# calculate rowsums
sum_loss$total <-rowSums(sum_loss[,5:19],na.rm = T) 
# remove singletons
sum_loss <- sum_loss %>% filter(total!=N & total!=P & total!=Q & total!=I)


