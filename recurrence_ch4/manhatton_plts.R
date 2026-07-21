# plot manhatton from cox data
#################
# plot manhatton
#################
source('../concordance_ch3/functions_new.R')
source('functions_feature_selection.R')

surv.lst <- c('OS','TTFT','PT')

df.cox.out <- list()
gr.cox.out <- GRangesList()
gr.top <- GRangesList()

for (sv in surv.lst){
  print(sv)
  #cox_all_features.Rds
  #load cox results
  df.cox.gain <- readRDS(paste0('CLL_filtered_disj_cox_Broad_WTscore_3Gain_',sv,'.Rds')) %>% mutate(cnv_type='Gain')
  df.cox.loss <- readRDS(paste0('CLL_filtered_disj_cox_Broad_WTscore_3Loss_',sv,'.Rds'))  %>% mutate(cnv_type='Loss')
  # comb
  df.cox.out[[sv]] <- rbind(df.cox.gain, df.cox.loss) 
  #names(df.cox.out[[sv]])
  print(table( df.cox.out[[sv]]$type))
  
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
  ## Get annot- top hits
  gr.top[[sv]] <- top_chr_hits(gr = gr.cox.out[[sv]])
  
}

# flatten gr list
df.top <- GenomicRanges::as.data.frame(gr.top)
table(df.top$seqnames)

#save annots for manual curation
write.csv(df.top, paste0('MH_annot_chr.csv'), quote = F, row.names = F)

#input manual  
df.top <- read.csv('MH_annot_chr_manual.csv')
gr.top <- GenomicRanges::makeGRangesFromDataFrame(df.top, keep.extra.columns = T)
table(gr.cox.out$OS$HR)
# auto manhatton
sv <- 'TTFT'
pdf(file = paste0('plots/',sv,'_MHplot.pdf'), width = 10, height = 7)
plot_cox_mh.auto(gr.cox.in  = gr.cox.out, gr.lab = gr.top, surv = sv)
dev.off()

# manual manhatton

annotations <- data.frame(
  search = c("LINC02732;FDX1",
             "LOC105374510", "SCAT8;PHF3;EYS",
             "SETD2","INTS6", "RB1"),
  label = c("ATM", "LOC105374510","EYS",
            "SETD2","DLEU1", "RB1"),stringsAsFactors = FALSE)

annotations <- data.frame(search = c("LOC105374618;LOC124901165",
                                     "LY6S","BAIAP2", "LOC105372698;MIR646HG",
                                     "CARD11", "SCAP"),  label = c("LOC105374618","LY6S","BAIAP2",
                                                                   "MIR646HG", "CARD11", "SCAP"),
                          stringsAsFactors = FALSE)

plot_cox_mh.manual(gr.cox.in = gr.cox.out,
                   annotations = annotations, p.col = 'plot_y',
                   sig.cutoff = 0.05)


####################################
### LRR BAF plots of features
##################################






