# All functions are here

# Count samples function will create disjoints from a all cnvs when given as a 
# gr object. Then it will count how many samples overlap each disjoint and return 
# the counts as a table
count_samples <- function(){
  # return disjs as granges list
  gr.cnv.disjoints <- disjoin(gr.cll.cnv)
  # create out columns from the disjoints
  #out_table <- out_table[,-5]
  #out_table$seqnames <- as.numeric(out_table$seqnames)
  #losses first
  gr_cnv_nex <- gr.cll.cnv %>% filter(method=='Nexus') %>% toGRanges()
  gr_cnv_pen <- gr.cll.cnv %>% filter(method=='PennCNV') %>% toGRanges()
  gr_cnv_qsnp <- gr.cll.cnv %>% filter(method=='QuantiSNP') %>% toGRanges()
  gr_cnv_ipn <- gr.cll.cnv %>% filter(method=='iPattern') %>% toGRanges()
  #calculate overlaps
  nex <- count_overlaps(gr.cnv.disjoints,gr_cnv_nex)
  pen <- count_overlaps(gr.cnv.disjoints,gr_cnv_pen)
  qsnp <- count_overlaps(gr.cnv.disjoints,gr_cnv_qsnp)
  ipn <- count_overlaps(gr.cnv.disjoints,gr_cnv_ipn)
  # add ovlps to disjoints object
  # add disjoint id
  gr.cnv.disjoints$pen <- pen
  gr.cnv.disjoints$qsnp <- qsnp
  gr.cnv.disjoints$ipn <- ipn
  gr.cnv.disjoints$nex <- nex

  return(gr.cnv.disjoints)
}
# Sample from range function to get sample list for a given range
# Return full sample list with 1s and 0s
samples_from_range <- function(i) {
  outcomes <- cll.outcome
  ovlps <- find_overlaps(gr.cll.cnv,range(ranges.out[i]))#filter_by_overlaps
  dat.ovlps <- data.frame(Sample=unique(ovlps$new_sam_id),alteration_all=1)
  outcomes <- left_join(outcomes,dat.ovlps,by='Sample')
  outcomes$alteration_all[is.na(outcomes$alteration_all)] <- 0
  # adding nexus
  ovlps_nex <- ovlps %>% filter(method=='Nexus')
  if(length(ovlps_nex) > 1) {
    dat.ovlps.nex <- data.frame(Sample=unique(ovlps_nex$new_sam_id),alteration_nex=1)
    outcomes <- left_join(outcomes,dat.ovlps.nex,by='Sample')
    outcomes$alteration_nex[is.na(outcomes$alteration_nex)] <- 0  } 
  return(outcomes)
}


# fun to get sample ids from particular algo
samples_from_algo <- function(i) {
  ovlps <- find_overlaps(gr.cll.cnv,range(ranges.out[i]))#filter_by_overlaps
  test <- data.frame(ovlps) %>% group_by(method) %>% reframe(new_sam_id)
  
}

# what we have as outputs now :- 1) sample counts for each algo
# 2) sample ids for all algos
# 3) sample ids for nexus
# 4) outcome data

# calculate survival for a given range - output a list with 6 vals 
#(os_nex,ttfs_nex,os_all,ttfs_all etc.)
# Same funs as lapply compatible ##
# assume lapply inputs ranges.out[1]
cal_survival <- function(i) {
  outcomes <<- samples_from_range(i)
  fitos <- survfit(Surv(OSdays/365.25, OSstatus) ~alteration_all, data = outcomes)
  hr.os <- summary(coxph(Surv(OSdays/365.25, OSstatus) ~alteration_all, data = outcomes))
  
  fitttfs <- survfit(Surv(TTFTdays/365.25, TTFTstatus) ~alteration_all, data = outcomes)
  fitpt <- survfit(Surv(OS_PTdays/365.25, OS_PTstatus) ~alteration_all, data = outcomes)
  # Get a gene list for range
  range.genes <- filter_by_overlaps(genes.data$genes,ranges.out[i])
  genes.list <- (range.genes$name)
  ran <- data.frame(ranges.out[i])
  ran <- ran[-5]
  ran$total <- table(outcomes$alteration_all)[2]
  #ran$nexus <- table(outcomes$alteration_nex)[2]
  ran$os.p <- surv_pvalue(fitos)$pval
  ran$ttfs.p <- surv_pvalue(fitttfs)$pval
  ran$pt.p <- surv_pvalue(fitpt)$pval
  #for nex
  ran$os.p.nex <- NA
  ran$ttfs.p.nex <- NA
  ran$pt.p.nex <- NA
  if (ncol(outcomes) > 49) {
    fitos.n <- survfit(Surv(OSdays/365.25, OSstatus) ~alteration_nex, data = outcomes)
    fitttfs.n <- survfit(Surv(TTFTdays/365.25, TTFTstatus) ~alteration_nex, data = outcomes)
    fitpt.n <- survfit(Surv(OS_PTdays/365.25, OS_PTstatus) ~alteration_nex, data = outcomes)
    ran$os.p.nex <- surv_pvalue(fitos.n)$pval
    ran$ttfs.p.nex <- surv_pvalue(fitttfs.n)$pval
    ran$pt.p.nex <- surv_pvalue(fitpt.n)$pval
    }
  ran$genes <- paste(genes.list,collapse = " ")
  coxph(Surv(OSdays/365.25, OSstatus) ~alteration_all, data = outcomes)
  return(ran)
}

# get survival curves for a disj
plot_survival <- function(i,input) {
  outcomes <<- samples_from_range(i)
  fitos <- survfit(Surv(OSdays/365.25, OSstatus) ~alteration_all, data = outcomes)
  fitttfs <- survfit(Surv(TTFTdays/365.25, TTFTstatus) ~alteration_all, data = outcomes)
  fitpt <- survfit(Surv(OS_PTdays/365.25, OS_PTstatus) ~alteration_all, data = outcomes)
  fit <- list(All_os=fitos,All_ttfs=fitttfs,All_pt=fitpt)
  # for nex
  fitos.n <- survfit(Surv(OSdays/365.25, OSstatus) ~alteration_nex, data = outcomes)
  fitttfs.n <- survfit(Surv(TTFTdays/365.25, TTFTstatus) ~alteration_nex, data = outcomes)
  fitpt.n <- survfit(Surv(OS_PTdays/365.25, OS_PTstatus) ~alteration_nex, data = outcomes)
  fit.n <- list(Nexus_os=fitos.n,Nexus_ttfs=fitttfs.n,Nexus_pt=fitpt.n)
  
  plot <- ggsurvplot(fit, pval = TRUE, conf.int = F,pval.coord = c(30, .95),
             risk.table = TRUE, # Add risk table
             risk.table.col = "strata", # Change risk table color by groups
             linetype = "strata",# Change line type by groups
             ggtheme = theme_bw(), # Change ggplot2 theme
             palette = c("#E7B800", "#2E9FDF"))
  plot_n <- ggsurvplot(fit.n, pval = TRUE, conf.int = F,pval.coord = c(20, .95),
             risk.table = TRUE, # Add risk table
             risk.table.col = "strata", # Change risk table color by groups
             linetype = "strata", # Change line type by groups
             ggtheme = theme_bw(), # Change ggplot2 theme
             palette = c("#E7B800", "#2E9FDF"))
  if (input=='nex') {
    return(plot_n)
  } else {return(plot)}
}
# function to plot disj when given a disjoint id

# IDEA2 - Reduce adjacent ranges and take the lowest p-val out 
# go into each reduced range and look for overlaps from del_pvals.ttfs and get back the sig
reduce_pval.list <- function(df,pval) { # this gives one p val per reduced range, then filter
  out.list <- list()
  ranges.reduced <- GenomicRanges::reduce(makeGRangesFromDataFrame(df))
  
  for (i in seq_along(ranges.reduced)) {
    t <- filter_by_overlaps(toGRanges(df),ranges.reduced[i])
    t <- t %>% filter(nex > 15)
    if (length(t) > 1) {
    # choose the lowest p after ordering but need to pick events with more than a certain sample no
    # or else they will be removed in the next filtering step
    if (pval=='os') {
      t <- t[order(t$os.p.nex),] # order and pick the lowest p
      out.list <- append(out.list,list(data.frame(t[1,]))) 
      
    } else if (pval=='ttfs') {
      t <- t[order(t$ttfs.p.nex),] # order and pick the lowest p
      out.list <- append(out.list,list(data.frame(t[1,])))
      
    } else {
      t <- t[order(t$pt.p.nex),] # order and pick the lowest p
      out.list <- append(out.list,list(data.frame(t[1,])))
    }
      } else {next}
  }
  out <- as.data.frame(do.call(rbind,out.list))
  out <- out[c(-5,-6)]
  return(out)
}




