# function to load rawCNV or gapmerged
load.raw_gap <- function(snp.ct, len, use.merged.cnv=F, only.common.sam=F, qc.pass.only=T) {
  library(GenomicRanges)
  library(plyranges)
  ##########################################################
  ### load raw CNVs or gap.merged
  #######################################################
  if (isTRUE(use.merged.cnv)) {
    # use merged cnvs
    gr.cll.cnv <- readRDS('/home/dean91/cnv_cll/concordance_ch3/Gap_0.2_merged.rds')
    #gr.cll.cnv <- unlist(GRangesList(gr.cll.cnv), use.names = F)
    #names(gr.cll.cnv) <- NULL
    
    # change names
   # names(mcols(gr.cll.cnv)) <- c("CN", "old_sample_id", "conf", "no_of_probes", "avgConf", 
  #                                "length", "loss_gain",  "caller", "study_centre", "sample_id")
  } else {
    
    # load all raw CNVs regardless of length or prob counts
    gr.cll.cnv <- read.csv('/home/dean91/cnv_cll/survival_mdr/All_cnvs_with_ox.csv') %>%  #& !study %in% c('Oxford','Oxford-ARC','Oxford-ADM')
      filter(!chr %in% c('X','Y','23')) %>%  toGRanges()
    
    # change names
    names(mcols(gr.cll.cnv)) <- c("CN", "old_sample_id", "conf",
                                  "no_of_probes", "avgConf", "length", "loss_gain",  "caller", "study_centre", "sample_id")
    
    seqlevelsStyle(gr.cll.cnv) <- 'UCSC'
    if (length(gr.cll.cnv) < 1) {stop("No CNAs left after algorithm filtering, check algo argument..")}
    
  }
  
  # filter
  gr.cll.cnv <- gr.cll.cnv %>% filter(no_of_probes > snp.ct)
  gr.cll.cnv <- gr.cll.cnv %>% filter(width(gr.cll.cnv) > len)
  
  ## harmonize sutdy names and sample ids
  gr.cll.cnv$study_centre <- gsub('^Bourn.*','Bournemouth', gr.cll.cnv$study_centre)
  gr.cll.cnv$study_centre <- gsub('^Hull.*','Hull', gr.cll.cnv$study_centre)
  gr.cll.cnv$study_centre <- gsub('^Newcastle.*','Newcastle', gr.cll.cnv$study_centre)  
  gr.cll.cnv$study_centre <- gsub('^Oxford.*','Oxford', gr.cll.cnv$study_centre) 
  gr.cll.cnv$study_centre <- gsub('^South.*','Southampton', gr.cll.cnv$study_centre) 
  gr.cll.cnv$sample_id <- gsub('1536_HRH-HOE', '1536_HRH_HOE', gr.cll.cnv$sample_id)
  
  if (isTRUE(only.common.sam)) {
  # only keep overlapping samples
  # Identify samples present in all 4 callers
  common_samples <- data.frame(gr.cll.cnv) %>% distinct(sample_id, caller) %>%
    dplyr::count(sample_id, name = "n_callers") %>%
    filter(n_callers == 4) %>% pull(sample_id)
  
  # Filter CNV table to only those samples
  gr.cll.cnv <- gr.cll.cnv %>% filter(sample_id %in% common_samples)
  }
  if (isTRUE(qc.pass.only)) {
    #load qc failed list
    qc.fail <-readLines('~/cnv_cll/concordance_ch3/sam_qc.fail')
    gr.cll.cnv <- gr.cll.cnv %>% filter(!sample_id %in% qc.fail)
  }
  return(gr.cll.cnv)
}

# new gap merging func
library(GenomicRanges)
library(data.table)

merge_by_gap_fraction_dt <- function(gr, max_gap_fraction = 0.2,
                                     sample_col = "sample_id", caller_col = "caller",
                                     type_col = "loss_gain", cn_col = "CN", 
                                     probe_col = "no_of_probes") {
  
  stopifnot(is(gr, "GRanges"))
  
  dt <- as.data.table(as.data.frame(gr))
  
  # Check required columns
  required_cols <- c(sample_col, "seqnames", "start", "end",
                     caller_col, type_col, cn_col)
  
  missing_cols <- setdiff(required_cols, names(dt))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ",
         paste(missing_cols, collapse = ", "))
  }
  
  group_cols <- c(sample_col, "seqnames", caller_col, type_col, cn_col)
  
  setorderv(dt, c(group_cols, "start", "end"))
  
  merged_dt <- dt[, {
    
    n <- .N
    merge_id <- integer(n)
    merge_id[1] <- 1L
    
    if (n > 1) {
      
      current_start <- start[1]
      current_end   <- end[1]
      current_id    <- 1L
      
      for (i in 2:n) {
        
        gap <- start[i] - current_end - 1L
        if (gap < 0) gap <- 0L
        
        current_width <- current_end - current_start + 1L
        next_width    <- end[i] - start[i] + 1L
        
        combined_length <- current_width + next_width + gap
        gap_fraction <- gap / combined_length
        
        if (gap_fraction < max_gap_fraction) {
          
          merge_id[i] <- current_id
          current_end <- max(current_end, end[i])
          
        } else {
          
          current_id <- current_id + 1L
          merge_id[i] <- current_id
          current_start <- start[i]
          current_end <- end[i]
        }
      }
    }
    
    cbind(.SD, merge_id = merge_id)
    
  }, by = group_cols]
  
  out_dt <- merged_dt[, {
    
    ans <- .SD[1]
    
    ans$start <- min(start)
    ans$end   <- max(end)
    ans$width <- ans$end - ans$start + 1L
    
    if (probe_col %in% names(.SD)) {
      ans[[probe_col]] <- sum(get(probe_col), na.rm = TRUE)
    }
    
    ans$n_merged_segments <- .N
    ans$merged_width_bp <- ans$width
    
    ans
    
  }, by = c(group_cols, "merge_id")]
  
  out_dt[, merge_id := NULL]
  
  gr_out <- GRanges(
    seqnames = out_dt$seqnames,
    ranges = IRanges(start = out_dt$start, end = out_dt$end),
    strand = if ("strand" %in% names(out_dt)) out_dt$strand else "*"
  )
  
  metadata_cols <- setdiff(
    names(out_dt),
    c("seqnames", "start", "end", "width", "strand")
  )
  
  mcols(gr_out) <- out_dt[, ..metadata_cols]
  
  gr_out
}

# get RO pairs
get_ro_pairs <- function(gr) {
  
  hits <- findOverlaps(gr, gr, ignore.strand = TRUE)
  
  hit_df <- data.frame(query = queryHits(hits),
    subject = subjectHits(hits)) %>%
    filter(query < subject)
  
  hit_df <- hit_df %>%
    mutate(min_probes = pmin(gr$no_of_probes[query], gr$no_of_probes[subject], na.rm = TRUE),
      max_probes = pmax(gr$no_of_probes[query], gr$no_of_probes[subject], na.rm = TRUE),
      mean_probes = rowMeans(cbind(gr$no_of_probes[query], gr$no_of_probes[subject]), na.rm = TRUE))
  
  hit_df <- hit_df %>%
    mutate(query_uid = gr$cnv_uid[query],
      subject_uid = gr$cnv_uid[subject],
      caller_query = gr$caller[query],
      caller_subject = gr$caller[subject],
      width_query = width(gr)[query],
      width_subject = width(gr)[subject]) %>%    
    filter(caller_query != caller_subject)
  
  ov <- pintersect(gr[hit_df$query],
    gr[hit_df$subject], ignore.strand = TRUE)
  
  hit_df %>%  mutate(overlap_bp = width(ov),
      reciprocal_overlap = pmin(
        overlap_bp / width_query,
        overlap_bp / width_subject),
      sample_id = unique(gr$sample_id),
      loss_gain = unique(gr$loss_gain))
}

make_consensus_events <- function(gr_sample_type) {
  
  # Merge overlapping CNVs into candidate consensus events
  consensus <- reduce(gr_sample_type, ignore.strand = TRUE)
  
  # Find which original CNV calls overlap each consensus event
  hits <- findOverlaps(consensus, gr_sample_type, ignore.strand = TRUE)
  
  # Caller support per consensus event
  caller_support <- tapply(
    subjectHits(hits),
    queryHits(hits),
    function(idx) {sort(unique(gr_sample_type$caller[idx]))
    }
  )
  
  # Add metadata
  consensus$n_callers <- lengths(caller_support)
  consensus$sample_id <- unique(gr_sample_type$sample_id)
  consensus$loss_gain <- unique(gr_sample_type$loss_gain)
  
  consensus$callers <- rep(NA_character_, length(consensus))
  consensus$n_callers <- rep(0L, length(consensus))
  consensus$callers[as.integer(names(caller_support))] <- 
    vapply(caller_support, paste, collapse = ";", FUN.VALUE = character(1))
  consensus$n_callers[as.integer(names(caller_support))] <- 
    lengths(caller_support)
  
  return(consensus)
}

plot_event_raw_cnvs <- function(event_id_to_plot, event_raw_df) {
  
  df <- event_raw_df %>%
    filter(event_id == event_id_to_plot) %>%
    mutate(
      caller = factor(caller),
      y = as.numeric(caller)
    )
  
  event_start <- unique(df$start_event)
  event_end   <- unique(df$end_event)
  
  ggplot(df) +
    geom_segment(
      aes(
        x = start_raw,
        xend = end_raw,
        y = y,
        yend = y
      ),
      linewidth = 3, col='gray',
    ) +
    geom_segment(
      aes(
        x = event_start,
        xend = event_end,
        y = max(y) + 1,
        yend = max(y) + 1
      ),
      linewidth = 4, col='red',
    ) +
    scale_y_continuous(
      breaks = c(sort(unique(df$y)), max(df$y) + 1),
      labels = c(levels(df$caller), "reduced_event")
    ) +
    theme_bw() +
    labs(
      x = "Genomic position",
      y = "Caller",
      title = paste("Reduced event", event_id_to_plot)
    )
}

# get RO clusters
get_ro_clusters <- function(df) {
  
  # If no reciprocal-overlap pairs exist
  if (nrow(df) == 0) return(NULL)
  
  # Build graph: CNVs are nodes, reciprocal-overlap matches are edges
  edges <- df %>%
    transmute(
      from = paste0(caller_query, "_", query),
      to   = paste0(caller_subject, "_", subject)
    )
  
  g <- graph_from_data_frame(edges, directed = FALSE)
  comp <- components(g)$membership
  
  cluster_df <- data.frame(
    node = names(comp),
    cluster_id = as.integer(comp)
  ) %>%
    mutate(
      caller = sub("_[^_]+$", "", node)
    )
  
  cluster_df %>%
    group_by(cluster_id) %>%
    summarise(
      n_callers = n_distinct(caller),
      callers = paste(sort(unique(caller)), collapse = ";"),
      n_cnvs =  dplyr::n(),
      .groups = "drop"
    )
}

# RO clusters for breakpoints
make_ro_clusters <- function(ro_pairs_filt) {
  
  edges <- ro_pairs_filt %>%
    dplyr::select(
      sample_id,
      loss_gain,
      from = query_uid,
      to = subject_uid
    )
  
  cluster_list <- edges %>%
    group_by(sample_id, loss_gain) %>%
    group_modify(~{
      
      g <- graph_from_data_frame(
        d = .x %>% dplyr::select(from, to),
        directed = FALSE
      )
      
      comp <- components(g)
      
      data.frame(
        cnv_uid = as.integer(names(comp$membership)),
        ro_cluster = as.integer(comp$membership)
      )
      
    }) %>%
    ungroup()
  
  cluster_list %>%
    mutate(
      ro_cluster_id = paste(sample_id, loss_gain, ro_cluster, sep = "_")
    )
}

# graphical params
theme_thesis <- function() {
  theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 11),
      axis.title = element_text(face = "bold"),
      strip.background = element_rect(fill = "grey95", colour = "black"),
      strip.text = element_text(face = "bold"),
      legend.position = "top",
      legend.title = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.2, colour = "grey90")
    )
}

cnv_cols <- c("Gain" = "#2C7BB6",
              "Loss" = "#D7191C")

caller_cols <- c("PennCNV"   = "#1B9E77",
                 "QuantiSNP" = "#D95F02",
                 "iPattern"  = "#7570B3",
                 "Nexus"     = "#E7298A")

### plot LRR/BAF
plot_LRR.BAF <- function(sample, chr, loss_gain){
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
  
  sam <- sample
  chr.plt <- paste0('chr',chr)
  type <- loss_gain
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
  #dir.create('thesis_out/LrrBaf_plts_missed4')
  png(filename = paste0('thesis_out/LrrBaf_plts_missed4/',sam,chr.plt,'.png'), width = 8, height = 5, res = 300, units = 'in')
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
  
}

plot_LRR.BAF.RO <- function(plt.dat, span=1, p.size=0.7, group){
  ################################
  ##kp plot coverage all sam
  ################################
  pp <- getDefaultPlotParams(plot.type= 1)
  pp$data1height <- 400
  pp$topmargin <- 20
  pp$ideogramheight <- 20

  ### change params for kp plot ####
  gr.cll.cnv <- load.raw_gap(snp.ct = 5, len = 5, use.merged.cnv = T)
  # other params
  point.size <- p.size
  # margins
  top_lrr_r1 <- 1
  top_lrr_r <- 0.72
  top_baf_r1 <- 0.66
  top_baf_r <- 0.38
  
  
  size <- plt.dat$median_call_width_bp
  n.clr <- plt.dat$n_callers
  clrs <- plt.dat$callers 
  sam <- plt.dat$sample_id
  chr.plt <- plt.dat$chr
  type <- plt.dat$loss_gain
  zoom.coord <- GRanges(paste0(chr.plt,':',plt.dat$min_start,'-',plt.dat$max_end))
  zoom.coord <- resize(zoom.coord, width = width(zoom.coord) + ifelse(size < 5e6, span*1e6, span*1e6*3), fix = "center")
  ### load raw CNVs or gap.merged
  gr.cll.cnv.sam <- gr.cll.cnv[gr.cll.cnv$sample_id==sam & gr.cll.cnv$loss_gain==type]
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
  dir.create(paste0('thesis_out/LrrBaf_plts/',group))
  png(filename = paste0('thesis_out/LrrBaf_plts/',group,'/',n.clr,'_',sam,'_',chr.plt,'_',clrs,'_.png'), width = 6, height = 4, res = 300, units = 'in')
  kp <- plotKaryotype(plot.type = 1, chromosomes = chr.plt, zoom = zoom.coord, plot.params = pp)#, zoom = target_region) #chromosomes = 'chr6'
  kpAddBaseNumbers(kp, tick.dist = ifelse(size < 5e6, 1e6, 10e6), tick.col="black", cex=0.7, add.units = T)
  try(kpPlotRegions(kp, data=gr.cll.cnv.caller[[1]], col='#E7298A', r0 =0, r1 = 0.08, avoid.overlapping = F))
  kpAddLabels(kp,labels = 'N',r0 =0, r1 = 0.08, cex=0.8, pos = 1, label.margin = 0.08)
  try(kpPlotRegions(kp, data=gr.cll.cnv.caller[[2]], col='#1B9E77', r0 =0.09, r1 = 0.17, avoid.overlapping = F))
  kpAddLabels(kp,labels = 'P', r0 =0.09, r1 = 0.17, cex=0.8, pos = 1,label.margin = 0.08)
  try(kpPlotRegions(kp, data=gr.cll.cnv.caller[[3]], col='#D95F02', r0 =0.18, r1 = 0.26, avoid.overlapping = F))
  kpAddLabels(kp,labels = 'Q', r0 =0.18, r1 = 0.26, cex=0.8, pos = 1,label.margin = 0.08)
  try(kpPlotRegions(kp, data=gr.cll.cnv.caller[[4]], col='#7570B3', r0 =0.27, r1 = 0.35, avoid.overlapping = F))
  kpAddLabels(kp,labels = 'I', r0 =0.27, r1 = 0.35, cex=0.8, pos = 1,label.margin = 0.08)
  
  # BAF
  kpPoints(kp, data=prb.file, y=prb.file$baf, r0=top_baf_r, r1=top_baf_r1, cex=point.size, data.panel = 1,col='#CDC0B0')#slategray3
  kpAxis(kp, ymax = 1, ymin = 0, r0 = top_baf_r, r1 = top_baf_r1, cex=0.7, data.panel = 1)
  kpAddLabels(kp,labels = 'BAF', r0 = top_baf_r, r1 = top_baf_r1, cex=0.8, srt=90, pos = 1, label.margin = 0.08)
  # LRR
  range(prb.file$lrr, na.rm = T)
  prb.zoom <- prb.file %>% filter_by_overlaps(zoom.coord)
  y.max <- ceiling(max(prb.zoom$lrr, na.rm = T))
  y.min <- floor(min(prb.zoom$lrr, na.rm = T))
  kpPoints(kp, data=prb.file, y=prb.file$lrr, r0=top_lrr_r, r1=top_lrr_r1,cex=point.size , col='#CDC0B0', 
           ymin = y.min, ymax = y.max)
  #kpPoints(kp, data=cnv, y=cnv$lrr, r0=top_lrr_r, r1=top_lrr_r1,cex=point.size , col='blue', ymin = y.min, ymax = y.max)
  kpAxis(kp,ymin = y.min, ymax = y.max, r0=top_lrr_r, r1=top_lrr_r1, cex=0.7)
  kpAddLabels(kp,labels = 'LRR',r0=top_lrr_r, r1=top_lrr_r1, cex=0.8,srt=90,pos = 1,label.margin = 0.08)
  dev.off()
  
}

plot_LRR.BAF.gene <- function(plt.dat, start.add=1, end.add=1, p.size=0.7, tick=1){
  ################################
  ##kp plot coverage all sam
  ################################
  pp <- getDefaultPlotParams(plot.type= 1)
  pp$data1height <- 400
  pp$topmargin <- 20
  pp$ideogramheight <- 20
  
  ### change params for kp plot ####
  gr.cll.cnv <- load.raw_gap(snp.ct = 5, len = 5, use.merged.cnv = T)
  # other params
  point.size <- p.size
  # margins
  top_lrr_r1 <- 1
  top_lrr_r <- 0.72
  top_baf_r1 <- 0.66
  top_baf_r <- 0.38
  
  size <- plt.dat$gene_width_bp
  clrs <- plt.dat$callers_overlapping_gene 
  sam <- plt.dat$sample_id
  chr.plt <- plt.dat$gene_chr
  type <- plt.dat$cnv_type
  zoom.coord <- GRanges(paste0(chr.plt,':',plt.dat$gene_start,'-',plt.dat$gene_end))
  start(zoom.coord) <- start(zoom.coord)-ifelse(size < 5e6, start.add*1e6, start.add*1e6*3) 
  end(zoom.coord) <-  end(zoom.coord)+ifelse(size < 5e6, end.add*1e6, end.add*1e6*3) 
  #zoom.coord <- resize(zoom.coord, width = width(zoom.coord) + ifelse(size < 5e6, span*1e6, span*1e6*3), fix = "center")
  ### load raw CNVs or gap.merged
  gr.cll.cnv.sam <- gr.cll.cnv[gr.cll.cnv$sample_id==sam & gr.cll.cnv$loss_gain==type]
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
  #dir.create('thesis_out/LrrBaf_plts_gene')
  png(filename = paste0('thesis_out/LrrBaf_plts_gene/',sam,'_',chr.plt,'_callers_',clrs,'.png'), width = 6, height = 4, res = 300, units = 'in')
  kp <- plotKaryotype(plot.type = 1, chromosomes = chr.plt, zoom = zoom.coord, plot.params = pp)#, zoom = target_region) #chromosomes = 'chr6'
  kpAddBaseNumbers(kp, tick.dist = tick*1e6, tick.col="black", cex=0.7, add.units = T)
  try(kpPlotRegions(kp, data=gr.cll.cnv.caller[[1]], col='#E7298A', r0 =0, r1 = 0.08, avoid.overlapping = F))
  kpAddLabels(kp,labels = 'N',r0 =0, r1 = 0.08, cex=0.8, pos = 1, label.margin = 0.08)
  try(kpPlotRegions(kp, data=gr.cll.cnv.caller[[2]], col='#1B9E77', r0 =0.09, r1 = 0.17, avoid.overlapping = F))
  kpAddLabels(kp,labels = 'P', r0 =0.09, r1 = 0.17, cex=0.8, pos = 1,label.margin = 0.08)
  try(kpPlotRegions(kp, data=gr.cll.cnv.caller[[3]], col='#D95F02', r0 =0.18, r1 = 0.26, avoid.overlapping = F))
  kpAddLabels(kp,labels = 'Q', r0 =0.18, r1 = 0.26, cex=0.8, pos = 1,label.margin = 0.08)
  try(kpPlotRegions(kp, data=gr.cll.cnv.caller[[4]], col='#7570B3', r0 =0.27, r1 = 0.35, avoid.overlapping = F))
  kpAddLabels(kp,labels = 'I', r0 =0.27, r1 = 0.35, cex=0.8, pos = 1,label.margin = 0.08)
  
  # BAF
  kpPoints(kp, data=prb.file, y=prb.file$baf, r0=top_baf_r, r1=top_baf_r1, cex=point.size, data.panel = 1,col='#CDC0B0')#slategray3
  kpAxis(kp, ymax = 1, ymin = 0, r0 = top_baf_r, r1 = top_baf_r1, cex=0.7, data.panel = 1)
  kpAddLabels(kp,labels = 'BAF', r0 = top_baf_r, r1 = top_baf_r1, cex=0.8, srt=90, pos = 1, label.margin = 0.08)
  # LRR
  range(prb.file$lrr, na.rm = T)
  prb.zoom <- prb.file %>% filter_by_overlaps(zoom.coord)
  y.max <- ceiling(max(prb.zoom$lrr, na.rm = T))
  y.min <- floor(min(prb.zoom$lrr, na.rm = T))
  kpPoints(kp, data=prb.file, y=prb.file$lrr, r0=top_lrr_r, r1=top_lrr_r1,cex=point.size , col='#CDC0B0', 
           ymin = y.min, ymax = y.max)
  #kpPoints(kp, data=cnv, y=cnv$lrr, r0=top_lrr_r, r1=top_lrr_r1,cex=point.size , col='blue', ymin = y.min, ymax = y.max)
  kpAxis(kp,ymin = y.min, ymax = y.max, r0=top_lrr_r, r1=top_lrr_r1, cex=0.7)
  kpAddLabels(kp,labels = 'LRR',r0=top_lrr_r, r1=top_lrr_r1, cex=0.8,srt=90,pos = 1,label.margin = 0.08)
  dev.off()
  
}

# new cox func to accept a matrix
cox_from_matrix <- function(cna_mat,
                            feature_metadata = NULL,
                            surv_type = "OS",
                            n.cores = 6) {
  
  library(survival)
  library(parallel)
  
  ## Filter rare features
 #keep <- Matrix::colSums(cna_mat) >= 5
 #cna_mat <- cna_mat[, keep, drop = FALSE]
  
 # if (!is.null(feature_metadata))
  #feature_metadata <- feature_metadata[keep, , drop = FALSE]
  
  ## Load outcome
  cll.outcome <- read.csv("../survival_mdr/All CLL outcome_spss.csv")
  #refine centre names
  cll.outcome$Study <- gsub('^Hull.*','Hull',cll.outcome$Study)
  cll.outcome$Study <- gsub('^Newcastle.*','Newcastle',cll.outcome$Study)
  cll.outcome$Study <- gsub('^Oxford.*','Oxford',cll.outcome$Study)
  
  ## Match samples
  sam.keep <- intersect(rownames(cna_mat), cll.outcome$Sample)
  
  cna_mat <- cna_mat[sam.keep, , drop = FALSE]
  cll.outcome <- cll.outcome[match(sam.keep, cll.outcome$Sample), ]
  
  stopifnot(all(rownames(cna_mat) == cll.outcome$Sample))
  
  ## Calculate sample support
  sample_support <- Matrix::colSums(cna_mat)
  
  ## Run Cox
  cox_results <- mclapply(seq_len(ncol(cna_mat)), function(i){
    
    x <- as.numeric(cna_mat[, i])
    
    if(length(unique(x)) < 2)
      return(NULL)
    
    fit <- tryCatch({
      
      if(surv_type == "OS"){
        coxph(Surv(OSdays, OSstatus) ~
                x +  Study + Array,
              data = cll.outcome)
        
      } else if(surv_type == "TTFT"){
        
        coxph(Surv(TTFTdays, TTFTstatus) ~
                x +  Study + Array,
              data = cll.outcome)
        
      } else {
        
        coxph(Surv(OS_PTdays, OS_PTstatus) ~
                x +  Study + Array,
              data = cll.outcome)
        
      }
      
    }, error = function(e) NULL)
    
    if(is.null(fit))
      return(NULL)
    
    s <- summary(fit)
    
    data.frame(
      feature = colnames(cna_mat)[i],
      HR = s$coefficients[1, "exp(coef)"],
      p = s$coefficients[1, "Pr(>|z|)"],
      cox_n_sam = sample_support[i]
    )
    
  }, mc.cores = n.cores)
  
  cox_results <- do.call(rbind, cox_results)
  cox_results$FDR <- p.adjust(cox_results$p, "fdr")
  
  # avoided merge
  idx <- match(cox_results$feature, rownames(feature_metadata))
  
  cox_results <- cbind(cox_results,
    feature_metadata[idx, , drop = FALSE]
  )
  
  ## Merge metadata
  # if(!is.null(feature_metadata)){
  #   
  #   feature_metadata$feature <- rownames(feature_metadata)
  #   
  #   cox_results <- merge(cox_results,
  #                        feature_metadata,
  #                        by = "feature",
  #                        all.x = TRUE,
  #                        sort = FALSE)
    
 # }
  
  return(cox_results)
  
}

# create mat from disjs
disjoint_matrix <- function(dj.filtered, harmonized.gr){
  
  samples <- unique(as.character(harmonized.gr$sample_id))
  
  hits <- findOverlaps(dj.filtered, harmonized.gr)
  
  # One row per sample-disjoint pair
  pairs <- unique(data.frame(
    sample = harmonized.gr$sample_id[subjectHits(hits)],
    disjoint = queryHits(hits)
  ))
  
  mat <- sparseMatrix(
    i = match(pairs$sample, samples),
    j = pairs$disjoint,
    x = 1,
    dims = c(length(samples), length(dj.filtered))
  )
  
  rownames(mat) <- samples
  
  colnames(mat) <- paste0(seqnames(dj.filtered), ":",
                          start(dj.filtered), "-",
                          end(dj.filtered))
  
  metadata <- data.frame(
  #  n_samples = dj.filtered$sample_support,
    mean_callers = dj.filtered$mean_callers,
    max_callers = dj.filtered$max_callers,
    n_samples = dj.filtered$sample_support,
    #mean_confidence = dj.filtered$mean_confidence,
    weighted_score = dj.filtered$weighted_score,
    #sum_confidence = dj.filtered$sum_confidence,
    row.names = colnames(mat)
  )
  
  list(matrix = mat,
       metadata = metadata)
  
}

# same for whole chr mat
whole_chr_matrix <- function(whole.gr, cnv_type) {
  
  library(Matrix)
  library(dplyr)
  
  ## One binary feature per chromosome gain/loss
  df <- data.frame(
    sample = whole.gr$sample_id,
    feature = paste0(as.character(seqnames(whole.gr)),
                     "_", whole.gr$loss_gain),
    loss_gain = whole.gr$loss_gain
  ) |>
    distinct()  %>% filter(loss_gain == cnv_type) %>% dplyr::select(-loss_gain)# remove duplicate caller records
#  feature <- paste0(as.character(seqnames(whole.gr)), "_", whole.gr$loss_gain)
  
  ## Sparse matrix
  mat <- sparseMatrix(
    i = match(df$sample, unique(df$sample)),
    j = match(df$feature, unique(df$feature)),
    x = 1,
    dims = c(length(unique(df$sample)),
             length(unique(df$feature)))
  )
  
  rownames(mat) <- unique(df$sample)
  colnames(mat) <- unique(df$feature)
  
  ## Metadata
  metadata <- df |>
    dplyr::count(feature, name = "n_samples")
  
  rownames(metadata) <- metadata$feature
  metadata$feature <- NULL
  
  list(matrix = mat,
    metadata = metadata
  )
}

# plot cox manhatton
##################
## plot mahnatton
##################
plot_cox_mh.manual <- function(gr.cox.in,
                        annotations = NULL,
                        p.col = "FDR",
                        sig.cutoff = 0.05){
  
  library(karyoploteR)
  library(GenomicRanges)
  
  ## -----------------------------
  ## Plot settings
  ## -----------------------------
  
  chr.num <- suppressWarnings(
    as.numeric(gsub("chr", "", as.character(seqnames(gr.cox.in))))
  )
  
  gr.cox.in$fill <- ifelse(chr.num %% 2 == 0,
                           "#3C8DAD",
                           "#F4A261")
  
  ## -----------------------------
  ## Create plot
  ## -----------------------------
  gr.cox.loss <- gr.cox.in %>% filter(cnv_type=='Loss')
  gr.cox.gain <- gr.cox.in %>% filter(cnv_type=='Gain')
  
  ymax.loss <- ceiling(max(abs(gr.cox.loss$plot_y)))
  ymax.gain <- ceiling(max(abs(gr.cox.gain$plot_y)))
  
  kp <- plotKaryotype(plot.type = 3, chromosomes = "autosomal", labels.plotter = NULL)
  kpAddChromosomeNames(kp, cex = 0.6, srt = 45)
  
  ## y axis
  kpAxis(kp, ymin = 0, data.panel = 2,
    ymax = ymax.loss, tick.pos = 0:ymax.loss, cex = 0.8)
  kpAxis(kp, ymin = 0, data.panel = 1,
         ymax = ymax.gain, tick.pos = 0:ymax.gain, cex = 0.8)
  
  ## significance line
  kpAbline(kp, h = -log10(sig.cutoff),
    ymin = 0,  ymax = ymax.gain, data.panel = 1, col = "red", lty = 3, lwd = 1.5)
  kpAbline(kp, h = -log10(sig.cutoff),
           ymin = 0,  ymax = ymax.loss, data.panel = 2, col = "red", lty = 3, lwd = 1.5)
  
  ## -----------------------------
  ## Manhattan points
  ## -----------------------------
  
  kpPoints(kp, data = gr.cox.loss,
    y = mcols(gr.cox.loss)[[p.col]],
    ymax = ymax.loss, pch = 21, data.panel = 2,
    bg = cnv_cols[[2]], #gr.cox.loss$fill
    col = "black", cex = 0.8, lwd = 0.2)
  
  kpPoints(kp, data = gr.cox.gain,
           y = mcols(gr.cox.gain)[[p.col]],
           ymax = ymax.gain, data.panel = 1,
           pch = 21, bg = cnv_cols[[1]],#gr.cox.gain$fill,
           col = "black",
           cex = 0.8, lwd = 0.2)
  
  ## -----------------------------
  ## Manual annotations
  ## -----------------------------
  
  if(!is.null(annotations)){
    
    annot.gr <- GRanges()
    annot.labels <- character()
    
    for(i in seq_len(nrow(annotations))){
      
      search.term <- annotations$search[i]
      label <- annotations$label[i]
      
      idx <- which(grepl(paste0("\\b", search.term, "\\b"),
          gr.cox.in$genes
        )
      )
      
      if(length(idx)==0)
        next
      
      ## choose most significant interval
      idx <- idx[which.min(mcols(gr.cox.in)[[p.col]][idx])]
      
      annot.gr <- c(annot.gr, gr.cox.in[idx])
      annot.labels <- c(annot.labels, label)
    }
    
    ## remove duplicated intervals
    keep <- !duplicated(paste0(
      seqnames(annot.gr), ":",
      start(annot.gr), "-",
      end(annot.gr)
    ))
    
    annot.gr <- annot.gr[keep]
    annot.labels <- annot.labels[keep]
    
    ## highlight points
    kpPoints(
      kp,
      data = annot.gr,
      y = -log10(mcols(annot.gr)[[p.col]]),
      ymax = ymax.loss,
      pch = 21,
      bg = "red",
      col = "black",
      cex = 1.4
    )
    
    ## labels
    kpText(kp,
      data = annot.gr,
      labels = annot.labels,
      y = -log10(mcols(annot.gr)[[p.col]]),
      ymax = ymax.loss,
      pos = 3,
      cex = 0.75
    )
  }
  
  invisible(kp)
}

plot_cox_mh.auto <- function(gr.cox.in, gr.lab, surv='OS'){
  
  library(karyoploteR)
  library(GenomicRanges) 
  
  gr.cox.out <- gr.cox.in[[surv]]
  df.lab <- gr.lab %>% filter(group_name == surv)
  # kp plot
  chr.num <- suppressWarnings(as.numeric(gsub("chr", "", as.character(seqnames(gr.cox.out)))))
  
  gr.cox.out$fill <- ifelse(chr.num %% 2 == 0,
                            "#3C8DAD", "#F4A261")

  ## -----------------------------
  ## Create plot
  ## -----------------------------
  gr.cox.loss <- gr.cox.out %>% filter(cnv_type=='Loss')
  gr.cox.gain <- gr.cox.out %>% filter(cnv_type=='Gain')
  
  ymax.loss <- ceiling(max(-log10(gr.cox.loss$FDR)))
  ymax.gain <- ceiling(max(-log10(gr.cox.gain$FDR)))+1
  
  kp <- plotKaryotype(plot.type = 3, chromosomes = "autosomal", labels.plotter = NULL)
  kpAddChromosomeNames(kp, cex = 0.6, srt = 45)
  
  ## y axis
  kpAxis(kp, ymin = 0, data.panel = 2,
         ymax = ymax.loss, tick.pos = 0:ymax.loss, cex = 0.8)
  kpAxis(kp, ymin = 0, data.panel = 1, 
         ymax = ymax.gain, tick.pos = 0:ymax.gain, cex = 0.8)
  
  ## significance line
  kpAbline(kp, h = -log10(0.05),
           ymin = 0,  ymax = ymax.gain, data.panel = 1, col = "red", lty = 3, lwd = 1.5)
  kpAbline(kp, h = -log10(0.05),
           ymin = 0,  ymax = ymax.loss, data.panel = 2, col = "red", lty = 3, lwd = 1.5)
  # mh plt
  kpPoints(kp, data = gr.cox.loss,
           y = -log10(gr.cox.loss$FDR),
           ymax = ymax.loss, pch = 21, data.panel = 2,
           bg = cnv_cols[[2]], #gr.cox.loss$fill
           col = "black", cex = 0.9, lwd = 0.08)
  

  # remove dups
  gr.lab.loss <- df.lab[df.lab$cnv_type=='Loss',]
 # df.lab.loss <- df.lab.loss[!duplicated(df.lab.loss$cytoband)]
  
  ## labels
  kpText(kp, data = gr.lab.loss,pos=1,
         labels = gr.lab.loss$cytoband, data.panel = 2,
         y = -log10(mcols(gr.lab.loss)$FDR),
         ymax = ymax.loss, 
         cex = 0.75)
  
  kpPoints(kp, data = gr.cox.gain, 
           y = -log10(gr.cox.gain$FDR),
           ymax = ymax.gain, data.panel = 1,
           pch = 21, bg = cnv_cols[[1]],#gr.cox.gain$fill, cnv_cols[[1]]
           col = "black",
           cex = 0.9, lwd = 0.08)
  
  gr.lab.gain <- df.lab[df.lab$cnv_type=='Gain',]
 # gr.lab.gain <- gr.lab.gain[!duplicated(gr.lab.gain$cytoband)]
  
  ## labels
  kpText(kp, data = gr.lab.gain, pos=3,
         labels = gr.lab.gain$cytoband, data.panel = 1,
         y = -log10(mcols(gr.lab.gain)$FDR),
         ymax = ymax.gain,
         cex = 0.75)
  
}


plot_cox_mh <- function(gr.cox.in) {
  #add chr col
  chr.col <- as.numeric(gsub('chr','', seqnames(gr.cox.in)))
  #chr.col <- seq(1:22)
  gr.cox.in$chr.col <- chr.col%%2
  table(gr.cox.in$chr.col)
  gr.cox.in$chr.col <- gsub('^1','#FFBD07AA', gr.cox.in$chr.col)
  gr.cox.in$chr.col <- gsub('^0','#00A6EDAA', gr.cox.in$chr.col)
  gr.cox.in$pos <- (end(gr.cox.in) + start(gr.cox.in))/2
  
  
  #png(paste(title,'.png'),width = 12,height = 7,units = 'in',res = 180)
  kp <- plotKaryotype(plot.type = 4, chromosomes = 'autosomal', labels.plotter = NULL)
  kpAddChromosomeNames(kp,srt=45,cex=0.6)
  kpPoints(kp, data=gr.cox.in, y = -log10(gr.cox.in$FDR), pch=21, col='black', cex = 0.8,
           lwd = 0.2, ymax = max(-log10(gr.cox.in$FDR),na.rm = T), bg=gr.cox.in$chr.col)
  
  #kpText(kp, data = gr.rep, labels = gr.rep$genes, y = -log10(gr.rep$FDR), ymax = max(-log10(gr.rep$FDR)),
  #       cex=0.55, col='black')

}

create.filter.disjs <- function(disjoints, cnv.gr) {
  #scoes
  caller_scores <- c(
    Nexus = 3,
    PennCNV = 2,
    QuantiSNP = 2,
    iPattern = 1 )
  # create one off disjoints
  cll_dj <- disjoints
  
  # create vectors
  mean_callers      <- numeric(length(cll_dj))
  max_callers       <- integer(length(cll_dj))
  n_samples_ge1     <- integer(length(cll_dj))
  n_samples_ge2     <- integer(length(cll_dj))
  n_samples_ge3     <- integer(length(cll_dj))
  n_samples_ge4     <- integer(length(cll_dj))
  
  mean_confidence   <- numeric(length(cll_dj))
  median_confidence <- numeric(length(cll_dj))
  max_confidence    <- numeric(length(cll_dj))
  sum_confidence    <- numeric(length(cll_dj))
  
  caller_combo_list <- vector("list", length(cll_dj))
  sample_conf_list  <- vector("list", length(cll_dj))
  
  # count overlaps
  hits <- findOverlaps(cll_dj, cnv.gr)
  
  ########################################
  ### caller support - per sample
  ######################################
  hit_split <- split(subjectHits(hits), queryHits(hits))
  idx <- as.integer(names(hit_split))
  
  for(i in seq_along(hit_split)) {
    
    cnv_idx <- hit_split[[i]]
    
    ## callers grouped by sample
    sample_callers <- split(
      as.character(gr.cll.cnv$caller[cnv_idx]),
      as.character(gr.cll.cnv$sample_id[cnv_idx])
    )
    
    ####################################################
    ## number of callers
    ####################################################
    caller_counts <- sapply(sample_callers,
                            function(x) length(unique(x)))
    
    mean_callers[idx[i]] <- mean(caller_counts)
    max_callers[idx[i]]  <- max(caller_counts)
    n_samples_ge1[idx[i]] <- sum(caller_counts == 1)
    n_samples_ge2[idx[i]] <- sum(caller_counts == 2)
    n_samples_ge3[idx[i]] <- sum(caller_counts == 3)
    n_samples_ge4[idx[i]] <- sum(caller_counts == 4)
    
    ####################################################
    ## caller combinations
    ####################################################
    
    sample_combos <- sapply(sample_callers, function(x)
      paste(sort(unique(x)), collapse=";"))
    
    caller_combo_list[[idx[i]]] <- sample_combos
    
    ####################################################
    ## empirical confidence score of each sample
    ####################################################
    
    sample_confidence <- sapply(sample_callers, function(x){
      
      sum(caller_scores[unique(x)])
      
    })
    
    sample_conf_list[[idx[i]]] <- sample_confidence
    
    ####################################################
    ## disjoint summaries
    ####################################################
    
    mean_confidence[idx[i]]   <- mean(sample_confidence)
    median_confidence[idx[i]] <- median(sample_confidence)
    max_confidence[idx[i]]    <- max(sample_confidence)
    
    ## recurrence-weighted confidence
    sum_confidence[idx[i]]    <- sum(sample_confidence)
    
  }
  
  ########################################
  ### add to disjoint object
  ########################################
  mcols(cll_dj)$mean_callers      <- mean_callers
  mcols(cll_dj)$max_callers       <- max_callers
  
  mcols(cll_dj)$n_samples_ge1     <- n_samples_ge1
  mcols(cll_dj)$n_samples_ge2     <- n_samples_ge2
  mcols(cll_dj)$n_samples_ge3     <- n_samples_ge3
  mcols(cll_dj)$n_samples_ge4     <- n_samples_ge4
  
  mcols(cll_dj)$mean_confidence   <- mean_confidence
  mcols(cll_dj)$median_confidence <- median_confidence
  mcols(cll_dj)$max_confidence    <- max_confidence
  mcols(cll_dj)$sum_confidence    <- sum_confidence
  
  mcols(cll_dj)$caller_combos     <- caller_combo_list
  mcols(cll_dj)$sample_confidence <- sample_conf_list
  
  ####################
  # count recurrence
  ####################
  sample_support <- tapply(subjectHits(hits),
                           queryHits(hits), function(x) {
                             length(unique(as.character(gr.cll.cnv$sample_id[x])))
                           })
  
  
  # attach to disjs
  mcols(cll_dj)$sample_support <- 0
  mcols(cll_dj)$sample_support[as.integer(names(sample_support))] <- sample_support
  # add weighted conf
  mcols(cll_dj)$weighted_score <- mcols(cll_dj)$mean_confidence * log2(mcols(cll_dj)$sample_support + 1)
  
  
  ########################
  ## rank disjoints
  #####################
  
  cll_dj <- cll_dj[order(mcols(cll_dj)$sum_confidence,
                         decreasing = TRUE)]
  return(cll_dj)
  
}
