library(GenomicRanges)

##########################################
# Function:
# Merge adjacent CNV calls
##########################################
merge_by_gap_fraction <- function(gr, max_gap_fraction = 0.2,
                                  sample_col = "sample_id", ncores=1) {
  library(parallel)
  
  stopifnot(is(gr, "GRanges"))
  
  # split by sample
  gr_list <- split(gr, mcols(gr)[[sample_col]])
  
  merged_all <- mclapply(gr_list, function(gr_sample) {
    
    gr_sample <- sort(gr_sample)
    
    # split by chromosome
    chr_list <- split(gr_sample, seqnames(gr_sample))
    
    merged_chr <- lapply(chr_list, function(gr_chr) {
      
      if(length(gr_chr) <= 1)
        return(gr_chr)
      
      current <- gr_chr[1]
      out <- GRanges()
      
      for(i in 2:length(gr_chr)) {
        
        next_range <- gr_chr[i]
        
        gap <- start(next_range) - end(current) - 1
        
        if(gap < 0)
          gap <- 0
        
        combined_length <-
          width(current) +
          width(next_range) + gap
        
        gap_fraction <- gap / combined_length
        
        same_cn <- mcols(current)$CN == mcols(next_range)$CN
        
        if(gap_fraction < max_gap_fraction && same_cn) {
          
          probes = mcols(current)$numSNP + mcols(next_range)$numSNP
    
          current <- GRanges(seqnames = seqnames(current),
            ranges = IRanges(start = min(start(current), start(next_range)),
              end   = max(end(current), end(next_range))))
          
          
          mcols(current)[[sample_col]] <- mcols(gr_chr)[[sample_col]][1]
          #get numSNP
          mcols(current)$numSNP <- probes
          #mcols(current)$sample <- mcols(next_range)$sample
          mcols(current)$CN <- mcols(next_range)$CN
          mcols(current)$CNV_type <- mcols(next_range)$CNV_type
          mcols(current)$study <- mcols(next_range)$study
          mcols(current)$method <- mcols(next_range)$method
          mcols(current)$length <- paste(combined_length)
          
          
        } else {
          
          out <- c(out, current)
          current <- next_range
        }
      }
      
      c(out, current)
    })
    
    # IMPORTANT FIX
    unlist(GRangesList(merged_chr), use.names = FALSE)
  }, mc.cores=ncores)
  
  # FINAL COMBINED GRanges
  unlist(GRangesList(merged_all), use.names = FALSE)
}
# what to do if there are biallelic del and mono deletion? do not merge them together?

##########################################
# Function:
# Create consensus CNVs within sample
# Using reduced ranges
##########################################
harmonize_sample <- function(gr_sample) {
  # Reduce overlapping ranges
  red <- disjoin(gr_sample)
  #red <- disjoin(gr.cll.cnv)#use cohort disj instead of sample level disj
  
  # Count supporting callers
  overlaps <- findOverlaps(red, gr_sample)
  
  caller_support <- tapply(subjectHits(overlaps), queryHits(overlaps),
    function(x) {length(unique(gr_sample$method[x]))
      }
    )
  # Add support metadata
  mcols(red)$n_callers <- as.numeric(caller_support)
  mcols(red)$sample <- unique(gr_sample$new_sam_id)
  mcols(red)$cohort <- unique(mcols(gr_sample)$study)
  mcols(red)$type <- unique(mcols(gr_sample)$CNV_type)
  
  return(red)
}

############################################################
# CREATE SAMPLE x DISJOINT CNA MATRIX
############################################################
# Goal:
# rows    = samples
# columns = disjoint CNA regions
#
# values:
#   1 = sample has CNV overlapping region
#   0 = no overlap
############################################################
# ASSUMPTIONS
############################################################
#
# harmonized_gr :
#   harmonized CNVs with metadata column:
#       sample
#
# dj_filtered :
#   final disjoint CNA regions
#
############################################################
cna.mat.cox <- function(dj.filtered, harmonized.gr, surv_type='OS', min_samples=5, n.cores=8) {
  library(SparseArray)
  library(survival)
  
  # CREATE UNIQUE REGION IDS
  region_ids <- paste0(seqnames(dj.filtered), ":",
                       start(dj.filtered), "-",
                       end(dj.filtered))
  
  samples <- unique(as.character(harmonized.gr$sample_id))
  
  # FIND OVERLAPS
  hits <- findOverlaps(dj.filtered, harmonized.gr)
  if (length(hits) == 0)
    stop("No overlaps found between filtered disjoints and CNVs.")
  
  cna_sparse <- sparseMatrix(
    i = match(harmonized.gr$sample_id[subjectHits(hits)], samples),
    j = queryHits(hits),
    x = 1,
    dims = c(length(samples), length(dj.filtered))
  )
  
  rownames(cna_sparse) <- samples
  
  colnames(cna_sparse) <- paste0(
    seqnames(dj.filtered), ":",
    start(dj.filtered), "-",
    end(dj.filtered)
  )
  #########################
  ## Filter sparse matrix 
  ## and do cox regression
  ########################
  
  keep <- Matrix::colSums(cna_sparse) >= min_samples
  table(keep)
  cna_sparse <- cna_sparse[, keep]
  
  #load outcome data
  cll.outcome <- read.csv('../survival_mdr/All CLL outcome_spss.csv')
  #remove oxford
  #cll.outcome <- cll.outcome[!cll.outcome$Study %in% c('Oxford-ADM', 'Oxford-ARC'), ]
  table(cll.outcome$Study)
  setdiff(rownames(cna_sparse), cll.outcome$Sample)
  
  # keep only intersect
  sam.keep <- intersect(cll.outcome$Sample, rownames(cna_sparse))
  cna_sparse <- cna_sparse[rownames(cna_sparse) %in% sam.keep, ]
  cll.outcome <- cll.outcome[cll.outcome$Sample %in% sam.keep, ]
  # check
  setdiff(rownames(cna_sparse), cll.outcome$Sample)
  
  # align with outcome data
  cll.outcome <- cll.outcome[match(rownames(cna_sparse), cll.outcome$Sample),]
  all(rownames(cna_sparse) == cll.outcome$Sample) #verify
  
  ####################################
  ### Do survival analysis
  ####################################
  names(cll.outcome)
  cox_results <- mclapply(seq_len(ncol(cna_sparse)),
                          function(i) {x <- as.numeric(cna_sparse[, i])
                          
                          if(length(unique(x)) < 2) {
                            return(NULL) }
                          
                          if (surv_type=='OS') {
                          fit <- tryCatch(coxph(Surv(OSdays, OSstatus) ~ x + TreatmentCenter + Study + Array,
                                                data = cll.outcome),
                                          error = function(e) NULL) 
                          }
                          else if (surv_type=='TTFT') {
                            fit <- tryCatch(coxph(Surv(TTFTdays, TTFTstatus) ~ x + TreatmentCenter + Study + Array,
                                                  data = cll.outcome),
                                            error = function(e) NULL)
                            
                          } else {
                            fit <- tryCatch(coxph(Surv(OS_PTdays, OS_PTstatus) ~ x + TreatmentCenter + Study + Array,
                                                  data = cll.outcome),
                                            error = function(e) NULL)
                          }
                          
                          if(is.null(fit)) {
                            return(NULL) }
                          
                          s <- summary(fit)
                          
                          data.frame(region = colnames(cna_sparse)[i],
                                     HR = s$coefficients[1,"exp(coef)"],
                                     p = s$coefficients[1,"Pr(>|z|)"])
                          }, 
                          mc.cores = n.cores)
  
  
  cox_results <- do.call(rbind, cox_results)
  cox_results$FDR <- p.adjust(cox_results$p, method = "fdr")
  
  gr.cox_results <- GRanges(cox_results$region)
  seqlevelsStyle(gr.cox_results) <- "UCSC"
  
  idx <- match(cox_results$region,
               paste0(seqnames(dj.filtered), ":",
                      start(dj.filtered), "-",
                      end(dj.filtered)))
  
  mcols(gr.cox_results)$HR <- cox_results$HR
  mcols(gr.cox_results)$pval <- cox_results$p
  mcols(gr.cox_results)$FDR <- cox_results$FDR
  
  ## recurrence
  mcols(gr.cox_results)$n_samples <- dj.filtered$n_samples[idx]
  
  ## caller and sample support
  mcols(gr.cox_results)$mean_callers <- dj.filtered$mean_callers[idx]
  mcols(gr.cox_results)$max_callers  <- dj.filtered$max_callers[idx]
  mcols(gr.cox_results)$sample_support  <- dj.filtered$sample_support[idx]
  
  ## confidence metrics
  mcols(gr.cox_results)$mean_confidence   <- dj.filtered$mean_confidence[idx]
  mcols(gr.cox_results)$weighted_score  <- dj.filtered$weighted_score[idx]
  mcols(gr.cox_results)$sum_confidence    <- dj.filtered$sum_confidence[idx]
  
  return(gr.cox_results)
  
}

# easy summary plots
# create summary table
length.plt <- function(gr, bin.size=50) {
raw.cnv_df <- data.frame(width = width(gr),
                         log10_width = log10(width(gr)),
                         CN = mcols(gr)$CN,
                         method = mcols(gr)$method,
                         study = mcols(gr)$study)

# create length bins
raw.cnv_df <- raw.cnv_df %>%
  mutate(size_bin = cut(width, breaks = c(1, 1e3, 1e4, 1e5, 1e6, 1e7, Inf),
                        labels = c("<1kb","1-10kb", "10-100kb",
                                   "100kb-1Mb", "1-10Mb", ">10 Mb"),  right = FALSE))

#raw.cnv_df %>% ggplot(aes(log10_width, fill=method)) + geom_histogram(bins = bin.size, alpha=0.7, position='identity')
#size bins
raw.cnv_df %>% ggplot(aes(size_bin, fill=method)) + geom_bar(position='dodge')
}




