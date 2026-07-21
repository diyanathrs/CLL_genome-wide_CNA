
library(GenomicRanges)
library(IRanges)
library(dplyr)

library(GenomicRanges)

reduce_cox_hits.old <- function(gr, max_gap = 0,
                            priority = c("priority","FDR","support")){
  
  priority <- match.arg(priority)
  
  stopifnot(is(gr, "GRanges"))
  stopifnot("cnv_type" %in% names(mcols(gr)))
  
  mcols(gr)$priority_score <-
    abs(log(mcols(gr)$HR)) *
    log2(mcols(gr)$n_samples + 1) *
    -log10(mcols(gr)$FDR + 1e-300)
  
  ## build clusters
  cl <- reduce(gr, ignore.strand = TRUE,
               min.gapwidth = max_gap + 1)
  
  hits <- findOverlaps(gr, cl)
  
  cluster.id <- subjectHits(hits)
  
  reps <- vector("list", length(cl))
  
  for(i in seq_along(cl)){
    
    idx <- queryHits(hits)[cluster.id == i]
    
    x <- gr[idx]
    
    if(priority=="priority"){
      
      ord <- order(
        -mcols(x)$priority_score,
        mcols(x)$FDR,
        -mcols(x)$n_samples
      )
      
    } else if(priority=="FDR"){
      
      ord <- order(
        mcols(x)$FDR
      #  -mcols(x)$n_samples,
      #  -mcols(x)$priority_score
      )
      
    } else{
      
      ord <- order(
        -mcols(x)$n_samples,
        mcols(x)$FDR,
        -mcols(x)$priority_score
      )
      
    }
    
    best <- x[ord[1]]
    
    mcols(best)$cluster_start <- start(cl[i])
    mcols(best)$cluster_end <- end(cl[i])
    mcols(best)$cluster_width <- width(cl[i])
    mcols(best)$cluster_size <- length(idx)
    
    reps[[i]] <- best
    
  }
  
  do.call(c, reps)
  
}

reduce_cox_hits <- function(gr,  max_gap = 0,
                            priority = c("priority", "FDR", "support")) {
  
  priority <- match.arg(priority)
  
  stopifnot(is(gr, "GRanges"))
  stopifnot("cnv_type" %in% names(mcols(gr)))
  
  ## Priority score
  mcols(gr)$priority_score <-
    abs(log(mcols(gr)$HR)) *
    log2(mcols(gr)$n_samples + 1) *
    -log10(mcols(gr)$FDR + 1e-300)
  
  ## Split by chromosome and CNV type
  groups <- split(gr, paste0(seqnames(gr), "_", mcols(gr)$cnv_type))
  
  reps <- lapply(groups, function(g) {
    
    ## Reduce only within this chromosome/CNV type
    cl <- reduce(g,  ignore.strand = TRUE,
                 min.gapwidth = max_gap + 1)
    
    hits <- findOverlaps(g, cl)
    
    cluster.id <- subjectHits(hits)
    
    out <- vector("list", length(cl))
    
    for(i in seq_along(cl)) {
      
      idx <- queryHits(hits)[cluster.id == i]
      
      x <- g[idx]
      
      if(priority == "priority") {
        
        ord <- order(
          -mcols(x)$priority_score,
          mcols(x)$FDR,
          -mcols(x)$n_samples
        )
        
      } else if(priority == "FDR") {
        
        ord <- order(
          mcols(x)$FDR,
          -mcols(x)$n_samples,
          -mcols(x)$priority_score
        )
        
      } else {
        
        ord <- order(
          -mcols(x)$n_samples,
          mcols(x)$FDR,
          -mcols(x)$priority_score
        )
        
      }
      
      best <- x[ord[1]]
      
      mcols(best)$cluster_start <- start(cl[i])
      mcols(best)$cluster_end <- end(cl[i])
      mcols(best)$cluster_width <- width(cl[i])
      mcols(best)$cluster_size <- length(idx)
      
      out[[i]] <- best
    }
    
    #do.call(c, out)
    
    unlist(GRangesList(out), use.names = F)
    
  })
  
  #do.call(c, reps)
  unlist(GRangesList(reps), use.names = F)
}


plot_km_feature <- function(df,
                            feature,
                            time = "OSdays",
                            status = "OSstatus",
                            feature_labels = c("No CNA", "CNA"),
                            palette = c("#1B9E77", "#D95F02"),
                            title = NULL, y.lab.cust = 'Overall survival',
                            show_hr = TRUE,
                            conf.int = FALSE) {
  
  library(survival)
  library(survminer)
  library(ggplot2)
  
  ## Check inputs
  stopifnot(feature %in% names(df))
  stopifnot(time %in% names(df))
  stopifnot(status %in% names(df))
  
  ## Extract required columns
  dat <- df[, c(time, status, feature, "Study")]
  dat <- dat[complete.cases(dat), ]
  #change time to years
  dat[time] <- dat[time]/365.25
  
  ## Create a temporary feature column
  dat$Feature <- dat[[feature]]
  
  ## Convert to factor if binary
  if (length(unique(dat$Feature)) == 2) {
    
    dat$Feature <- factor(
      dat$Feature,
      levels = c(0, 1),
      labels = feature_labels
    )
    
  } else {
    
    dat$Feature <- factor(dat$Feature)
    
  }
  
  ## Kaplan-Meier fit
  surv.formula <- as.formula(paste0("Surv(", time, ", ", status, ") ~ Feature"))
  
  fit <- do.call(survfit, list(formula = surv.formula, data = dat))
  
  fit$call$formula
  
  #fit <- survfit(Surv(dat[[time]], dat[[status]]) ~ Feature, data = dat)
  
  ## Cox model
  hr_text <- NULL
  
  if (show_hr) {
    
    rhs <- paste(c("Feature", "Study"),
                 collapse = " + ")
    
    cox.formula <- as.formula(
      paste0(
        "Surv(",
        time,
        ",",
        status,
        ") ~ ",
        rhs))
    
    cox.fit <- coxph(cox.formula, data = dat)
    
    # coxfit <- coxph(
    #   Surv(dat[[time]], dat[[status]]), ~ Feature,
    #   data = dat)
    # 
    s <- summary(cox.fit)
    
    HR <- round(s$coefficients[1, "exp(coef)"], 2)
    LCL <- round(s$conf.int[1, "lower .95"], 2)
    UCL <- round(s$conf.int[1, "upper .95"], 2)
    P <- signif(s$coefficients[1, "Pr(>|z|)"], 3)
    
    hr_text <- sprintf("HR = %.2f\n95%% CI %.2f–%.2f\np = %s",
      HR, LCL, UCL, P
    )
  }
  
  if (is.null(title))
    title <- feature
  
  p <- ggsurvplot(
    fit,
    data = dat,
    risk.table = TRUE,
    pval = TRUE,
    conf.int = conf.int,
    legend.title = "",
    legend.labs = levels(dat$Feature),
    palette = palette,
    title = title,
    xlab = "Time (years)",
    ylab = y.lab.cust,
    ggtheme = theme_bw(),
    risk.table.height = 0.25)
  
  if (show_hr) {
    
    p$plot <- p$plot +
      annotate(
        "text",
        x = Inf,
        y = 0.15,
        hjust = 1.05,
        label = hr_text,
        size = 4
      )
    
  }
  
  return(p)
}

# get top CHR hits for annot
top_chr_hits <- function(gr,
                         by = c("FDR", "p"),
                         ignore_whole_chr = TRUE,
                         keep_ties = FALSE) {
  
  stopifnot(inherits(gr, "GRanges"))
  by <- match.arg(by)
  
  # only sig
  gr <- gr %>% filter(FDR < 0.01)
  
  ## Remove whole chromosome events if requested
  if (ignore_whole_chr && "type" %in% names(mcols(gr))) {
    gr <- gr[mcols(gr)$type != "whole_chr"]
  }
  
  ## Remove missing scores
  score <- mcols(gr)[[by]]
  gr <- gr[!is.na(score)]
  score <- mcols(gr)[[by]]
  
  ## Group by chromosome and CNV type
  group <- paste(seqnames(gr), mcols(gr)$cnv_type, sep = "_")
  
  idx <- unlist(
    lapply(split(seq_along(gr), group), function(i) {
      
      best <- which.min(score[i])
      
      if (keep_ties) {
        i[score[i] == score[i][best]]
      } else {
        i[best]
      }
      
    }),
    use.names = FALSE
  )
  
  gr.out <- gr[idx]
  
  gr.out
}


elastic_net_cox <- function(cna.mat,
                            surv.dat,
                            alpha = 0.5,
                            nfolds = 10,
                            seed = 1234) {
  
  # harmonize samples
  missing.sam <- surv.dat$Sample[is.na(surv.dat$OSdays)]
  surv.dat <- surv.dat[!surv.dat$Sample %in% missing.sam,]
  sam <- intersect(rownames(cna.mat), surv.dat$Sample)
  # keep
  cna.mat <- as.matrix(cna.mat[sam,])
  surv.dat <- surv.dat[surv.dat$Sample %in% sam,]
  
  stopifnot(nrow(cna.mat) == length(surv.dat$OSstatus))
  #stopifnot(length(surv.time) == length(surv.status))
  
  sum(is.na(surv.dat$OSdays))
  
  y <- Surv(event = surv.dat$OSstatus, time = surv.dat$OSdays)
  
  set.seed(seed)
  
  cv.fit <- cv.glmnet(
    x = cna.mat,
    y = y,
    family = "cox",
    alpha = alpha,
    nfolds = nfolds,
    standardize = FALSE
  )
  
  fit <- glmnet(
    x = cna.mat,
    y = y,
    family = "cox",
    alpha = alpha,
    lambda = cv.fit$lambda.min,
    standardize = FALSE
  )
  
  coef.df <- as.matrix(coef(fit))
  
  coef.df <- data.frame(
    feature = rownames(coef.df),
    beta = coef.df[,1],
    stringsAsFactors = FALSE
  )
  
  coef.df <- subset(coef.df, beta != 0)
  
  coef.df$HR <- exp(coef.df$beta)
  
  coef.df <- coef.df[
    order(abs(coef.df$beta), decreasing = TRUE),
  ]
  
  list(
    cv = cv.fit,
    fit = fit,
    features = coef.df
  )
}