# This is the latest and cleanest method
# Just use disjoints for both concordance and recurrence
library(dplyr)
library(GenomicRanges)
library(plyranges)
library(parallel)
library(karyoploteR)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)

#####################################################
## notes
# Keep
# All broad CNVs (>1 Mb)
# All whole chromosome events - separate and add later to cox
# Focal CNVs detected by all 4 callers - scoring 
# Focal CNVs detected by PennCNV + QuantiSNP + Nexus - scoring
#########################################################

# hg19 chromosome lengths
hg19_chr_lengths <- c(
  chr1 = 249250621, chr2 = 243199373,
  chr3 = 198022430, chr4 = 191154276,
  chr5 = 180915260, chr6 = 171115067,
  chr7 = 159138663, chr8 = 146364022,
  chr9 = 141213431, chr10 = 135534747,
  chr11 = 135006516, chr12 = 133851895,
  chr13 = 115169878, chr14 = 107349540,
  chr15 = 102531392, chr16 = 90354753,
  chr17 = 81195210, chr18 = 78077248,
  chr19 = 59128983, chr20 = 63025520,
  chr21 = 48129895, chr22 = 51304566,
  chrX = 155270560, chrY = 59373566
)

source('../concordance_ch3/functions_new.R')
source('../concordance_ch3/functions.R')
# Process one type (del or gain) at a time

# get genes
genes.gr <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene)

genes.gr$gene_name <- mapIds(org.Hs.eg.db, keys = genes.gr$gene_id,
  column = "SYMBOL",  keytype = "ENTREZID",
  multiVals = "first")

######################################
## score based on the combo weights ##
######################################
type <- 'Gain'
size <- 'Broad'
snp <- 50
len.ct <- 50

# load all raw CNVs regardless of length or prob counts
gr.cll.cnv <- load.raw_gap(snp.ct = snp, len = len.ct, use.merged.cnv = T, only.common.sam = F, qc.pass.only = F)

gr.cll.cnv$length <- width(gr.cll.cnv)
table(gr.cll.cnv$study_centre)

#############################
## separate whole chr events - nex is enough
#############################
# Fraction of chromosome covered
cov_frac <- width(gr.cll.cnv) / hg19_chr_lengths[as.character(seqnames(gr.cll.cnv))]
# Whole chromosome CNVs (≥90% of chromosome)
gr.whole_chr <- gr.cll.cnv[cov_frac >= 0.9]
table(seqnames(gr.whole_chr))
length(gr.whole_chr)
table(gr.whole_chr$caller)

# Remaining CNVs
gr.non.whole_chr <- gr.cll.cnv[cov_frac < 0.9]
#red.non.whole_chr.loss <- GenomicRanges::reduce(gr.non.whole_chr[gr.non.whole_chr$loss_gain=='Loss'])

#################################
## back to non.whole chr
################################
#filter
gr.cll.cnv <- gr.non.whole_chr %>% filter(loss_gain == type) # filter on broad vs focal & length >= 1e6
if (size=='Focal') {
  gr.cll.cnv.fil <- gr.cll.cnv %>% filter(length <= 1e6)
} else if (size =='Broad') {
  gr.cll.cnv.fil <- gr.cll.cnv %>% filter(length >= 1e6)
} else {
  gr.cll.cnv.fil <- gr.cll.cnv
}

range(gr.cll.cnv.fil$length)
range(width(gr.cll.cnv.fil))

##########################
# create and filter disjs
###########################
#disjs <- disjoin(gr.cll.cnv.fil[gr.cll.cnv.fil$caller=='Nexus'])
disjs <- disjoin(gr.cll.cnv.fil)
length(disjs)
disjs <- disjs %>% filter(width(disjs) > 10)
length(disjs)

# run disj
cll_dj <- create.filter.disjs(disjoints = disjs, cnv.gr = gr.cll.cnv)

# plot
hist(cll_dj$mean_confidence)
hist(cll_dj$weighted_score)
wgt_cutoff <- 3
cll_dj.fil <- cll_dj %>% filter(mean_callers >= 1.5) %>% # final cutoff = mean_conf 3 , mean callers = 1.5
  filter(sample_support > 10) %>% filter(mean_confidence >= wgt_cutoff)# mean callers > 2.1 lowest djs
length(cll_dj.fil)
#cll_dj.fil <- cll_dj %>% filter(weighted_score >= 15) %>% filter(sample_support > 10) 
length(cll_dj.fil)
length(reduce(cll_dj.fil))
range(cll_dj.fil$sample_support)

# plot
# crude kp plot
kp <- plotKaryotype(plot.type = 1, chromosomes = 'autosomal') #chromosomes = 'chr6'
kpPlotCoverage(kp, data=cll_dj.fil, show.0.cov = T, col=ifelse(type=='Loss','red','blue'))

################### 
## cox regression
###################
surv <- 'PT'
#create sparse matrix from whole chr only - add back to cox mat later
wc <- whole_chr_matrix(gr.whole_chr, cnv_type = type)
stopifnot(identical(sort(colnames(wc$matrix)),
                    sort(rownames(wc$metadata))))

# create disj matrix before cox
if (!surv == 'OS') {
  gr.non.whole_chr.cox <- gr.cll.cnv #%>% filter(!study_centre=='Oxford')
} else {
  gr.non.whole_chr.cox <- gr.cll.cnv
}

#gr.non.whole_chr.cox <- gr.cll.cnv %>% filter(!study_centre=='Oxford')

table(gr.non.whole_chr.cox$study_centre)
table(gr.cll.cnv$study_centre)
disj.mat <- disjoint_matrix(dj.filtered = cll_dj.fil, harmonized.gr = gr.non.whole_chr.cox)
range(disj.mat$metadata$n_samples)
# combine with whole chr mat
disj.mat$metadata$type <- "disjoint"
wc$metadata$type <- "whole_chr"

# align before comb
## Union of all samples
all_samples <- union(rownames(disj.mat$matrix), rownames(wc$matrix))

## Expand matrices (missing rows become all zeros)
dj_mat <- Matrix::Matrix(0, nrow = length(all_samples), ncol = ncol(disj.mat$matrix), sparse = TRUE)
wc_mat <- Matrix::Matrix(0, nrow = length(all_samples), ncol = ncol(wc$matrix), sparse = TRUE)

rownames(dj_mat) <- rownames(wc_mat) <- all_samples
colnames(dj_mat) <- colnames(disj.mat$matrix)
colnames(wc_mat) <- colnames(wc$matrix)

dj_mat[rownames(disj.mat$matrix), ] <- disj.mat$matrix
wc_mat[rownames(wc$matrix), ] <- wc$matrix

cna_mat <- cbind(dj_mat, wc_mat)

#combine metadata
common_cols <- intersect(names(disj.mat$metadata), names(wc$metadata))

feature_metadata <- rbind(disj.mat$metadata[, common_cols, drop = FALSE],
  wc$metadata[, common_cols, drop = FALSE])
stopifnot(all(colnames(cna_mat) %in% rownames(feature_metadata)))
table(feature_metadata$type, useNA="always")
range(feature_metadata$n_samples)
#save cna mat
saveRDS(cna_mat, paste0('CNA_mat_',size,'_WTscore_',wgt_cutoff,type,'.Rds'))

####################################
# Run cox on disj and whole.chr mat
###################################
cox.out <- cox_from_matrix(cna_mat = cna_mat, feature_metadata = feature_metadata, surv_type = surv, n.cores = 8)
#view(cox.out)
#view(data.frame(gr.cox.out))
range(cox.out$HR, na.rm = T)
range(cox.out$FDR, na.rm = T)

###################
# annotate genes 
#################
## Identify disjoint rows
table(cox.out$type, useNA = "always")
is_disjoint <- cox.out$type == "disjoint"

## Convert only disjoint features to GRanges
gr <- GRanges(cox.out$feature[is_disjoint])
## Gene annotation
gene_hits <- findOverlaps(gr, genes.gr)

gene_list <- split(mcols(genes.gr)$gene_name[subjectHits(gene_hits)],
  queryHits(gene_hits))

genes <- rep(NA_character_, nrow(cox.out))
genes[is_disjoint][as.integer(names(gene_list))] <-
  vapply(gene_list,
    function(x) paste(unique(x), collapse = ";"),
    character(1))
cox.out$genes <- genes

## Cytoband annotation
cyto_hits <- findOverlaps(gr, kp$cytoband)
cyto_list <- split(paste0(as.character(seqnames(kp$cytoband))[subjectHits(cyto_hits)],
    mcols(kp$cytoband)$name[subjectHits(cyto_hits)]),queryHits(cyto_hits))
cytobands <- rep(NA_character_, nrow(cox.out))
cytobands[is_disjoint][as.integer(names(cyto_list))] <-
  vapply(cyto_list,
    function(x) paste(unique(x), collapse = ";"),
    character(1))

cox.out$cytoband <- cytobands
#add cyoband to whole chr << complete
cox.out$cytoband[cox.out$type == "whole_chr"] <- cox.out$feature[cox.out$type == "whole_chr"]
View(cox.out)

##################
# save disj dat
################
nrow(cox.out)
saveRDS(cox.out, file = paste0('CLL_filtered_disj_cox_',size,'_WTscore_',wgt_cutoff,type,'_', surv,'.Rds'))

# quick mh
gr.cox.out <- cox.out[cox.out$type == "disjoint",] %>%
  tidyr::separate(feature, into = c("chr", "coords"), sep = ":") %>%
  tidyr::separate(coords, into = c("start", "end"), sep = "-") %>%
  mutate(start = as.numeric(start),
         end = as.numeric(end)
  ) %>%
  makeGRangesFromDataFrame(
    keep.extra.columns = TRUE,
    seqnames.field = "chr",
    start.field = "start",
    end.field = "end"
  )

# initial mh plot
plot_cox_mh(gr.cox.in = gr.cox.out)

