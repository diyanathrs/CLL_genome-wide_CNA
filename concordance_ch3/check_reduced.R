##################################################################
# check if reduced ranges create consensus that are too larger than its candidates
# if so there's over merging issue - one caller calling large segment
##################################################################
library(GenomicRanges)
library(dplyr)
library(tidyr)
library(ggplot2)


source('functions_new.R')

# Process one type (del or gain) at a time
type <- 'Loss'
use.merged.cnv <- T

if (isTRUE(use.merged.cnv)) {
  # use merged cnvs
  gr.cll.cnv <- readRDS('Gap_0.2_merged_Loss.rds')
  gr.cll.cnv <- unlist(GRangesList(gr.cll.cnv), use.names = F)
  names(gr.cll.cnv) <- NULL
  
  # change names
  names(mcols(gr.cll.cnv)) <- c("copy_number", "old_sample_id", "conf", "no_of_probes", "avgConf", 
                                "length", "loss_gain",  "caller", "study_centre", "sample_id")
} else {
    

# load all raw CNVs regardless of length or prob counts
gr.cll.cnv <- read.csv('../survival_mdr/All_cnvs_with_ox.csv') %>%  #& !study %in% c('Oxford','Oxford-ARC','Oxford-ADM')
  filter(!chr %in% c('X','Y','23')) %>% filter(CNV_type==type & numSNP > 5) %>% toGRanges()

# change names
names(mcols(gr.cll.cnv)) <- c("copy_number", "old_sample_id", "conf",
                              "no_of_probes", "avgConf", "length", "loss_gain",  "caller", "study_centre", "sample_id")

seqlevelsStyle(gr.cll.cnv) <- 'UCSC'
if (length(gr.cll.cnv) < 1) {stop("No CNAs left after algorithm filtering, check algo argument..")}

}

# Split by sample and CNV type
gr.sam.split <- split(gr.cll.cnv, gr.cll.cnv$sample_id)

# Build consensus events
consensus_list <- mclapply(gr.sam.split, make_consensus_events, mc.cores = 10)

# Combine back into one GRanges object
consensus_gr <- unlist(GRangesList(consensus_list), use.names = F)
range(consensus_gr$n_callers)

hits <- findOverlaps(consensus_gr, gr.cll.cnv, ignore.strand = TRUE)

event_raw_links <- data.frame(event_id = queryHits(hits),
                              raw_id   = subjectHits(hits))

event_info <- as.data.frame(consensus_gr) %>%
  mutate(event_id = row_number(),
         event_width_bp = width,
         event_width_kb = width / 1000) %>% 
  dplyr::select(event_id,
         seqnames,
         start,
         end,
         event_width_bp,
         event_width_kb,
         sample_id,
         loss_gain,
         n_callers,
         callers)


raw_info <- as.data.frame(gr.cll.cnv) %>%
  mutate(raw_id = row_number(),
    raw_width_bp = width,
    raw_width_kb = width / 1000) %>%
  dplyr::select(raw_id, seqnames,
    start,
    end,
    raw_width_bp,
    raw_width_kb,
    sample_id,
    loss_gain,
    caller)

# collapse samples
event_raw_df <- event_raw_links %>%
  left_join(event_info, by = "event_id") %>%
  left_join(raw_info, by = "raw_id", suffix = c("_event", "_raw")) %>%
  filter(sample_id_event == sample_id_raw,
         loss_gain_event == loss_gain_raw)

# summarise
event_size_qc <- event_raw_df %>%
  dplyr::group_by(event_id, seqnames_event,
    start_event, end_event, sample_id_event, 
    loss_gain_event, event_width_bp, event_width_kb, 
    n_callers, callers) %>% 
  dplyr::summarise(n_raw_cnvs = length(raw_id),
    n_unique_callers = dplyr::n_distinct(caller),
    min_raw_width_kb = min(raw_width_kb, na.rm = TRUE),
    median_raw_width_kb = median(raw_width_kb, na.rm = TRUE),
    max_raw_width_kb = max(raw_width_kb, na.rm = TRUE),
    event_to_median_raw_ratio =
      dplyr::first(event_width_kb) / median(raw_width_kb, na.rm = TRUE),
    event_to_max_raw_ratio =
      dplyr::first(event_width_kb) / max(raw_width_kb, na.rm = TRUE),
    raw_start_min = min(start_raw, na.rm = TRUE),
    raw_start_max = max(start_raw, na.rm = TRUE),
    raw_end_min = min(end_raw, na.rm = TRUE),
    raw_end_max = max(end_raw, na.rm = TRUE), .groups = "drop")

#save outputs
saveRDS(event_size_qc, 'CNV_size_summary.Rds')
saveRDS(event_raw_df, 'CNV_event_raw_df.Rds')

# free memory
rm(event_raw_links, hits, raw_info, event_info)
gc()
# inspect 
head(event_size_qc, 20)

suspicious_events <- event_size_qc %>%
  filter(n_raw_cnvs >= 3 | event_to_median_raw_ratio >= 5 |
      event_width_kb >= 5000) %>%
  arrange(desc(event_to_median_raw_ratio))

good_events <- event_size_qc %>%
  filter(n_raw_cnvs >= 4 & event_to_median_raw_ratio <= 5) %>%
  arrange(desc(event_to_median_raw_ratio))

# crude view
one_event <- suspicious_events$event_id[2]
event_raw_df %>%
  filter(event_id == one_event) %>% arrange(start_raw) %>%
  dplyr::select(event_id,
    sample_id_event,
    loss_gain_event,
    caller,
    seqnames_raw,
    start_raw,
    end_raw,
    raw_width_kb,
    event_width_kb)

# plot
plot_event_raw_cnvs(one_event, event_raw_df)

# red size vs median contributing cnv size
ggplot(event_size_qc, aes(x = median_raw_width_kb,
           y = event_width_kb)) +
  geom_point(alpha = 0.4) +
  scale_x_log10() +
  scale_y_log10() +
  theme_bw() +
  labs(x = "Median raw CNV size within reduced event (kb)",
    y = "Reduced event size (kb)")

ggplot(event_size_qc,
       aes(x = event_to_median_raw_ratio)) +
  geom_histogram(bins = 50) +
  scale_x_log10() +
  theme_bw() +
  labs(x = "Reduced event size / median raw CNV size",
    y = "Number of reduced events")
