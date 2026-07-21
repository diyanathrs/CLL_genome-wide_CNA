# Do initial QC including PCA
# load sig intensities

library(data.table)
library(dplyr)
library(ggplot2)
library(parallel)

sig.files <- list.files('../sig_intensities', pattern = '.txt', full.names = T)
sam.ids <- gsub('.txt','',sig.files)
sam.ids <- gsub('../sig_intensities/','',sam.ids)
other.sam <- data.frame(old.name=sam.ids, V2=sam.ids)
#remove ox
other.sam <- other.sam[-c(1:27,637:870),]

# check samIds against dataset
samples <- read.csv('../survival_mdr/All CLL outcome_spss.csv')

# separate ox ARC
ox.arc.adm <- data.frame(old.name=sam.ids[c(1:27,637:870)])
ox.arc.adm$V2 <- sub("(_[^_]*){2}$", "", ox.arc.adm$old.name)
ox.arc.adm$V2 <- gsub("_R$", "", ox.arc.adm$V2)
ox.arc.adm$V2[grep("^[0-9]", ox.arc.adm$V2)] <- paste0('ADM',ox.arc.adm$V2[grep("^[0-9]", ox.arc.adm$V2)])
ox.arc.adm$V2 <- gsub('_','',ox.arc.adm$V2)

#check
setdiff(ox.arc.adm$V2, samples$Sample)
# final metadata
metadata <- rbind(other.sam, ox.arc.adm)
names(metadata) <- c('Sample_ID', 'new_sampleID')
metadata <- right_join(metadata, samples, by= c('new_sampleID' = 'Sample'))
metadata <- metadata %>% mutate(centre = case_when(TreatmentCenter == "L" ~ "Leicester",
                                                   TreatmentCenter == "B" ~ "Bournemouth",
                                                   TreatmentCenter == "C" ~ "Cardiff",
                                                   TreatmentCenter == "H" ~ "Hull",
                                                   TreatmentCenter == "N" ~ "Newcastle",
                                                   TreatmentCenter == "Oc" ~ "Oxford",
                                                   TreatmentCenter == "Od" ~ "Oxford",
                                                   TreatmentCenter == "S" ~ "Southampton"))
# add genotyping array
metadata <- metadata %>% mutate(array = case_when(Array == "A" ~ "InfiniumOEE-8v1.3",
                                                   Array == "B" ~ "HumanOEE-8v1",
                                                   Array == "C" ~ "GSA-24v1-0",
                                                   Array == "D" ~ "InfiniumOEE-8v1.4",
                                                   Array == "E" ~ "HumanOEE-8v1.2",
                                                   Array == "F" ~ "GSA-24v1-0"))

#check
setdiff(metadata$new_sampleID, samples$Sample)
setdiff(samples$Sample, metadata$new_sampleID)
setdiff(metadata$Sample_ID, sam.ids)

####################
### load sig intensitite
##################
ref <- fread(sig.files[1])

marker_set <- ref %>%
  filter(Chr %in% as.character(1:22)) %>%
  arrange(as.numeric(Chr), Position) %>%
  group_by(Chr) %>%
  slice(seq(1, n(), by = 50)) %>%   # keep every 100th marker
  ungroup() %>%
  pull(Name)

length(marker_set)

# extract pre-defined markers
read_lrr_sample <- function(file, marker_set) {
  
  x <- fread(file)
  
  x <- x[Name %in% marker_set]
  
  # Ensure same marker order for every sample
  x <- x[match(marker_set, x$Name)]
  
  sample_id <- tools::file_path_sans_ext(basename(file))
  
  out <- x[[2]]
  names(out) <- marker_set
  
  data.table(Sample_ID = sample_id,
    t(out))
}

lrr_list <- mclapply(sig.files, read_lrr_sample, marker_set = marker_set, mc.cores = 10 )

lrr_dt <- rbindlist(lrr_list, fill = TRUE)

sample_ids <- lrr_dt$Sample_ID

lrr_mat <- as.matrix(lrr_dt[, -1])
rownames(lrr_mat) <- sample_ids

dim(lrr_mat)

# Remove bad markers and impute missing values
# Remove markers missing in >5% of samples
keep <- colMeans(is.na(lrr_mat)) < 0.05
lrr_mat <- lrr_mat[, keep]

# Median imputation per marker
for (j in seq_len(ncol(lrr_mat))) {
  lrr_mat[is.na(lrr_mat[, j]), j] <- median(lrr_mat[, j], na.rm = TRUE)
}

# run PCA
pca_lrr <- prcomp(lrr_mat,
  center = TRUE, scale. = TRUE)

pca_df <- as.data.frame(pca_lrr$x[, 1:10])
pca_df$Sample_ID <- rownames(pca_lrr$x)

# plot PCA
pca_df <- inner_join(pca_df, metadata, by="Sample_ID")

names(pca_df)

centre <- ggplot(pca_df, aes(PC1, PC2, colour = centre)) +
  geom_point(size=1, alpha = 0.8) +
  theme_bw(base_size = 12) + 
  labs(colour = "Centre") + scale_color_brewer(palette = 'Set2')

centre2 <- ggplot(pca_df, aes(PC3, PC4, colour = centre)) +
  geom_point(size=1, alpha = 0.8) +
  theme_bw(base_size = 12) +
  labs(colour = "Centre") + scale_color_brewer(palette = 'Set2')

centre3 <- ggplot(pca_df, aes(PC1, PC2, col = array)) +
  geom_point(alpha = 0.8, size=1) + theme_bw(base_size = 12) + 
  labs(colour = "Array") + scale_color_brewer(palette = 'Paired') 

centre4 <- ggplot(pca_df, aes(PC3, PC4, col = array)) +
  geom_point(size=1, alpha = 0.8) + theme_bw(base_size = 12) + 
  labs(colour = "Array") + scale_color_brewer(palette = 'Paired') 

shared_legend.cen <- get_legend(centre + theme(legend.position = "right"))
shared_legend.arr <- get_legend(centre3 + theme(legend.position = "right"))
# 3. Remove legends from all plots
centre_nolgg <- centre + theme(legend.position = "none") 
centre2_nolgg <- centre2 + theme(legend.position = "none") 
centre3_nolgg <- centre3 + theme(legend.position = "none") 
centre4_nolgg <- centre4 + theme(legend.position = "none") 

library(cowplot)
# 4. Combine plots into 2x2 grid
plot_grid_1 <- plot_grid(centre_nolgg, centre2_nolgg, shared_legend.cen, ncol = 3,
                         rel_widths = c(1, 1, 0.5), axis = 'lr', align = 'vh')

plot_grid_2 <- plot_grid(centre3_nolgg, centre4_nolgg, shared_legend.arr, ncol = 3, 
                         rel_widths = c(1, 1, 0.5), axis = 'lr', align = 'vh')

final_plot <- plot_grid(plot_grid_1, plot_grid_2, nrow = 2, labels = c('A','B'))
final_plot

ggsave('thesis_out/CLL_lrr_PCA.png', plot = final_plot, width = 9,
       height = 6, units = "in", dpi = 300, bg = "white")
