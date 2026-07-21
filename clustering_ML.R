library(vegan)

#test with surv data and typed meatadata
cna.mat <- read.csv('All CLL outcome_spss.csv')

# use unsupervised clusrting using the minimum CNAs
# plot Survival
names(cna.mat)
cna_features <- cna.mat[c(22:ncol(cna.mat))]
rownames(cna_features) <- cna.mat$Sample
prevalence <- colMeans(cna_features != 0, na.rm = TRUE)
plot(prevalence)

cna_filtered <- cna_features#[, prevalence >= 0.03 & prevalence <= 0.97]
cna_filtered2 <- cna_filtered[rowSums(cna_filtered != 0) > 0, ]

dim(cna_filtered)

# create distance mat
jaccard_dist <- vegdist(cna_filtered2,
  method = "jaccard", binary = TRUE)

# hierarchical clustering
hc <- hclust(jaccard_dist, method = "average")

plot(hc, labels = FALSE,
  hang = -1,
  main = "Hierarchical clustering of recurrent CNAs",
  xlab = "Patients", ylab = "Jaccard distance")

clusters_k2 <- cutree(hc, k = 2)
clusters_k3 <- cutree(hc, k = 3)
clusters_k5 <- cutree(hc, k = 5)

table(clusters_k3)
table(clusters_k6)

# check clusters are stable
library(cluster)

sil_results <- lapply(2:10, function(k) {
  cl <- cutree(hc, k = k)
  sil <- silhouette(cl, jaccard_dist)
  
  data.frame(k = k,
    mean_silhouette = mean(sil[, "sil_width"]))
})

sil_results <- do.call(rbind, sil_results)

plot(sil_results$k,
  sil_results$mean_silhouette,
  type = "b",
  xlab = "Number of clusters",
  ylab = "Mean silhouette width")

# heatmap
library(ComplexHeatmap)
library(circlize)

cluster_factor <- factor(clusters_k6)

Heatmap(as.matrix(cna_filtered2),
  name = "CNA",
  show_row_names = FALSE,
  show_column_names = FALSE,
  cluster_rows = hc,
  cluster_columns = FALSE,
  row_split = 6,
  column_title = "Recurrent CNA regions",
  row_title = "Patients")

## test clusters in clinical measures

cna.mat.filtered <- cna.mat[match(rownames(cna_filtered2),
                                  cna.mat$Sample), ]

stopifnot(identical(cna.mat.filtered$Sample, rownames(cna_filtered2)))

cna.mat.filtered$cluster <- factor(clusters_k2)

# plot KM
km_fit <- survfit(Surv(OSdays, OSstatus) ~ cluster,
  data = cna.mat.filtered)

ggsurvplot(km_fit,
  data = cna.mat.filtered,
  pval = TRUE,
  risk.table = TRUE,
  conf.int = F,
  xlab = "Time",
  legend.title = "CNA cluster")


##################################
### survival forest
############################
library(randomForestSRC)
library(survival)
library(survminer)
library(timeROC)
library(ggplot2)
library(dplyr)

cna.mat <- read.csv('All CLL outcome_spss.csv')
# Select the 27 CNA variables
cna_cols <- names(cna.mat[c(22:ncol(cna.mat))])

length(cna_cols)
# Should return 27

rsf_data <- cna.mat %>%
  dplyr::select(OSdays, OSstatus, all_of(c(22:ncol(cna.mat)))) %>%
  filter(complete.cases(.))

# Ensure correct data types
rsf_data$OSdays  <- as.numeric(rsf_data$OSdays)
rsf_data$OSstatus <- as.numeric(rsf_data$OSstatus)

# CNA variables should usually be numeric/binary
rsf_data[cna_cols] <- lapply(rsf_data[cna_cols],
  function(x) as.numeric(as.character(x)))

table(rsf_data$OSstatus)
dim(rsf_data)

# fit random surv forest
set.seed(123)

rsf_model <- rfsrc(Surv(OSdays, OSstatus) ~ .,
  data       = rsf_data,
  ntree      = 1000,
  nodesize   = 15,
  mtry       = 5,
  importance = "permute",
  block.size = 1,
  na.action  = "na.impute")

rsf_model

# obtain c index
oob_error <- tail(rsf_model$err.rate, 1)
oob_cindex <- 1 - oob_error

oob_cindex

performance <- data.frame(Metric = "Out-of-bag C-index",
  Value = round(oob_cindex, 3))

performance

# model convergence
error_df <- data.frame(
  Trees = seq_along(rsf_model$err.rate),
  OOB_Cindex = 1 - rsf_model$err.rate)

p_error <- ggplot(error_df, aes(x = Trees, y = OOB_Cindex)) +
  geom_line(linewidth = 0.8) +
  labs(
    title = "Random Survival Forest performance",
    subtitle = paste0(
      "Final OOB C-index = ",
      round(oob_cindex, 3)
    ),
    x = "Number of trees",
    y = "Out-of-bag C-index"
  ) +
  theme_classic(base_size = 13)

p_error

#5. Plot CNA variable importance

vimp_df <- data.frame(CNA = names(rsf_model$importance),
  Importance = as.numeric(rsf_model$importance)) %>%
  arrange(desc(Importance))

vimp_df

top_vimp <- vimp_df %>%  slice_head(n = 10) %>%
  mutate(CNA = reorder(CNA, Importance))

p_vimp <- ggplot(top_vimp, aes(x = CNA, y = Importance)) +
  geom_col() +
  coord_flip() +
  labs(
    title = "Most prognostic recurrent CNAs",
    subtitle = "Random Survival Forest permutation importance",
    x = NULL,
    y = "Permutation importance") +
  theme_classic(base_size = 13)

p_vimp

rsf_risk <- rsf_model$predicted.oob

summary(rsf_risk)
