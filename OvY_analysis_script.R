

# From Frank Script
knitr::opts_chunk$set(echo = TRUE, message = FALSE, warning = FALSE, cache = TRUE)

setwd("C:/Users/Potakowskyi/BioMed X GmbH/TDA - General/Iris_Teams/Paper/1010_Iris Proteomics Analysis/Data/")

## Defining a working directory
# We defined the working directory for the analysis. This directory should contain the metadata file and the tab delimited text files containing the data.

library(locfdr)
library(vsn)
library(gplots)
library(limma)
library(fdrtool)
library(biobroom)
library(tidyverse)
library(pheatmap)
library(colorspace)
library(RColorBrewer)
library(readxl)
library(writexl)

# Qualitative Hues: "#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#b15928", "#ffed6f", "#7570b3", "#e7298a"
customPlot <- list(
  theme_bw(base_size = 12), 
  scale_shape_manual(values = c(16, 17, 15, 3, 7, 8)), 
  scale_fill_manual(values = c("#d0d1e6", "#a6bddb", "#67a9cf", "#3690c0", "#02818a", "#fee090", "#fdae61", "#f46d43", "#d73027", "#a50026", "#cab2d6", "#6a3d9a", "#7570b3", "#e7298a")), 
  scale_colour_manual(values = c("#d0d1e6", "#a6bddb", "#67a9cf", "#3690c0", "#02818a", "#fee090", "#fdae61", "#f46d43", "#d73027", "#a50026", "#cab2d6", "#6a3d9a", "#7570b3", "#e7298a"))
)
script.version <- "V1"
dir_save <- paste0("data_analysis_results_", script.version)
if (!dir.exists(dir_save)){
  dir.create(dir_save)
}

# Loading data and annotating experimental conditions
# Same list as we get the excel file for small columns
# sorting out all values that are = 0 
conditions <- read.csv("metadata.csv")
files <- file.path(unique(conditions$file))
files <- files[file.exists(files)]
data <- NULL
for (i in seq_along(files))
{
  print(files[i])
  x <- read_delim(file.path(files[i]), delim = "	")
  names(x) <- make.names(names(x))
  x$Gene[is.na(x$Gene)] <- x$Protein.ID[is.na(x$Gene)]
  #keep only proteins with at least two unique peptide matches
  x <- x %>%
    dplyr::filter(Razor.Peptides >= 2)
  #remove reverse hits
  x <- x %>%
    dplyr::filter(!grepl("rev_", Protein))
  #remove known contaminants
  x <- x %>%
    dplyr::filter(!grepl("contam_", Protein))
  x$file <- files[i]
  keepnames <- as.character(subset(conditions, file == files[i])$col.name)
  x <- x[, c("Gene", "Protein.ID", "Protein.Description", "Organism", "Unique.Peptides", "Razor.Peptides", "Total.Intensity", "file", keepnames)]
  x <- x %>%
    group_by(Gene, Protein.ID, Protein.Description, Organism, Unique.Peptides, Razor.Peptides, Total.Intensity, file) %>%
    gather(key = "col.name", value = "value", keepnames)
  x$Total.Intensity <- as.numeric(as.character(x$Total.Intensity))
  data <- bind_rows(data, x)
  rm(x, keepnames)
}
rm(i, files)
data <- left_join(data, conditions)
data <- subset(data, value > 0)

# Protein identification overview
## summarizes in different ways where peptides where found/ how many unique etc.
## sub contains no raw values 
sub <- unique(data[, c("Gene", "condition", "rep")])
sub$found <- 1
sub <- unique(sub)
ggplot(data = sub, aes(condition, fill = rep)) + 
  geom_bar(position = position_dodge()) + 
  customPlot + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
with(sub, table(condition, rep))
sub_i <- data %>%
  group_by(Gene) %>%
  dplyr::select(Gene, rep, condition) %>%
  unique() %>%
  summarise(found.in.conditions = n())
ggplot(data = sub_i, aes(found.in.conditions)) + 
  geom_bar(position = position_dodge()) + 
  customPlot + 
  xlab("no of identifications per condition")
table(sub_i$found.in.conditions)
sub_i <- data %>%
  group_by(Gene, condition) %>%
  dplyr::select(Gene, rep, condition) %>%
  unique() %>%
  summarise(found.in.conditions = n())
ggplot(data = sub_i, aes(found.in.conditions)) + 
  geom_bar(position = position_dodge()) + 
  customPlot + 
  facet_wrap(~ condition) +
  xlab("no of identifications per condition")
table(sub_i$condition, sub_i$found.in.conditions)
rm(sub_i)
rm(sub)

# Transforming long data to wide data
# visualization of log intensity of each condition as a boxplot
ggplot(data = data, aes(condition, log2(value), fill = rep)) + 
  geom_boxplot() + 
  ylab("log2(Intensity)") + 
  customPlot + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
## cdata includes the raw intensity of each protein per condition, has NA values
cdata <- data %>%
  group_by(Gene) %>%
  mutate(key = paste("Intensity", condition, sep = "_")) %>%
  dplyr::select(Gene, key, value) %>%
  group_by(Gene, key) %>%
  summarise(value = sum(value, na.rm = TRUE)) %>%
  spread(key = key, value = value)


#found.in.conditions
## summarizes all conditions in which a peptide is found
fdata <- data %>%
  group_by(Gene) %>%
  dplyr::select(Gene, condition) %>%
  unique() %>%
  summarise(found.in.conditions = n())

## Yes/ No > it the peptide found in a condition or not? (= 1 and NA)
fdata_i <- data %>%
  mutate(file2 = make.names(file)) %>%
  group_by(Gene, file2) %>%
  dplyr::select(Gene, condition, file2) %>%
  unique() %>%
  summarise(found.in.condition = n()) %>%
  mutate(file2 = paste0("found.in.condition_", file2)) %>%
  spread(key = file2, value = found.in.condition)

## adds fdata_i to fdata
fdata <- left_join(fdata, fdata_i)

# repeated using Gene grouping
fdata_i <- data %>%
  group_by(Gene, condition) %>%
  dplyr::select(Gene, condition, rep) %>%
  unique() %>%
  summarise(found.in.conditions = n()) %>%
  mutate(condition = paste0("found.in.conditions_", condition)) %>%
  spread(key = condition, value = found.in.conditions)
fdata <- left_join(fdata, fdata_i)


#found.in.files
fdata_i <- data %>%
  group_by(Gene) %>%
  dplyr::select(Gene, file) %>%
  unique() %>%
  summarise(found.in.files = n())
fdata <- left_join(fdata, fdata_i)
#found.in.reps
fdata_i <- data %>%
  group_by(Gene) %>%
  dplyr::select(Gene, rep) %>%
  unique() %>%
  summarise(found.in.reps = n())
fdata <- left_join(fdata, fdata_i)
#Unique.Peptides overview
fdata_i <- data %>%
  group_by(Gene) %>%
  dplyr::select(Gene, Unique.Peptides) %>%
  unique() %>%
  summarise(max.Unique.Peptides = max(Unique.Peptides, na.rm = TRUE))
fdata <- left_join(fdata, fdata_i)
#Razor.Peptides overview
fdata_i <- data %>%
  group_by(Gene) %>%
  dplyr::select(Gene, Razor.Peptides) %>%
  unique() %>%
  summarise(max.Razor.Peptides = max(Razor.Peptides, na.rm = TRUE))
fdata <- left_join(fdata, fdata_i)
fdata_i <- data %>%
  group_by(Gene, condition) %>%
  dplyr::select(Gene, condition, Razor.Peptides) %>%
  unique() %>%
  summarise(max.Razor.Peptides = max(Razor.Peptides, na.rm = TRUE)) %>%
  mutate(condition = paste0("max.Razor.Peptides_", condition)) %>%
  spread(key = condition, value = max.Razor.Peptides)
fdata <- left_join(fdata, fdata_i)
#abundance overview
fdata_i <- data %>%
  group_by(Gene) %>%
  dplyr::select(Gene, Total.Intensity) %>%
  unique() %>%
  summarise(average.Total.Intensity = mean(Total.Intensity, na.rm = TRUE))
fdata <- left_join(fdata, fdata_i)
id_data <- data %>%
  ungroup() %>%
  dplyr::select(Gene, Protein.ID, Protein.Description, Organism) %>%
  unique()
dups <- id_data %>%
  group_by(Gene) %>%
  summarize(n = n()) %>%
  dplyr::filter(n > 1) %>%
  dplyr::select(Gene)
if (nrow(dups) > 0)
{
  dups <- dups$Gene
  dup.data <- subset(id_data, Gene %in% dups)
  id_data <- subset(id_data, !Gene %in% dups)
  combine_ids <- function(ids)
  {
    ids <- as.character(ids)
    return(paste(
      sort(unique(unlist(strsplit(ids, split = "[;]")))), collapse = ";"))
  }
  dup.data <- dup.data %>%
    group_by(Gene) %>%
    summarise(Protein.ID = combine_ids(Protein.ID), 
              Protein.Description = combine_ids(Protein.Description), 
              Organism = combine_ids(Organism))
  id_data <- bind_rows(id_data, dup.data)
  rm(dup.data, combine_ids)
}
rm(dups)
fdata <- left_join(id_data, fdata)
rm(id_data)
cdata <- left_join(fdata, cdata)
rm(fdata_i, fdata)

# Filter data
#Only proteins that were quantified with two unique peptide matches are kept for the analysis. Proteins were filtered according to these condition already when loading the data. Moreover, only proteins were kept if they were quantified in at least 2/3 of the replicates.
dim(cdata)
min.found.in.conditions <- 1
min.Razor.Peptides <- 2
for (cond in unique(conditions$condition))
{
  cname <- paste0("keep_", cond)
  cdata[, cname] <- ifelse(
    cdata[, paste0("found.in.conditions_", cond)] >= 
      min.found.in.conditions & cdata[, paste0("max.Razor.Peptides_", cond)] >= min.Razor.Peptides, TRUE, FALSE)
  cdata[is.na(cdata[, cname]), cname] <- FALSE
  rm(cname)
}
rm(cond)
cdata$keep <- FALSE
for (i in seq_along(rownames(cdata)))
{
  cdata$keep[i] <- any(as.logical(cdata[i, grep("keep_", names(cdata), value = TRUE)]))
}
rm(i)
cdata <- subset(cdata, keep)
dim(cdata)
rm(min.found.in.conditions, min.Razor.Peptides)

# Building an expression set object
# Constructing assay data
raw_data <- cdata %>%
  dplyr::select(starts_with("Intensity")) %>%
  as.data.frame()
rownames(raw_data) <- cdata$Gene
names(raw_data) <- gsub("Intensity_", "", names(raw_data))

# Constructing metadata
conditions_i <- data.frame(ID = names(raw_data))
conditions_i <- left_join(conditions_i, conditions)

# Constructing fdata
fdata <- cdata %>%
  dplyr::select(-starts_with("Intensity")) %>%
  as.data.frame()

# Defining ID columns
rownames(fdata) <- fdata$Gene
rownames(conditions_i) <- conditions_i$ID
colnames(raw_data) <- conditions_i$ID

# Log2-transformation of raw_ Intensitys
raw_data_m <- log2(as.matrix(raw_data))
raw_data_m[is.infinite((raw_data_m))] <- NA
raw_data_m[is.na((raw_data_m))] <- NA

# Creating an expression set
raw_dataE <- ExpressionSet(assayData = raw_data_m, 
                           phenoData = AnnotatedDataFrame(conditions_i), 
                           featureData = AnnotatedDataFrame(fdata))
validObject(raw_dataE)
rm(raw_data, raw_data_m, conditions_i, fdata)

# Calculation of ctrl-fold changes
fc.data <- as.data.frame(2^exprs(raw_dataE))
names.orig <- names(fc.data)
for (i in seq_along(conditions$ID))
{
  fc.data[, paste0(conditions$ID[i], ".ctrl.ratio")] <-
    fc.data[, as.character(conditions$ID[i])] /
    apply(fc.data[, names.orig], 1, sum, na.rm = TRUE)
  # apply(fc.data[, grep(paste0("^", conditions$ctrl[i]), names.orig, value = TRUE)], 1, median, na.rm = TRUE)
  # fc.data[, as.character(conditions$ctrl.ID[i])]
}
rm(i)
fc.data <- fc.data %>%
  dplyr::select(ends_with(".ctrl.ratio"))
fc.data$Gene <- rownames(fc.data)
fc.data <- fc.data %>%
  dplyr::select(ends_with(".ctrl.ratio"))
ctrl.ratio_dataE <- raw_dataE
names(fc.data) <- gsub(".ctrl.ratio", "", names(fc.data))
fc.data <- fc.data[ , colnames(ctrl.ratio_dataE)]
exprs(ctrl.ratio_dataE) <- fc.data %>%
  as.matrix() 
rm(fc.data)

# PCA - principal component analysis
sets <- list("raw data" = raw_dataE, "ctrl.ratio data" = ctrl.ratio_dataE)
PCA_data <- NULL
for (i in seq_along(sets)) {
  set <- sets[[i]]
  set.name <- names(sets)[i]
  PCA_m <- t(na.omit(exprs(set)))
  if (length(PCA_m) > 50) {
    PCA <- prcomp(PCA_m, scale = FALSE)
    perc_var <-
      round(100 * PCA$sdev ^ 2 /
              sum(PCA$sdev ^ 2), 1)
    PCA_data_i <-
      data.frame(PC1 = PCA$x[, 1],
                 PC2 = PCA$x[, 2],
                 PC3 = PCA$x[, 3],
                 PC1.var = perc_var[1],
                 PC2.var = perc_var[2],
                 condition = pData(set)$condition,
                 rep = pData(set)$rep,
                 measurement = set.name)
    PCA_data <- bind_rows(PCA_data, PCA_data_i)
    rm(PCA_data_i)
  }
}
rm(i)
PCA_data$measurement <- factor(PCA_data$measurement, ordered = TRUE, levels = c("raw data", "ctrl.ratio data"))

PCA_plot <- ggplot(data = PCA_data, aes(PC1, PC2)) +
  geom_point(aes(colour = condition, shape = rep), size = 4) +
  guides(size = "none") +
  customPlot +
  ggtitle("PCA analysis") +
  xlab("PC1") +
  ylab("PC2") +
  facet_wrap(. ~ measurement + paste("PC1:", PC1.var, "% var - PC2:", PC2.var, "% var"), scale = "free", ncol = 2)

ggsave(file.path(dir_save, paste0("PCA_analysis_", script.version, ".pdf")),
       width = 8, height = 5)
write.csv(PCA_data, file = file.path(dir_save, paste0("PCA_analysis_data_", script.version, ".csv")), row.names = FALSE)

print(PCA_plot)

# Merge modified data back into 'cdata'
cdata_i <- as.data.frame(2^exprs(ctrl.ratio_dataE))
names(cdata_i) <- paste0("ctrl.ratio_", names(cdata_i))
cdata_i$Gene <- rownames(cdata_i)
cdata <- left_join(cdata, cdata_i)
rm(cdata_i)

write.csv(cdata, file = file.path(dir_save, paste0("Full_dataset_", script.version, ".csv")), row.names = FALSE)

# Create tidy data
mdata <- NULL

mdata_i <- tidy(raw_dataE)
mdata_i <- mdata_i %>%
  mutate(value = 2 ^ value)
names(mdata_i)[1] <- "Gene"
names(mdata_i)[2] <- "ID"
mdata_i <- left_join(mdata_i, conditions)
mdata_i$measurement <- "raw_Intensity"
mdata <- bind_rows(mdata, mdata_i)
rm(mdata_i)

mdata_i <- tidy(ctrl.ratio_dataE)
mdata_i <- mdata_i %>%
  mutate(value = value)
names(mdata_i)[1] <- "Gene"
names(mdata_i)[2] <- "ID"
mdata_i <- left_join(mdata_i, conditions)
mdata_i$measurement <- "ctrl.ratio"
mdata <- bind_rows(mdata, mdata_i)
rm(mdata_i)

mdata$condition <- factor(mdata$condition, ordered = TRUE, levels = c("SEC_35_fr1", "SEC_35_fr2", "SEC_35_fr3", "SEC_35_fr4", "SEC_35_fr5", "SEC_70_fr1", "SEC_70_fr2", "SEC_70_fr3", "SEC_70_fr4", "SEC_70_fr5"))
mdata$measurement <- factor(mdata$measurement, ordered = TRUE, levels = c("raw_Intensity", "ctrl.ratio"))

# Data transformation overview
ggplot(data = subset(mdata, !grepl("ratio", measurement)), aes(condition, log2(value))) +
  geom_boxplot(aes(fill = rep)) +
  customPlot +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  facet_grid(. ~ measurement)
ggsave(file.path(dir_save, paste0("Normalization_overview_", script.version, ".pdf")), width = 8.4, height = 5)

ggplot(data = subset(mdata, grepl("ratio", measurement)), aes(condition, (value))) +
  geom_boxplot(aes(fill = rep)) +
  customPlot +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  facet_grid(. ~ measurement, scale = "free")
ggsave(file.path(dir_save, paste0("Normalization_overview_ratios_", script.version, ".pdf")), width = 8.4, height = 5)

# Ratio analysis
#A protein is considered a 'hit', if the false discovery rate is smaller 0.05 and 
# a fold change of at least 2-fold is observed. 
#A protein is considered a 'candidate', if the false discovery rate is smaller 0.2 and 
# a fold change of at least 1.5-fold is observed.
fc.data <- as.data.frame(2^exprs(raw_dataE))
AveExpr.data <- NULL
log2mean <- function(x)
{
  return(mean(log2(x), na.rm = TRUE))
}
orig.ncol <- ncol(fc.data)
comparison.list <- c("SEC_70_fr1_vs_SEC_35_fr1" = "SEC_70_fr1 - SEC_35_fr1",
                     "SEC_70_fr2_vs_SEC_35_fr2" = "SEC_70_fr2 - SEC_35_fr2",
                     "SEC_70_fr3_vs_SEC_35_fr3" = "SEC_70_fr3 - SEC_35_fr3",
                     "SEC_70_fr4_vs_SEC_35_fr4" = "SEC_70_fr4 - SEC_35_fr4",
                     "SEC_70_fr5_vs_SEC_35_fr5" = "SEC_70_fr5 - SEC_35_fr5")
for (i in seq_along(comparison.list))
{
  comparison <- names(comparison.list)[i]
  comparison.cols <- unlist(strsplit(comparison.list[i], " - "))
  nominator <- comparison.cols[1]
  denominator <- comparison.cols[2]
  fc.data[, comparison] <-
    fc.data[, nominator] / fc.data[, denominator]
  AveExpr <- 
    base::apply(fc.data[, comparison.cols],
                MARGIN = 1, log2mean)
  AveExpr.data_i <- data.frame(Gene = rownames(fc.data),
                               comparison = comparison,
                               AveExpr = AveExpr)
  AveExpr.data <- bind_rows(AveExpr.data, AveExpr.data_i)
}
rm(AveExpr, AveExpr.data_i, log2mean)
fc.data <- fc.data[-c(1:orig.ncol)]

## Calculating p-values based on fold change distribution
ratio.treshold <- 2
ratio_results <- NULL
for (cn in names(fc.data))
{
  ratios <- log2(fc.data[, cn])
  ratio_results_i <- data.frame(Gene = rownames(fc.data),
                                comparison = cn,
                                logFC = ratios,
                                z.score = NA,
                                pvalue = NA,
                                fdr = NA)
  nas <- is.na(ratios)
  cleaned.ratio <- ratios[!nas]
  pratios <- cleaned.ratio >= 0
  nratios <- !pratios
  if (abs(median(cleaned.ratio)) > log2(ratio.treshold))
  {
    fdr.estimation.method <- "ratios -> student's t-distribution -> p.adjust -> fdr"
    hist(cleaned.ratio, main = cn, breaks = 50)
    pvalues <- NA
    pvalues[nratios] <- pt(cleaned.ratio[nratios], df = length(cleaned.ratio[nratios])) * 2
    pvalues[pratios] <- pt(cleaned.ratio[pratios], df = length(cleaned.ratio[pratios]), lower.tail = FALSE) * 2
    ratio_results_i$z.score <- ratios
    ratio_results_i$pvalue[!nas] <- pvalues
    ratio_results_i$fdr[!nas] <- p.adjust(pvalues, method = "fdr")
    rm(pvalues)
  }else
  {
    lfdr <- NULL
    try(lfdr <- locfdr(cleaned.ratio, main = cn))
    if (!is.null(lfdr))
    {
      est_mean <- lfdr$fp0[, "delta"]["mlest"]
      est_sd <- lfdr$fp0[, "sigma"]["mlest"]
      if (is.na(est_mean))
      {
        est_mean <- lfdr$fp0[, "delta"]["cmest"]
        est_sd <- lfdr$fp0[, "sigma"]["cmest"]
        fdr.estimation.method <- "locfdr(CME) -> z"
      }else
      {
        fdr.estimation.method <- "locfdr(MLE) -> z"
      }
      z <- (cleaned.ratio - est_mean) / est_sd
      rm(est_mean, est_sd)
    }else
    {
      fdr.estimation.method <- "quantile normalization -> z"
      hist(cleaned.ratio, main = cn, breaks = 50)
      quants <- quantile(cleaned.ratio, probs = c(0.1587, 0.5, 0.8413), na.rm = TRUE)
      z <- cleaned.ratio - quants[2]
      pratios <- z >= 0
      nratios <- !pratios
      z[nratios] <- z[nratios] / as.numeric(diff(quants)[1])
      z[pratios] <- z[pratios] / as.numeric(diff(quants)[2])
      rm(quants)
    }
    ratio_results_i$z.score[!nas] <- z
    rm(lfdr)
    fdrtool.results <- NULL
    fdrtool.results <- try(fdrtool(z))
    if (!is.null(fdrtool.results))
    {
      fdr.estimation.method <- paste(fdr.estimation.method, "-> fdrtool(qval) -> fdr")
      ratio_results_i$pvalue[!nas] <- fdrtool.results$pval
      ratio_results_i$fdr[!nas] <- fdrtool.results$qval
    }else
    {
      fdr.estimation.method <- paste(fdr.estimation.method, "-> student's t-distribution -> p.adjust -> fdr")
      pz <- z >= 0
      nz <- !pratios
      pvalues <- NA
      pvalues[nz] <- pt(z[nz], df = length(z[nz])) * 2
      pvalues[pz] <- pt(z[pz], df = length(z[pz]), lower.tail = FALSE) * 2
      ratio_results_i$pvalue[!nas] <- pvalues
      ratio_results_i$fdr[!nas] <- p.adjust(pvalues, method = "fdr")
      rm(pz, nz, pvalues)
    }
    rm(fdrtool.results, z)
  }
  ratio_results_i$fdr.estimation.method <- fdr.estimation.method
  ratio_results <- bind_rows(ratio_results, ratio_results_i)
  rm(ratio_results_i, fdr.estimation.method, pratios, nratios, nas)
}
rm(fc.data, cn, ratio.treshold)

## Merge ratio results with feature Data
ratio_results <- left_join(AveExpr.data, ratio_results)
ratio_results <- left_join(fData(raw_dataE), ratio_results)

## Hit annotation
fdr_hit_threshold <- 0.05
fdr_candidate_threshold = 0.2
fc_hit_threshold <- 2
fc_candidate_threshold <- 1.5
ratio_results$hit <- 
  with(ratio_results, ifelse(fdr <= fdr_hit_threshold & abs(logFC) >= 
                               log2(fc_hit_threshold), TRUE, FALSE))
ratio_results$hit_annotation <- 
  with(ratio_results, 
       ifelse(fdr <= fdr_hit_threshold & abs(logFC) >= 
                log2(fc_hit_threshold), "hit",
              ifelse(fdr <= fdr_candidate_threshold & abs(logFC) >= 
                       log2(fc_candidate_threshold), "candidate", "no hit")))
ratio_results$hit_annotation <- 
  factor(ratio_results$hit_annotation, 
         ordered = TRUE, 
         levels = c("hit", "candidate", "no hit"))
with(ratio_results, table(comparison, hit_annotation))

## Volcano-plot
ggplot(data = ratio_results, aes(logFC, -log10(pvalue), colour = hit_annotation)) +
  geom_vline(aes(xintercept = 0)) + 
  geom_point() +
  geom_text(aes(label = gsub("[|].+", "", Gene)), 
            data = subset(ratio_results, hit_annotation != "no hit"), 
            vjust = 0, nudge_y = 0.1, size = 2, check_overlap = TRUE) + 
  facet_wrap( ~ comparison) +
  xlab("log2(fold change)") +
  customPlot
ggsave(file.path(dir_save, paste0("Volcano_plot_", script.version, ".pdf")), width = 16, height = 9)

# Clustering of data
## Preparation of cluster data
# hits <- unique(subset(ratio_results, hit_annotation %in% c("hit", "candidate"))$Gene)
hits <- cdata$Gene

# The median is calculated if there would be Gene-Condition duplicates?
m_clust.data <- subset(mdata, measurement == "ctrl.ratio" & Gene %in% hits) %>%
  group_by(Gene, condition) %>%
  summarise(median.value = median(value, na.rm = TRUE))

clust.data_m <- m_clust.data %>%
  mutate(key = paste(condition, sep = "_")) %>%
  ungroup() %>%
  dplyr::select(Gene, key, median.value) %>%
  group_by(Gene) %>%
  spread(key = key, value = median.value) %>%
  as.data.frame()

write.csv(clust.data_m, file.path(dir_save, paste0("Clustering_data_", nrow(clust.data_m), "_proteins_", script.version, ".csv")), row.names = FALSE)
rownames(clust.data_m) <- clust.data_m$Gene

m_clust.data <- m_clust.data %>%
  dplyr::filter(Gene %in% clust.data_m$Gene)
clust.data_m$Gene <- NULL
clust.data_m <- clust.data_m %>%
  as.matrix() 

hca_1 <- clust.data_m
hca_1[is.na(hca_1)] <- 0
clust.data_m[is.na(clust.data_m)] <- 0

## Set number of clusters
cluster.number <- ceiling(sqrt(sqrt(nrow(clust.data_m)) * ncol(clust.data_m)))

## K-means clustering
no <- ceiling(nrow(clust.data_m) / 4)
wss <- (no - 1) * sum(apply(clust.data_m, 2, var))
for (i in 2:(no - 1)) wss[i] <- sum(kmeans(clust.data_m, 
                                           centers = i)$withinss)
plot(1:(no - 1), wss, type = "b", xlab = "number of Clusters",
     ylab = "within groups sum of squares")
abline(v = cluster.number, col = "red")
kmeans.fit <- kmeans(clust.data_m, cluster.number)
cluster.groups <- data.frame(Gene = names(kmeans.fit$cluster),
                             kmeans.cluster.group = kmeans.fit$cluster)
rm(wss, kmeans.fit, no)

## Hierarchical clustering
d <- dist(clust.data_m, method = "euclidean")
hclust.fit <- hclust(d, method = "ward.D2")
cluster.groups_hclust <- 
  data.frame(Gene = names(cutree(hclust.fit, k = cluster.number)),
             hclust.cluster.group = cutree(hclust.fit, k = cluster.number))
cluster.groups <- left_join(cluster.groups, cluster.groups_hclust)
rm(cluster.groups_hclust)
m_clust.data <- left_join(m_clust.data, cluster.groups)
m_clust.data$Gene <- factor(m_clust.data$Gene, ordered = TRUE, 
                            levels = hclust.fit$labels[hclust.fit$order])
rm(d, hclust.fit)

# Plotting all genes on a heatmap, without clustering
gr.width <- 2.5 + ncol(clust.data_m) * 0.1 + max(nchar(rownames(clust.data_m))) * 0.008
gr.height <- 2.5 + nrow(clust.data_m) * 0.01
ggplot(data = m_clust.data, aes(condition, Gene)) +
  geom_tile(aes(fill = median.value)) +
  scale_fill_gradientn(colours = c("#377eb8", "#984ea3", "#e41a1c", "#ff7f00", "#ffff33"), name = "proportion") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
        axis.text.y = element_text(size = 1)) +
  ylab("gene")

#Plotting hierarchical clustering (EMBL)
cl_graph_2 <- ggplot(data = m_clust.data, aes(condition, Gene)) +
  geom_tile(aes(fill = median.value)) +
  scale_fill_gradientn(colours = c("#377eb8", "#984ea3", "#e41a1c", "#ff7f00", "#ffff33"), name = "proportion") +
  facet_grid(hclust.cluster.group ~ ., scale = "free", space = "free") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
        axis.text.y = element_text(size = 1)) +
  ylab("gene") +
  ggtitle("hierarchical clustering")


# Plotting k-means clustering (EMBL)
ggplot(data = m_clust.data, aes(condition, Gene)) +
  geom_tile(aes(fill = median.value)) +
  scale_fill_gradientn(colours = c("#377eb8", "#984ea3", "#e41a1c", "#ff7f00", "#ffff33"), name = "proportion") +
  facet_grid(kmeans.cluster.group ~ ., scale = "free", space = "free") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
        axis.text.y = element_text(size = 1)) +
  ylab("gene") +
  ggtitle("kmeans clustering")

#### Heatmap from Veerman et al.
colourscheme <- (c("skyblue4", "skyblue3", "orangered2", "lightgoldenrod2", "darkseagreen4", "darkseagreen3", "sienna1", "plum3", "plum1"))
pal1 <- colourscheme
pal2 <- sequential_hcl(10, palette = "Mint")
pal3 <- (c("aquamarine3",  "tomato3", "cadetblue", "azure3", "darkolivegreen3", "lightpink2", "navajowhite2", "steelblue4"))
pal4 <- c("#feb24c", "#fd8d3c", "#fc4e2a", "#e31a1c", "#b10026")
#pal5 <- c("#d0d1e6", "#a6bddb", "#67a9cf", "#3690c0", "#02818a", "#fee090", "#fdae61", "#f46d43", "#d73027", "#a50026", "#cab2d6", "#6a3d9a", "#7570b3", "#e7298a")
pal5 <- colorRampPalette(brewer.pal(5, "Set3"))(15)

# Generate heatmap with cluster numbers in row labels (veerman) - use Frank input
cl_graph_1 <- pheatmap(
  clust.data_m,
  cellheight = 0.6,
  #scale = "row",
  color = c("#fff7ec", colorRampPalette(pal4)(99)),
  border_color = NA,
  cluster_cols = F,
  cutree_rows = 15,
  cex = 1,
  fontsize_col = 10,
  show_rownames = F,
  treeheight_col = 30,
)


# Generate heatmap with cluster numbers in row labels (veerman)
cl_graph_1 <- pheatmap(
  hca_1,
  clustering_method = "ward.D2",
  cellheight = 0.6,
  #scale = "row",
  color = c("#fff7ec", colorRampPalette(pal4)(99)),
  border_color = NA,
  cluster_cols = F,
  cutree_rows = 15,
  cex = 1,
  fontsize_col = 10,
  show_rownames = F,
  treeheight_col = 30,
)

clusters <- cutree(cl_graph_1$tree_row, k = 15)[cl_graph_1$tree_row[["order"]]]
annot_row <- data.frame(row.names = names(clusters), cluster = as.factor(clusters))

# Assign the palette to heatmap_ann
heatmap_ann <- pal5


# Generate heatmap with cluster annotations
cl_graph_1 <- pheatmap(
  hca_1,
  clustering_method = "ward.D2",
  cellheight = 0.6,
  #scale = "row",
  color = c("#fff7ec", colorRampPalette(pal4)(99)),
  border_color = NA,
  cluster_cols = F,
  cutree_rows = 15,
  cex = 1,
  fontsize_col = 10,
  angle_col = 45,
  show_rownames = FALSE,  # Hide row names
  annotation_row = annot_row,
  treeheight_col = 30
)

pdf(file.path(dir_save, paste0("Heatmap_", script.version, ".pdf")), width = 9, height = 9)
print(cl_graph_1)
dev.off()

#make annotation for the rows
#cluster_n <- list()
#cluster_n[["clust.data_m"]]<- data.frame(cluster_plasma=factor(gene_clusters_plasma) , row.names = names(gene_clusters_plasma) )

####################################################### Veerman et al. GO analysis
#Get the proteins in the clusters
gene_clusters_plasma <- cutree(cl_graph_1$tree_row, k = 15)

gene_clusters_plasma_list <- as.data.frame(gene_clusters_plasma)
gene_clusters_plasma_list <- rownames_to_column(gene_clusters_plasma_list, var = "Protein")
           
cl_plasma_sep <- list()

for (k in levels(factor(gene_clusters_plasma))) {
  cl_plasma_sep[[paste("cluster", k, sep = "_")]]<- names(gene_clusters_plasma[gene_clusters_plasma==k])
}

clusters_plasma_proteins <- do.call("rbind", lapply(cl_plasma_sep, FUN = as.data.frame))
write.csv(clusters_plasma_proteins, file.path(dir_save, paste0("Proteins_cluster", script.version, ".csv")), row.names = FALSE)

# GO analysis
library(devtools)
library("enrichR")
library(rafalib)


listEnrichrDbs()
go_cluster_plasma <- sapply(cl_plasma_sep, function(x) { enrichr(genes = x, databases = "GO_Cellular_Component_2018")})

# Set the width and height of the plotting device
pdf(file.path(dir_save, paste0("GO_Veerman", script.version, ".pdf")), width = 10, height = 6)

for (q in 1:length(go_cluster_plasma)) {
  # Check if the cluster is 12 or 14, if so, skip to the next iteration
  if (q %in% c(12, 14)) {
    next
  }
  
  # Check if there are any rows in the cluster
  if (sum(!is.na(hca_1[clusters == q])) > 0) {
    barplot(-log10(go_cluster_plasma[[q]][["Adjusted.P.value"]][10:1]), 
            horiz = TRUE, 
            names.arg = sub("\\(GO.*", "", go_cluster_plasma[[q]][["Term"]])[10:1], 
            las = 1, 
            cex.names = 0.9, 
            col = rev(pal2), 
            main = sub(".G.*", "", names(go_cluster_plasma[q])), 
            border = NA, 
            xlab = "-log10 adjusted P value", 
            space = 0.4)
  }
}

# Close the pdf device
dev.off()


#### Get cluster numbers for EV associated proteins

database <- read_excel("C:/Users/Potakowskyi/BioMed X GmbH/TDA - General/Iris_Teams/Paper/1010_Iris Proteomics Analysis/Data/template anaylsis.xlsx")
database <- as.data.frame(database)

#### Get cluster numbers for EV associated proteins
gene_clusters_plasma_list <- read_excel("Protein_clusters_VeermanV1.xlsx")
database <- read_excel("C:/Users/Potakowskyi/BioMed X GmbH/TDA - General/Iris_Teams/Paper/1010_Iris Proteomics Analysis/Data/template anaylsis.xlsx")
database <- as.data.frame(database)

# Initialize an empty vector to store the column names
matching_columns <- character(nrow(gene_clusters_plasma_list))

# Loop over each row in gene_clusters_plasma_list
for (i in 1:nrow(gene_clusters_plasma_list)) {
  # Get the protein name from gene_clusters_plasma_list
  protein <- gene_clusters_plasma_list$Protein[i]
  
  # Find the matching rows in the database
  matching_rows <- which(database == protein, arr.ind = TRUE)
  
  # If there are matches, record the column names
  if (nrow(matching_rows) > 0) {
    column_names <- colnames(database)[matching_rows[, "col"]]
    matching_columns[i] <- paste(column_names, collapse = ", ")
  } else {
    matching_columns[i] <- NA
  }
}

# Add the matching column names to gene_clusters_plasma_list
gene_clusters_plasma_list$matching_columns <- matching_columns

write_xlsx(gene_clusters_plasma_list, file.path(dir_save, paste0("Protein_clusters_Veerman", script.version, ".xlsx")))

setwd("C:/Users/Potakowskyi/BioMed X GmbH/TDA - General/Iris_Teams/Paper/1010_Iris Proteomics Analysis/Small Columns/data_analysis_results_V1/")

library(readxl)
library(tidyverse)
library(writexl)

data <- read_excel("Protein_clusters_VeermanV1.xlsx")
str(data)

# Define the categories
categories <- c("Plasma proteins", "GO-EV proteins", "Top100 EV proteins")

# Create a function to filter the data based on categories
filter_category <- function(category) {
  data %>%
    filter(str_detect(matching_columns, category))
}

# Create separate data frames for each category
plasma_proteins_df <- filter_category("Plasma proteins")
go_ev_proteins_df <- filter_category("GO-EV proteins")
top100_ev_proteins_df <- filter_category("Top100 EV proteins")

# Create a list of data frames
data_frames <- list(
  "Plasma Proteins" = plasma_proteins_df,
  "GO-EV Proteins" = go_ev_proteins_df,
  "Top100 EV Proteins" = top100_ev_proteins_df
)

write_xlsx(data_frames, "protein_categories.xlsx")

# Function to calculate frequency for a given data frame
calculate_frequency <- function(df, category_name) {
  df %>%
    count(gene_clusters_plasma, name = paste0("frequency_", category_name)) %>%
    rename(cluster_number = gene_clusters_plasma)
}

# Apply the function to each data frame with appropriate category names
plasma_proteins_freq <- calculate_frequency(plasma_proteins_df, "plasma")
go_ev_proteins_freq <- calculate_frequency(go_ev_proteins_df, "go_ev")
top100_ev_proteins_freq <- calculate_frequency(top100_ev_proteins_df, "top100_ev")

# Create a list of data frames
frequency_dfs <- list(
  "Plasma Proteins Frequencies" = plasma_proteins_freq,
  "GO-EV Proteins Frequencies" = go_ev_proteins_freq,
  "Top100 EV Proteins Frequencies" = top100_ev_proteins_freq
)

# Save the list of data frames to an Excel file
write_xlsx(frequency_dfs, "Number_of_proteins_per_cluster.xlsx")
