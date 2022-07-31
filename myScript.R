# ---
# name: dana hughes
# title: DGE analysis in Ewing sarcoma cells
# date: saturday july 30 2022
# ---

# libraries
library(tidyverse)
library(DESeq2)
library(EnhancedVolcano)
library(pheatmap)
library(enrichR)
library(msigdbr)
library(clusterProfiler)

# load data
rse <- readRDS("EwS.rds")
rowRanges(rse)$symbol <- sapply(rowRanges(rse)$symbol, "[", 1)
rse$condition <- c(rep("shCTR", 3), rep("shEF1", 4))

# create dds object and analyse
dds <- DESeqDataSet(rse, design = ~ condition)

# pre-filtering
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds$condition <- droplevels(dds$condition)

# run analysis
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "shEF1", "shCTR"))

# plot PCA
rld <- rlog(dds)
plotPCA(rld)

# MA plot
plotMA(res, main = "1. MA plot")
resNorm <- lfcShrink(dds = dds, res = res, type = "normal", coef = 2)
plotMA(resNorm, main = "2. MA Plot")

# write table as csv file
resdf <- as.data.frame(res)
write_csv(resdf, file = "DEGs.csv")

# volcano plot
EnhancedVolcano(resdf,
                lab = rownames(resdf),
                x = 'log2FoldChange', y = 'padj',
                title = "Volcano Plot", pCutoff = 1e-100, FCcutoff = 3,
                pointSize = 1.5)
# heat map
pheatmap(resdf)

# enrichment analysis
resdf <- resdf %>% # add GENEID to column
  rownames_to_column() %>%
  mutate(GENEID = gsub(rowname, pattern = "\\..+", replacement = "")) %>%
  dplyr::select(-rowname)

resdf %>% # over expressed genes
  dplyr::filter(padj < .01 & log2FoldChange > 2) %>%
  write_csv(file = "over_expressed_genes.csv")

resdf %>% # under expressed genes
  dplyr::filter(padj < .01 & log2FoldChange < -2) %>%
  write_csv(file = "under_expressed_genes.csv")
