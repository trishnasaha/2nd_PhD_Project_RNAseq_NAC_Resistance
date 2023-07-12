# Reference
# http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

# Installation
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()

# sudo apt-get install libxml2-dev
# sudo apt-get install r-cran-xml
# sudo apt-get install libcurl4-openssl-dev
# install.packages("https://cran.r-project.org/src/contrib/Archive/locfit/locfit_1.5-9.4.tar.gz", repos = NULL, type="source")

BiocManager::install("DESeq2")

# Load libraries
library(DESeq2)
library(EnhancedVolcano)

# Read a merged_counts file
cts = read.table("cts.tsv", sep = '\t')

# Read phenotype sample data
coldata = read.csv("coldata.tsv", sep = '\t')

# Validate that samples in `cts` and `coldata` are the same
all(colnames(cts) %in% rownames(coldata))

# Validate that sample order in `cts` and `coldata` are the same
all(colnames(cts) == rownames(coldata))

# DESeq2 Differential Expression Analysis
# Create DESeq2 object 
dds = DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~ condition)

# Pre-filter
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Setting reference level
dds$condition <- relevel(dds$condition, ref = "Control")

# Run Differential expression analysis
dds <- DESeq(dds)
res <- results(dds)
res

summary(res)

# identify genes with FDR < 0.05
res_sig = subset(res, res$padj < 0.05)
summary(res_sig)

# up-regulated genes
up <- subset(res_sig, log2FoldChange >= 2)
summary(up)
write.table(up, file='up_regulated.tsv', quote=FALSE, sep='\t')

# down-regulated genes
down <- subset(res_sig, log2FoldChange <= -2)
summary(down)
write.table(down, file='down_regulated.tsv', quote=FALSE, sep='\t')

# Volcano Plot
png("volcano-plot.png", units="px", width=4500, height=4000, res=600)
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoff = 1e-05,
                FCcutoff = 2,
                selectLab = c("NA"))
dev.off()