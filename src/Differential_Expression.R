library(gplots)
library(dplyr)
library("biomaRt")

ESET_FPKM = readRDS("./data_proc/ESET_FPKM.rds")

ESET_FPKM_mat = exprs(ESET_FPKM)
dim(ESET_FPKM_mat)
pd = pData(ESET_FPKM)

# ttest = apply(ESET_FPKM_mat, 1, function(x) {
#     tt = t.test(x ~ pd$V8)
#     return (c (tt$p.value, tt$statistic))
# })

# saveRDS(ttest, "./data_proc/ttest.rds")
ttest = readRDS("./data_proc/ttest.rds")

ttest_df = as.data.frame(t(ttest))
colnames(ttest_df) = c("pval", "tval")
head(ttest_df)

hist(ttest_df$pval, breaks = 1000)
ttest_df$pval_adj = p.adjust(ttest_df$pval, method = "bonferroni")
hist(ttest_df$pval_adj, breaks = 1000)

padj_sort = ttest_df[order(ttest_df$pval_adj)[1:200], ]
head(padj_sort)

ESET_FPKM_mat_sel = ESET_FPKM_mat[rownames(padj_sort), ]
ROWV = hclust(dist(ESET_FPKM_mat_sel), method = "ward.D")
COLV = hclust(dist(t(ESET_FPKM_mat_sel)), method = "ward.D")
heatmap.2(ESET_FPKM_mat_sel, trace = "none", Colv = as.dendrogram(COLV), Rowv = as.dendrogram(ROWV), scale = "row", col = bluered(24), ColSideColors = c("white", "red", "blue")[pd$V8])

# Get gene annotation -----------------------------------------------------

ensembl = useMart("ensembl")
ld = listDatasets(ensembl)
hg38 = useDataset("hsapiens_gene_ensembl", mart = ensembl)
hg38 = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
attr = listAttributes(hg38)

## up-regulated genes in cancer
ttest_df_pos = ttest_df[order(ttest_df$tval, decreasing = T), ][1:200, ]

ensemblID_pos = sapply(rownames(ttest_df_pos), function(x) strsplit(x, "\\.")[[1]][1])
anno_df_NM_pos = getBM(values = ensemblID_pos, 
                       attributes = c("hgnc_symbol",  "entrezgene", "description", "ensembl_gene_id"),  
                       filters = "ensembl_gene_id", 
                       mart = hg38)
head (anno_df_NM_pos)
write.csv(anno_df_NM_pos, "./output/P1RCC_up_genes.csv")

## down-regulated genes in cancer
ttest_df_neg = ttest_df[order(ttest_df$tval), ][1:200, ]

ensemblID_neg = sapply(rownames(ttest_df_neg), function(x) strsplit(x, "\\.")[[1]][1])
anno_df_NM_neg = getBM(values = ensemblID_neg, 
                       attributes = c("hgnc_symbol",  "entrezgene", "description", "ensembl_gene_id"), 
                       filters = "ensembl_gene_id", 
                       mart = hg38)
head (anno_df_NM_neg)
write.csv(anno_df_NM_neg, "./output/P1RCC_down_genes.csv")

ESET_FPKM_mat_sel = ESET_FPKM_mat[rownames(ttest_df_neg), ]
ROWV = hclust(dist(ESET_FPKM_mat_sel), method = "ward.D")
COLV = hclust(dist(t(ESET_FPKM_mat_sel)), method = "ward.D")
heatmap.2(ESET_FPKM_mat_sel, trace = "none", Colv = as.dendrogram(COLV), Rowv = as.dendrogram(ROWV), scale = "row", col = bluered(24), ColSideColors = c("white", "red", "blue")[pd$V8])
