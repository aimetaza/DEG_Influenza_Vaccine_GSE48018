if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c("GEOquery", "limma"), ask = FALSE, update = FALSE)
BiocManager::install("illuminaHumanv4.db", ask = FALSE, update = FALSE)

install.packages(c("pheatmap", "ggplot2", "dplyr"))

if (!requireNamespace("umap", quietly = TRUE)) {
  install.packages("umap")
}
install.packages("umap")
library(GEOquery)
library(limma)
library(pheatmap)
library(ggplot2)
library(dplyr)
library(illuminaHumanv4.db)
library(AnnotationDbi)
library(umap)
gset <- getGEO(filename = file.choose())
class(gset)
dim(exprs(gset))
ex <- exprs(gset)
dim(ex)
qx <- as.numeric(quantile(ex, c(0, 0.25, 0.5, 0.75, 0.99, 1), na.rm = TRUE))
qx
pheno <- pData(gset)
head(pheno)
titles <- pData(gset)$title
head(titles)

# Ambil bagian setelah underscore terakhir
day_info <- sub(".*_", "", titles)

# Lihat hasilnya
head(day_info)
gset$Day <- factor(day_info)
levels(gset$Day)

design <- model.matrix(~0 + gset$Day)
colnames(design) <- levels(gset$Day)
design[1:5, ]

contrast_matrix <- makeContrasts(
  Day1 - Day0,
  levels = design
)

contrast_matrix

fit <- lmFit(ex, design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

results <- topTable(
  fit2,
  adjust = "fdr",
  number = Inf
)

head(results)

results$status <- "NO"

results$status[results$logFC > 0 & results$adj.P.Val < 0.05] <- "UP"
results$status[results$logFC < 0 & results$adj.P.Val < 0.05] <- "DOWN"

table(results$status)

volcano_data <- data.frame(
  logFC = results$logFC,
  adj.P.Val = results$adj.P.Val,
  status = results$status
)

library(ggplot2)

ggplot(volcano_data, aes(x = logFC, y = -log10(adj.P.Val), color = status)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("DOWN" = "blue", "NO" = "grey", "UP" = "red")) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal() +
  ggtitle("Volcano Plot: Day1 vs Day0 (Influenza Vaccine)")

results_sorted <- results[order(results$adj.P.Val), ]
top50 <- head(results_sorted, 50)

head(top50)

mat_heatmap <- ex[rownames(top50), ]
dim(mat_heatmap)

mat_heatmap <- mat_heatmap[rowSums(is.na(mat_heatmap)) == 0, ]

gene_variance <- apply(mat_heatmap, 1, var)
mat_heatmap <- mat_heatmap[gene_variance > 0, ]

annotation_col <- data.frame(
  Day = gset$Day
)

rownames(annotation_col) <- colnames(mat_heatmap)

library(pheatmap)

pheatmap(
  mat_heatmap,
  scale = "row",
  annotation_col = annotation_col,
  show_colnames = FALSE,
  fontsize_row = 6,
  main = "Top 50 DEGs - Day1 vs Day0"
)

if (!require("clusterProfiler")) {
  BiocManager::install("clusterProfiler")
}

if (!require("org.Hs.eg.db")) {
  BiocManager::install("org.Hs.eg.db")
}

library(clusterProfiler)
library(org.Hs.eg.db)

up_genes <- rownames(results[results$status == "UP", ])
length(up_genes)

library(AnnotationDbi)

gene_symbols <- AnnotationDbi::select(
  illuminaHumanv4.db,
  keys = up_genes,
  columns = "SYMBOL",
  keytype = "PROBEID"
)

gene_symbols <- na.omit(gene_symbols$SYMBOL)
head(gene_symbols)

entrez_ids <- bitr(
  gene_symbols,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

head(entrez_ids)

ego <- enrichGO(
  gene          = entrez_ids$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  ont           = "BP",      # Biological Process
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

head(ego)

dotplot(ego, showCategory = 15) +
  ggtitle("GO Enrichment - Day1 vs Day0 (UP Genes)")


ekegg <- enrichKEGG(
  gene         = entrez_ids$ENTREZID,
  organism     = "hsa",
  pvalueCutoff = 0.05
)

head(ekegg)


dotplot(ekegg, showCategory = 15) +
  ggtitle("KEGG Pathway Enrichment - Day1 vs Day0")

umap_input <- t(ex)

dim(umap_input)

set.seed(123)  # biar hasil stabil

umap_result <- umap(umap_input)

umap_df <- data.frame(
  UMAP1 = umap_result$layout[,1],
  UMAP2 = umap_result$layout[,2],
  Day = gset$Day
)

head(umap_df)

library(ggplot2)

ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Day)) +
  geom_point(size = 2.5, alpha = 0.8) +
  theme_minimal() +
  ggtitle("UMAP Plot - Influenza Vaccine Response") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

gene_var <- apply(ex, 1, var)


top_var_genes <- names(sort(gene_var, decreasing = TRUE))[1:1000]

ex_var <- ex[top_var_genes, ]


umap_input2 <- t(ex_var)

set.seed(123)
umap_result2 <- umap(umap_input2)

umap_df2 <- data.frame(
  UMAP1 = umap_result2$layout[,1],
  UMAP2 = umap_result2$layout[,2],
  Day = gset$Day
)


ggplot(umap_df2, aes(x = UMAP1, y = UMAP2, color = Day)) +
  geom_point(size = 2.5, alpha = 0.8) +
  theme_minimal() +
  ggtitle("UMAP - Top 1000 Most Variable Genes") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

deg_genes <- rownames(results[results$adj.P.Val < 0.05, ])
length(deg_genes)

ex_deg <- ex[deg_genes, ]
umap_input3 <- t(ex_deg)

set.seed(123)
umap_result3 <- umap(umap_input3)

umap_df3 <- data.frame(
  UMAP1 = umap_result3$layout[,1],
  UMAP2 = umap_result3$layout[,2],
  Day = gset$Day
)

ggplot(umap_df3, aes(x = UMAP1, y = UMAP2, color = Day)) +
  geom_point(size = 2.5, alpha = 0.8) +
  theme_minimal() +
  ggtitle("UMAP - Significant DEGs Only")


