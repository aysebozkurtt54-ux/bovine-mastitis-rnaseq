# DESeq2: enfekte (E. coli) ve kontrol meme çeyreği, 24. saat, eşleşmiş tasarım (~ inek + durum)
# Kullanım: Rscript 02_deseq2.R <klasör>   (klasörde counts.csv ve coldata.csv olmalı)
args <- commandArgs(trailingOnly = TRUE)
outdir <- args[1]
suppressPackageStartupMessages({
  library(DESeq2); library(ggplot2); library(pheatmap)
})

counts  <- as.matrix(read.csv(file.path(outdir, "counts.csv"), row.names = 1, check.names = FALSE))
coldata <- read.csv(file.path(outdir, "coldata.csv"), row.names = 1, stringsAsFactors = FALSE)
counts  <- counts[, rownames(coldata)]
mode(counts) <- "integer"
coldata$cow       <- factor(coldata$cow)
coldata$condition <- factor(coldata$condition, levels = c("control", "infected"))

dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ cow + condition)

# Çok düşük sayımlı genleri at (en az 3 örnekte >= 10 okuma)
keep <- rowSums(counts(dds) >= 10) >= 3
dds  <- dds[keep, ]
cat("Analize giren gen sayısı:", nrow(dds), "\n")

dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "infected", "control"), alpha = 0.05)
cat("\n")
summary(res)

tab <- as.data.frame(res)
tab$gene <- rownames(tab)
tab <- tab[order(tab$padj), c("gene", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
write.csv(tab, file.path(outdir, "deseq2_results_all.csv"), row.names = FALSE)
sig <- subset(tab, !is.na(padj) & padj < 0.05 & abs(log2FoldChange) >= 1)
write.csv(sig, file.path(outdir, "deseq2_significant.csv"), row.names = FALSE)
write.csv(as.data.frame(counts(dds, normalized = TRUE)), file.path(outdir, "normalized_counts.csv"))
cat("\nAnlamlı gen (padj < 0.05 ve |log2FC| >= 1):", nrow(sig),
    " (artan:", sum(sig$log2FoldChange > 0), " azalan:", sum(sig$log2FoldChange < 0), ")\n")

# --- Grafikler ---
vsd <- vst(dds, blind = FALSE)

pca <- plotPCA(vsd, intgroup = c("condition", "cow"), returnData = TRUE)
pv  <- round(100 * attr(pca, "percentVar"))
p1 <- ggplot(pca, aes(PC1, PC2, color = condition, shape = cow)) +
  geom_point(size = 4) +
  labs(x = paste0("PC1: %", pv[1]), y = paste0("PC2: %", pv[2]),
       title = "PCA (vst)", color = "Condition", shape = "Cow") +
  theme_bw()
ggsave(file.path(outdir, "pca.png"), p1, width = 6.5, height = 5, dpi = 150)

vd <- tab[!is.na(tab$padj), ]
vd$grup <- "not significant"
vd$grup[vd$padj < 0.05 & vd$log2FoldChange >=  1] <- "up"
vd$grup[vd$padj < 0.05 & vd$log2FoldChange <= -1] <- "down"
p2 <- ggplot(vd, aes(log2FoldChange, -log10(padj), color = grup)) +
  geom_point(alpha = 0.6, size = 1) +
  scale_color_manual(values = c("not significant" = "grey70", "up" = "firebrick", "down" = "steelblue")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  labs(title = "Infected vs control (24 h)", x = "log2 fold change", y = "-log10(padj)", color = NULL) +
  theme_bw()
ggsave(file.path(outdir, "volcano.png"), p2, width = 6.5, height = 5, dpi = 150)

top <- head(tab$gene[!is.na(tab$padj)], 30)
mat <- assay(vsd)[top, , drop = FALSE]
ann <- data.frame(condition = coldata$condition, cow = coldata$cow, row.names = rownames(coldata))
pheatmap(mat, scale = "row", annotation_col = ann, show_rownames = TRUE, fontsize_row = 7,
         filename = file.path(outdir, "heatmap_top30.png"), width = 7, height = 8)

writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo.txt"))
cat("\nBitti. Çıktılar:", outdir, "\n")
