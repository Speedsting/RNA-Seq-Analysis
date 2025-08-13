# Set up: install only if needed
if (!requireNamespace("BiocManager", quitely=TRUE))
  install.packages("BiocManager")
BiocManager::install(c("DESeq2", "edgeR", "limma", "cmapR", "data.table", "EnchancedVolcano"))

library(data.table)
library(cmapR)
library(DESeq2)
library(edgeR)
library(limma)
library(ggplot2)
library(pheatmap)
library(EnhancedVolcano)

gct_path <- "data/samples_raw.gct"

# parse_gctx from cmapR
gct <- parse_gctx(gct_path)
expr <- gct@mat
pheno <- as.data.frame(gct@cdesc)
pheno$disease <- factor(pheno$disease)
stopifnot("disease" %in% colnames(pheno))

dim(expr)   # briefly show dimensions
table(pheno$disease)

dds <- DESeqDataSetFromMatrix(
  countData = expr,
  colData = pheno,
  design = ~ disease
)
dds <- DESeq(dds)
res <- results(dds, contrast=c("disease", "LN", "healthy"))
res_ordered <- res[order(res$padj), ]
head(res_ordered, 20)

# Save DESeq2 results
write.csv(as.data.frame(res_ordered), file="results/deseq2_results.csv")

EnhancedVolcano(
  res,
  lab = rownames(res),
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.05,
  FCcutoff = 1,
  title = "Volcano: diseased vs healthy"
)
ggsave("results/volcano_disease_vs_healthy.png", width=8, height=6)

rld <- rlog(dds, blind=FALSE)
p <- plotPCA(rld, intgroup="disease") +
  ggtitle("PCA on rlog-transformed counts")
ggsave("results/pca_rlog.png", width=8, height=6)

topVarGenes <- head(order(rowVars(assay(rld)), decreasing=TRUE), 50)
mat <- assay(rld)[topVarGenes,]
annotation <- as.data.frame(colData(rld)["disease"])
gt <- pheatmap(mat, annotation_col=annotation, show_rownames=FALSE,
               main="Top 50 variable genes", filename="results/heatmap.png")$gtable

dge <- DGEList(counts=expr, samples=pheno)
keep <- filterByExpr(dge, design=model.matrix(~ disease, data=pheno))
dge <- dge[keep, , keep.lib.sizes=FALSE]
dge <- calcNormFactors(dge)
design <- model.matrix(~ disease, data=pheno)
dge <- estimateDisp(dge, design)

fit <- glmQLFit(dge, design)
qlf <- glmQLFTest(fit, coef=2)
edgeR_table <- topTags(qlf, n=Inf)$table
write.csv(edgeR_table, "results/edgeR_results.csv")

v <- voom(dge, design, plot=TRUE, save.plot=TRUE)
voom_df <- data.frame(
  gene=rownames(v$E),
  x=v$voom.xy$x,
  y=v$voom.xy$y,
  stringsAsFactors=FALSE
)

top50_idx <- order(voom_df$y, decreasing=TRUE)[1:50]
top50 <- voom_df[top50_idx, ]
ggplot(voom_df, aes(x=x, y=y)) +
  geom_point(alpha=0.3) +
  geom_point(data=top50, aes(x=x, y=y), color="red") +
  geom_text(data=top50, aes(x=x, y=y, label=gene), vjust=-0.5, size=2) +
  labs(x="Average log-CPM",
       y="Sqrt residual variance",
       title="Mean-variance trend with top 50 genes") +
  theme_minimal()
ggsave("results/top50_idx.png", width=8, height=6)


fit2 <- lmFit(v, design)
fit2 <- eBayes(fit2)
limma_table <- topTable(fit2, coef=2, number=Inf)
write.csv(limma_table, "results/limma_voom_results.csv")

deseq_genes <- rownames(subset(res, padj < 0.05))
edgeR_genes <- rownames(topTags(qlf, n=Inf)$table)[topTags(qlf, n=Inf)$table$FDR < 0.05]
limma_genes <- rownames(topTable(fit2, coef=2, p.value=0.05))

intersect(deseq_genes, edgeR_genes)
intersect(deseq_genes, limma_genes)
intersect(edgeR_genes, limma_genes)

# Training and testing data w/Machine Learning
Col_mat = as.matrix(t(v$E)) # x data
Col_mat = varFilter(v$E, var.func=IQR, var.cutoff=0.95)
Col_mat = t(Col_mat)

Col_resp = as.factor(v$targets$disease) # y data
Col_mat = as.matrix(t(v$E))
Col_resp = as.factor(v$targets$disease)

# Split the data into train and testing sets
set.seed(42)
train_data = createDataPartition(Col_resp, p=0.75, list=FALSE)

x_train_data = Col_mat[train_data, ]
x_test_data = Col_mat[-train_data, ]
y_train_data = Col_resp[train_data]
y_test_data = Col_resp[-train_data]

prop.table(table(y_train_data))
prop.table(table(y_test_data))

res_data = cv.glmnet(
  x = x_train_data,
  y = y_train_data,
  alpha = 0.5,
  family = "binomial"
)

y_pred_data = predict(
  res_data,
  newx = x_test_data,
  type = "class",
  s = "lambda.min"
)
confusion_matrix = table(y_pred_data, y_test_data)

print(confusion_matrix)
print(paste0("Sensitivity: ", sensitivity(confusion_matrix)))
print(paste0("Specificity: ", specificity(confusion_matrix)))
print(paste0("Precision: ", precision(onfusion_matrix)))
print(paste0("kappa: ", kappa(confusion_matrix)))
print(paste0("negPredValue: ", negPredValue(confusion_matrix)))
print(paste0("posPredValue: ", posPredValue(confusion_matrix)))

sessionInfo()
writeLines(capture.output(sessionInfo()), "results/session_info.txt")