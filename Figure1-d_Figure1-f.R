library(tidyverse)
library(pheatmap)
library(DESeq2)
library(RColorBrewer)
library(PoiClaClu)
library(genefilter)
library(gplots)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(apeglm)
library(fdrtool) 

##################################################################################################
# processing: from count to deseq results
##################################################################################################
ddsMat <- DESeqDataSetFromMatrix(countData = countdata,
                                 colData = coldata,
                                 design = ~ condition )
ddsMat <- ddsMat[rowSums(counts(ddsMat)) > 1, ]
ddsMat <- estimateSizeFactors(ddsMat)
ddsMat$condition <- relevel(ddsMat$condition, ref="untreated")
ddsMat <- DESeq(ddsMat,betaPrior=T)

# rld 
rld <- rlog(ddsMat, blind=FALSE)

# deseq results
res <- results(ddsMat)
# add fdr padj to results
FDR.res <- fdrtool(res$stat, statistic= "normal", plot = T)
res[,"padj"]  <- p.adjust(FDR.res$pval, method = "BH")
res <- res[which(!is.na(res$pvalue) &!is.na(res$padj) & &res$baseMean > 0),] 

## add symbol and entrez names to the deseq2 results
res$symbol <- mapIds(org.Mm.eg.db,
                     keys=row.names(res),
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")
res$entrez <- mapIds(org.Mm.eg.db,
                     keys=row.names(res),
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")

##################################################################################################
# Figure1d: PCA plots
##################################################################################################
pcaData <- plotPCA(rld, intgroup = c( "condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color=condition)) + geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed()

##################################################################################################
# Figure1f: Visualization using heatmap plots: highlight top significantly DEG genes
##################################################################################################

resSig <- res[which(res$padj<0.01),]
resSigUp30 <- head(resSig[order(resSig$log2FoldChange, decreasing = TRUE),],30)
resSigDw30 <- tail(resSig[order(resSig$log2FoldChange, decreasing = TRUE),],30)
resSigTOP30 <- rbind(resSigUp30,resSigDw30)
sigGenes <- rownames(resSigTOP30[order(-resSigTOP30$log2FoldChange),])

mat.sig <- assay(rld)[sigGenes, ]
mat.sig <- mat.sig - rowMeans(mat.sig)

## annotation of heatmap
df <- data.frame(Group=factor(colData(rld)[, c("condition")]))
rownames(df) = colnames(mat.sig)

## plotting
pheatmap(mat.sig, 
         annotation=df,
         main="DEG: untreated vs TLI/ATS", 
         fontsize= 8, fontsize_row = 7, 
         cluster_rows=F, 
         cluster_cols=F, border_color = NA, 
         annotation_colors = list(Group = c(Neg="grey",Pos="black"))))

##################################################################################################
##################################################################################################
