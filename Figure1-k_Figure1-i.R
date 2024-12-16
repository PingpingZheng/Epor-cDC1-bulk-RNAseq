library(tidyverse)
library(pheatmap)
library(DESeq2)
library(RColorBrewer)
library(PoiClaClu)
library(genefilter)
library(gplots)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(EnhancedVolcano)
library(apeglm)
library(fdrtool) 

##################################################################################################
# processing: from count to deseq results
##################################################################################################
ddsMat <- DESeqDataSetFromMatrix(countData = countdata,
                                 colData = coldata,
                                 design = ~ GroupID)
ddsMat <- estimateSizeFactors(ddsMat)
ddsMat <- ddsMat[rowSums( counts(ddsMat, normalized=TRUE) ) >= 1,]
ddsMat <- DESeq(ddsMat)

# rld 
rld <- rlog(ddsMat, blind=FALSE)

# deseq results
res <- results(ddsMat, contrast = c("GroupID","Pos","Neg"), pAdjustMethod = "BH")
# add fdr padj to results
FDR.res <- fdrtool(res$stat, statistic= "normal", plot = T)
res[,"padjFdr"]  <- p.adjust(FDR.res$pval, method = "BH")
res <- res[which(!is.na(res$pvalue) &!is.na(res$padjFdr) & &res$baseMean > 0),] 

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
# Figure 1K: Visualization using volcano plots: highlight epor related genes
##################################################################################################

selected.upgenes <- c("Epor","Siglece","Vsig4","Tgfb2","Arg1","Aldh1a1","Axl","Gas6","Pros1","Lrp1","Timd4","Cd5l","Fabp4","Lpl","Ptgs1","Pparg","Nr1h3","Apoe","Abca1","Hif1a","Epas1","Hmox1","Slc40a1")
selected.downgenes <- c("Havcr2","Notch4","Fndc7","Pdia5","Ppt1","Plcb4","Clnk","Casp6","Naaa","Fam149a","Hepacam2","Rab7b","Wdfy4","Ciita","Tlr3","Tlr11")

keyvals <- ifelse(
  res$log2FoldChange < 0 & res$symbol %in%  c(selected.upgenes,selected.downgenes), 'royalblue',
  ifelse(res$log2FoldChange > 0 & res$symbol %in%  c(selected.upgenes,selected.downgenes), 'red4',
         'grey60'))
keyvals[is.na(keyvals)] <- 'grey60'
names(keyvals)[keyvals == 'red4'] <- 'high'
names(keyvals)[keyvals == 'grey60'] <- 'others'
names(keyvals)[keyvals == 'royalblue'] <- 'low'
EnhancedVolcano(res,
                lab = res$symbol,
                x = 'log2FoldChange',
                y = 'padjFdr',
                selectLab = c(selected.upgenes,selected.downgenes),
                boxedLabels = TRUE,
                drawConnectors = TRUE,
                max.overlaps = Inf,
                pointSize = c(ifelse(res$symbol %in% c(selected.upgenes,selected.downgenes), 4, 2)),
                labSize = 3.0,
                ylim = c(0,15),
                borderWidth = 1.5,
                borderColour = 'black',
                #shape = 21,
                pCutoff = 0.0,
                FCcutoff = 12,
                colAlpha = 0.8,
                ylab = bquote(~-Log[10] ~ italic(padj)),
                colCustom = keyvals,
                title = 'Selected Differentially Expressed Genes',
                subtitle = "Epor-tdT+ vs Epor-tdT- XCR1+CD8a+cDC1s"
)
##################################################################################################
# Figure1I: Visualization using heatmap plots: highlight epor related genes
##################################################################################################

## selected genes
selectedGenes <- c(selected.upgenes,selected.downgenes)
resSelected <-  res[res$symbol %in% selectedGenes,]
## selected matrix
mat.sig <- assay(rld)[rownames(resSelected), ]
mat.sig <- mat.sig - rowMeans(mat.sig)

## change rowname of matrix data from esembl to symb
all.equal(rownames(resSelected),rownames(mat.sig))
rownames(mat.sig) <- resSelected$symbol
mat.sig <- mat.sig[match(selectedGenes, rownames(mat.sig)),]

## annotation of heatmap
df <- data.frame(Group=factor(colData(rld)[, c("GroupID")]))
rownames(df) = colnames(mat.sig)

## plotting
pheatmap(mat.sig, 
         annotation=df,
         main="Selected DEG: Epor-tdT+ vs Epor-tdT- XCR1+CD8a+cDC1s", 
         fontsize= 8, fontsize_row = 7, 
         cluster_rows=F, 
         cluster_cols=F, border_color = NA, 
         annotation_colors = list(Group = c(Neg="grey",Pos="black"))))
##################################################################################################
##
##################################################################################################
