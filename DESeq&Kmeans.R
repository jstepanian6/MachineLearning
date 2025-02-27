if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
BiocManager::install("EnhancedVolcano")
install.packages("Hmisc")

library(EnhancedVolcano)
library(DESeq2)
library(ggplot2)
library(rgl)
library(plotly)
setwd("/home/jstepanian/Documents/BreastCancer/Step2_DESeq_filteredMatrix-20241216T140227Z-001/Step2_DESeq_filteredMatrix")
counts <- read.csv("geneCountMatrix_280Samples_antisense_guide_RNA_lncRNA_miRNA_pseudogenes_rRNA_sequencefeature_snoRNA_tRNA_uncharacterized.csv", sep = ",", row.names = 1)
head(counts)
metadata <- read.csv("SuplementaryMaterial_Metadata.tsv", sep = "\t")
# head(metadata) 
all(colnames(counts)==rownames(metadata))
ncol(counts) == nrow(colData) 

dds <- DESeqDataSetFromMatrix(countData=counts, 
                              colData=metadata, 
                                design=~ClusterJAN23_v2)
dds$Group <- relevel(dds$Cluster, ref = "Control")
vsd <- vst(dds, blind=FALSE)
vstNormalized <-(assay(vsd))

write.table(cbind(metadata, t(vstNormalized)), 'inputKmeans.tsv', sep="\t")
#Export PCA coordinates
rv <- rowVars(assay(vsd))
select <- order(rv, decreasing = TRUE)[seq_len(min(500,length(rv)))]
pca <- prcomp(t(assay(vsd)[select, ]))
percentVar <- pca$sdev^2/sum(pca$sdev^2)
d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], PC3=pca$x[,3])
ggplot(data = d, aes_string(x = "PC1", y = "PC2")) +  geom_point(size = 3) + xlab(paste0("PC1: ", round(percentVar[1] *100), "% variance")) + ylab(paste0("PC2: ", round(percentVar[2] *100), "% variance"))

plot_ly(x=pca$x[, 1], y=pca$x[, 2], z=pca$x[, 3], type="scatter3d", mode="markers")

k <-  3 # Number of clusters, you can adjust this number as needed
kmeans_result <- kmeans(d, centers = k)
d$cluster <- as.factor(kmeans_result$cluster)
plot_ly(x = pca$x[, 1], y = pca$x[, 2], z = pca$x[, 3], type = "scatter3d", mode = "markers", color = d$cluster)


plot_ly(data = d, x = ~PC1, y = ~PC2, z = ~PC3, color = ~cluster, type = "scatter3d", mode = "markers") %>%
  layout(scene = list(xaxis = list(title = paste0("PC1: ", round(percentVar[1] * 100), "% variance")),
                      yaxis = list(title = paste0("PC2: ", round(percentVar[2] * 100), "% variance")),
                      zaxis = list(title = paste0("PC3: ", round(percentVar[3] * 100), "% variance"))))
ggplot(data = d, aes(x = PC1, y = PC2, color = cluster)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) +
  ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) +
  scale_color_manual(values = rainbow(k)) +  # Change colors if needed
  theme_minimal()

k <- 3  # Number of clusters, you can adjust this number as needed
kmeans_result <- kmeans(d, centers = k)

# Add cluster labels to the data frame
d$cluster <- as.factor(kmeans_result$cluster)
metadata$kmeans_cluster <- as.factor(kmeans_result$cluster)
write.table(metadata, file ="SuplementaryMaterial_Metadata_withKmeansResults.tsv", sep="\t")


# Plot the PCA results with K-means colors and variance information

ggplot(data = d, aes(x = PC1, y = PC2, color = cluster)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) +
  ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) +
  scale_color_manual(values = rainbow(k)) +  # Change colors if needed
  theme_minimal()

write.table(d, "PCACoordinatesDESeq2.tsv", sep ="\t")

### HERE New Deseq using k-means as trait 
counts <- counts
head(counts)
metadata_k <- read.csv("SuplementaryMaterial_Metadata_withKmeansResults.tsv", sep = "\t")

head(metadata_k)
all(colnames(data)==rownames(metadata_k))
ncol(counts) == nrow(colData) 

dds <- DESeqDataSetFromMatrix(countData=counts, 
                              colData=metadata_k, 
                              design=~kmeans_cluster)
dds <- dds[ rowSums(counts(dds)) > 30, ]
#dds <- dds[apply(counts(dds), 1, function(x) all(x > 10)), ]
##Prefiltering more agressive - reduces from 53309 genes to 47345
smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]
sizeFactors(dds)

vsd <- vst(dds, blind=FALSE)
vstNormalized <-(assay(vsd))
svg("PCA_Kmeans_controls.svg")
plotPCA(vsd, intgroup=c("kmeans_cluster", "kmeans_cluster_control"))
dev.off()

dds <- DESeq(dds)
resultsNames(dds)

plotDispEsts(dds)


#Create constrast
C1C2 <- results(dds, contrast=c("kmeans_cluster","C1","C2"))
C1C3 <- results(dds, contrast=c("kmeans_cluster","C1","C3"))
C3C2 <- results(dds, contrast=c("kmeans_cluster","C3","C2"))

#PlotMA
svg("MA_C1_C2.svg")
plotMA(C1C2, ylim=c(-30,30), cex=.8)
dev.off()


svg("MA_C1_C3.svg")
plotMA(C1C3, ylim=c(-30,30), cex=.8)
dev.off()

svg("MA_C2_C3.svg")
plotMA(C3C2, ylim=c(-30,30), cex=.8)
dev.off()

#Volcano Plot
svg("volcano_C1_C2.svg")
pl <- EnhancedVolcano(C1C2,
                lab = NA,
                x = 'log2FoldChange',
                y = 'pvalue',
                xlab = bquote(~Log[2]~ 'fold change'),
                #pCutoff = 10e-32,
                FCcutoff = 15.0,
                pointSize = 1.0,
                labSize = 5.0,
                colAlpha = 1,
                legendPosition = 'bottom',title = "C1 vs C2",
                subtitle = "")

pl+ ggplot2::coord_cartesian(xlim=c(-30, 30)) 
pl + coord_flip()
dev.off()


svg("volcano_C1_C3.svg")
pl <- EnhancedVolcano(C1C3,
                      lab = NA,
                      x = 'log2FoldChange',
                      y = 'pvalue',
                      xlab = bquote(~Log[2]~ 'fold change'),
                      #pCutoff = 10e-32,
                      FCcutoff = 10.0,
                      pointSize = 1.5,
                      labSize = 5.0,
                      colAlpha = 1,
                      legendPosition = 'bottom',title = "C1 vs C3",
                      subtitle = "")
pl+ ggplot2::coord_cartesian(xlim=c(-25, 25)) 
pl + coord_flip()
dev.off()

svg("volcano_C3_C2.svg")
pl <- EnhancedVolcano(C3C2,
                      lab = NA,
                      x = 'log2FoldChange',
                      y = 'pvalue',
                      xlab = bquote(~Log[2]~ 'fold change'),
                      #pCutoff = 10e-32,
                      FCcutoff = 15.0,
                      pointSize = 1.5,
                      labSize = 5.0,
                      colAlpha = 1,
                      legendPosition = 'bottom',title = "C3 vs C2",
                      subtitle = "")
pl+ ggplot2::coord_cartesian(xlim=c(-25, 25)) 
pl + coord_flip()
dev.off()


#Export tables with results 
write.table(C1C2, "Contrast_C1_C2.tsv", sep="\t")
write.table(C1C3, "Contrast_C1_C3.tsv", sep="\t")
write.table(C2C3, "Contrast_C2_C3.tsv", sep="\t")

######################################################################
### HERE New Deseq creating contrast CX vs others  
counts <- counts
head(counts)
metadata_k <- read.csv("SuplementaryMaterial_Metadata_withKmeansResults_ClusterComparisson.csv", sep = "\t")

head(metadata_k)
all(colnames(data)==rownames(metadata_k))
ncol(counts) == nrow(colData) 

dds <- DESeqDataSetFromMatrix(countData=counts, 
                              colData=metadata_k, 
                              design=~C1vsOthers)
dds <- dds[ rowSums(counts(dds)) > 30, ]
#dds <- dds[apply(counts(dds), 1, function(x) all(x > 10)), ]
##Prefiltering more agressive - reduces from 53309 genes to 47345
smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]
sizeFactors(dds)

vsd <- vst(dds, blind=FALSE)
vstNormalized <-(assay(vsd))
svg("PCA_Kmeans_controls.svg")
plotPCA(vsd, intgroup=c("C1vsOthers"))
dev.off()
design <- ~ relevel(condition, ref = "Other")  

dds <- DESeq(dds)
resultsNames(dds)

plotDispEsts(dds)

res <- results(dds)
plotMA(res)

plotMA(res, ylim=c(-15,15), cex=.8)
#Volcano Plot
svg("volcano_C1_C2.svg")
pl <- EnhancedVolcano(res,
                      lab = NA,
                      x = 'log2FoldChange',
                      y = 'pvalue',
                      xlab = bquote(~Log[2]~ 'fold change'),
                      #pCutoff = 10e-32,
                      FCcutoff = 15.0,
                      pointSize = 1.0,
                      labSize = 5.0,
                      colAlpha = 1,
                      legendPosition = 'bottom',title = "C1 vs C2",
                      subtitle = "")

pl
#Export tables with results 
write.table(res, "Contrast_C1_others.tsv", sep="\t")

##Cluster2
counts <- counts
head(counts)
metadata_k <- read.csv("SuplementaryMaterial_Metadata_withKmeansResults_ClusterComparisson.csv", sep = "\t")

head(metadata_k)
all(colnames(data)==rownames(metadata_k))
ncol(counts) == nrow(colData) 

dds <- DESeqDataSetFromMatrix(countData=counts, 
                              colData=metadata_k, 
                              design=~C2vsOthers)
dds <- dds[ rowSums(counts(dds)) > 30, ]
#dds <- dds[apply(counts(dds), 1, function(x) all(x > 10)), ]
##Prefiltering more agressive - reduces from 53309 genes to 47345
smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]
sizeFactors(dds)

vsd <- vst(dds, blind=FALSE)
vstNormalized <-(assay(vsd))
plotPCA(vsd, intgroup=c("C2vsOthers"))

dds <- DESeq(dds)
resultsNames(dds)

plotDispEsts(dds)

res <- results(dds)
plotMA(res)

plotMA(res, ylim=c(-20,30), cex=.8)
#Volcano Plot
pl <- EnhancedVolcano(res,
                      lab = NA,
                      x = 'log2FoldChange',
                      y = 'pvalue',
                      xlab = bquote(~Log[2]~ 'fold change'),
                      #pCutoff = 10e-32,
                      FCcutoff = 15.0,
                      pointSize = 1.0,
                      labSize = 5.0,
                      colAlpha = 1,
                      legendPosition = 'bottom',title = "C1 vs C2",
                      subtitle = "")

pl
#Export tables with results 
write.table(res, "Contrast_C2_others.tsv", sep="\t")

##Cluster3
counts <- counts
head(counts)
metadata_k <- read.csv("SuplementaryMaterial_Metadata_withKmeansResults_ClusterComparisson.csv", sep = "\t")

head(metadata_k)
all(colnames(data)==rownames(metadata_k))
ncol(counts) == nrow(colData) 

dds <- DESeqDataSetFromMatrix(countData=counts, 
                              colData=metadata_k, 
                              design=~C3vsOthers)
dds <- dds[ rowSums(counts(dds)) > 30, ]
#dds <- dds[apply(counts(dds), 1, function(x) all(x > 10)), ]
##Prefiltering more agressive - reduces from 53309 genes to 47345
smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]
sizeFactors(dds)

vsd <- vst(dds, blind=FALSE)
vstNormalized <-(assay(vsd))
plotPCA(vsd, intgroup=c("C3vsOthers"))

dds <- DESeq(dds)
resultsNames(dds)

plotDispEsts(dds)

res <- results(dds)
plotMA(res)

plotMA(res, ylim=c(-30,30), cex=.8)
#Volcano Plot
pl <- EnhancedVolcano(res,
                      lab = NA,
                      x = 'log2FoldChange',
                      y = 'pvalue',
                      xlab = bquote(~Log[2]~ 'fold change'),
                      #pCutoff = 10e-32,
                      FCcutoff = 15.0,
                      pointSize = 1.0,
                      labSize = 5.0,
                      colAlpha = 1,
                      legendPosition = 'bottom',title = "C1 vs C2",
                      subtitle = "")

pl
#Export tables with results 
write.table(res, "Contrast_C3_others.tsv", sep="\t")




##Validate results

plotCounts(dds, gene="gene-TRPS1|TRPS1", intgroup="kmeans_cluster")
#DEFB130B
#KRTAP9-2

##Contrast of Cluster 1 vs Controls - in groups of 5
allSamples<- rownames(metadata_k)

controlSamples <- rownames(metadata_k[metadata_k[,"kmeans_cluster_control"]=="Control",])
c1samples<- rownames(metadata_k[metadata_k[,"kmeans_cluster_control"]=="C1",])

n=25

datalistUpC1Control =list ()
datalistDownC1Control =list ()

for (i in 1:n)
{
  sc1 <- sample(c1samples,size = 6, replace = F)
  selSamples <- sort(c(sc1,controlSamples))
  print(selSamples)
  selIdx <- allSamples %in% selSamples
  
  selCounts <- counts[,selIdx]
  selMetadata <- metadata_k[selIdx,]
  all(colnames(selCounts)==selMetadata$Sample)
  dds <- DESeqDataSetFromMatrix(countData=selCounts, 
                                colData=selMetadata, 
                                design=~kmeans_cluster_control)
  dds <- dds[ rowSums(counts(dds)) > 20, ]
  #Run DeSeq2
  dds <- DESeq(dds)
  resultsNames(dds)
  
  #Export results
  resControlvsC1 <- lfcShrink(dds, coef="kmeans_cluster_control_Control_vs_C1", type="normal", lfcThreshold=1)
  resT <- as.data.frame(resControlvsC1)
  indUpT <- resT[,2]>1 & !is.na(resT[,6]) & resT[,6]<0.1
  datalistUpC1Control[[i]] <- resT[indUpT,]
  indDownT <- resT[,2]< -1 & !is.na(resT[,6]) & resT[,6]<0.1
  datalistDownC1Control[[i]] <- resT[indDownT,]
  
}
consolidatedUpC1Control <- do.call(rbind, datalistUpC1Control)
consolidatedDownC1Control <- do.call(rbind, datalistDownC1Control)
write.table(consolidatedUpC1Control, file="Genes_Control_vs_C1_25Iter_Up.tsv", sep="\t", quote=F)
write.table(consolidatedDownC1Control, file="Genes_Control_vs_C1_25Iter_Down.tsv", sep="\t", quote=F)


plup <- EnhancedVolcano(consolidatedUpC1Control,
                      lab = NA,
                      x = 'log2FoldChange',
                      y = 'pvalue',
                      xlab = bquote(~Log[2]~ 'fold change'),
                      #pCutoff = 10e-32,
                      FCcutoff = 15.0,
                      pointSize = 1.5,
                      labSize = 5.0,
                      colAlpha = 1,
                      legendPosition = 'bottom',title = "consolidatedUpC1Control",
                      subtitle = "")
pldown <- EnhancedVolcano(consolidatedDownC1Control,
                        lab = NA,
                        x = 'log2FoldChange',
                        y = 'pvalue',
                        xlab = bquote(~Log[2]~ 'fold change'),
                        #pCutoff = 10e-32,
                        FCcutoff = 15.0,
                        pointSize = 1.5,
                        labSize = 5.0,
                        colAlpha = 1,
                        legendPosition = 'bottom',title = "consolidatedUpC1Control",
                        subtitle = "")
plup + pldown


dds <- dds[ rowSums(counts(dds)) > 30, ]
#dds <- dds[apply(counts(dds), 1, function(x) all(x > 10)), ]
##Prefiltering more agressive - reduces from 53309 genes to 47345
smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]
sizeFactors(dds)

vsd <- vst(dds, blind=FALSE)
vstNormalized <-(assay(vsd))
dds <- DESeq(dds)
resultsNames(dds)

plotDispEsts(dds)


#Create constrast
ControlC1 <- results(dds, contrast=c("kmeans_cluster_control","Control","C1"))
ControlC2 <- results(dds, contrast=c("kmeans_cluster_control","Control","C2"))
ControlC3 <- results(dds, contrast=c("kmeans_cluster_control","Control","C3"))

write.table(ControlC1, "Contrast_ControlC1.tsv", sep="\t")
write.table(ControlC2, "Contrast_ControlC2.tsv", sep="\t")
write.table(ControlC3, "Contrast_ControlC3.tsv", sep="\t")

EnhancedVolcano(ControlC1,
                lab = NA,
                x = 'log2FoldChange',
                y = 'pvalue',
                xlab = bquote(~Log[2]~ 'fold change'),
                #pCutoff = 10e-32,
                FCcutoff = 15.0,
                pointSize = 1.5,
                labSize = 5.0,
                colAlpha = 1,
                legendPosition = 'bottom',title = "Control C1",
                subtitle = "")

##Plot counts by ancestry 
counts <- counts
metadata_k <- metadata_k
head(metadata_k)
all(colnames(data)==rownames(metadata_k))
ncol(counts) == nrow(colData) 

dds <- DESeqDataSetFromMatrix(countData=counts, 
                              colData=metadata_k, 
                              design=~POP_NoMIX)
dds <- dds[ rowSums(counts(dds)) > 30, ]
#dds <- dds[apply(counts(dds), 1, function(x) all(x > 10)), ]
##Prefiltering more agressive - reduces from 53309 genes to 47345
smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]
sizeFactors(dds)

vsd <- vst(dds, blind=FALSE)
vstNormalized <-(assay(vsd))
plotPCA(vsd, intgroup=c("kmeans_cluster", "kmeans_cluster_control", "POP_NoMIX"))

dds <- DESeq(dds)
resultsNames(dds)




plotCounts(dds, gene="gene-ERBB2|ERBB2", intgroup="POP_NoMIX")
plotCounts(dds, gene="gene-NBPF12|NBPF12", intgroup="Status")
plotCounts(dds, gene="gene-SMG1P3|SMG1P3", intgroup="Status")
plotCounts(dds, gene="gene-RGPD6|RGPD6", intgroup="Status")
plotCounts(dds, gene="gene-NBPF10|NBPF10", intgroup="Status")
plotCounts(dds, gene="gene-RGPD8|RGPD8", intgroup="Status")
plotCounts(dds, gene="gene-SPDYE14P|SPDYE14P", intgroup="Status")
plotCounts(dds, gene="gene-SPDYE17|SPDYE17", intgroup="Status")
plotCounts(dds, gene="gene-SPDYE11|SPDYE11", intgroup="Status")
plotCounts(dds, gene="gene-RANBP2|RANBP2", intgroup="Status")
#p-value e-3**
plotCounts(dds, gene="gene-PMS2P9|PMS2P9", intgroup="Status")
plotCounts(dds, gene="gene-PMS2P2|PMS2P2", intgroup="Status")
plotCounts(dds, gene="gene-WASHC2A|WASHC2A", intgroup="Status")
plotCounts(dds, gene="gene-NBPF19|NBPF19", intgroup="Status")
#p-value e-2**
plotCounts(dds, gene="gene-NORAD|NORAD", intgroup="Cluster")
plotCounts(dds, gene="gene-NPIPB2|NPIPB2", intgroup="Cluster")
plotCounts(dds, gene="gene-CTAGE6|CTAGE6", intgroup="Cluster")
plotCounts(dds, gene="gene-LRRC37A4P-2|LRRC37A4P", intgroup="Cluster")
plotCounts(dds, gene="gene-WASHC2C|WASHC2C", intgroup="Cluster")
plotCounts(dds, gene="gene-LRRC37A2-3|LRRC37A2", intgroup="Cluster")



#RunDeseq
dds <- DESeq(dds)
res <- results(dds)
order(res$pvalue)
resOrdered <- res[order(res$pvalue),]
resOrdered <- na.omit(resOrdered)
resOrdered$pvalue
head(resOrdered)
#genes 34, 35, 36, 

write.csv(resOrdered[1:50, ], "DEG_raw.csv")
#p-value 0
plotCounts(dds, gene="gene-EIF4A1|EIF4A1", intgroup="Cluster")
plotCounts(dds, gene="gene-EEF1A1|EEF1A1", intgroup="Cluster")
plotCounts(dds, gene="gene-C18orf32|C18orf32", intgroup="Cluster")
plotCounts(dds, gene="gene-HNRNPA2B1|HNRNPA2B1", intgroup="Cluster")
#Results
resultsNames(dds)
#Contrast results Clusters 
Cluster0_vs_Cluster1 <- results(dds, contrast = c("Cluster", "Cluster_1", "Cluster_0"))

C0C1sorted <- Cluster0_vs_Cluster1[order(Cluster0_vs_Cluster1$pvalue),]
C0C1sorted <- na.omit(C0C1sorted)
C0C1sorted$pvalue
write.csv(resOrdered[1:250, ],file = "Cluster0_vs_Cluster1.csv")

 #los genes con p values pequeños son RNA o AS, el primer gen con p-valor pequeño y codificante es el 109 (ADAM50_8.3e-186)
plotCounts(dds, gene="gene-ADAM20|ADAM20", intgroup="Cluster")
plotCounts(dds, gene="gene-CAMKMT|CAMKMT", intgroup="Cluster")
plotCounts(dds, gene="gene-DEFB108B|DEFB108B", intgroup="Cluster")
plotCounts(dds, gene="gene-NOTCH2NLR|NOTCH2NLR", intgroup = "Cluster") #1.27315107495635e-149
#Contrast other clusters 
Cluster0_vs_Cluster2 <- results(dds, contrast = c("Cluster", "Cluster_0", "Cluster_2"))
Cluster1_vs_Cluster2 <- results(dds, contrast = c("Cluster", "Cluster_1", "Cluster_2"))
  #Sort Cluster0vsCluster1 
resOrderedC0_C1 <- Cluster0_vs_Cluster1[order(Cluster0_vs_Cluster1$pvalue),]
resOrderedC0_C1[(resOrderedC0_C1$pvalue<0.05),]
any(is.na(resOrderedC0_C1$pvalue))
resOrderedC0_C1 <- na.omit(resOrderedC0_C1)
resOrderedC0_C1_p0.05 <-resOrderedC0_C1[(resOrderedC0_C1$pvalue<0.05),]
resOrderedC0_C1_p0.05 <- resOrderedC0_C1_p0.05[order(resOrderedC0_C1_p0.05$pvalue),]
head(resOrderedC0_C1_p0.05, n=10)



#Plot
plotMA(Cluster0_vs_Cluster1, ylim=c(-12,30), main="Cluster0 vs Cluster1")
plotMA(Cluster0_vs_Cluster2, ylim=c(-30,30), main="Cluster 0 vs Cluster 2")
plotMA(Cluster1_vs_Cluster2, ylim=c(-20,20), main="Cluster 1 vs Cluster 2")

#write results 
write.csv(as.data.frame(), 
          file="condition_treated_results.csv")

#Volvano plot
res <- results(dds, contrast = c("Cluster", "Cluster_1", "Cluster_0"))
svglite("Cluster1_vs_cluster0.svg")
EnhancedVolcano(res,  lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Cluster 1 versus Cluster 0',
                pCutoff = 10e-16,
                FCcutoff = 15,
                pointSize = 1.0,
                labSize = 6.0,
                col=c('black', 'black', 'black', 'red3'),
                colAlpha = 1)
dev.off()
res2 <- results(dds, contrast = c("Cluster", "Cluster_1", "Cluster_2"))

svg("Cluster1_vs_Cluster2.svg")
EnhancedVolcano(res2,  lab = rownames(res2),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Cluster 1 versus Cluster 2',
                pCutoff = 10e-16,
                FCcutoff = 15,
                pointSize = 1.0,
                labSize = 6.0,
                col=c('black', 'black', 'black', 'red3'),
                colAlpha = 1)
dev.off()
res3 <- results(dds, contrast = c("Cluster", "Cluster_0", "Cluster_2"))
svg('Cluster0_vs_Cluster2.svg')
EnhancedVolcano(res3,  lab = rownames(res3),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Cluster 0 versus Cluster 2',
                pCutoff = 10e-16,
                FCcutoff = 15,
                pointSize = 1.0,
                labSize = 6.0,
                col=c('black', 'black', 'black', 'red3'),
                colAlpha = 1)
dev.off()



#reset par
par(mfrow=c(1,1))
# Make a basic volcano plot
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-20,20)))

# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<4.081299e-10 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

vsdata <- vst(dds, blind=FALSE)
svg("PCA_fondogris.svg")
plotPCA(vsdata, intgroup="cell.type")
dev.off()

#Plot MA
svg("PlotMA_fondogris.svg")
plotMA(res, ylim=c(-26,18))
dev.off()
##Normalización Shrinkage

resNorm <- lfcShrink(dds, coef=2, type="normal")
resAsh <- lfcShrink(dds, coef=2, type="ashr")
resApeglm <- lfcShrink(dds, coef=2, type="apeglm")
svg("ShrinkagePlots.svg")
par(mfrow=c(1,3), mar=c(4,4,2,1))
xlim <- c(200,1e6); ylim <- c(-20,20)
plotMA(resNorm, xlim=xlim, ylim=ylim, main="apeglm")
plotMA(resAsh, xlim=xlim, ylim=ylim, main="normal")
plotMA(resApeglm, xlim=xlim, ylim=ylim, main="ashr")
dev.off()