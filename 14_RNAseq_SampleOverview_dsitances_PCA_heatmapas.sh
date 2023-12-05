##14. Further exploration of the RA-seq data; you should start (continue) from #13. (DifferentialGeneExpressionWithDEseq2.sh)

library("vsn")
library("pheatmap")
library("RColorBrewer")


vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)





# Principal component plot of the samples
pdf("/home/yaaroba/ENCODEproject/figures/PCA.pdf", height=6, width=6)
plotPCA(vsd, intgroup=c("condition", "type"))
dev.off()






# Heatmap of the sample-to-sample distances:
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$type, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

pdf("/home/yaaroba/ENCODEproject/figures/AllSampleClustering.pdf", height=6, width=8)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.off()







# Heatmap of the count matrix - 20 top expressed genes
# this gives log2(n + 1)
ntd <- normTransform(dds)

select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("condition","type")])

pdf("/home/yaaroba/ENCODEproject/figures/MatrixTop20highlyExpressedAcrossSamples.pdf", height=6, width=6)
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)
dev.off()







