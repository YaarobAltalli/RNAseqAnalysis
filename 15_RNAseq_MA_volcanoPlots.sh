##15 MA & volcano plots:


# for E11.5 vs E14.5 of Forebrain:
# creating "DESeqDataSet" object for estimating sizing factors and further differential gee expression analysis
ddsFBe11_FBe14 = DESeqDataSetFromMatrix(countData = countData_FBe11_FBe14, colData = coldata_FBe11_FBe14, design = ~ condition)
ddsFBe11_FBe14<-estimateSizeFactors(ddsFBe11_FBe14)
SF_FBe11_FBe14<-sizeFactors(ddsFBe11_FBe14)
ddsFBe11_FBe14 <- DESeq(ddsFBe11_FBe14)
res_FBe11_FBe14 <- results(ddsFBe11_FBe14)
resFBe11_FBe14_Df<-as.data.frame(res_FBe11_FBe14, stringsAsFactors=FALSE)
resFBe11_FBe14_Df<- cbind(annotatedNormCounts,resFBe11_FBe14_Df)   

res_FBe11_FBe14_05 <- results(ddsFBe11_FBe14, alpha=0.05) # to inlude only differentially expressed genes with padj<0.05
summary(res_FBe11_FBe14_05) # summary; outputs pasted below workflow

# generating shrunken log2 fold changes for better visualization:
resultsNames(ddsFBe11_FBe14)
> resultsNames(ddsFBe11_FBe14)
[1] "Intercept"                    "condition_FBe14.5_vs_FBe11.5"

resLFC_FBe11_FBe14 <- lfcShrink(ddsFBe11_FBe14, coef="condition_FBe14.5_vs_FBe11.5", type="apeglm")


head(resFBe11_FBe14_Df[, 14:17], 4)


resFBe11_FBe14_Df$FBe11mean<-rowMeans(resFBe11_FBe14_Df[,10:13])
resFBe11_FBe14_Df$FBe14mean<-rowMeans(resFBe11_FBe14_Df[,14:17])
resFBe11_FBe14_Df$FBe11_FBe14mean<-rowMeans(resFBe11_FBe14_Df[,10:17])
resFBe11_FBe14_Df$LFC_log2FoldChange<-resLFC_FBe11_FBe14$log2FoldChange
resFBe11_FBe14_Df$LFC_padj<-resLFC_FBe11_FBe14$padj
resFBe11_FBe14_Df_fig<-resFBe11_FBe14_Df[which(is.na(resFBe11_FBe14_Df$log2FoldChange)==FALSE),]
resFBe11_FBe14_Df_fig$padj[which(is.na(resFBe11_FBe14_Df_fig$padj)==TRUE)]<-1 

resFBe11_FBe14_Df_fig$col<- "gray20"

colFun<- function(i,tab){
if (tab$padj[i]<=0.05 ) a<- "red"
else if (tab$padj[i]>0.05 & tab$padj[i]<=0.1) a<- "green"
else if (tab$padj[i]>0.1 & tab$padj[i]<1 ) a<- "gray20"
else if (tab$padj[i]==1) a<- "white" # making data points with padj=NA (1) "invisible"
else a<-"gray20"
return(a)
}

resFBe11_FBe14_Df_fig$col<-unlist(lapply(1:length(resFBe11_FBe14_Df_fig[,1]), colFun, resFBe11_FBe14_Df_fig))

## MA plots

# regular log2 fold change
pdf("/home/yaaroba/ENCODEproject/figures/FBe11_FBe14_RNAseq_MA_plot.pdf", height=6, width=6)
par(mar=c(3,3,3,3)+0.1,mgp=c(1.7,0.6,0))
plot(resFBe11_FBe14_Df_fig$FBe11_FBe14mean, resFBe11_FBe14_Df_fig$log2FoldChange, pch=20, cex=0.2, col=resFBe11_FBe14_Df_fig$col, cex.main=1.2, xlab="mean normalized counts", ylab="log2FoldChange", xlim=c(0,10000), ylim=c(-5,5),
main="Gene expression in forebrain of E11.5 and E14.5 mouse embryos",bty="n", cex.main=1)
abline(h = 0, col = "black")
dev.off()

# shrunken regular log2 fold change
pdf("/home/yaaroba/ENCODEproject/figures/FBe11_FBe14_RNAseq_MA_plot_shrunken_log2foldChanges.pdf", height=6, width=6)
par(mar=c(3,3,3,3)+0.1,mgp=c(1.7,0.6,0))
plot(resFBe11_FBe14_Df_fig$FBe11_FBe14mean, resFBe11_FBe14_Df_fig$LFC_log2FoldChange, pch=20, cex=0.2, col=resFBe11_FBe14_Df_fig$col, cex.main=1.2, xlab="mean normalized counts", ylab="log2FoldChange", xlim=c(0,10000), ylim=c(-5,5),
main="Gene expression in forebrain of E11.5 and E14.5 mouse embryos",bty="n", cex.main=1)
abline(h = 0, col = "black")
dev.off()

## Vulcano plots

# regular log2 fold change
pdf("/home/yaaroba/ENCODEproject/figures/FBe11_FBe14_RNAseq_volcano_plot.pdf", height=6, width=6)
# regular log2 fold change
par(mar=c(3,3,3,3)+0.1,mgp=c(1.7,0.6,0))
plot(resFBe11_FBe14_Df_fig$log2FoldChange, -log10(resFBe11_FBe14_Df_fig$padj), pch=20, cex=0.2, col=resFBe11_FBe14_Df_fig$col, cex.main=1.2, xlab="Fold Change: log2(E14.5/E11.5)", ylab="-log10(padj)", xlim=c(-10,10), ylim=c(0.5,15), asp=2,
main="Gene expression in forebrain of E11.5 and E14.5  mouse embryos",bty="n", cex.main=1)
dev.off()

pdf("/home/yaaroba/ENCODEproject/figures/FBe11_FBe14_RNAseq_volcano_plot_shrunken_log2foldChanges.pdf", height=6, width=6)
# shrunken log2 fold change
par(mar=c(3,3,3,3)+0.1,mgp=c(1.7,0.6,0))
plot(resFBe11_FBe14_Df_fig$LFC_log2FoldChange, -log10(resFBe11_FBe14_Df_fig$padj), pch=20, cex=0.2, col=resFBe11_FBe14_Df_fig$col, cex.main=1.2, xlab="Fold Change: log2(E14.5/E11.5)", ylab="-log10(padj)", xlim=c(-10,10), ylim=c(0.1,15),
main="Gene expression in forebrain of E11.5 and E14.5 mouse embryos",bty="n", cex.main=1)
dev.off()






# for E14.5 vs P0 of Forebrain:
# creating "DESeqDataSet" object for estimating sizing factors and further differential gee expression analysis
ddsFBe14_FBp0 = DESeqDataSetFromMatrix(countData = countData_FBe14_FBp0, colData = coldata_FBe14_FBp0, design = ~ condition)
ddsFBe14_FBp0<-estimateSizeFactors(ddsFBe14_FBp0)
SF_FBe14_FBp0<-sizeFactors(ddsFBe14_FBp0)
ddsFBe14_FBp0 <- DESeq(ddsFBe14_FBp0)
res_FBe14_FBp0 <- results(ddsFBe14_FBp0)
resFBe14_FBp0_Df<-as.data.frame(res_FBe14_FBp0, stringsAsFactors=FALSE)
resFBe14_FBp0_Df<- cbind(annotatedNormCounts,resFBe14_FBp0_Df)   

res_FBe14_FBp0_05 <- results(ddsFBe14_FBp0, alpha=0.05) # to inlude only differentially expressed genes with padj<0.05
summary(res_FBe14_FBp0_05) # summary; outputs pasted below workflow

# generating shrunken log2 fold changes for better visualization:
resultsNames(ddsFBe14_FBp0)
> resultsNames(ddsFBe14_FBp0)
[1] "Intercept"                    "condition_FBp0_vs_FBe14.5"

resLFC_FBe14_FBp0 <- lfcShrink(ddsFBe14_FBp0, coef="condition_FBp0_vs_FBe14.5", type="apeglm")


head(resFBe14_FBp0_Df[, c(14:19, 24:25)], 2)


resFBe14_FBp0_Df$FBe14mean<-rowMeans(resFBe14_FBp0_Df[,14:17])
resFBe14_FBp0_Df$FBp0mean<-rowMeans(resFBe14_FBp0_Df[,c(18:19, 24:25)])
resFBe14_FBp0_Df$FBe14_FBp0mean<-rowMeans(resFBe14_FBp0_Df[,c(14:19, 24:25)])
resFBe14_FBp0_Df$LFC_log2FoldChange<-resLFC_FBe14_FBp0$log2FoldChange
resFBe14_FBp0_Df$LFC_padj<-resLFC_FBe14_FBp0$padj
resFBe14_FBp0_Df_fig<-resFBe14_FBp0_Df[which(is.na(resFBe14_FBp0_Df$log2FoldChange)==FALSE),]
resFBe14_FBp0_Df_fig$padj[which(is.na(resFBe14_FBp0_Df_fig$padj)==TRUE)]<-1 

resFBe14_FBp0_Df_fig$col<- "gray20"

colFun<- function(i,tab){
if (tab$padj[i]<=0.05 ) a<- "red"
else if (tab$padj[i]>0.05 & tab$padj[i]<=0.1) a<- "green"
else if (tab$padj[i]>0.1 & tab$padj[i]<1 ) a<- "gray20"
else if (tab$padj[i]==1) a<- "white" # making data points with padj=NA (1) "invisible"
else a<-"gray20"
return(a)
}

resFBe14_FBp0_Df_fig$col<-unlist(lapply(1:length(resFBe14_FBp0_Df_fig[,1]), colFun, resFBe14_FBp0_Df_fig))

## MA plots

# regular log2 fold change
pdf("/home/yaaroba/ENCODEproject/figures/FBe14_FBp0_RNAseq_MA_plot.pdf", height=6, width=6)
par(mar=c(3,3,3,3)+0.1,mgp=c(1.7,0.6,0))
plot(resFBe14_FBp0_Df_fig$FBe14_FBp0mean, resFBe14_FBp0_Df_fig$log2FoldChange, pch=20, cex=0.2, col=resFBe14_FBp0_Df_fig$col, cex.main=1.2, xlab="mean normalized counts", ylab="log2FoldChange", xlim=c(0,10000), ylim=c(-5,5),
main="Gene expression in forebrain of e14.5 and p0 mouse embryos",bty="n", cex.main=1)
abline(h = 0, col = "black")
dev.off()

# shrunken regular log2 fold change
pdf("/home/yaaroba/ENCODEproject/figures/FBe14_FBp0_RNAseq_MA_plot_shrunken_log2foldChanges.pdf", height=6, width=6)
par(mar=c(3,3,3,3)+0.1,mgp=c(1.7,0.6,0))
plot(resFBe14_FBp0_Df_fig$FBe14_FBp0mean, resFBe14_FBp0_Df_fig$LFC_log2FoldChange, pch=20, cex=0.2, col=resFBe14_FBp0_Df_fig$col, cex.main=1.2, xlab="mean normalized counts", ylab="log2FoldChange", xlim=c(0,10000), ylim=c(-5,5),
main="Gene expression in forebrain of e14.5 and p0 mouse embryos",bty="n", cex.main=1)
abline(h = 0, col = "black")
dev.off()

## Vulcano plots

# regular log2 fold change
pdf("/home/yaaroba/ENCODEproject/figures/FBe14_FBp0_RNAseq_volcano_plot.pdf", height=6, width=6)
# regular log2 fold change
par(mar=c(3,3,3,3)+0.1,mgp=c(1.7,0.6,0))
plot(resFBe14_FBp0_Df_fig$log2FoldChange, -log10(resFBe14_FBp0_Df_fig$padj), pch=20, cex=0.2, col=resFBe14_FBp0_Df_fig$col, cex.main=1.2, xlab="Fold Change: log2(p0.5/e14.5)", ylab="-log10(padj)", xlim=c(-10,10), ylim=c(0.5,15), asp=2,
main="Gene expression in forebrain of e14.5 and p0 mouse embryos",bty="n", cex.main=1)
dev.off()

pdf("/home/yaaroba/ENCODEproject/figures/FBe14_FBp0_RNAseq_volcano_plot_shrunken_log2foldChanges.pdf", height=6, width=6)
# shrunken log2 fold change
par(mar=c(3,3,3,3)+0.1,mgp=c(1.7,0.6,0))
plot(resFBe14_FBp0_Df_fig$LFC_log2FoldChange, -log10(resFBe14_FBp0_Df_fig$padj), pch=20, cex=0.2, col=resFBe14_FBp0_Df_fig$col, cex.main=1.2, xlab="Fold Change: log2(p0.5/e14.5)", ylab="-log10(padj)", xlim=c(-10,10), ylim=c(0.1,15),
main="Gene expression in forebrain of e14.5 and p0 mouse embryos",bty="n", cex.main=1)
dev.off()







# for e11.5 vs P0 of Forebrain:
# creating "DESeqDataSet" object for estimating sizing factors and further differential gee expression analysis
ddsFBe11_FBp0 = DESeqDataSetFromMatrix(countData = countData_FBe11_FBp0, colData = coldata_FBe11_FBp0, design = ~ condition)
ddsFBe11_FBp0<-estimateSizeFactors(ddsFBe11_FBp0)
SF_FBe11_FBp0<-sizeFactors(ddsFBe11_FBp0)
ddsFBe11_FBp0 <- DESeq(ddsFBe11_FBp0)
res_FBe11_FBp0 <- results(ddsFBe11_FBp0)
resFBe11_FBp0_Df<-as.data.frame(res_FBe11_FBp0, stringsAsFactors=FALSE)
resFBe11_FBp0_Df<- cbind(annotatedNormCounts,resFBe11_FBp0_Df)   

res_FBe11_FBp0_05 <- results(ddsFBe11_FBp0, alpha=0.05) # to inlude only differentially expressed genes with padj<0.05
summary(res_FBe11_FBp0_05) # summary; outputs pasted below workflow

# generating shrunken log2 fold changes for better visualization:
resultsNames(ddsFBe11_FBp0)
> resultsNames(ddsFBe11_FBp0)
[1] "Intercept"                    "condition_FBp0_vs_FBe11.5"

resLFC_FBe11_FBp0 <- lfcShrink(ddsFBe11_FBp0, coef="condition_FBp0_vs_FBe11.5", type="apeglm")


head(resFBe11_FBp0_Df[, c(10:13 ,18:19, 24:25)], 2)


resFBe11_FBp0_Df$FBe11mean<-rowMeans(resFBe11_FBp0_Df[,10:13])
resFBe11_FBp0_Df$FBp0mean<-rowMeans(resFBe11_FBp0_Df[,c(18:19, 24:25)])
resFBe11_FBp0_Df$FBe11_FBp0mean<-rowMeans(resFBe11_FBp0_Df[,c(10:13, 18:19, 24:25)])
resFBe11_FBp0_Df$LFC_log2FoldChange<-resLFC_FBe11_FBp0$log2FoldChange
resFBe11_FBp0_Df$LFC_padj<-resLFC_FBe11_FBp0$padj
resFBe11_FBp0_Df_fig<-resFBe11_FBp0_Df[which(is.na(resFBe11_FBp0_Df$log2FoldChange)==FALSE),]
resFBe11_FBp0_Df_fig$padj[which(is.na(resFBe11_FBp0_Df_fig$padj)==TRUE)]<-1 

resFBe11_FBp0_Df_fig$col<- "gray20"

colFun<- function(i,tab){
if (tab$padj[i]<=0.05 ) a<- "red"
else if (tab$padj[i]>0.05 & tab$padj[i]<=0.1) a<- "green"
else if (tab$padj[i]>0.1 & tab$padj[i]<1 ) a<- "gray20"
else if (tab$padj[i]==1) a<- "white" # making data points with padj=NA (1) "invisible"
else a<-"gray20"
return(a)
}

resFBe11_FBp0_Df_fig$col<-unlist(lapply(1:length(resFBe11_FBp0_Df_fig[,1]), colFun, resFBe11_FBp0_Df_fig))

## MA plots

# regular log2 fold change
pdf("/home/yaaroba/ENCODEproject/figures/FBe11_FBp0_RNAseq_MA_plot.pdf", height=6, width=6)
par(mar=c(3,3,3,3)+0.1,mgp=c(1.7,0.6,0))
plot(resFBe11_FBp0_Df_fig$FBe11_FBp0mean, resFBe11_FBp0_Df_fig$log2FoldChange, pch=20, cex=0.2, col=resFBe11_FBp0_Df_fig$col, cex.main=1.2, xlab="mean normalized counts", ylab="log2FoldChange", xlim=c(0,10000), ylim=c(-5,5),
main="Gene expression in forebrain of e11.5 and p0 mouse embryos",bty="n", cex.main=1)
abline(h = 0, col = "black")
dev.off()

# shrunken regular log2 fold change
pdf("/home/yaaroba/ENCODEproject/figures/FBe11_FBp0_RNAseq_MA_plot_shrunken_log2foldChanges.pdf", height=6, width=6)
par(mar=c(3,3,3,3)+0.1,mgp=c(1.7,0.6,0))
plot(resFBe11_FBp0_Df_fig$FBe11_FBp0mean, resFBe11_FBp0_Df_fig$LFC_log2FoldChange, pch=20, cex=0.2, col=resFBe11_FBp0_Df_fig$col, cex.main=1.2, xlab="mean normalized counts", ylab="log2FoldChange", xlim=c(0,10000), ylim=c(-5,5),
main="Gene expression in forebrain of e11.5 and p0 mouse embryos",bty="n", cex.main=1)
abline(h = 0, col = "black")
dev.off()

## Vulcano plots

# regular log2 fold change
pdf("/home/yaaroba/ENCODEproject/figures/FBe11_FBp0_RNAseq_volcano_plot.pdf", height=6, width=6)
# regular log2 fold change
par(mar=c(3,3,3,3)+0.1,mgp=c(1.7,0.6,0))
plot(resFBe11_FBp0_Df_fig$log2FoldChange, -log10(resFBe11_FBp0_Df_fig$padj), pch=20, cex=0.2, col=resFBe11_FBp0_Df_fig$col, cex.main=1.2, xlab="Fold Change: log2(p0.5/e11.5)", ylab="-log10(padj)", xlim=c(-10,10), ylim=c(0.5,15), asp=2,
main="Gene expression in forebrain of e11.5 and p0 mouse embryos",bty="n", cex.main=1)
dev.off()

pdf("/home/yaaroba/ENCODEproject/figures/FBe11_FBp0_RNAseq_volcano_plot_shrunken_log2foldChanges.pdf", height=6, width=6)
# shrunken log2 fold change
par(mar=c(3,3,3,3)+0.1,mgp=c(1.7,0.6,0))
plot(resFBe11_FBp0_Df_fig$LFC_log2FoldChange, -log10(resFBe11_FBp0_Df_fig$padj), pch=20, cex=0.2, col=resFBe11_FBp0_Df_fig$col, cex.main=1.2, xlab="Fold Change: log2(p0.5/e11.5)", ylab="-log10(padj)", xlim=c(-10,10), ylim=c(0.1,15),
main="Gene expression in forebrain of e11.5 and p0 mouse embryos",bty="n", cex.main=1)
dev.off()

#### !!!!!!!!  The rest of the samples (for Liver and Heart) were done the same way but with changing the columns to those of each case !!!!!!!! ####













########################### HEART:


# for E11.5 vs E14.5 of Heart:
# creating "DESeqDataSet" object for estimating sizing factors and further differential gee expression analysis
ddsHTe11_HTe14 = DESeqDataSetFromMatrix(countData = countData_HTe11_HTe14, colData = coldata_HTe11_HTe14, design = ~ condition)
ddsHTe11_HTe14<-estimateSizeFactors(ddsHTe11_HTe14)
SF_HTe11_HTe14<-sizeFactors(ddsHTe11_HTe14)
ddsHTe11_HTe14 <- DESeq(ddsHTe11_HTe14)
res_HTe11_HTe14 <- results(ddsHTe11_HTe14)
resHTe11_HTe14_Df<-as.data.frame(res_HTe11_HTe14, stringsAsFactors=FALSE)
resHTe11_HTe14_Df<- cbind(annotatedNormCounts,resHTe11_HTe14_Df)   

res_HTe11_HTe14_05 <- results(ddsHTe11_HTe14, alpha=0.05) # to inlude only differentially expressed genes with padj<0.05
summary(res_HTe11_HTe14_05) # summary; outputs pasted below workflow

# generating shrunken log2 fold changes for better visualization:
resultsNames(ddsHTe11_HTe14)
> resultsNames(ddsHTe11_HTe14)
[1] "Intercept"                    "condition_HTe14.5_vs_HTe11.5"

resLFC_HTe11_HTe14 <- lfcShrink(ddsHTe11_HTe14, coef="condition_HTe14.5_vs_HTe11.5", type="apeglm")


head(resHTe11_HTe14_Df[, 30:37], 4)


resHTe11_HTe14_Df$HTe11mean<-rowMeans(resHTe11_HTe14_Df[,30:33])
resHTe11_HTe14_Df$HTe14mean<-rowMeans(resHTe11_HTe14_Df[,34:37])
resHTe11_HTe14_Df$HTe11_HTe14mean<-rowMeans(resHTe11_HTe14_Df[,30:37])
resHTe11_HTe14_Df$LFC_log2FoldChange<-resLFC_HTe11_HTe14$log2FoldChange
resHTe11_HTe14_Df$LFC_padj<-resLFC_HTe11_HTe14$padj
resHTe11_HTe14_Df_fig<-resHTe11_HTe14_Df[which(is.na(resHTe11_HTe14_Df$log2FoldChange)==FALSE),]
resHTe11_HTe14_Df_fig$padj[which(is.na(resHTe11_HTe14_Df_fig$padj)==TRUE)]<-1 

resHTe11_HTe14_Df_fig$col<- "gray20"

colFun<- function(i,tab){
if (tab$padj[i]<=0.05 ) a<- "red"
else if (tab$padj[i]>0.05 & tab$padj[i]<=0.1) a<- "green"
else if (tab$padj[i]>0.1 & tab$padj[i]<1 ) a<- "gray20"
else if (tab$padj[i]==1) a<- "white" # making data points with padj=NA (1) "invisible"
else a<-"gray20"
return(a)
}

resHTe11_HTe14_Df_fig$col<-unlist(lapply(1:length(resHTe11_HTe14_Df_fig[,1]), colFun, resHTe11_HTe14_Df_fig))

## MA plots

# regular log2 fold change
pdf("/home/yaaroba/ENCODEproject/figures/HTe11_HTe14_RNAseq_MA_plot.pdf", height=6, width=6)
par(mar=c(3,3,3,3)+0.1,mgp=c(1.7,0.6,0))
plot(resHTe11_HTe14_Df_fig$HTe11_HTe14mean, resHTe11_HTe14_Df_fig$log2FoldChange, pch=20, cex=0.2, col=resHTe11_HTe14_Df_fig$col, cex.main=1.2, xlab="mean normalized counts", ylab="log2FoldChange", xlim=c(0,10000), ylim=c(-5,5),
main="Gene expression in Heart of E11.5 and E14.5 mouse embryos",bty="n", cex.main=1)
abline(h = 0, col = "black")
dev.off()

# shrunken regular log2 fold change
pdf("/home/yaaroba/ENCODEproject/figures/HTe11_HTe14_RNAseq_MA_plot_shrunken_log2foldChanges.pdf", height=6, width=6)
par(mar=c(3,3,3,3)+0.1,mgp=c(1.7,0.6,0))
plot(resHTe11_HTe14_Df_fig$HTe11_HTe14mean, resHTe11_HTe14_Df_fig$LFC_log2FoldChange, pch=20, cex=0.2, col=resHTe11_HTe14_Df_fig$col, cex.main=1.2, xlab="mean normalized counts", ylab="log2FoldChange", xlim=c(0,10000), ylim=c(-5,5),
main="Gene expression in Heart of E11.5 and E14.5 mouse embryos",bty="n", cex.main=1)
abline(h = 0, col = "black")
dev.off()

## Vulcano plots

# regular log2 fold change
pdf("/home/yaaroba/ENCODEproject/figures/HTe11_HTe14_RNAseq_volcano_plot.pdf", height=6, width=6)
# regular log2 fold change
par(mar=c(3,3,3,3)+0.1,mgp=c(1.7,0.6,0))
plot(resHTe11_HTe14_Df_fig$log2FoldChange, -log10(resHTe11_HTe14_Df_fig$padj), pch=20, cex=0.2, col=resHTe11_HTe14_Df_fig$col, cex.main=1.2, xlab="Fold Change: log2(E14.5/E11.5)", ylab="-log10(padj)", xlim=c(-10,10), ylim=c(0.5,15), asp=2,
main="Gene expression in Heart of E11.5 and E14.5  mouse embryos",bty="n", cex.main=1)
dev.off()

pdf("/home/yaaroba/ENCODEproject/figures/HTe11_HTe14_RNAseq_volcano_plot_shrunken_log2foldChanges.pdf", height=6, width=6)
# shrunken log2 fold change
par(mar=c(3,3,3,3)+0.1,mgp=c(1.7,0.6,0))
plot(resHTe11_HTe14_Df_fig$LFC_log2FoldChange, -log10(resHTe11_HTe14_Df_fig$padj), pch=20, cex=0.2, col=resHTe11_HTe14_Df_fig$col, cex.main=1.2, xlab="Fold Change: log2(E14.5/E11.5)", ylab="-log10(padj)", xlim=c(-10,10), ylim=c(0.1,15),
main="Gene expression in Heart of E11.5 and E14.5 mouse embryos",bty="n", cex.main=1)
dev.off()






# for E14.5 vs P0 of Heart:
# creating "DESeqDataSet" object for estimating sizing factors and further differential gee expression analysis
ddsHTe14_HTp0 = DESeqDataSetFromMatrix(countData = countData_HTe14_HTp0, colData = coldata_HTe14_HTp0, design = ~ condition)
ddsHTe14_HTp0<-estimateSizeFactors(ddsHTe14_HTp0)
SF_HTe14_HTp0<-sizeFactors(ddsHTe14_HTp0)
ddsHTe14_HTp0 <- DESeq(ddsHTe14_HTp0)
res_HTe14_HTp0 <- results(ddsHTe14_HTp0)
resHTe14_HTp0_Df<-as.data.frame(res_HTe14_HTp0, stringsAsFactors=FALSE)
resHTe14_HTp0_Df<- cbind(annotatedNormCounts,resHTe14_HTp0_Df)   

res_HTe14_HTp0_05 <- results(ddsHTe14_HTp0, alpha=0.05) # to inlude only differentially expressed genes with padj<0.05
summary(res_HTe14_HTp0_05) # summary; outputs pasted below workflow

# generating shrunken log2 fold changes for better visualization:
resultsNames(ddsHTe14_HTp0)
> resultsNames(ddsHTe14_HTp0)
[1] "Intercept"                    "condition_HTp0_vs_HTe14.5"

resLFC_HTe14_HTp0 <- lfcShrink(ddsHTe14_HTp0, coef="condition_HTp0_vs_HTe14.5", type="apeglm")


head(resHTe14_HTp0_Df[, c(20:23, 34:37)], 2)


resHTe14_HTp0_Df$HTe14mean<-rowMeans(resHTe14_HTp0_Df[,34:37])
resHTe14_HTp0_Df$HTp0mean<-rowMeans(resHTe14_HTp0_Df[,20:23])
resHTe14_HTp0_Df$HTe14_HTp0mean<-rowMeans(resHTe14_HTp0_Df[,c(20:23, 34:37)])
resHTe14_HTp0_Df$LFC_log2FoldChange<-resLFC_HTe14_HTp0$log2FoldChange
resHTe14_HTp0_Df$LFC_padj<-resLFC_HTe14_HTp0$padj
resHTe14_HTp0_Df_fig<-resHTe14_HTp0_Df[which(is.na(resHTe14_HTp0_Df$log2FoldChange)==FALSE),]
resHTe14_HTp0_Df_fig$padj[which(is.na(resHTe14_HTp0_Df_fig$padj)==TRUE)]<-1 

resHTe14_HTp0_Df_fig$col<- "gray20"

colFun<- function(i,tab){
if (tab$padj[i]<=0.05 ) a<- "red"
else if (tab$padj[i]>0.05 & tab$padj[i]<=0.1) a<- "green"
else if (tab$padj[i]>0.1 & tab$padj[i]<1 ) a<- "gray20"
else if (tab$padj[i]==1) a<- "white" # making data points with padj=NA (1) "invisible"
else a<-"gray20"
return(a)
}

resHTe14_HTp0_Df_fig$col<-unlist(lapply(1:length(resHTe14_HTp0_Df_fig[,1]), colFun, resHTe14_HTp0_Df_fig))

## MA plots

# regular log2 fold change
pdf("/home/yaaroba/ENCODEproject/figures/HTe14_HTp0_RNAseq_MA_plot.pdf", height=6, width=6)
par(mar=c(3,3,3,3)+0.1,mgp=c(1.7,0.6,0))
plot(resHTe14_HTp0_Df_fig$HTe14_HTp0mean, resHTe14_HTp0_Df_fig$log2FoldChange, pch=20, cex=0.2, col=resHTe14_HTp0_Df_fig$col, cex.main=1.2, xlab="mean normalized counts", ylab="log2FoldChange", xlim=c(0,10000), ylim=c(-5,5),
main="Gene expression in Heart of e14.5 and p0 mouse embryos",bty="n", cex.main=1)
abline(h = 0, col = "black")
dev.off()

# shrunken regular log2 fold change
pdf("/home/yaaroba/ENCODEproject/figures/HTe14_HTp0_RNAseq_MA_plot_shrunken_log2foldChanges.pdf", height=6, width=6)
par(mar=c(3,3,3,3)+0.1,mgp=c(1.7,0.6,0))
plot(resHTe14_HTp0_Df_fig$HTe14_HTp0mean, resHTe14_HTp0_Df_fig$LFC_log2FoldChange, pch=20, cex=0.2, col=resHTe14_HTp0_Df_fig$col, cex.main=1.2, xlab="mean normalized counts", ylab="log2FoldChange", xlim=c(0,10000), ylim=c(-5,5),
main="Gene expression in Heart of e14.5 and p0 mouse embryos",bty="n", cex.main=1)
abline(h = 0, col = "black")
dev.off()

## Vulcano plots

# regular log2 fold change
pdf("/home/yaaroba/ENCODEproject/figures/HTe14_HTp0_RNAseq_volcano_plot.pdf", height=6, width=6)
# regular log2 fold change
par(mar=c(3,3,3,3)+0.1,mgp=c(1.7,0.6,0))
plot(resHTe14_HTp0_Df_fig$log2FoldChange, -log10(resHTe14_HTp0_Df_fig$padj), pch=20, cex=0.2, col=resHTe14_HTp0_Df_fig$col, cex.main=1.2, xlab="Fold Change: log2(p0.5/e14.5)", ylab="-log10(padj)", xlim=c(-10,10), ylim=c(0.5,15), asp=2,
main="Gene expression in Heart of e14.5 and p0 mouse embryos",bty="n", cex.main=1)
dev.off()

pdf("/home/yaaroba/ENCODEproject/figures/HTe14_HTp0_RNAseq_volcano_plot_shrunken_log2foldChanges.pdf", height=6, width=6)
# shrunken log2 fold change
par(mar=c(3,3,3,3)+0.1,mgp=c(1.7,0.6,0))
plot(resHTe14_HTp0_Df_fig$LFC_log2FoldChange, -log10(resHTe14_HTp0_Df_fig$padj), pch=20, cex=0.2, col=resHTe14_HTp0_Df_fig$col, cex.main=1.2, xlab="Fold Change: log2(p0.5/e14.5)", ylab="-log10(padj)", xlim=c(-10,10), ylim=c(0.1,15),
main="Gene expression in Heart of e14.5 and p0 mouse embryos",bty="n", cex.main=1)
dev.off()







# for e11.5 vs P0 of Heart:
# creating "DESeqDataSet" object for estimating sizing factors and further differential gee expression analysis
ddsHTe11_HTp0 = DESeqDataSetFromMatrix(countData = countData_HTe11_HTp0, colData = coldata_HTe11_HTp0, design = ~ condition)
ddsHTe11_HTp0<-estimateSizeFactors(ddsHTe11_HTp0)
SF_HTe11_HTp0<-sizeFactors(ddsHTe11_HTp0)
ddsHTe11_HTp0 <- DESeq(ddsHTe11_HTp0)
res_HTe11_HTp0 <- results(ddsHTe11_HTp0)
resHTe11_HTp0_Df<-as.data.frame(res_HTe11_HTp0, stringsAsFactors=FALSE)
resHTe11_HTp0_Df<- cbind(annotatedNormCounts,resHTe11_HTp0_Df)   

res_HTe11_HTp0_05 <- results(ddsHTe11_HTp0, alpha=0.05) # to inlude only differentially expressed genes with padj<0.05
summary(res_HTe11_HTp0_05) # summary; outputs pasted below workflow

# generating shrunken log2 fold changes for better visualization:
resultsNames(ddsHTe11_HTp0)
> resultsNames(ddsHTe11_HTp0)
[1] "Intercept"                    "condition_HTp0_vs_HTe11.5"

resLFC_HTe11_HTp0 <- lfcShrink(ddsHTe11_HTp0, coef="condition_HTp0_vs_HTe11.5", type="apeglm")


head(resHTe11_HTp0_Df[, c(20:23, 30:33)], 2)


resHTe11_HTp0_Df$HTe11mean<-rowMeans(resHTe11_HTp0_Df[,30:33])
resHTe11_HTp0_Df$HTp0mean<-rowMeans(resHTe11_HTp0_Df[,20:23])
resHTe11_HTp0_Df$HTe11_HTp0mean<-rowMeans(resHTe11_HTp0_Df[,c(20:23, 30:33)])
resHTe11_HTp0_Df$LFC_log2FoldChange<-resLFC_HTe11_HTp0$log2FoldChange
resHTe11_HTp0_Df$LFC_padj<-resLFC_HTe11_HTp0$padj
resHTe11_HTp0_Df_fig<-resHTe11_HTp0_Df[which(is.na(resHTe11_HTp0_Df$log2FoldChange)==FALSE),]
resHTe11_HTp0_Df_fig$padj[which(is.na(resHTe11_HTp0_Df_fig$padj)==TRUE)]<-1 

resHTe11_HTp0_Df_fig$col<- "gray20"

colFun<- function(i,tab){
if (tab$padj[i]<=0.05 ) a<- "red"
else if (tab$padj[i]>0.05 & tab$padj[i]<=0.1) a<- "green"
else if (tab$padj[i]>0.1 & tab$padj[i]<1 ) a<- "gray20"
else if (tab$padj[i]==1) a<- "white" # making data points with padj=NA (1) "invisible"
else a<-"gray20"
return(a)
}

resHTe11_HTp0_Df_fig$col<-unlist(lapply(1:length(resHTe11_HTp0_Df_fig[,1]), colFun, resHTe11_HTp0_Df_fig))

## MA plots

# regular log2 fold change
pdf("/home/yaaroba/ENCODEproject/figures/HTe11_HTp0_RNAseq_MA_plot.pdf", height=6, width=6)
par(mar=c(3,3,3,3)+0.1,mgp=c(1.7,0.6,0))
plot(resHTe11_HTp0_Df_fig$HTe11_HTp0mean, resHTe11_HTp0_Df_fig$log2FoldChange, pch=20, cex=0.2, col=resHTe11_HTp0_Df_fig$col, cex.main=1.2, xlab="mean normalized counts", ylab="log2FoldChange", xlim=c(0,10000), ylim=c(-5,5),f
main="Gene expression in Heart of e11.5 and p0 mouse embryos",bty="n", cex.main=1)
abline(h = 0, col = "black")
dev.off()

# shrunken regular log2 fold change
pdf("/home/yaaroba/ENCODEproject/figures/HTe11_HTp0_RNAseq_MA_plot_shrunken_log2foldChanges.pdf", height=6, width=6)
par(mar=c(3,3,3,3)+0.1,mgp=c(1.7,0.6,0))
plot(resHTe11_HTp0_Df_fig$HTe11_HTp0mean, resHTe11_HTp0_Df_fig$LFC_log2FoldChange, pch=20, cex=0.2, col=resHTe11_HTp0_Df_fig$col, cex.main=1.2, xlab="mean normalized counts", ylab="log2FoldChange", xlim=c(0,10000), ylim=c(-5,5),
main="Gene expression in Heart of e11.5 and p0 mouse embryos",bty="n", cex.main=1)
abline(h = 0, col = "black")
dev.off()

## Vulcano plots

# regular log2 fold change
pdf("/home/yaaroba/ENCODEproject/figures/HTe11_HTp0_RNAseq_volcano_plot.pdf", height=6, width=6)
# regular log2 fold change
par(mar=c(3,3,3,3)+0.1,mgp=c(1.7,0.6,0))
plot(resHTe11_HTp0_Df_fig$log2FoldChange, -log10(resHTe11_HTp0_Df_fig$padj), pch=20, cex=0.2, col=resHTe11_HTp0_Df_fig$col, cex.main=1.2, xlab="Fold Change: log2(p0.5/e11.5)", ylab="-log10(padj)", xlim=c(-10,10), ylim=c(0.5,15), asp=2,
main="Gene expression in Heart of e11.5 and p0 mouse embryos",bty="n", cex.main=1)
dev.off()

pdf("/home/yaaroba/ENCODEproject/figures/HTe11_HTp0_RNAseq_volcano_plot_shrunken_log2foldChanges.pdf", height=6, width=6)
# shrunken log2 fold change
par(mar=c(3,3,3,3)+0.1,mgp=c(1.7,0.6,0))
plot(resHTe11_HTp0_Df_fig$LFC_log2FoldChange, -log10(resHTe11_HTp0_Df_fig$padj), pch=20, cex=0.2, col=resHTe11_HTp0_Df_fig$col, cex.main=1.2, xlab="Fold Change: log2(p0.5/e11.5)", ylab="-log10(padj)", xlim=c(-10,10), ylim=c(0.1,15),
main="Gene expression in Heart of e11.5 and p0 mouse embryos",bty="n", cex.main=1)
dev.off()














################### LIVER:


# for E11.5 vs E14.5 of Liver:
# creating "DESeqDataSet" object for estimating sizing factors and further differential gee expression analysis
ddsLVe11_LVe14 = DESeqDataSetFromMatrix(countData = countData_LVe11_LVe14, colData = coldata_LVe11_LVe14, design = ~ condition)
ddsLVe11_LVe14<-estimateSizeFactors(ddsLVe11_LVe14)
SF_LVe11_LVe14<-sizeFactors(ddsLVe11_LVe14)
ddsLVe11_LVe14 <- DESeq(ddsLVe11_LVe14)
res_LVe11_LVe14 <- results(ddsLVe11_LVe14)
resLVe11_LVe14_Df<-as.data.frame(res_LVe11_LVe14, stringsAsFactors=FALSE)
resLVe11_LVe14_Df<- cbind(annotatedNormCounts,resLVe11_LVe14_Df)   

res_LVe11_LVe14_05 <- results(ddsLVe11_LVe14, alpha=0.05) # to inlude only differentially expressed genes with padj<0.05
summary(res_LVe11_LVe14_05) # summary; outputs pasted below workflow

# generating shrunken log2 fold changes for better visualization:
resultsNames(ddsLVe11_LVe14)
> resultsNames(ddsLVe11_LVe14)
[1] "Intercept"                    "condition_LVe14.5_vs_LVe11.5"

resLFC_LVe11_LVe14 <- lfcShrink(ddsLVe11_LVe14, coef="condition_LVe14.5_vs_LVe11.5", type="apeglm")


head(resLVe11_LVe14_Df[, c(26:29, 38:41)], 4)


resLVe11_LVe14_Df$LVe11mean<-rowMeans(resLVe11_LVe14_Df[,26:29])
resLVe11_LVe14_Df$LVe14mean<-rowMeans(resLVe11_LVe14_Df[,38:41])
resLVe11_LVe14_Df$LVe11_LVe14mean<-rowMeans(resLVe11_LVe14_Df[,c(26:29, 38:41)])
resLVe11_LVe14_Df$LFC_log2FoldChange<-resLFC_LVe11_LVe14$log2FoldChange
resLVe11_LVe14_Df$LFC_padj<-resLFC_LVe11_LVe14$padj
resLVe11_LVe14_Df_fig<-resLVe11_LVe14_Df[which(is.na(resLVe11_LVe14_Df$log2FoldChange)==FALSE),]
resLVe11_LVe14_Df_fig$padj[which(is.na(resLVe11_LVe14_Df_fig$padj)==TRUE)]<-1 

resLVe11_LVe14_Df_fig$col<- "gray20"

colFun<- function(i,tab){
if (tab$padj[i]<=0.05 ) a<- "red"
else if (tab$padj[i]>0.05 & tab$padj[i]<=0.1) a<- "green"
else if (tab$padj[i]>0.1 & tab$padj[i]<1 ) a<- "gray20"
else if (tab$padj[i]==1) a<- "white" # making data points with padj=NA (1) "invisible"
else a<-"gray20"
return(a)
}

resLVe11_LVe14_Df_fig$col<-unlist(lapply(1:length(resLVe11_LVe14_Df_fig[,1]), colFun, resLVe11_LVe14_Df_fig))

## MA plots

# regular log2 fold change
pdf("/home/yaaroba/ENCODEproject/figures/LVe11_LVe14_RNAseq_MA_plot.pdf", height=6, width=6)
par(mar=c(3,3,3,3)+0.1,mgp=c(1.7,0.6,0))
plot(resLVe11_LVe14_Df_fig$LVe11_LVe14mean, resLVe11_LVe14_Df_fig$log2FoldChange, pch=20, cex=0.2, col=resLVe11_LVe14_Df_fig$col, cex.main=1.2, xlab="mean normalized counts", ylab="log2FoldChange", xlim=c(0,10000), ylim=c(-5,5),
main="Gene expression in Liver of E11.5 and E14.5 mouse embryos",bty="n", cex.main=1)
abline(h = 0, col = "black")
dev.off()

# shrunken regular log2 fold change
pdf("/home/yaaroba/ENCODEproject/figures/LVe11_LVe14_RNAseq_MA_plot_shrunken_log2foldChanges.pdf", height=6, width=6)
par(mar=c(3,3,3,3)+0.1,mgp=c(1.7,0.6,0))
plot(resLVe11_LVe14_Df_fig$LVe11_LVe14mean, resLVe11_LVe14_Df_fig$LFC_log2FoldChange, pch=20, cex=0.2, col=resLVe11_LVe14_Df_fig$col, cex.main=1.2, xlab="mean normalized counts", ylab="log2FoldChange", xlim=c(0,10000), ylim=c(-5,5),
main="Gene expression in Liver of E11.5 and E14.5 mouse embryos",bty="n", cex.main=1)
abline(h = 0, col = "black")
dev.off()

## Vulcano plots

# regular log2 fold change
pdf("/home/yaaroba/ENCODEproject/figures/LVe11_LVe14_RNAseq_volcano_plot.pdf", height=6, width=6)
# regular log2 fold change
par(mar=c(3,3,3,3)+0.1,mgp=c(1.7,0.6,0))
plot(resLVe11_LVe14_Df_fig$log2FoldChange, -log10(resLVe11_LVe14_Df_fig$padj), pch=20, cex=0.2, col=resLVe11_LVe14_Df_fig$col, cex.main=1.2, xlab="Fold Change: log2(E14.5/E11.5)", ylab="-log10(padj)", xlim=c(-10,10), ylim=c(0.5,15), asp=2,
main="Gene expression in Liver of E11.5 and E14.5  mouse embryos",bty="n", cex.main=1)
dev.off()

pdf("/home/yaaroba/ENCODEproject/figures/LVe11_LVe14_RNAseq_volcano_plot_shrunken_log2foldChanges.pdf", height=6, width=6)
# shrunken log2 fold change
par(mar=c(3,3,3,3)+0.1,mgp=c(1.7,0.6,0))
plot(resLVe11_LVe14_Df_fig$LFC_log2FoldChange, -log10(resLVe11_LVe14_Df_fig$padj), pch=20, cex=0.2, col=resLVe11_LVe14_Df_fig$col, cex.main=1.2, xlab="Fold Change: log2(E14.5/E11.5)", ylab="-log10(padj)", xlim=c(-10,10), ylim=c(0.1,15),
main="Gene expression in Liver of E11.5 and E14.5 mouse embryos",bty="n", cex.main=1)
dev.off()






# for E14.5 vs P0 of Liver:
# creating "DESeqDataSet" object for estimating sizing factors and further differential gee expression analysis
ddsLVe14_LVp0 = DESeqDataSetFromMatrix(countData = countData_LVe14_LVp0, colData = coldata_LVe14_LVp0, design = ~ condition)
ddsLVe14_LVp0<-estimateSizeFactors(ddsLVe14_LVp0)
SF_LVe14_LVp0<-sizeFactors(ddsLVe14_LVp0)
ddsLVe14_LVp0 <- DESeq(ddsLVe14_LVp0)
res_LVe14_LVp0 <- results(ddsLVe14_LVp0)
resLVe14_LVp0_Df<-as.data.frame(res_LVe14_LVp0, stringsAsFactors=FALSE)
resLVe14_LVp0_Df<- cbind(annotatedNormCounts,resLVe14_LVp0_Df)   

res_LVe14_LVp0_05 <- results(ddsLVe14_LVp0, alpha=0.05) # to inlude only differentially expressed genes with padj<0.05
summary(res_LVe14_LVp0_05) # summary; outputs pasted below workflow

# generating shrunken log2 fold changes for better visualization:
resultsNames(ddsLVe14_LVp0)
> resultsNames(ddsLVe14_LVp0)
[1] "Intercept"                    "condition_LVp0_vs_LVe14.5"

resLFC_LVe14_LVp0 <- lfcShrink(ddsLVe14_LVp0, coef="condition_LVp0_vs_LVe14.5", type="apeglm")


head(resLVe14_LVp0_Df[, c(38:45)], 2)


resLVe14_LVp0_Df$LVe14mean<-rowMeans(resLVe14_LVp0_Df[,38:41])
resLVe14_LVp0_Df$LVp0mean<-rowMeans(resLVe14_LVp0_Df[,42:45])
resLVe14_LVp0_Df$LVe14_LVp0mean<-rowMeans(resLVe14_LVp0_Df[,38:45])
resLVe14_LVp0_Df$LFC_log2FoldChange<-resLFC_LVe14_LVp0$log2FoldChange
resLVe14_LVp0_Df$LFC_padj<-resLFC_LVe14_LVp0$padj
resLVe14_LVp0_Df_fig<-resLVe14_LVp0_Df[which(is.na(resLVe14_LVp0_Df$log2FoldChange)==FALSE),]
resLVe14_LVp0_Df_fig$padj[which(is.na(resLVe14_LVp0_Df_fig$padj)==TRUE)]<-1 

resLVe14_LVp0_Df_fig$col<- "gray20"

colFun<- function(i,tab){
if (tab$padj[i]<=0.05 ) a<- "red"
else if (tab$padj[i]>0.05 & tab$padj[i]<=0.1) a<- "green"
else if (tab$padj[i]>0.1 & tab$padj[i]<1 ) a<- "gray20"
else if (tab$padj[i]==1) a<- "white" # making data points with padj=NA (1) "invisible"
else a<-"gray20"
return(a)
}

resLVe14_LVp0_Df_fig$col<-unlist(lapply(1:length(resLVe14_LVp0_Df_fig[,1]), colFun, resLVe14_LVp0_Df_fig))

## MA plots

# regular log2 fold change
pdf("/home/yaaroba/ENCODEproject/figures/LVe14_LVp0_RNAseq_MA_plot.pdf", height=6, width=6)
par(mar=c(3,3,3,3)+0.1,mgp=c(1.7,0.6,0))
plot(resLVe14_LVp0_Df_fig$LVe14_LVp0mean, resLVe14_LVp0_Df_fig$log2FoldChange, pch=20, cex=0.2, col=resLVe14_LVp0_Df_fig$col, cex.main=1.2, xlab="mean normalized counts", ylab="log2FoldChange", xlim=c(0,10000), ylim=c(-5,5),
main="Gene expression in Liver of e14.5 and p0 mouse embryos",bty="n", cex.main=1)
abline(h = 0, col = "black")
dev.off()

# shrunken regular log2 fold change
pdf("/home/yaaroba/ENCODEproject/figures/LVe14_LVp0_RNAseq_MA_plot_shrunken_log2foldChanges.pdf", height=6, width=6)
par(mar=c(3,3,3,3)+0.1,mgp=c(1.7,0.6,0))
plot(resLVe14_LVp0_Df_fig$LVe14_LVp0mean, resLVe14_LVp0_Df_fig$LFC_log2FoldChange, pch=20, cex=0.2, col=resLVe14_LVp0_Df_fig$col, cex.main=1.2, xlab="mean normalized counts", ylab="log2FoldChange", xlim=c(0,10000), ylim=c(-5,5),
main="Gene expression in Liver of e14.5 and p0 mouse embryos",bty="n", cex.main=1)
abline(h = 0, col = "black")
dev.off()

## Vulcano plots

# regular log2 fold change
pdf("/home/yaaroba/ENCODEproject/figures/LVe14_LVp0_RNAseq_volcano_plot.pdf", height=6, width=6)
# regular log2 fold change
par(mar=c(3,3,3,3)+0.1,mgp=c(1.7,0.6,0))
plot(resLVe14_LVp0_Df_fig$log2FoldChange, -log10(resLVe14_LVp0_Df_fig$padj), pch=20, cex=0.2, col=resLVe14_LVp0_Df_fig$col, cex.main=1.2, xlab="Fold Change: log2(p0.5/e14.5)", ylab="-log10(padj)", xlim=c(-10,10), ylim=c(0.5,15), asp=2,
main="Gene expression in Liver of e14.5 and p0 mouse embryos",bty="n", cex.main=1)
dev.off()

pdf("/home/yaaroba/ENCODEproject/figures/LVe14_LVp0_RNAseq_volcano_plot_shrunken_log2foldChanges.pdf", height=6, width=6)
# shrunken log2 fold change
par(mar=c(3,3,3,3)+0.1,mgp=c(1.7,0.6,0))
plot(resLVe14_LVp0_Df_fig$LFC_log2FoldChange, -log10(resLVe14_LVp0_Df_fig$padj), pch=20, cex=0.2, col=resLVe14_LVp0_Df_fig$col, cex.main=1.2, xlab="Fold Change: log2(p0.5/e14.5)", ylab="-log10(padj)", xlim=c(-10,10), ylim=c(0.1,15),
main="Gene expression in Liver of e14.5 and p0 mouse embryos",bty="n", cex.main=1)
dev.off()







# for e11.5 vs P0 of Liver:
# creating "DESeqDataSet" object for estimating sizing factors and further differential gee expression analysis
ddsLVe11_LVp0 = DESeqDataSetFromMatrix(countData = countData_LVe11_LVp0, colData = coldata_LVe11_LVp0, design = ~ condition)
ddsLVe11_LVp0<-estimateSizeFactors(ddsLVe11_LVp0)
SF_LVe11_LVp0<-sizeFactors(ddsLVe11_LVp0)
ddsLVe11_LVp0 <- DESeq(ddsLVe11_LVp0)
res_LVe11_LVp0 <- results(ddsLVe11_LVp0)
resLVe11_LVp0_Df<-as.data.frame(res_LVe11_LVp0, stringsAsFactors=FALSE)
resLVe11_LVp0_Df<- cbind(annotatedNormCounts,resLVe11_LVp0_Df)   

res_LVe11_LVp0_05 <- results(ddsLVe11_LVp0, alpha=0.05) # to inlude only differentially expressed genes with padj<0.05
summary(res_LVe11_LVp0_05) # summary; outputs pasted below workflow

# generating shrunken log2 fold changes for better visualization:
resultsNames(ddsLVe11_LVp0)
> resultsNames(ddsLVe11_LVp0)
[1] "Intercept"                    "condition_LVp0_vs_LVe11.5"

resLFC_LVe11_LVp0 <- lfcShrink(ddsLVe11_LVp0, coef="condition_LVp0_vs_LVe11.5", type="apeglm")


head(resLVe11_LVp0_Df[, c(26:29, 42:45)], 2)



resLVe11_LVp0_Df$LVe11mean<-rowMeans(resLVe11_LVp0_Df[,26:29])
resLVe11_LVp0_Df$LVp0mean<-rowMeans(resLVe11_LVp0_Df[,42:45])
resLVe11_LVp0_Df$LVe11_LVp0mean<-rowMeans(resLVe11_LVp0_Df[,c(26:29, 42:45)])
resLVe11_LVp0_Df$LFC_log2FoldChange<-resLFC_LVe11_LVp0$log2FoldChange
resLVe11_LVp0_Df$LFC_padj<-resLFC_LVe11_LVp0$padj
resLVe11_LVp0_Df_fig<-resLVe11_LVp0_Df[which(is.na(resLVe11_LVp0_Df$log2FoldChange)==FALSE),]
resLVe11_LVp0_Df_fig$padj[which(is.na(resLVe11_LVp0_Df_fig$padj)==TRUE)]<-1 

resLVe11_LVp0_Df_fig$col<- "gray20"

colFun<- function(i,tab){
if (tab$padj[i]<=0.05 ) a<- "red"
else if (tab$padj[i]>0.05 & tab$padj[i]<=0.1) a<- "green"
else if (tab$padj[i]>0.1 & tab$padj[i]<1 ) a<- "gray20"
else if (tab$padj[i]==1) a<- "white" # making data points with padj=NA (1) "invisible"
else a<-"gray20"
return(a)
}

resLVe11_LVp0_Df_fig$col<-unlist(lapply(1:length(resLVe11_LVp0_Df_fig[,1]), colFun, resLVe11_LVp0_Df_fig))

## MA plots

# regular log2 fold change
pdf("/home/yaaroba/ENCODEproject/figures/LVe11_LVp0_RNAseq_MA_plot.pdf", height=6, width=6)
par(mar=c(3,3,3,3)+0.1,mgp=c(1.7,0.6,0))
plot(resLVe11_LVp0_Df_fig$LVe11_LVp0mean, resLVe11_LVp0_Df_fig$log2FoldChange, pch=20, cex=0.2, col=resLVe11_LVp0_Df_fig$col, cex.main=1.2, xlab="mean normalized counts", ylab="log2FoldChange", xlim=c(0,10000), ylim=c(-5,5),
main="Gene expression in Liver of e11.5 and p0 mouse embryos",bty="n", cex.main=1)
abline(h = 0, col = "black")
dev.off()

# shrunken regular log2 fold change
pdf("/home/yaaroba/ENCODEproject/figures/LVe11_LVp0_RNAseq_MA_plot_shrunken_log2foldChanges.pdf", height=6, width=6)
par(mar=c(3,3,3,3)+0.1,mgp=c(1.7,0.6,0))
plot(resLVe11_LVp0_Df_fig$LVe11_LVp0mean, resLVe11_LVp0_Df_fig$LFC_log2FoldChange, pch=20, cex=0.2, col=resLVe11_LVp0_Df_fig$col, cex.main=1.2, xlab="mean normalized counts", ylab="log2FoldChange", xlim=c(0,10000), ylim=c(-5,5),
main="Gene expression in Liver of e11.5 and p0 mouse embryos",bty="n", cex.main=1)
abline(h = 0, col = "black")
dev.off()

## Vulcano plots

# regular log2 fold change
pdf("/home/yaaroba/ENCODEproject/figures/LVe11_LVp0_RNAseq_volcano_plot.pdf", height=6, width=6)
# regular log2 fold change
par(mar=c(3,3,3,3)+0.1,mgp=c(1.7,0.6,0))
plot(resLVe11_LVp0_Df_fig$log2FoldChange, -log10(resLVe11_LVp0_Df_fig$padj), pch=20, cex=0.2, col=resLVe11_LVp0_Df_fig$col, cex.main=1.2, xlab="Fold Change: log2(p0.5/e11.5)", ylab="-log10(padj)", xlim=c(-10,10), ylim=c(0.5,15), asp=2,
main="Gene expression in Liver of e11.5 and p0 mouse embryos",bty="n", cex.main=1)
dev.off()

pdf("/home/yaaroba/ENCODEproject/figures/LVe11_LVp0_RNAseq_volcano_plot_shrunken_log2foldChanges.pdf", height=6, width=6)
# shrunken log2 fold change
par(mar=c(3,3,3,3)+0.1,mgp=c(1.7,0.6,0))
plot(resLVe11_LVp0_Df_fig$LFC_log2FoldChange, -log10(resLVe11_LVp0_Df_fig$padj), pch=20, cex=0.2, col=resLVe11_LVp0_Df_fig$col, cex.main=1.2, xlab="Fold Change: log2(p0.5/e11.5)", ylab="-log10(padj)", xlim=c(-10,10), ylim=c(0.1,15),
main="Gene expression in Liver of e11.5 and p0 mouse embryos",bty="n", cex.main=1)
dev.off()