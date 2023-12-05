##3. Differential Gene Expression Analysis

mkdir /home/yaaroba/ENCODEproject/figures

R
library(DESeq2)
library(data.table)
library(apeglm)

# loading "annotatedNormCountsFPKM" object genertaed within "GeneratingNormalizedBedgraphsDEseq2.sh" workflow
load("/home/yaaroba/ENCODEproject/Tables/annotatedNormCountsFPKM.RData") 
load("/home/yaaroba/ENCODEproject/Tables/annotatedNormCounts.RData")
#If you already did part of this code and you came back to proceed with it you can load the saved workspace and proceed further without doing everything again:
load("/home/yaaroba/ENCODEproject/STAR_geneCounts/.RData")

### building sample table for all samples analysed
where<-"/home/yaaroba/ENCODEproject/STAR_geneCounts/"
setwd(where)

# buliding sample table necessary for DEseq2 analysis
sampleFiles <- grep("ReadsPerGene.out.tab",list.files(where),value=TRUE)

### in analyzed data three different types of RNA-seq libraries: files 1-16 unstranded => 2nd column; files 17-25 first strand library => 3rd column; files 26-28 second strand library  => 4th column

countData = data.frame(fread(sampleFiles[1]))[c(1,2)]

# Loop and read the 4th column remaining files
for(i in 2:36) {
        countData = cbind(countData, data.frame(fread(sampleFiles[i]))[2])
}

# get sample names from the column of files printed out after typing "sampleFiles"; use sublime funcinalities: mark file column with mouse + shift, get the coursor at the end of all rows at once with ctrl +s hift + L and tailor the file names unless you remove all redundant parts of teh file name; add "", go to teh beginning of teh line with home and press backspace

sampleNames<-c("geneID","forebrain_E11.5_rep1","forebrain_E11.5_rep2","forebrain_E11.5_rep3","forebrain_E11.5_rep4","forebrain_E14.5_rep1","forebrain_E14.5_rep2","forebrain_E14.5_rep3","forebrain_E14.5_rep4","forebrain_P0_rep1","forebrain_P0_rep2","heart_P0_rep1","heart_P0_rep2","heart_P0_rep3","heart_P0_rep4","forebrain_P0_rep3","forebrain_P0_rep4","liver_E11.5_rep1","liver_E11.5_rep2","liver_E11.5_rep3","liver_E11.5_rep4","heart_E11.5_rep1","heart_E11.5_rep2","heart_E11.5_rep3","heart_E11.5_rep4","heart_E14.5_rep1","heart_E14.5_rep2","heart_E14.5_rep3","heart_E14.5_rep4","liver_E14.5_rep1","liver_E14.5_rep2","liver_E14.5_rep3","liver_E14.5_rep4","liver_P0_rep1","liver_P0_rep2","liver_P0_rep3","liver_P0_rep4")

colnames(countData) = sampleNames

# retrieving original head of STAR gene reads count data with sample stats for further summary
countDataHead<- countData[1:40,]   
rownames(countDataHead) = countDataHead[,1]
countDataHead<-countDataHead[,2:ncol(countDataHead)]

# Skip first 4 lines, count data starts on the 5th line
countData = countData[c(5:nrow(countData)),]
rownames(countData) = countData[,1]  ### moving geneIDs to rownames
countData = countData[,c(2:ncol(countData))] ### removing column with geneIDs

condition<-c("FBe11.5","FBe11.5","FBe11.5","FBe11.5","FBe14.5","FBe14.5","FBe14.5","FBe14.5","FBp0", "FBp0", "HTp0", "HTp0","HTp0", "HTp0", "FBp0", "FBp0","LVe11.5", "LVe11.5", "LVe11.5", "LVe11.5", "HTe11.5", "HTe11.5", "HTe11.5", "HTe11.5", "HTe14.5", "HTe14.5", "HTe14.5", "HTe14.5","LVe14.5","LVe14.5","LVe14.5","LVe14.5","LVp0","LVp0","LVp0","LVp0")
replicate<-c(1,2,3,4,1,2,3,4,1,2,1,2,3,4,3,4,1,2,3,4,1,2,3,4,1,2,3,4,1,2,3,4,1,2,3,4)

# generating table with sample information
coldata = data.frame(cbind(condition, replicate))
coldata$type<-"single-read"
rownames(coldata) = sampleNames[2:length(sampleNames)]   # first name refers to geneIDs
coldata$condition <- factor(coldata$condition)
coldata$replicate <- factor(coldata$replicate)
coldata$replicate <- factor(coldata$replicate)



# building table specific for FBe11_FBe14 analysis (Forebrain E11.5 vs Forebrain E14.5)
#to check whether the selected set is true for the conditions:
head(countData[, c(1:8)], 4)

countData_FBe11_FBe14<-countData[,1:8]
coldata_FBe11_FBe14<-coldata[1:8,]

# creating "DESeqDataSet" object for estimating sizing factors and further differential gee expression analysis
ddsFBe11_FBe14 = DESeqDataSetFromMatrix(countData = countData_FBe11_FBe14, colData = coldata_FBe11_FBe14, design = ~ condition)

ddsFBe11_FBe14<-estimateSizeFactors(ddsFBe11_FBe14)
SF_FBe11_FBe14<-sizeFactors(ddsFBe11_FBe14)

ddsFBe11_FBe14 <- DESeq(ddsFBe11_FBe14)

res_FBe11_FBe14 <- results(ddsFBe11_FBe14)

resFBe11_FBe14_Df<-as.data.frame(res_FBe11_FBe14, stringsAsFactors=FALSE)
# annotatedNormCounts cerated with GeneratingNormalizedBedgraphsDEseq2.sh workflow
# may be dowloaded with:
# annotatedNormCounts<- read.table("/home/micgdu/Analysis/Epidermis/geneExpression_Fan_Nayat/analysis/RNAseq/RNAseqNormalizedReadsPerGene_FPKM_Nayak2023_Fan2018_5June2023.txt", sep="\t", header=F, stringsAsFactors = FALSE)

resFBe11_FBe14_Df<- cbind(annotatedNormCounts,resFBe11_FBe14_Df)   

res_FBe11_FBe14_05 <- results(ddsFBe11_FBe14, alpha=0.05) # to inlude only differentially expressed genes with padj<0.05
summary(res_FBe11_FBe14_05) # summary; outputs pasted below workflow

write.table(resFBe11_FBe14_Df, file="/home/yaaroba/ENCODEproject/Tables/RNAseq_differentialExpression_FBe11_FBe14_HeP2020_26Aug23.txt",quote=FALSE,sep="\t",row.names=TRUE, col.names=TRUE)

# generating shrunken log2 fold changes for better visualization:

resultsNames(ddsFBe11_FBe14)
> resultsNames(ddsFBe11_FBe14)
[1] "Intercept"                    "condition_FBe14.5_vs_FBe11.5"

resLFC_FBe11_FBe14 <- lfcShrink(ddsFBe11_FBe14, coef="condition_FBe14.5_vs_FBe11.5", type="apeglm")





# building table specific for FBe14_FBp0 analysis (Forebrain E14.5 vs Forebrain P0)
#to check whether the selected set is true for the conditions:
head(countData[, c(5:10, 15:16)], 4)

countData_FBe14_FBp0<-countData[, c(5:10, 15:16)]
coldata_FBe14_FBp0<-coldata[c(5:10, 15:16),]

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

write.table(resFBe14_FBp0_Df, file="/home/yaaroba/ENCODEproject/Tables/RNAseq_differentialExpression_FBe14_FBp0_HeP2020_26Aug2023.txt",quote=FALSE,sep="\t",row.names=TRUE, col.names=TRUE)




# building table specific for FBe11_FBp0 analysis (Forebrain E11.5 vs Forebrain P0)
#to check whether the selected set is true for the conditions:
head(countData[, c(1:4, 9:10, 15:16)], 4)

countData_FBe11_FBp0<-countData[, c(1:4, 9:10, 15:16)]
coldata_FBe11_FBp0<-coldata[c(1:4, 9:10, 15:16),]

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

write.table(resFBe11_FBp0_Df, file="/home/yaaroba/ENCODEproject/Tables/RNAseq_differentialExpression_FBe11_FBp0_HeP2020_26Aug2023.txt",quote=FALSE,sep="\t",row.names=TRUE, col.names=TRUE)





# building table specific for HTe11_HTe14 analysis (Heart E11.5 vs Heart E14.5)
#to check whether the selected set is true for the conditions:
head(countData[, c(21:28)], 4)

countData_HTe11_HTe14<-countData[,21:28] 
coldata_HTe11_HTe14<-coldata[21:28,]

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

write.table(resHTe11_HTe14_Df, file="/home/yaaroba/ENCODEproject/Tables/RNAseq_differentialExpression_HTe11_HTe14_Hep2020_26Aug2023.txt",quote=FALSE,sep="\t",row.names=TRUE, col.names=TRUE)




# building table specific for HTe14_HTp0 analysis (Heart E14.5 vs Heart P0)
#to check whether the selected set is true for the conditions:
head(countData[, c(11:14,25:28)], 4)

countData_HTe14_HTp0<-countData[,c(11:14,25:28)]
coldata_HTe14_HTp0<-coldata[c(11:14,25:28),]

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

write.table(resHTe14_HTp0_Df, file="/home/yaaroba/ENCODEproject/Tables/RNAseq_differentialExpression_HTe14_HTp0_HeP2020_27Aug2023.txt",quote=FALSE,sep="\t",row.names=TRUE, col.names=TRUE)




# building table specific for HTe11_HTp0 analysis (Heart E11.5 vs Heart P0)
#to check whether the selected set is true for the conditions:
head(countData[, c(11:14,21:24)], 4)

countData_HTe11_HTp0<-countData[,c(11:14,21:24)]
coldata_HTe11_HTp0<-coldata[c(11:14,21:24),]

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

write.table(resHTe11_HTp0_Df, file="/home/yaaroba/ENCODEproject/Tables/RNAseq_differentialExpression_HTe11_HTp0_HeP2020_27Aug2023.txt",quote=FALSE,sep="\t",row.names=TRUE, col.names=TRUE)




# building table specific for LVe11_LVe14 analysis (Liver E11.5 vs Liver E14.5)
#to check whether the selected set is true for the conditions:
head(countData[, c(17:20,29:32)], 4)

countData_LVe11_LVe14<-countData[,c(17:20,29:32)] 
coldata_LVe11_LVe14<-coldata[c(17:20,29:32),]

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

write.table(resLVe11_LVe14_Df, file="/home/yaaroba/ENCODEproject/Tables/RNAseq_differentialExpression_LVe11_LVe14_Hep2020_27Aug2023.txt",quote=FALSE,sep="\t",row.names=TRUE, col.names=TRUE)




# building table specific for LVe14_LVp0 analysis (Liver E14.5 vs Heart P0)
#to check whether the selected set is true for the conditions:
head(countData[, c(29:36)], 4)

countData_LVe14_LVp0<-countData[,29:36]
coldata_LVe14_LVp0<-coldata[29:36,]

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

write.table(resLVe14_LVp0_Df, file="/home/yaaroba/ENCODEproject/Tables/RNAseq_differentialExpression_LVe14_LVp0_HeP2020_27Aug2023.txt",quote=FALSE,sep="\t",row.names=TRUE, col.names=TRUE)




# building table specific for LVe11_LVp0 analysis (Liver E11.5 vs Heart P0)
#to check whether the selected set is true for the conditions:
head(countData[, c(17:20,33:36)], 4)

countData_LVe11_LVp0<-countData[,c(17:20,33:36)]
coldata_LVe11_LVp0<-coldata[c(17:20,33:36),]

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

write.table(resLVe11_LVp0_Df, file="/home/yaaroba/ENCODEproject/Tables/RNAseq_differentialExpression_LVe11_LVp0_HeP2020_27Aug2023.txt",quote=FALSE,sep="\t",row.names=TRUE, col.names=TRUE)




### R outputs, processing intermidaite results:
> summary(res_FBe11_FBe14_05)
out of 34782 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)       : 11243, 32%
LFC < 0 (down)     : 3077, 8.8%
outliers [1]       : 0, 0%
low counts [2]     : 9772, 28%
(mean count < 1)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results


> summary(res_FBe14_FBp0_05)
out of 35493 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)       : 1155, 3.3%
LFC < 0 (down)     : 1808, 5.1%
outliers [1]       : 0, 0%
low counts [2]     : 11325, 32%
(mean count < 1)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results



> summary(res_FBe11_FBp0_05)
out of 32715 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)       : 3083, 9.4%
LFC < 0 (down)     : 935, 2.9%
outliers [1]       : 0, 0%
low counts [2]     : 15847, 48%
(mean count < 4)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results



> summary(res_HTe11_HTe14_05) 
out of 36050 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)       : 4283, 12%
LFC < 0 (down)     : 4732, 13%
outliers [1]       : 0, 0%
low counts [2]     : 10841, 30%
(mean count < 1)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results


> summary(res_HTe14_HTp0_05)
out of 35713 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)       : 6353, 18%
LFC < 0 (down)     : 6346, 18%
outliers [1]       : 0, 0%
low counts [2]     : 10731, 30%
(mean count < 1)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results


> summary(res_HTe11_HTp0_05)
out of 35880 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)       : 7234, 20%
LFC < 0 (down)     : 7531, 21%
outliers [1]       : 0, 0%
low counts [2]     : 8763, 24%
(mean count < 1)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results



> summary(res_LVe11_LVe14_05)
out of 31751 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)       : 8717, 27%
LFC < 0 (down)     : 2705, 8.5%
outliers [1]       : 0, 0%
low counts [2]     : 8847, 28%
(mean count < 1)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results


> summary(res_LVe14_LVp0_05)
out of 30581 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)       : 3164, 10%
LFC < 0 (down)     : 14023, 46%
outliers [1]       : 0, 0%
low counts [2]     : 10186, 33%
(mean count < 1)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results



> summary(res_LVe11_LVp0_05)
out of 23678 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)       : 2816, 12%
LFC < 0 (down)     : 6970, 29%
outliers [1]       : 0, 0%
low counts [2]     : 8073, 34%
(mean count < 0)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results
