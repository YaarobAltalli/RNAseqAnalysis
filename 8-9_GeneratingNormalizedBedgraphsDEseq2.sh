### 8. Estimating sizing factors Deseq2 per milon factors, generating tables with raw and normalized gene counts and stats for the libraries:

# move ReadsPerGene.out.tab files generateed by STAR to previously created "STAR_geneCounts" folder
# check in STAR manual what is their structure: https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf
# Use the first column of ReadsPerGene.out.tab as the first column of the input to DEseq2
#Use the 2nd, 3rd or 4th column of ReadsPerGene.out.tab as the 2nd column depending on the strandedness of your library:
# 2nd column - for unstranded data
# 3rd column  - for 1st read agreeing with RNA strand
# 4th column  - for 2nd read agreeing with RNA strand (typical for Illumina stranded Tru-seq)
#  Cut out the first 4 lines of the ReadsPerGene.out.tab that contain counts for non-genic read (unmapped/multimappers/ambiguous/noFeature).
# variant to do this before launching R: for instance for Illumina stranded Tru-seq: 
# awk 'NR>4 {print $1 "\t" $4}' ReadsPerGene.out.tab DEseq.input

mv /home/yaaroba/ENCODEproject/STARalingments/*_ReadsPerGene.out.tab /home/yaaroba/ENCODEproject/STAR_geneCounts
cd /home/yaaroba/ENCODEproject/STAR_geneCounts/

R
library(DESeq2)
library(data.table)

where<-"/home/yaaroba/ENCODEproject/STAR_geneCounts/"
setwd(where)

# buliding sample table necessary for DEseq2 analysis
sampleFiles <- grep("ReadsPerGene.out.tab",list.files(where),value=TRUE)

### in analyzed data all RNA-seq libraries: files 1-36 unstranded => 2nd column

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

# creating "DESeqDataSet" object for estimating sizing factors and further differential gee expression analysis
dds = DESeqDataSetFromMatrix(countData = countData, colData = coldata, design = ~ condition)

readCounts<-counts(dds, normalized = FALSE, replaced = FALSE)
readCountsDF<-as.data.frame(readCounts, stringsAsFactors = FALSE)
readCountsDF$id<-rownames(readCountsDF)

dds<-estimateSizeFactors(dds)
SF<-sizeFactors(dds)

readCountsNorm<-counts(dds, normalized = TRUE, replaced = FALSE)
readCountsNormDF<-as.data.frame(readCountsNorm, stringsAsFactors = FALSE)
readCountsNormDF$id<-rownames(readCountsNormDF)

# creating table with gene counts for all samples; adding gene info: coordinates, gene name, gene type
library(rtracklayer)
gtf<- import("/home/micgdu/GenomicData/genomesDec2022/gencode.vM31.primary_assembly.annotation.gtf.gz")
gtf_genes<-gtf[which(gtf$type=="gene")]
geneInfo<-as.data.frame(gtf_genes[,c(5,6,7)], stringsAsFactors = FALSE)
rownames(geneInfo)<-geneInfo$gene_id

# retrieving transcript length from gtf file 
library(GenomicFeatures)
gtf <- "/home/micgdu/GenomicData/genomesDec2022/gencode.vM31.primary_assembly.annotation.gtf.gz"
txdb.filename <- "vM31.primary_assembly.mm39.txdb"
txdb <- makeTxDbFromGFF(gtf)
saveDb(txdb, txdb.filename)

mm39_txlens <- transcriptLengths(txdb)
head(mm39_txlens)

gene_transcript_length<- function(i,db) {
        median(db$tx_len[which(db$gene_id==gtf_genes$gene_id[i])])
}

Median_Transcript_length<-unlist(lapply(1:length(gtf_genes), gene_transcript_length, mm39_txlens))
geneInfo$txLens<-Median_Transcript_length

annotatedCounts<-cbind(geneInfo,readCountsDF)

write.table(annotatedCounts, file="/home/yaaroba/ENCODEproject/Tables/RNAseqRawReadsPerGene_HeP_Nature_2020_8aug2023.txt",quote=FALSE,sep="\t",row.names=TRUE, col.names=TRUE)

# analogous table with gene counts normalized with sizing factors
annotatedNormCounts<-cbind(geneInfo,readCountsNormDF)

save(annotatedNormCounts, file="/home/yaaroba/ENCODEproject/Tables/annotatedNormCounts.RData")

write.table(annotatedNormCounts, file="/home/yaaroba/ENCODEproject/Tables/RNAseqNormalizedReadsPerGene_HeP_Nature_2020_8aug2023.txt",quote=FALSE,sep="\t",row.names=TRUE, col.names=TRUE)

# analogous table with fpkm (fragment counts normalized per kilobase of feature length per million
# mapped fragments based on annotatedNormCounts and retrieved median transcript length

readCountsNormDF_FPKM<-readCountsNormDF
readCountsNormDF_FPKM$id<-NULL
readCountsNormDF_FPKM<-readCountsNormDF_FPKM/(geneInfo$txLens/1000)
annotatedNormCountsFPKM<-cbind(geneInfo,readCountsNormDF_FPKM)

save(annotatedNormCountsFPKM, file="/home/yaaroba/ENCODEproject/Tables/annotatedNormCountsFPKM.RData")

write.table(annotatedNormCountsFPKM, file="/home/yaaroba/ENCODEproject/Tables/RNAseqNormalizedReadsPerGene_FPKM_HeP_Nature_2020_8aug2023.txt",quote=FALSE,sep="\t",row.names=TRUE, col.names=TRUE)

# generating stats for the libraries:

sampleStats<-rbind(countDataHead[1:4,]/1e6,colSums(readCounts)/1e6)
featuredNorm<- sampleStats[5,]/SF
SF_rpm<-SF*mean(as.numeric(featuredNorm)) #sizing factors adjusted by single factor: mean level of rpm within exons in analyzed samples for convienience of track comparisons between experiment sets 
sampleStats<-rbind(colSums(sampleStats),sampleStats,featuredNorm, SF, SF_rpm)
rownames(sampleStats)<- c("all_reads", "N_unmapped","N_multimapping","N_noFeature","N_ambiguous","exons","exonsNorm", "SF", "SF_rpm")
sampleStats<-t(sampleStats)

write.table(sampleStats, file="/home/yaaroba/ENCODEproject/summaries/RNAseqReadContStats_HeP_Nature_2020_8aug2023.txt",quote=FALSE,sep="\t",row.names=TRUE, col.names=TRUE)


##9. Creating normalized bedGraph files with estimated rpm normalized sizzing factors

# continuew with:
where<-"/home/yaaroba/ENCODEproject/RNAseqBedGraphs"
setwd(where)
sampleFiles <- grep(".bedGraph.gz",list.files(where),value=TRUE)

sampleNames<-c("forebrain_E11.5_rep1_Fwd","forebrain_E11.5_rep1_Rev","forebrain_E11.5_rep2_Fwd","forebrain_E11.5_rep2_Rev","forebrain_E11.5_rep3_Fwd","forebrain_E11.5_rep3_Rev","forebrain_E11.5_rep4_Fwd","forebrain_E11.5_rep4_Rev","forebrain_E14.5_rep1_Fwd","forebrain_E14.5_rep1_Rev","forebrain_E14.5_rep2_Fwd","forebrain_E14.5_rep2_Rev","forebrain_E14.5_rep3_Fwd","forebrain_E14.5_rep3_Rev","forebrain_E14.5_rep4_Fwd","forebrain_E14.5_rep4_Rev","forebrain_P0_rep1_Fwd","forebrain_P0_rep1_Rev","forebrain_P0_rep2_Fwd","forebrain_P0_rep2_Rev","heart_P0_rep1_Fwd","heart_P0_rep1_Rev","heart_P0_rep2_Fwd","heart_P0_rep2_Rev","heart_P0_rep3_Fwd","heart_P0_rep3_Rev","heart_P0_rep4_Fwd","heart_P0_rep4_Rev","forebrain_P0_rep3_Fwd","forebrain_P0_rep3_Rev","forebrain_P0_rep4_Fwd","forebrain_P0_rep4_Rev","liver_E11.5_rep1_Fwd","liver_E11.5_rep1_Rev","liver_E11.5_rep2_Fwd","liver_E11.5_rep2_Rev","liver_E11.5_rep3_Fwd","liver_E11.5_rep3_Rev","liver_E11.5_rep4_Fwd","liver_E11.5_rep4_Rev","heart_E11.5_rep1_Fwd","heart_E11.5_rep1_Rev","heart_E11.5_rep2_Fwd","heart_E11.5_rep2_Rev","heart_E11.5_rep3_Fwd","heart_E11.5_rep3_Rev","heart_E11.5_rep4_Fwd","heart_E11.5_rep4_Rev","heart_E14.5_rep1_Fwd","heart_E14.5_rep1_Rev","heart_E14.5_rep2_Fwd","heart_E14.5_rep2_Rev","heart_E14.5_rep3_Fwd","heart_E14.5_rep3_Rev","heart_E14.5_rep4_Fwd","heart_E14.5_rep4_Rev","liver_E14.5_rep1_Fwd","liver_E14.5_rep1_Rev","liver_E14.5_rep2_Fwd","liver_E14.5_rep2_Rev","liver_E14.5_rep3_Fwd","liver_E14.5_rep3_Rev","liver_E14.5_rep4_Fwd","liver_E14.5_rep4_Rev","liver_P0_rep1_Fwd","liver_P0_rep1_Rev","liver_P0_rep2_Fwd","liver_P0_rep2_Rev","liver_P0_rep3_Fwd","liver_P0_rep3_Rev","liver_P0_rep4_Fwd","liver_P0_rep4_Rev")

strand<- rep(a<-c("Fwd","Rev"), times=36)

SampleFactor<-rep(SF_rpm, each = 2)

sampleTable <- data.frame(fileName = sampleFiles, name=sampleNames, SampleFactor, strand)

bedtoolsBdgNorm<-function(x){
  # function based on sampleTable and DESeq2 sizing factors (SF) & stranded bedGraphs created before
  # all SF multiplied by the same factor: mean ~rpm of uniq reads in exons for analysed batch of samples
  # to faciliate vizual comparison between sample batches from different experiments   
  f<-sampleTable$SampleFactor[x]
  a<-read.table(sampleTable[x,1], header=FALSE, stringsAsFactors=FALSE)
  a[,2]<-format(a[,2], scientific=FALSE)
  a[,3]<-format(a[,3], scientific=FALSE)
  if (sampleTable$strand[x]=="Fwd") {
  a[,4]<-format(a[,4]/f,scientific=FALSE)} else { 
  a[,4]<-format(-a[,4]/f, scientific=FALSE)}
  a<-a[which(a[,1]%in% chr),]
  write.table(a ,file=paste(sampleTable$name[x], "_norm_",sampleTable$strand[x],".bedGraph", sep=""),quote=FALSE,sep="\t",row.names=FALSE, col.names=FALSE)
} 

chr<-c(paste("chr",1:22,sep=""),"chrX", "chrY", "chrM")      # limiting output to standard chromosomes

for (i in 1:dim(sampleTable)[1]) bedtoolsBdgNorm(i)

q()

mv /home/yaaroba/ENCODEproject/RNAseqBedGraphs/*.bedGraph /home/yaaroba/ENCODEproject/RNAseqNormBedGraphs





################################################################################################################################################

> saveDb(txdb, txdb.filename)
TxDb object:
# Db type: TxDb
# Supporting package: GenomicFeatures
# Data source: /home/micgdu/GenomicData/genomesDec2022/gencode.vM31.primary_assembly.annotation.gtf.gz
# Organism: NA
# Taxonomy ID: NA
# miRBase build ID: NA
# Genome: NA
# Nb of transcripts: 149482
# Db created by: GenomicFeatures package from Bioconductor
# Creation time: 2023-08-08 14:09:48 +0000 (Tue, 08 Aug 2023)
# GenomicFeatures version at creation time: 1.52.0
# RSQLite version at creation time: 2.3.1
# DBSCHEMAVERSION: 1.2
> head(mm39_txlens)
  tx_id              tx_name              gene_id nexon tx_len
1     1 ENSMUST00000193812.2 ENSMUSG00000102693.2     1   1070
2     2 ENSMUST00000082908.3 ENSMUSG00000064842.3     1    110
3     3 ENSMUST00000192857.2 ENSMUSG00000102851.2     1    480
4     4 ENSMUST00000161581.2 ENSMUSG00000089699.2     2    250
5     5 ENSMUST00000192183.2 ENSMUSG00000103147.2     1    926
6     6 ENSMUST00000193244.2 ENSMUSG00000102348.2     1   1634




> sampleFiles


 [1] "SRR3191782_GSM2071304_RNA-seq_forebrain_E11.5_ENCLB917PKP_Mus_musculus_RNA-Seq_1_trimmed.fq.gz_Aligned.out.bam_Fwd.bedGraph.gz"
 [2] "SRR3191782_GSM2071304_RNA-seq_forebrain_E11.5_ENCLB917PKP_Mus_musculus_RNA-Seq_1_trimmed.fq.gz_Aligned.out.bam_Rev.bedGraph.gz"
 [3] "SRR3191783_GSM2071304_RNA-seq_forebrain_E11.5_ENCLB917PKP_Mus_musculus_RNA-Seq_1_trimmed.fq.gz_Aligned.out.bam_Fwd.bedGraph.gz"
 [4] "SRR3191783_GSM2071304_RNA-seq_forebrain_E11.5_ENCLB917PKP_Mus_musculus_RNA-Seq_1_trimmed.fq.gz_Aligned.out.bam_Rev.bedGraph.gz"
 [5] "SRR3191784_GSM2071305_RNA-seq_forebrain_E11.5_ENCLB026BHP_Mus_musculus_RNA-Seq_1_trimmed.fq.gz_Aligned.out.bam_Fwd.bedGraph.gz"
 [6] "SRR3191784_GSM2071305_RNA-seq_forebrain_E11.5_ENCLB026BHP_Mus_musculus_RNA-Seq_1_trimmed.fq.gz_Aligned.out.bam_Rev.bedGraph.gz"
 [7] "SRR3191785_GSM2071305_RNA-seq_forebrain_E11.5_ENCLB026BHP_Mus_musculus_RNA-Seq_1_trimmed.fq.gz_Aligned.out.bam_Fwd.bedGraph.gz"
 [8] "SRR3191785_GSM2071305_RNA-seq_forebrain_E11.5_ENCLB026BHP_Mus_musculus_RNA-Seq_1_trimmed.fq.gz_Aligned.out.bam_Rev.bedGraph.gz"
 [9] "SRR3191859_GSM2071345_RNA-seq_forebrain_E14.5_ENCLB964APA_Mus_musculus_RNA-Seq_1_trimmed.fq.gz_Aligned.out.bam_Fwd.bedGraph.gz"
[10] "SRR3191859_GSM2071345_RNA-seq_forebrain_E14.5_ENCLB964APA_Mus_musculus_RNA-Seq_1_trimmed.fq.gz_Aligned.out.bam_Rev.bedGraph.gz"
[11] "SRR3191860_GSM2071345_RNA-seq_forebrain_E14.5_ENCLB964APA_Mus_musculus_RNA-Seq_1_trimmed.fq.gz_Aligned.out.bam_Fwd.bedGraph.gz"
[12] "SRR3191860_GSM2071345_RNA-seq_forebrain_E14.5_ENCLB964APA_Mus_musculus_RNA-Seq_1_trimmed.fq.gz_Aligned.out.bam_Rev.bedGraph.gz"
[13] "SRR3191861_GSM2071346_RNA-seq_forebrain_E14.5_ENCLB336GOL_Mus_musculus_RNA-Seq_1_trimmed.fq.gz_Aligned.out.bam_Fwd.bedGraph.gz"
[14] "SRR3191861_GSM2071346_RNA-seq_forebrain_E14.5_ENCLB336GOL_Mus_musculus_RNA-Seq_1_trimmed.fq.gz_Aligned.out.bam_Rev.bedGraph.gz"
[15] "SRR3191862_GSM2071346_RNA-seq_forebrain_E14.5_ENCLB336GOL_Mus_musculus_RNA-Seq_1_trimmed.fq.gz_Aligned.out.bam_Fwd.bedGraph.gz"
[16] "SRR3191862_GSM2071346_RNA-seq_forebrain_E14.5_ENCLB336GOL_Mus_musculus_RNA-Seq_1_trimmed.fq.gz_Aligned.out.bam_Rev.bedGraph.gz"
[17] "SRR3191960_GSM2071398_RNA-seq_forebrain_P0_ENCLB055JUC_Mus_musculus_RNA-Seq_1_trimmed.fq.gz_Aligned.out.bam_Fwd.bedGraph.gz"
[18] "SRR3191960_GSM2071398_RNA-seq_forebrain_P0_ENCLB055JUC_Mus_musculus_RNA-Seq_1_trimmed.fq.gz_Aligned.out.bam_Rev.bedGraph.gz"
[19] "SRR3191961_GSM2071398_RNA-seq_forebrain_P0_ENCLB055JUC_Mus_musculus_RNA-Seq_1_trimmed.fq.gz_Aligned.out.bam_Fwd.bedGraph.gz"
[20] "SRR3191961_GSM2071398_RNA-seq_forebrain_P0_ENCLB055JUC_Mus_musculus_RNA-Seq_1_trimmed.fq.gz_Aligned.out.bam_Rev.bedGraph.gz"
[21] "SRR3192004_GSM2071437_RNA-seq_heart_P0_ENCLB658ICO_Mus_musculus_RNA-Seq_1_trimmed.fq.gz_Aligned.out.bam_Fwd.bedGraph.gz"
[22] "SRR3192004_GSM2071437_RNA-seq_heart_P0_ENCLB658ICO_Mus_musculus_RNA-Seq_1_trimmed.fq.gz_Aligned.out.bam_Rev.bedGraph.gz"
[23] "SRR3192005_GSM2071437_RNA-seq_heart_P0_ENCLB658ICO_Mus_musculus_RNA-Seq_1_trimmed.fq.gz_Aligned.out.bam_Fwd.bedGraph.gz"
[24] "SRR3192005_GSM2071437_RNA-seq_heart_P0_ENCLB658ICO_Mus_musculus_RNA-Seq_1_trimmed.fq.gz_Aligned.out.bam_Rev.bedGraph.gz"
[25] "SRR3192006_GSM2071438_RNA-seq_heart_P0_ENCLB741KQB_Mus_musculus_RNA-Seq_1_trimmed.fq.gz_Aligned.out.bam_Fwd.bedGraph.gz"
[26] "SRR3192006_GSM2071438_RNA-seq_heart_P0_ENCLB741KQB_Mus_musculus_RNA-Seq_1_trimmed.fq.gz_Aligned.out.bam_Rev.bedGraph.gz"
[27] "SRR3192007_GSM2071438_RNA-seq_heart_P0_ENCLB741KQB_Mus_musculus_RNA-Seq_1_trimmed.fq.gz_Aligned.out.bam_Fwd.bedGraph.gz"
[28] "SRR3192007_GSM2071438_RNA-seq_heart_P0_ENCLB741KQB_Mus_musculus_RNA-Seq_1_trimmed.fq.gz_Aligned.out.bam_Rev.bedGraph.gz"
[29] "SRR3192012_GSM2071440_RNA-seq_forebrain_P0_ENCLB181TCJ_Mus_musculus_RNA-Seq_1_trimmed.fq.gz_Aligned.out.bam_Fwd.bedGraph.gz"
[30] "SRR3192012_GSM2071440_RNA-seq_forebrain_P0_ENCLB181TCJ_Mus_musculus_RNA-Seq_1_trimmed.fq.gz_Aligned.out.bam_Rev.bedGraph.gz"
[31] "SRR3192013_GSM2071440_RNA-seq_forebrain_P0_ENCLB181TCJ_Mus_musculus_RNA-Seq_1_trimmed.fq.gz_Aligned.out.bam_Fwd.bedGraph.gz"
[32] "SRR3192013_GSM2071440_RNA-seq_forebrain_P0_ENCLB181TCJ_Mus_musculus_RNA-Seq_1_trimmed.fq.gz_Aligned.out.bam_Rev.bedGraph.gz"
[33] "SRR3192046_GSM2071461_RNA-seq_liver_E11.5_ENCLB905LVV_Mus_musculus_RNA-Seq_1_trimmed.fq.gz_Aligned.out.bam_Fwd.bedGraph.gz"
[34] "SRR3192046_GSM2071461_RNA-seq_liver_E11.5_ENCLB905LVV_Mus_musculus_RNA-Seq_1_trimmed.fq.gz_Aligned.out.bam_Rev.bedGraph.gz"
[35] "SRR3192047_GSM2071461_RNA-seq_liver_E11.5_ENCLB905LVV_Mus_musculus_RNA-Seq_1_trimmed.fq.gz_Aligned.out.bam_Fwd.bedGraph.gz"
[36] "SRR3192047_GSM2071461_RNA-seq_liver_E11.5_ENCLB905LVV_Mus_musculus_RNA-Seq_1_trimmed.fq.gz_Aligned.out.bam_Rev.bedGraph.gz"
[37] "SRR3192048_GSM2071462_RNA-seq_liver_E11.5_ENCLB449LBZ_Mus_musculus_RNA-Seq_1_trimmed.fq.gz_Aligned.out.bam_Fwd.bedGraph.gz"
[38] "SRR3192048_GSM2071462_RNA-seq_liver_E11.5_ENCLB449LBZ_Mus_musculus_RNA-Seq_1_trimmed.fq.gz_Aligned.out.bam_Rev.bedGraph.gz"
[39] "SRR3192049_GSM2071462_RNA-seq_liver_E11.5_ENCLB449LBZ_Mus_musculus_RNA-Seq_1_trimmed.fq.gz_Aligned.out.bam_Fwd.bedGraph.gz"
[40] "SRR3192049_GSM2071462_RNA-seq_liver_E11.5_ENCLB449LBZ_Mus_musculus_RNA-Seq_1_trimmed.fq.gz_Aligned.out.bam_Rev.bedGraph.gz"
[41] "SRR3192141_GSM2071523_RNA-seq_heart_E11.5_ENCLB347FRI_Mus_musculus_RNA-Seq_1_trimmed.fq.gz_Aligned.out.bam_Fwd.bedGraph.gz"
[42] "SRR3192141_GSM2071523_RNA-seq_heart_E11.5_ENCLB347FRI_Mus_musculus_RNA-Seq_1_trimmed.fq.gz_Aligned.out.bam_Rev.bedGraph.gz"
[43] "SRR3192142_GSM2071523_RNA-seq_heart_E11.5_ENCLB347FRI_Mus_musculus_RNA-Seq_1_trimmed.fq.gz_Aligned.out.bam_Fwd.bedGraph.gz"
[44] "SRR3192142_GSM2071523_RNA-seq_heart_E11.5_ENCLB347FRI_Mus_musculus_RNA-Seq_1_trimmed.fq.gz_Aligned.out.bam_Rev.bedGraph.gz"
[45] "SRR3192143_GSM2071524_RNA-seq_heart_E11.5_ENCLB601XLL_Mus_musculus_RNA-Seq_1_trimmed.fq.gz_Aligned.out.bam_Fwd.bedGraph.gz"
[46] "SRR3192143_GSM2071524_RNA-seq_heart_E11.5_ENCLB601XLL_Mus_musculus_RNA-Seq_1_trimmed.fq.gz_Aligned.out.bam_Rev.bedGraph.gz"
[47] "SRR3192144_GSM2071524_RNA-seq_heart_E11.5_ENCLB601XLL_Mus_musculus_RNA-Seq_1_trimmed.fq.gz_Aligned.out.bam_Fwd.bedGraph.gz"
[48] "SRR3192144_GSM2071524_RNA-seq_heart_E11.5_ENCLB601XLL_Mus_musculus_RNA-Seq_1_trimmed.fq.gz_Aligned.out.bam_Rev.bedGraph.gz"
[49] "SRR3192159_GSM2071534_RNA-seq_heart_E14.5_ENCLB931ORG_Mus_musculus_RNA-Seq_1_trimmed.fq.gz_Aligned.out.bam_Fwd.bedGraph.gz"
[50] "SRR3192159_GSM2071534_RNA-seq_heart_E14.5_ENCLB931ORG_Mus_musculus_RNA-Seq_1_trimmed.fq.gz_Aligned.out.bam_Rev.bedGraph.gz"
[51] "SRR3192160_GSM2071534_RNA-seq_heart_E14.5_ENCLB931ORG_Mus_musculus_RNA-Seq_1_trimmed.fq.gz_Aligned.out.bam_Fwd.bedGraph.gz"
[52] "SRR3192160_GSM2071534_RNA-seq_heart_E14.5_ENCLB931ORG_Mus_musculus_RNA-Seq_1_trimmed.fq.gz_Aligned.out.bam_Rev.bedGraph.gz"
[53] "SRR3192161_GSM2071535_RNA-seq_heart_E14.5_ENCLB584VFZ_Mus_musculus_RNA-Seq_1_trimmed.fq.gz_Aligned.out.bam_Fwd.bedGraph.gz"
[54] "SRR3192161_GSM2071535_RNA-seq_heart_E14.5_ENCLB584VFZ_Mus_musculus_RNA-Seq_1_trimmed.fq.gz_Aligned.out.bam_Rev.bedGraph.gz"
[55] "SRR3192162_GSM2071535_RNA-seq_heart_E14.5_ENCLB584VFZ_Mus_musculus_RNA-Seq_1_trimmed.fq.gz_Aligned.out.bam_Fwd.bedGraph.gz"
[56] "SRR3192162_GSM2071535_RNA-seq_heart_E14.5_ENCLB584VFZ_Mus_musculus_RNA-Seq_1_trimmed.fq.gz_Aligned.out.bam_Rev.bedGraph.gz"
[57] "SRR3192259_GSM2071597_RNA-seq_liver_E14.5_ENCLB278MMD_Mus_musculus_RNA-Seq_1_trimmed.fq.gz_Aligned.out.bam_Fwd.bedGraph.gz"
[58] "SRR3192259_GSM2071597_RNA-seq_liver_E14.5_ENCLB278MMD_Mus_musculus_RNA-Seq_1_trimmed.fq.gz_Aligned.out.bam_Rev.bedGraph.gz"
[59] "SRR3192260_GSM2071597_RNA-seq_liver_E14.5_ENCLB278MMD_Mus_musculus_RNA-Seq_1_trimmed.fq.gz_Aligned.out.bam_Fwd.bedGraph.gz"
[60] "SRR3192260_GSM2071597_RNA-seq_liver_E14.5_ENCLB278MMD_Mus_musculus_RNA-Seq_1_trimmed.fq.gz_Aligned.out.bam_Rev.bedGraph.gz"
[61] "SRR3192261_GSM2071598_RNA-seq_liver_E14.5_ENCLB200EFP_Mus_musculus_RNA-Seq_1_trimmed.fq.gz_Aligned.out.bam_Fwd.bedGraph.gz"
[62] "SRR3192261_GSM2071598_RNA-seq_liver_E14.5_ENCLB200EFP_Mus_musculus_RNA-Seq_1_trimmed.fq.gz_Aligned.out.bam_Rev.bedGraph.gz"
[63] "SRR3192262_GSM2071598_RNA-seq_liver_E14.5_ENCLB200EFP_Mus_musculus_RNA-Seq_1_trimmed.fq.gz_Aligned.out.bam_Fwd.bedGraph.gz"
[64] "SRR3192262_GSM2071598_RNA-seq_liver_E14.5_ENCLB200EFP_Mus_musculus_RNA-Seq_1_trimmed.fq.gz_Aligned.out.bam_Rev.bedGraph.gz"
[65] "SRR3192328_GSM2071642_RNA-seq_liver_P0_ENCLB356IIP_Mus_musculus_RNA-Seq_1_trimmed.fq.gz_Aligned.out.bam_Fwd.bedGraph.gz"
[66] "SRR3192328_GSM2071642_RNA-seq_liver_P0_ENCLB356IIP_Mus_musculus_RNA-Seq_1_trimmed.fq.gz_Aligned.out.bam_Rev.bedGraph.gz"
[67] "SRR3192329_GSM2071642_RNA-seq_liver_P0_ENCLB356IIP_Mus_musculus_RNA-Seq_1_trimmed.fq.gz_Aligned.out.bam_Fwd.bedGraph.gz"
[68] "SRR3192329_GSM2071642_RNA-seq_liver_P0_ENCLB356IIP_Mus_musculus_RNA-Seq_1_trimmed.fq.gz_Aligned.out.bam_Rev.bedGraph.gz"
[69] "SRR3192330_GSM2071642_RNA-seq_liver_P0_ENCLB356IIP_Mus_musculus_RNA-Seq_1_trimmed.fq.gz_Aligned.out.bam_Fwd.bedGraph.gz"
[70] "SRR3192330_GSM2071642_RNA-seq_liver_P0_ENCLB356IIP_Mus_musculus_RNA-Seq_1_trimmed.fq.gz_Aligned.out.bam_Rev.bedGraph.gz"
[71] "SRR3192331_GSM2071642_RNA-seq_liver_P0_ENCLB356IIP_Mus_musculus_RNA-Seq_1_trimmed.fq.gz_Aligned.out.bam_Fwd.bedGraph.gz"
[72] "SRR3192331_GSM2071642_RNA-seq_liver_P0_ENCLB356IIP_Mus_musculus_RNA-Seq_1_trimmed.fq.gz_Aligned.out.bam_Rev.bedGraph.gz"

