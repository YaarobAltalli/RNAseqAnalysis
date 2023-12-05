### Gorkin et al. Nature 2020 - Bing Ren and Len A. Pennacchio - "An atlas of dynamic chromatin landscapes in mouse fetal development" 
# RAW data: https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA63471&o=acc_s%3Aa
# Type of data analyzed: ATAC-seq 35-51 nucleotide
#						 RNA-seq 100 nucleotide


##0. Setting up folders and preparing structure of the workflow:
#Workflow Folder Structure:
mkdir /home/yaaroba/ENCODEproject
mkdir /home/yaaroba/ENCODEproject/fastqcCheck/
mkdir /home/yaaroba/ENCODEproject/multiQC
mkdir /home/yaaroba/ENCODEproject/fastqcCheckTrimmed/
mkdir -p /home/yaaroba/ENCODEproject/trimmedFastq/RNAseq  
mkdir /home/yaaroba/ENCODEproject/trimmedFastq/ATACseq
mkdir /home/yaaroba/ENCODEproject/STARalingments
mkdir /home/yaaroba/ENCODEproject/STARFinalOut
mkdir /home/yaaroba/ENCODEproject/STAR_geneCounts
mkdir /home/yaaroba/ENCODEproject/bamFiles
mkdir /home/yaaroba/ENCODEproject/RNAseqBamFiles
mkdir /home/yaaroba/ENCODEproject/strandedBamFiles
mkdir /home/yaaroba/ENCODEproject/RNAseqBedGraphs
mkdir /home/yaaroba/ENCODEproject/RNAseqNormBedGraphs
mkdir /home/yaaroba/ENCODEproject/RNAseqBigWigs
mkdir -p /home/yaaroba/ENCODEproject/analysis/HeP_Nature_2020_RNAseq
mkdir /home/yaaroba/ENCODEproject/analysis/GorkinD_Nature_2020_ATACseq
mkdir /home/yaaroba/ENCODEproject/summaries

#Folders and files needed to download files using SRAdownloader:
#make folder for this project within sraRunTables folder:
mkdir /home/yaaroba/sraRunTables/ENCODEproject #this is the folder where I will be saving all SRAtables of ENCODE project following certain nomiclature criteria: SraRunTable_authorInitialsOrName_journalInitials_yearPublished_typeOfDataIncluded_theLabProducingTheData
mkdir /home/yaaroba/sraRunTables/ENCODEproject/part1 #this folder will include the first part of analysis: 3 tissues (Forebrain FB, Heart HT, Liver LV) at 3 different timepoints (E11.5, E14.5 and P0) 
#upload the SraRunTables to this folder maunally (using upload bottom on the terminal)
#tables:
/home/yaaroba/sraRunTables/ENCODEproject/part1/SraRunTable_GorkinD_Nature_2020_ATAC_LiverHeartForebrain_LPennacchio.txt
/home/yaaroba/sraRunTables/ENCODEproject/part1/SraRunTable_HeP_Nature_2020_RNA_LiverHeartForebrain_LPennacchio_BWold.txt

#Preparing folders within Dysk2/SRA_repository folder to store the rawdata of each scientific article:
mkdir -p /dysk2/SRA_repository/ENCODEproject/GorkinD_Nature_2020
mkdir /dysk2/SRA_repository/ENCODEproject/HeP_Nature_2020


##1. Downloading the data using SRAdownloader: 
#!!!!!(please dont do this step unless you are downloading new data to this repository, the provided data are already downloaded)!!!!!#

#sradownloader --outdir (where to save the output) (direction to the saved SraRunTable.txt in your home directory)
sradownloader --outdir /dysk2/SRA_repository/ENCODEproject/GorkinD_Nature_2020 /home/yaaroba/sraRunTables/ENCODEproject/part1/SraRunTable_GorkinD_Nature_2020_ATAC_LiverHeartForebrain_LPennacchio.txt

sradownloader --outdir /dysk2/SRA_repository/ENCODEproject/HeP_Nature_2020 /home/yaaroba/sraRunTables/ENCODEproject/part1/SraRunTable_HeP_Nature_2020_RNA_LiverHeartForebrain_LPennacchio_BWold.txt

#In the naming of these files the age is not included and there is "from" word which is not needed so to remove this word I applied this code:
cd /dysk2/SRA_repository/ENCODEproject/GorkinD_Nature_2020/
for file in *.fastq.gz; do
    # Extract all words by splitting the filename using '_'
    words=($(echo "$file" | tr '_' ' '))
    # Remove the forth word "from" from the array
    unset "words[3]"
    # Reconstruct the filename with the remaining words separated by '_'
    new_filename=$(IFS=_ ; echo "${words[*]}")
    # Rename the file
    mv "$file" "$new_filename"
done #Repeat the same for datasets here: /dysk2/SRA_repository/ENCODEproject/HeP_Nature_2020
#I added the age manually to the name of each fastq file


##2. Analyse read properties/ quality with fastqc (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/); 
#command line; all files in one go:
fastqc --outdir /home/yaaroba/ENCODEproject/fastqcCheck/ --format fastq --threads 40 /dysk2/SRA_repository/ENCODEproject/GorkinD_Nature_2020/*fastq.gz
fastqc --outdir /home/yaaroba/ENCODEproject/fastqcCheck/ --format fastq --threads 40 /dysk2/SRA_repository/ENCODEproject/HeP_Nature_2020/*fastq.gz


##3. Perform Adapter trimming with TrimGalore (cutadapt wrppper; https://github.com/FelixKrueger/TrimGalore):
# basic syntax: trim_galore [options] <filename(s)>

### trimming single-end reads
for file in /dysk2/SRA_repository/ENCODEproject/HeP_Nature_2020/*.fastq.gz; do trim_galore $file -o /home/yaaroba/ENCODEproject/trimmedFastq/HeP_Nature_2020_RNAseq & done &

### trimming paired-end reads
for R1 in /dysk2/SRA_repository/ENCODEproject/GorkinD_Nature_2020/*_ATAC-seq_1.fastq.gz
do
   R2=${R1//_ATAC-seq_1.fastq.gz/_ATAC-seq_2.fastq.gz}
   trim_galore --paired $R1 $R2 -o /home/yaaroba/ENCODEproject/trimmedFastq/GorkinD_Nature_2020_ATACseq & done 
done


##4. Fatstqc quality check after adapter trimming:
fastqc --outdir /home/yaaroba/ENCODEproject/fastqcCheckTrimmed --format fastq --threads 40 /home/yaaroba/ENCODEproject/trimmedFastq/HeP_Nature_2020_RNAseq/*fq.gz
fastqc --outdir /home/yaaroba/ENCODEproject/fastqcCheckTrimmed --format fastq --threads 40 /home/yaaroba/ENCODEproject/trimmedFastq/GorkinD_Nature_2020_ATACseq/*fq.gz


##5. MultiQC: summarize fastQC
multiqc /home/yaaroba/ENCODEproject/fastqcCheckTrimmed -n multiQC_trimmed -o /home/yaaroba/ENCODEproject/multiQC