##6. STAR Alignment for RNA-seq datasets: (seRNA-seq, single file) (HeP et al., Nature, 2020)
cd /home/yaaroba/ENCODEproject/trimmedFastq/HeP_Nature_2020_RNAseq
samples=($(ls -f *.fq.gz))
index="/home/micgdu/GenomicData/genomicIndices/Mmusculus/STAR/STAR-2.7.6a_GRCm39_99bp/"
outFolder="/home/yaaroba/ENCODEproject/STARalingments"

for ((i=0;i<=${#samples[@]};i++)); do
   /home/yaaroba/scripts/STAR_SE_singleFile_3may.sh $index $outFolder "${samples[i]}" &
done

mv /home/yaaroba/ENCODEproject/STARalingments/*.bam /home/yaaroba/ENCODEproject/RNAseqBamFiles

mv /home/yaaroba/ENCODEproject/STARalingments/*Log.final.out /home/yaaroba/ENCODEproject/STARFinalOut