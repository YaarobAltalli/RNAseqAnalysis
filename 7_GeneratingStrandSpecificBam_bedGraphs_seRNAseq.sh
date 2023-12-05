##7. generating strand specific bam & bedGraphs files for ssRNAseq reads

cd /home/yaaroba/ENCODEproject/RNAseqBamFiles
samples=($(ls -f *.bam))

for ((i=0;i<=${#samples[@]};i++)); do
   /home/yaaroba/scripts/ssRNAbamStrandSplit19june2023.sh "${samples[i]}" &
done

mv *_Fwd*bam *_Rev*bam /home/yaaroba/ENCODEproject/strandedBamFiles
mv *bedGraph.gz /home/yaaroba/ENCODEproject/RNAseqBedGraphs

# note this will split generate "starnd-specific" RNAseq bam files but its only valid for Fan et al. data as Nayak et al. RNAseq is not strand-specific  