##10. generating normalized bigwig files

cd /home/yaaroba/ENCODEproject/RNAseqNormBedGraphs/
samples=($(ls -f *.bedGraph))

for ((i=0;i<=${#samples[@]};i++)); do
   /home/micgdu/workflows/RNAseq/scripts/normBedGrpaphToBigWigM.sh "${samples[i]}" mm39 &
done


mv /home/yaaroba/ENCODEproject/RNAseqNormBedGraphs/*.bw /home/yaaroba/ENCODEproject/RNAseqBigWigs


##11. preparing bigwig files for viewing in UCSC
# currently bigwig files are placed at private google cloud account; files in big_wigs/devEpiPapers/# Fan2018_Nayak2023_RNAseq1june2023
# 2. For displaying of the bigwig files in UCSC or other genomic browsers they need to be placed in web-accessible folder. I have set up private google cloud account with web accessible storage folder which works fine with UCSC and is relatively cheap. To # explore your data in the UCSC genomic browser:
# - deposit your new bigwigs in your subfolder within: https://uam-my.sharepoint.com/:f:/r/personal/micgdu_amu_edu_pl/Documents/# Lab_Organization_Files/DataAnalysis/bigWig_temp?csf=1&web=1&e=eqOD6r   
# - Let me know that your new files are there
# - I will transfer the files and send you back the web address which you could use for the UCSC or other browser, on most # occasions I should be able to do this within several hours
# 3.    For browsing genomic files you can also try the IGV browser which is simpler, have less display options and is not linked # with so many online databases than UCSC but works faster: https://software.broadinstitute.org/software/igv/download


### 12. constructing UCSC bigwig file track descriptors 
# example of url address
https://storage.googleapis.com/big_wigs/yaarob/Gorkin_chrLandscapes2019/RNA-seq/forebrain_E11.5_rep1_Rev_norm_Rev.bedGraphS.bw 

#UCSC tracks
# HeP_Nature_2020 not strand specific library prep => fwd=rev

#Forebrain
track type=bigWig name="FB_E11.5_1" description="forebrain_E11.5_rep1" bigDataUrl=https://storage.googleapis.com/big_wigs/yaarob/Gorkin_chrLandscapes2019/RNA-seq/forebrain_E11.5_rep1_Rev_norm_Rev.bedGraphS.bw  visibility=full color=0,0,0 autoScale=off alwaysZero=on gridDefault=on  graphType=bar windowingFunction=mean negateValues=on  viewLimits=0:3 maxHeightPixels=200:50:10

track type=bigWig name="FB_E11.5_2" description="forebrain_E11.5_rep2" bigDataUrl=https://storage.googleapis.com/big_wigs/yaarob/Gorkin_chrLandscapes2019/RNA-seq/forebrain_E11.5_rep2_Rev_norm_Rev.bedGraphS.bw  visibility=full color=0,0,0 autoScale=off alwaysZero=on gridDefault=on  graphType=bar windowingFunction=mean negateValues=on  viewLimits=0:3 maxHeightPixels=200:50:10

track type=bigWig name="FB_E11.5_3" description="forebrain_E11.5_rep3" bigDataUrl=https://storage.googleapis.com/big_wigs/yaarob/Gorkin_chrLandscapes2019/RNA-seq/forebrain_E11.5_rep3_Rev_norm_Rev.bedGraphS.bw  visibility=full color=0,0,0 autoScale=off alwaysZero=on gridDefault=on  graphType=bar windowingFunction=mean negateValues=on  viewLimits=0:3 maxHeightPixels=200:50:10

track type=bigWig name="FB_E11.5_4" description="forebrain_E11.5_rep4" bigDataUrl=https://storage.googleapis.com/big_wigs/yaarob/Gorkin_chrLandscapes2019/RNA-seq/forebrain_E11.5_rep4_Rev_norm_Rev.bedGraphS.bw  visibility=full color=0,0,0 autoScale=off alwaysZero=on gridDefault=on  graphType=bar windowingFunction=mean negateValues=on  viewLimits=0:3 maxHeightPixels=200:50:10


track type=bigWig name="FB_E14.5_1" description="forebrain_E14.5_rep1" bigDataUrl=https://storage.googleapis.com/big_wigs/yaarob/Gorkin_chrLandscapes2019/RNA-seq/forebrain_E14.5_rep1_Rev_norm_Rev.bedGraphS.bw  visibility=full color=0,0,0 autoScale=off alwaysZero=on gridDefault=on  graphType=bar windowingFunction=mean negateValues=on  viewLimits=0:3 maxHeightPixels=200:50:10

track type=bigWig name="FB_E14.5_2" description="forebrain_E14.5_rep2" bigDataUrl=https://storage.googleapis.com/big_wigs/yaarob/Gorkin_chrLandscapes2019/RNA-seq/forebrain_E14.5_rep2_Rev_norm_Rev.bedGraphS.bw  visibility=full color=0,0,0 autoScale=off alwaysZero=on gridDefault=on  graphType=bar windowingFunction=mean negateValues=on  viewLimits=0:3 maxHeightPixels=200:50:10

track type=bigWig name="FB_E14.5_3" description="forebrain_E14.5_rep3" bigDataUrl=https://storage.googleapis.com/big_wigs/yaarob/Gorkin_chrLandscapes2019/RNA-seq/forebrain_E14.5_rep3_Rev_norm_Rev.bedGraphS.bw  visibility=full color=0,0,0 autoScale=off alwaysZero=on gridDefault=on  graphType=bar windowingFunction=mean negateValues=on  viewLimits=0:3 maxHeightPixels=200:50:10

track type=bigWig name="FB_E14.5_4" description="forebrain_E14.5_rep4" bigDataUrl=https://storage.googleapis.com/big_wigs/yaarob/Gorkin_chrLandscapes2019/RNA-seq/forebrain_E14.5_rep4_Rev_norm_Rev.bedGraphS.bw  visibility=full color=0,0,0 autoScale=off alwaysZero=on gridDefault=on  graphType=bar windowingFunction=mean negateValues=on  viewLimits=0:3 maxHeightPixels=200:50:10


track type=bigWig name="FB_P0_1" description="forebrain_P0_rep1" bigDataUrl=https://storage.googleapis.com/big_wigs/yaarob/Gorkin_chrLandscapes2019/RNA-seq/forebrain_P0_rep1_Rev_norm_Rev.bedGraphS.bw  visibility=full color=0,0,0 autoScale=off alwaysZero=on gridDefault=on  graphType=bar windowingFunction=mean negateValues=on  viewLimits=0:3 maxHeightPixels=200:50:10

track type=bigWig name="FB_P0_2" description="forebrain_P0_rep2" bigDataUrl=https://storage.googleapis.com/big_wigs/yaarob/Gorkin_chrLandscapes2019/RNA-seq/forebrain_P0_rep2_Rev_norm_Rev.bedGraphS.bw  visibility=full color=0,0,0 autoScale=off alwaysZero=on gridDefault=on  graphType=bar windowingFunction=mean negateValues=on  viewLimits=0:3 maxHeightPixels=200:50:10

track type=bigWig name="FB_P0_3" description="forebrain_P0_rep3" bigDataUrl=https://storage.googleapis.com/big_wigs/yaarob/Gorkin_chrLandscapes2019/RNA-seq/forebrain_P0_rep3_Rev_norm_Rev.bedGraphS.bw  visibility=full color=0,0,0 autoScale=off alwaysZero=on gridDefault=on  graphType=bar windowingFunction=mean negateValues=on  viewLimits=0:3 maxHeightPixels=200:50:10

track type=bigWig name="FB_P0_4" description="forebrain_P0_rep4" bigDataUrl=https://storage.googleapis.com/big_wigs/yaarob/Gorkin_chrLandscapes2019/RNA-seq/forebrain_P0_rep4_Rev_norm_Rev.bedGraphS.bw  visibility=full color=0,0,0 autoScale=off alwaysZero=on gridDefault=on  graphType=bar windowingFunction=mean negateValues=on  viewLimits=0:3 maxHeightPixels=200:50:10

#Heart:
track type=bigWig name="HT_E11.5_1" description="Heart_E11.5_rep1" bigDataUrl=https://storage.googleapis.com/big_wigs/yaarob/Gorkin_chrLandscapes2019/RNA-seq/heart_E11.5_rep1_Rev_norm_Rev.bedGraphS.bw  visibility=full color=0,0,0 autoScale=off alwaysZero=on gridDefault=on  graphType=bar windowingFunction=mean negateValues=on  viewLimits=0:3 maxHeightPixels=200:50:10

track type=bigWig name="HT_E11.5_2" description="Heart_E11.5_rep2" bigDataUrl=https://storage.googleapis.com/big_wigs/yaarob/Gorkin_chrLandscapes2019/RNA-seq/heart_E11.5_rep2_Rev_norm_Rev.bedGraphS.bw  visibility=full color=0,0,0 autoScale=off alwaysZero=on gridDefault=on  graphType=bar windowingFunction=mean negateValues=on  viewLimits=0:3 maxHeightPixels=200:50:10

track type=bigWig name="HT_E11.5_3" description="Heart_E11.5_rep3" bigDataUrl=https://storage.googleapis.com/big_wigs/yaarob/Gorkin_chrLandscapes2019/RNA-seq/heart_E11.5_rep3_Rev_norm_Rev.bedGraphS.bw  visibility=full color=0,0,0 autoScale=off alwaysZero=on gridDefault=on  graphType=bar windowingFunction=mean negateValues=on  viewLimits=0:3 maxHeightPixels=200:50:10

track type=bigWig name="HT_E11.5_4" description="Heart_E11.5_rep4" bigDataUrl=https://storage.googleapis.com/big_wigs/yaarob/Gorkin_chrLandscapes2019/RNA-seq/heart_E11.5_rep4_Rev_norm_Rev.bedGraphS.bw  visibility=full color=0,0,0 autoScale=off alwaysZero=on gridDefault=on  graphType=bar windowingFunction=mean negateValues=on  viewLimits=0:3 maxHeightPixels=200:50:10


track type=bigWig name="HT_E14.5_1" description="Heart_E14.5_rep1" bigDataUrl=https://storage.googleapis.com/big_wigs/yaarob/Gorkin_chrLandscapes2019/RNA-seq/heart_E14.5_rep1_Rev_norm_Rev.bedGraphS.bw  visibility=full color=0,0,0 autoScale=off alwaysZero=on gridDefault=on  graphType=bar windowingFunction=mean negateValues=on  viewLimits=0:3 maxHeightPixels=200:50:10

track type=bigWig name="HT_E14.5_2" description="Heart_E14.5_rep2" bigDataUrl=https://storage.googleapis.com/big_wigs/yaarob/Gorkin_chrLandscapes2019/RNA-seq/heart_E14.5_rep2_Rev_norm_Rev.bedGraphS.bw  visibility=full color=0,0,0 autoScale=off alwaysZero=on gridDefault=on  graphType=bar windowingFunction=mean negateValues=on  viewLimits=0:3 maxHeightPixels=200:50:10

track type=bigWig name="HT_E14.5_3" description="Heart_E14.5_rep3" bigDataUrl=https://storage.googleapis.com/big_wigs/yaarob/Gorkin_chrLandscapes2019/RNA-seq/heart_E14.5_rep3_Rev_norm_Rev.bedGraphS.bw  visibility=full color=0,0,0 autoScale=off alwaysZero=on gridDefault=on  graphType=bar windowingFunction=mean negateValues=on  viewLimits=0:3 maxHeightPixels=200:50:10

track type=bigWig name="HT_E14.5_4" description="Heart_E14.5_rep4" bigDataUrl=https://storage.googleapis.com/big_wigs/yaarob/Gorkin_chrLandscapes2019/RNA-seq/heart_E14.5_rep4_Rev_norm_Rev.bedGraphS.bw  visibility=full color=0,0,0 autoScale=off alwaysZero=on gridDefault=on  graphType=bar windowingFunction=mean negateValues=on  viewLimits=0:3 maxHeightPixels=200:50:10


track type=bigWig name="HT_P0_1" description="Heart_P0_rep1" bigDataUrl=https://storage.googleapis.com/big_wigs/yaarob/Gorkin_chrLandscapes2019/RNA-seq/heart_P0_rep1_Rev_norm_Rev.bedGraphS.bw  visibility=full color=0,0,0 autoScale=off alwaysZero=on gridDefault=on  graphType=bar windowingFunction=mean negateValues=on  viewLimits=0:3 maxHeightPixels=200:50:10

track type=bigWig name="HT_P0_2" description="Heart_P0_rep2" bigDataUrl=https://storage.googleapis.com/big_wigs/yaarob/Gorkin_chrLandscapes2019/RNA-seq/heart_P0_rep2_Rev_norm_Rev.bedGraphS.bw  visibility=full color=0,0,0 autoScale=off alwaysZero=on gridDefault=on  graphType=bar windowingFunction=mean negateValues=on  viewLimits=0:3 maxHeightPixels=200:50:10

track type=bigWig name="HT_P0_3" description="Heart_P0_rep3" bigDataUrl=https://storage.googleapis.com/big_wigs/yaarob/Gorkin_chrLandscapes2019/RNA-seq/heart_P0_rep3_Rev_norm_Rev.bedGraphS.bw  visibility=full color=0,0,0 autoScale=off alwaysZero=on gridDefault=on  graphType=bar windowingFunction=mean negateValues=on  viewLimits=0:3 maxHeightPixels=200:50:10

track type=bigWig name="HT_P0_4" description="Heart_P0_rep4" bigDataUrl=https://storage.googleapis.com/big_wigs/yaarob/Gorkin_chrLandscapes2019/RNA-seq/heart_P0_rep4_Rev_norm_Rev.bedGraphS.bw  visibility=full color=0,0,0 autoScale=off alwaysZero=on gridDefault=on  graphType=bar windowingFunction=mean negateValues=on  viewLimits=0:3 maxHeightPixels=200:50:10


#Liver:
track type=bigWig name="LV_E11.5_1" description="Liver_E11.5_rep1" bigDataUrl=https://storage.googleapis.com/big_wigs/yaarob/Gorkin_chrLandscapes2019/RNA-seq/liver_E11.5_rep1_Rev_norm_Rev.bedGraphS.bw  visibility=full color=0,0,0 autoScale=off alwaysZero=on gridDefault=on  graphType=bar windowingFunction=mean negateValues=on  viewLimits=0:3 maxHeightPixels=200:50:10

track type=bigWig name="LV_E11.5_2" description="Liver_E11.5_rep2" bigDataUrl=https://storage.googleapis.com/big_wigs/yaarob/Gorkin_chrLandscapes2019/RNA-seq/liver_E11.5_rep2_Rev_norm_Rev.bedGraphS.bw  visibility=full color=0,0,0 autoScale=off alwaysZero=on gridDefault=on  graphType=bar windowingFunction=mean negateValues=on  viewLimits=0:3 maxHeightPixels=200:50:10

track type=bigWig name="LV_E11.5_3" description="Liver_E11.5_rep3" bigDataUrl=https://storage.googleapis.com/big_wigs/yaarob/Gorkin_chrLandscapes2019/RNA-seq/liver_E11.5_rep3_Rev_norm_Rev.bedGraphS.bw  visibility=full color=0,0,0 autoScale=off alwaysZero=on gridDefault=on  graphType=bar windowingFunction=mean negateValues=on  viewLimits=0:3 maxHeightPixels=200:50:10

track type=bigWig name="LV_E11.5_4" description="Liver_E11.5_rep4" bigDataUrl=https://storage.googleapis.com/big_wigs/yaarob/Gorkin_chrLandscapes2019/RNA-seq/liver_E11.5_rep4_Rev_norm_Rev.bedGraphS.bw  visibility=full color=0,0,0 autoScale=off alwaysZero=on gridDefault=on  graphType=bar windowingFunction=mean negateValues=on  viewLimits=0:3 maxHeightPixels=200:50:10


track type=bigWig name="LV_E14.5_1" description="Liver_E14.5_rep1" bigDataUrl=https://storage.googleapis.com/big_wigs/yaarob/Gorkin_chrLandscapes2019/RNA-seq/liver_E14.5_rep1_Rev_norm_Rev.bedGraphS.bw  visibility=full color=0,0,0 autoScale=off alwaysZero=on gridDefault=on  graphType=bar windowingFunction=mean negateValues=on  viewLimits=0:3 maxHeightPixels=200:50:10

track type=bigWig name="LV_E14.5_2" description="Liver_E14.5_rep2" bigDataUrl=https://storage.googleapis.com/big_wigs/yaarob/Gorkin_chrLandscapes2019/RNA-seq/liver_E14.5_rep2_Rev_norm_Rev.bedGraphS.bw  visibility=full color=0,0,0 autoScale=off alwaysZero=on gridDefault=on  graphType=bar windowingFunction=mean negateValues=on  viewLimits=0:3 maxHeightPixels=200:50:10

track type=bigWig name="LV_E14.5_3" description="Liver_E14.5_rep3" bigDataUrl=https://storage.googleapis.com/big_wigs/yaarob/Gorkin_chrLandscapes2019/RNA-seq/liver_E14.5_rep3_Rev_norm_Rev.bedGraphS.bw  visibility=full color=0,0,0 autoScale=off alwaysZero=on gridDefault=on  graphType=bar windowingFunction=mean negateValues=on  viewLimits=0:3 maxHeightPixels=200:50:10

track type=bigWig name="LV_E14.5_4" description="Liver_E14.5_rep4" bigDataUrl=https://storage.googleapis.com/big_wigs/yaarob/Gorkin_chrLandscapes2019/RNA-seq/liver_E14.5_rep4_Rev_norm_Rev.bedGraphS.bw  visibility=full color=0,0,0 autoScale=off alwaysZero=on gridDefault=on  graphType=bar windowingFunction=mean negateValues=on  viewLimits=0:3 maxHeightPixels=200:50:10


track type=bigWig name="LV_P0_1" description="Liver_P0_rep1" bigDataUrl=https://storage.googleapis.com/big_wigs/yaarob/Gorkin_chrLandscapes2019/RNA-seq/liver_P0_rep1_Rev_norm_Rev.bedGraphS.bw  visibility=full color=0,0,0 autoScale=off alwaysZero=on gridDefault=on  graphType=bar windowingFunction=mean negateValues=on  viewLimits=0:3 maxHeightPixels=200:50:10

track type=bigWig name="LV_P0_2" description="Liver_P0_rep2" bigDataUrl=https://storage.googleapis.com/big_wigs/yaarob/Gorkin_chrLandscapes2019/RNA-seq/liver_P0_rep2_Rev_norm_Rev.bedGraphS.bw  visibility=full color=0,0,0 autoScale=off alwaysZero=on gridDefault=on  graphType=bar windowingFunction=mean negateValues=on  viewLimits=0:3 maxHeightPixels=200:50:10

track type=bigWig name="LV_P0_3" description="Liver_P0_rep3" bigDataUrl=https://storage.googleapis.com/big_wigs/yaarob/Gorkin_chrLandscapes2019/RNA-seq/liver_P0_rep3_Rev_norm_Rev.bedGraphS.bw  visibility=full color=0,0,0 autoScale=off alwaysZero=on gridDefault=on  graphType=bar windowingFunction=mean negateValues=on  viewLimits=0:3 maxHeightPixels=200:50:10

track type=bigWig name="LV_P0_4" description="Liver_P0_rep4" bigDataUrl=https://storage.googleapis.com/big_wigs/yaarob/Gorkin_chrLandscapes2019/RNA-seq/liver_P0_rep4_Rev_norm_Rev.bedGraphS.bw  visibility=full color=0,0,0 autoScale=off alwaysZero=on gridDefault=on  graphType=bar windowingFunction=mean negateValues=on  viewLimits=0:3 maxHeightPixels=200:50:10