#!/bin/sh

#  ./04_Special-nucs_motifs.sh
#  Created by Sara Wernig Zorc on 07.04.21.

project_folder=/Volumes/Sara_RAID_20TB/Sara_data_2018/MNase_project/MACC_project
fasta=/Users/sara/data/D.Melanogaster_genome/Drosophila_melanogaster/UCSC/Sequence/Fasta/dm3.fa
chromSizes=/Users/sara/data/D.Melanogaster_genome/Drosophila_melanogaster/dm3/dm3.genome.sizes
out_folder=$project_folder/results/Hypo-Hyper_StabilityNucs/H4_data

# --- BED --- #
#H3 data
#awk '{FS=OFS="\t"} NR >1 {print $1,$2,$3,$4}' $project_folder/results/DANPOS/ChIP-H4/len50-139/MNase_ChIP-H4_1.5U/pooled/nucMACC_filt5/nucMACC_scores.csv > $project_folder/results/DANPOS/ChIP-H4/nucMACC_characterization/all_subNucs.bed

#H4 data
#awk '{FS=OFS="\t"} NR >1 {print $1,$2,$3,$4}' $project_folder/results/DANPOS/ChIP-H4/len50-139/1.5U/pooled/nucMACC/sub-nucMACC_scores.csv > $project_folder/results/DANPOS/ChIP-H4/nucMACC_characterization/all_subNucs.bed


#  --- Input files ---  #

altNucs=/Users/sara/data/MACC_project/nextflow_out/H4_ChIP_v2.1_featureCOUNTS/RUN/10_sub-nucMACC/nonCanonical_subNucs.bed
fragNucs=/Users/sara/data/MACC_project/nextflow_out/H4_ChIP_v2.1_featureCOUNTS/RUN/10_sub-nucMACC/unStable_subNucs.bed

###############################  INPUT  ###########################################

high_exp=$project_folder/results/RNA-seq/genes/Expression-HIGH_genes.bed
mid_exp=$project_folder/results/RNA-seq/genes/Expression-MODERATE_genes.bed
low_exp=$project_folder/results/RNA-seq/genes/Expression-LOW_genes.bed

not_exp=$project_folder/results/RNA-seq/genes/Expression-SILENT_genes.bed
exp_exp=$project_folder/results/RNA-seq/genes/Expressed_genes.bed

exp_TSS=$project_folder/results/RNA-seq/readCount/Expressed_genes_TSS.bed
not_exp_TSS=$project_folder/results/RNA-seq/readCount/Silent_genes_TSS.bed

################# FN at TSS #################

awk '{FS=OFS="\t"}''function abs(x){return ((x < 0.0) ? -x : x)} {if (abs($16) < 100) print $1,$2,$3,$14,$7,$5}' /Users/sara/data/MACC_project/nextflow_out/H4_ChIP_v2.1_featureCOUNTS/RUN/10_sub-nucMACC/fragNucs_annotated.csv > /Users/sara/data/MACC_project/nextflow_out/H4_ChIP_v2.1_featureCOUNTS/RUN/10_sub-nucMACC/fragNucs_-100_TSS_+100.bed

awk '{FS=OFS="\t"}''function abs(x){return ((x < 0.0) ? -x : x)} {if (abs($16) < 100) print $1,$2,$3,$14,$7,$5}' /Users/sara/data/MACC_project/nextflow_out/H3_ChIP_v2_featureCOUNTS/RUN/10_sub-nucMACC/fragNucs_annotated.csv > /Users/sara/data/MACC_project/nextflow_out/H3_ChIP_v2_featureCOUNTS/RUN/10_sub-nucMACC/fragNucs_-100_TSS_+100.bed

############## Expressed genes TSS coordinates ################

awk '{FS=OFS="\t"}{b= $2 + 200}{print $1,$2,b,$4,$5 }' $project_folder/results/RNA-seq/genes/Expressed_genes.bed > $project_folder/results/RNA-seq/genes/Expressed_genes_TSS.bed

############## Expressed genes gene body coordinates ################

awk '{FS=OFS="\t"}{a = $2 + 201}{print $1,a,$3,$4,$5 }' $project_folder/results/RNA-seq/genes/Expressed_genes.bed > $project_folder/results/RNA-seq/genes/Expressed_genes_geneBody.bed
#snoRNAs are smaller than 200bp and get excluded in the analysis!!

############## FN at promoters ##############

awk '{FS=OFS="\t"} {if ($8=="Promoter") print $$1,$2,$3,$14,$7,$5}' /Users/sara/data/MACC_project/nextflow_out/H4_ChIP_v2.1_featureCOUNTS/RUN/10_sub-nucMACC/fragNucs_annotated.csv > /Users/sara/data/MACC_project/nextflow_out/H4_ChIP_v2.1_featureCOUNTS/RUN/10_sub-nucMACC/fragNucs_promoter.bed

awk '{FS=OFS="\t"} {if ($8=="Promoter") print $$1,$2,$3,$14,$7,$5}' /Users/sara/data/MACC_project/nextflow_out/H3_ChIP_v2_featureCOUNTS/RUN/10_sub-nucMACC/fragNucs_annotated.csv > /Users/sara/data/MACC_project/nextflow_out/H3_ChIP_v2_featureCOUNTS/RUN/10_sub-nucMACC/fragNucs_promoter.bed

############## H4 ##############
FN_promoters=/Users/sara/data/MACC_project/nextflow_out/H4_ChIP_v2.1_featureCOUNTS/RUN/10_sub-nucMACC/fragNucs_promoter.bed
FN_TSS=/Users/sara/data/MACC_project/nextflow_out/H4_ChIP_v2.1_featureCOUNTS/RUN/10_sub-nucMACC/fragNucs_-100_TSS_+100.bed

altNucs=/Users/sara/data/MACC_project/nextflow_out/H4_ChIP_v2.1_featureCOUNTS/RUN/10_sub-nucMACC/nonCanonical_subNucs.bed
fragNucs=/Users/sara/data/MACC_project/nextflow_out/H4_ChIP_v2.1_featureCOUNTS/RUN/10_sub-nucMACC/unStable_subNucs.bed

############## H3 ##############
FN_promoters=/Users/sara/data/MACC_project/nextflow_out/H3_ChIP_v2_featureCOUNTS/RUN/10_sub-nucMACC/fragNucs_promoter.bed
FN_TSS=/Users/sara/data/MACC_project/nextflow_out/H3_ChIP_v2_featureCOUNTS/RUN/10_sub-nucMACC/fragNucs_-100_TSS_+100.bed

altNucs=/Users/sara/data/MACC_project/nextflow_out/H3_ChIP_v2_featureCOUNTS/RUN/10_sub-nucMACC/nonCanonical_subNucs.bed
fragNucs=/Users/sara/data/MACC_project/nextflow_out/H3_ChIP_v2_featureCOUNTS/RUN/10_sub-nucMACC/unStable_subNucs.bed

##### Separate expressed genes by FN at TSS ######
awk -F'\t' 'NR==FNR{c[$4]++;next};!c[$4] > 0' $FN_TSS $not_exp_TSS > /Users/sara/data/MACC_project/nextflow_out/H4_ChIP_v2.1_featureCOUNTS/RUN/10_sub-nucMACC/Silent_genes_TSS_withOut_FN.bed
awk -F'\t' 'NR==FNR{c[$4]++;next};c[$4] > 0' $FN_TSS $not_exp_TSS > /Users/sara/data/MACC_project/nextflow_out/H4_ChIP_v2.1_featureCOUNTS/RUN/10_sub-nucMACC/Silent_genes_TSS_with_FN.bed

awk -F'\t' 'NR==FNR{c[$4]++;next};!c[$4] > 0' $FN_TSS $exp_TSS > /Users/sara/data/MACC_project/nextflow_out/H4_ChIP_v2.1_featureCOUNTS/RUN/10_sub-nucMACC/Expressed_genes_TSS_withOut_FN.bed
awk -F'\t' 'NR==FNR{c[$4]++;next};c[$4] > 0' $FN_TSS $exp_TSS > /Users/sara/data/MACC_project/nextflow_out/H4_ChIP_v2.1_featureCOUNTS/RUN/10_sub-nucMACC/Expressed_genes_TSS_with_FN.bed

FN_exp=/Users/sara/data/MACC_project/nextflow_out/H4_ChIP_v2.1_featureCOUNTS/RUN/10_sub-nucMACC/Expressed_genes_TSS_with_FN.bed
noFN_exp=/Users/sara/data/MACC_project/nextflow_out/H4_ChIP_v2.1_featureCOUNTS/RUN/10_sub-nucMACC/Expressed_genes_TSS_withOut_FN.bed

FN_silent=/Users/sara/data/MACC_project/nextflow_out/H4_ChIP_v2.1_featureCOUNTS/RUN/10_sub-nucMACC/Silent_genes_TSS_with_FN.bed
noFN_silent=/Users/sara/data/MACC_project/nextflow_out/H4_ChIP_v2.1_featureCOUNTS/RUN/10_sub-nucMACC/Silent_genes_TSS_withOut_FN.bed

wc -l $FN_silent $noFN_silent $FN_exp $noFN_exp $fragNucs $FN_promoters

# --- FASTA --- #
bedtools getfasta -fi $fasta -bed $FN_exp -fo ${FN_exp/%.bed/.fasta}
bedtools getfasta -fi $fasta -bed -bed $noFN_exp -fo ${noFN_exp/%.bed/.fasta}

#   --- Motifs ---    #

###Promoters with FN - Over promoters without FN (+/- 100bp from TSS)
findMotifs.pl /Users/sara/data/MACC_project/nextflow_out/H4_ChIP_v2.1_featureCOUNTS/RUN/10_sub-nucMACC/FN-TSS_promoterIDs.txt fly /Users/sara/data/MACC_project/nextflow_out/H4_ChIP_v2.1_featureCOUNTS/RUN/11_motifs/motifs/un.stable/FN_over_notFN/ -start -100 -end 100 -p 20 -bg /Users/sara/data/MACC_project/nextflow_out/H4_ChIP_v2.1_featureCOUNTS/RUN/10_sub-nucMACC/expressedGenes-withoutFN_promoterIDs.txt

###Expressed Promoters with FN - Over all expressed promoters
findMotifs.pl /Users/sara/data/MACC_project/nextflow_out/H4_ChIP_v2.1_featureCOUNTS/RUN/10_sub-nucMACC/FN-TSS_promoterIDs.txt fly /Users/sara/data/MACC_project/nextflow_out/H4_ChIP_v2.1_featureCOUNTS/RUN/11_motifs/motifs/un.stable/FN_overExpressed/ -start -100 -end 100 -p 20 -bg /Users/sara/data/MACC_project/nextflow_out/H4_ChIP_v2.1_featureCOUNTS/RUN/10_sub-nucMACC/allExpressedGenes_promoterIDs.txt

###Expressed Promoters without FN - Over all expressed promoters
findMotifs.pl /Users/sara/data/MACC_project/nextflow_out/H4_ChIP_v2.1_featureCOUNTS/RUN/10_sub-nucMACC/expressedGenes-withoutFN_promoterIDs.txt fly /Users/sara/data/MACC_project/nextflow_out/H4_ChIP_v2.1_featureCOUNTS/RUN/11_motifs/motifs/un.stable/notFN_overExpressed/ -start -100 -end 100 -p 20 -bg /Users/sara/data/MACC_project/nextflow_out/H4_ChIP_v2.1_featureCOUNTS/RUN/10_sub-nucMACC/allExpressedGenes_promoterIDs.txt


##########################################  TF enrichment (ChIP-seq)  ##############################################################
#M1BP TF
#Expressed genes w/ or w/o FN
computeMatrix scale-regions -S /Users/sara/data/D.Melanogaster_genome/Public_datasets/M1BP_ChIP-seq/GSM1208162_M1BP.bigwig -R $noFN_exp $FN_exp -a 2000 -b 2000 -o /Users/sara/data/MACC_project/nextflow_out/H4_ChIP_v2.1_featureCOUNTS/RUN/11_motifs/M1BP/M1BP_expGenes_FNgroups.mat.gz -p 20 --missingDataAsZero --startLabel "Start" --endLabel "End"

plotProfile -m /Users/sara/data/MACC_project/nextflow_out/H4_ChIP_v2.1_featureCOUNTS/RUN/11_motifs/M1BP/M1BP_expGenes_FNgroups.mat.gz -out /Users/sara/data/MACC_project/nextflow_out/H4_ChIP_v2.1_featureCOUNTS/RUN/11_motifs/M1BP/M1BP_expGenes_FNgroups_profile.pdf --dpi 400 --colors "#D2ac1e" "#2166AC" --regionsLabel "Expressed genes w/o FN" "Expressed genes w/ FN"

#RNA pol-II pausing - 1. GRO-seq
computeMatrix reference-point -S $project_folder/public_data/GSM3274626_GRO-seq/GSM32746*.bw -R $noFN_exp $FN_exp -a 1500 -b 1500 -o /Users/sara/data/MACC_project/nextflow_out/H4_ChIP_v2.1_featureCOUNTS/RUN/11_motifs/GRO-seq/GROseq_enrich_groups.mat.gz -p 20 --missingDataAsZero --referencePoint center

plotProfile -m /Users/sara/data/MACC_project/nextflow_out/H4_ChIP_v2.1_featureCOUNTS/RUN/11_motifs/GRO-seq/GROseq_enrich_groups.mat.gz --outFileName /Users/sara/data/MACC_project/nextflow_out/H4_ChIP_v2.1_featureCOUNTS/RUN/11_motifs/GRO-seq/GROseq_expGenes_FNgroups.pdf --dpi 400 --colors "#D2ac1e" "#2166AC" --regionsLabel "Expressed genes w/o FN" "Expressed genes w/ FN"

#Pausing index
 computeMatrix scale-regions -S $project_folder/public_data/GSM3274626_GRO-seq/GSM32746*.bw -R $project_folder/results/RNA-seq/genes/Expressed_genes_TSS.bed $project_folder/results/RNA-seq/genes/Expressed_genes_geneBody.bed -a 500 -b 2000 --outFileNameMatrix /Users/sara/data/MACC_project/nextflow_out/H4_ChIP_v2.1_featureCOUNTS/RUN/11_motifs/GRO-seq/GRO-seq_pausing_values.tab --outFileSortedRegions /Users/sara/data/MACC_project/nextflow_out/H4_ChIP_v2.1_featureCOUNTS/RUN/11_motifs/GRO-seq/GRO-seq_pausing_sortedRegions.bed -o /Users/sara/data/MACC_project/nextflow_out/H4_ChIP_v2.1_featureCOUNTS/RUN/11_motifs/GRO-seq/GRO-seq_pausing.mat.gz -p 20 --missingDataAsZero --startLabel "Start" --endLabel "End" --sortUsing mean 2> log.out
 
plotHeatmap -m /Users/sara/data/MACC_project/nextflow_out/H4_ChIP_v2.1_featureCOUNTS/RUN/11_motifs/GRO-seq/GRO-seq_pausing.mat.gz --dpi 400 --outFileName /Users/sara/data/MACC_project/nextflow_out/H4_ChIP_v2.1_featureCOUNTS/RUN/11_motifs/GRO-seq/GRO-seq_pausing_heatmap.pdf
 
#RNA pol-II pausing - 2. PRO-seq
#Expressed genes w/ or w/o FN
computeMatrix reference-point -S $project_folder/public_data/GSM1032757_PRO-seq/GSE42117_RAW/*.bigwig -R $noFN_exp $FN_exp -a 1500 -b 1500 -o PROseq_expGenes_FNgroups.mat.gz -p 20 --missingDataAsZero --referencePoint center

plotProfile -m PROseq_expGenes_FNgroups.mat.gz --outFileName PROseq_expGenes_FNgroups.pdf --dpi 400

#RNA pol-II pausing - 2. RNA pol-II ChIP-seq
#Expressed genes w/ or w/o FN
computeMatrix reference-point -S /Users/sara/data/D.Melanogaster_genome/Public_datasets/modEncode_329_RNApol-II/signal_data_files/RNA-polymerase-II:Cell-Line=S2-DRSC#Developmental-Stage=Late-Embryonic-stage#Tissue=Embryo-derived-cell-line:ChIP-chip:Rep-1::Dmel_r5.32:modENCODE_329:repset.4162648.smoothedM.bigwig -R $noFN_exp $FN_exp -a 1500 -b 1500 -o RNApol-II/RNApol-II_expGenes_FNgroups.mat.gz -p 20 --missingDataAsZero --smartLabels --referencePoint center --samplesLabel RNApol-II

plotProfile -m RNApol-II/RNApol-II_expGenes_FNgroups.mat.gz --outFileName RNApol-II/RNApol-II_expGenes_FNgroups.pdf --dpi 400 --samplesLabel RNApol-II --colors "#D2ac1e" "#2166AC" --regionsLabel "Expressed genes w/o FN" "Expressed genes w/ FN"

################################################################################################
################################################################################################

#All DM3 promoters (https://epd.epfl.ch/drosophila/drosophila_database.php?db=drosophila)
#/Users/sara/data/D.Melanogaster_genome/Dm_EPDnew_promoters.bed

#GRO-seq
computeMatrix reference-point -S /Volumes/Sara_RAID_20TB/Sara_data_2018/MNase_project/GSM3274626_GRO-seq/GSM32746*.bw -R $project_folder/results/DANPOS/ChIP-H4/nucMACC_characterization/fragileNucs_annotated_80bp-TSS.bed $project_folder/results/DANPOS/ChIP-H4/nucMACC_characterization/fragileNucs_annotated_not-TSS.bed /Users/sara/data/D.Melanogaster_genome/Dm_EPDnew_promoters.bed -a 2000 -b 2000 -o /Volumes/Sara_RAID_20TB/Sara_data_2018/MNase_project/GSM3274626_GRO-seq/GROseq_enrich_promoters_refPoint.mat.gz -p 20 --missingDataAsZero --smartLabels

plotProfile -m /Volumes/Sara_RAID_20TB/Sara_data_2018/MNase_project/GSM3274626_GRO-seq/GROseq_enrich_promoters_refPoint.mat.gz --outFileName /Volumes/Sara_RAID_20TB/Sara_data_2018/MNase_project/GSM3274626_GRO-seq/GROseq_enrich_promoters_refPoint.pdf --dpi 400

bedtools intersect -a /Users/sara/data/D.Melanogaster_genome/Dm3_promoters.bed -b $project_folder/results/DANPOS/ChIP-H4/nucMACC_characterization/fragileNucs_80bp-TSS_promoters.bed


#########################################################################################
#########################################################################################
#########################################################################################


chromSizes="/Users/sara/data/D.Melanogaster_genome/Drosophila_melanogaster/BDGP5/Sequence/dm3.chrom.sizes"

wigtobigwig /Users/sara/data/D.Melanogaster_genome/Public_datasets/modEncode_296_H3K27ac/signal_data_files/H3K27ac:Cell-Line=S2-DRSC#Developmental-Stage=Late-Embryonic-stage#Tissue=Embryo-derived-cell-line:ChIP-chip:Rep-1::Dmel_r5.32:modENCODE_296:repset.3211339.smoothedM.wig $chromSizes /Users/sara/data/D.Melanogaster_genome/Public_datasets/modEncode_296_H3K27ac/signal_data_files/H3K27ac_S2_modEncode296.smooth.bigwig -clip


#H3K4me3
computeMatrix reference-point -S /Users/sara/data/D.Melanogaster_genome/Public_datasets/modEncode_305_H3K4me3/signal_data_files/H3K4me3:Developmental-Stage=Late-Embryonic-stage#Tissue=Embryo-derived-cell-line:ChIP-chip:Rep-1::Dmel_r5.32:modENCODE_305:repset.3211336.smoothedM.bigwig /Users/sara/data/D.Melanogaster_genome/Public_datasets/modEncode_3761_H3K4me3/signal_data_files/H3K4me3:Cell-Line=S2-DRSC#Developmental-Stage=Late-Embryonic-stage#Tissue=Embryo-derived-cell-line:ChIP-chip:Rep-1::Dmel_r5.32:modENCODE_3761:repset.15664090.smoothedM.bigwig -R $project_folder/results/DANPOS/ChIP-H4/nucMACC_characterization/fragileNucs_annotated_80bp-TSS.bed $project_folder/results/DANPOS/ChIP-H4/nucMACC_characterization/fragileNucs_annotated_promotersOnly.bed $project_folder/results/DANPOS/ChIP-H4/nucMACC_characterization/fragileNucs_annotated_not-TSS.bed  -a 1500 -b 1500 -o $project_folder/results/fragileNucs/H3K4me3/H3K4me3_FNgroups_refPoint.mat.gz -p 20 --missingDataAsZero --referencePoint center --samplesLabel H3K4me3_1 H3K4me3_2


plotProfile -m $project_folder/results/fragileNucs/H3K4me3/H3K4me3_FNgroups_refPoint.mat.gz --outFileName $project_folder/results/fragileNucs/H3K4me3/H3K4me3_FNgroups_refPoint.pdf --plotHeight 20 --plotWidth 30 --dpi 400 --regionsLabel TSS promoter other  --samplesLabel H3K4me3_1 H3K4me3_2

plotHeatmap -m $project_folder/results/fragileNucs/H3K4me3/H3K4me3_FNgroups_refPoint.mat.gz --outFileName $project_folder/results/fragileNucs/H3K4me3/H3K4me3_FNgroups_heatmap_refPoint.pdf --dpi 400 --regionsLabel TSS promoter other  --samplesLabel H3K4me3_1 H3K4me3_2

#scale regions
#computeMatrix scale-regions -S /Users/sara/data/D.Melanogaster_genome/Public_datasets/modEncode_305_H3K4me3/signal_data_files/H3K4me3:Developmental-Stage=Late-Embryonic-stage#Tissue=Embryo-derived-cell-line:ChIP-chip:Rep-1::Dmel_r5.32:modENCODE_305:repset.3211336.smoothedM.bigwig -R $project_folder/results/DANPOS/ChIP-H4/nucMACC_characterization/fragileNucs_annotated_80bp-TSS.bed $project_folder/results/DANPOS/ChIP-H4/nucMACC_characterization/fragileNucs_annotated_promotersOnly.bed $project_folder/results/DANPOS/ChIP-H4/nucMACC_characterization/fragileNucs_annotated_not-TSS.bed  -a 1500 -b 1500 -o $project_folder/results/fragileNucs/H3K4me3/H3K4me3_FNgroups_scaleRegions.mat.gz -p 20 --missingDataAsZero --smartLabels

#plotProfile -m $project_folder/results/fragileNucs/H3K4me3/H3K4me3_FNgroups_scaleRegions.mat.gz --outFileName $project_folder/results/fragileNucs/H3K4me3/H3K4me3_FNgroups_scaleRegions.pdf --plotHeight 20 --plotWidth 30 --dpi 400 --regionsLabel TSS promoter other  --samplesLabel H3K4me3

#H3K27ac
computeMatrix reference-point -S /Users/sara/data/D.Melanogaster_genome/Public_datasets/modEncode_296_H3K27ac/signal_data_files/H3K27ac_S2_modEncode296.smooth.bigwig -R $project_folder/results/DANPOS/ChIP-H4/nucMACC_characterization/fragileNucs_annotated_80bp-TSS.bed $project_folder/results/DANPOS/ChIP-H4/nucMACC_characterization/fragileNucs_annotated_promotersOnly.bed $project_folder/results/DANPOS/ChIP-H4/nucMACC_characterization/fragileNucs_annotated_not-TSS.bed  -a 1500 -b 1500 -o $project_folder/results/fragileNucs/H3K27ac/H3K27ac_FNgroups_refPoint.mat.gz -p 20 --missingDataAsZero --samplesLabel H3K27ac --referencePoint center

plotProfile -m $project_folder/results/fragileNucs/H3K27ac/H3K27ac_FNgroups_refPoint.mat.gz --outFileName $project_folder/results/fragileNucs/H3K27ac/H3K27ac_FNgroups_refPoint.pdf --plotHeight 20 --plotWidth 30 --dpi 400 --regionsLabel TSS promoter other --samplesLabel H3K27ac

plotHeatmap -m $project_folder/results/fragileNucs/H3K27ac/H3K27ac_FNgroups_refPoint.mat.gz --outFileName $project_folder/results/fragileNucs/H3K27ac/H3K27ac_FNgroups_heatmap_refPoint.pdf --dpi 400 --regionsLabel TSS promoter other --samplesLabel H3K27ac

computeMatrix reference-point -S /Users/sara/data/D.Melanogaster_genome/Public_datasets/modEncode_296_H3K27ac/signal_data_files/H3K27ac_S2_modEncode296.smooth.bigwig -R /Volumes/Sara_RAID_20TB/Sara_data_2018/MNase_project/GSM3274626_GRO-seq/FN_promoters.bed /Volumes/Sara_RAID_20TB/Sara_data_2018/MNase_project/GSM3274626_GRO-seq/other_promoters.bed -a 1500 -b 1500 -o $project_folder/results/fragileNucs/H3K27ac/H3K27ac_promoters_refPoint.mat.gz -p 20 --missingDataAsZero --samplesLabel H3K27ac --referencePoint center

plotProfile -m $project_folder/results/fragileNucs/H3K27ac/H3K27ac_promoters_refPoint.mat.gz --outFileName $project_folder/results/fragileNucs/H3K27ac/H3K27ac_promoters_refPoint.pdf --plotHeight 20 --plotWidth 30 --dpi 400 --regionsLabel "FN promoters" "other promoters" --samplesLabel H3K27ac

plotHeatmap -m $project_folder/results/fragileNucs/H3K27ac/H3K27ac_promoters_refPoint.mat.gz --outFileName $project_folder/results/fragileNucs/H3K27ac/H3K27ac_promoters_heatmap_refPoint.pdf --dpi 400 --regionsLabel "FN promoters" "other promoters" --samplesLabel H3K27ac

#H3K27me3
computeMatrix reference-point -S /Users/sara/data/D.Melanogaster_genome/Public_datasets/modEncode_298_H3K27me3/signal_data_files/H3K27me3:Cell-Line=S2-DRSC#Developmental-Stage=Late-Embryonic-stage#Tissue=Embryo-derived-cell-line:ChIP-chip:Rep-1::Dmel_r5.32:modENCODE_298:repset.4621697.smoothedM.bigwig -R $project_folder/results/DANPOS/ChIP-H4/nucMACC_characterization/fragileNucs_annotated_80bp-TSS.bed $project_folder/results/DANPOS/ChIP-H4/nucMACC_characterization/fragileNucs_annotated_promotersOnly.bed $project_folder/results/DANPOS/ChIP-H4/nucMACC_characterization/fragileNucs_annotated_not-TSS.bed -a 1500 -b 1500 -o $project_folder/results/fragileNucs/H3K27me3/H3K27me3_FNgroups_refPoint.mat.gz -p 20 --missingDataAsZero --smartLabels --referencePoint center --samplesLabel H3K27me3

plotProfile -m $project_folder/results/fragileNucs/H3K27me3/H3K27me3_FNgroups_refPoint.mat.gz --outFileName $project_folder/results/fragileNucs/H3K27me3/H3K27me3_FNgroups_refPoint.pdf --plotHeight 20 --plotWidth 30 --dpi 400 --regionsLabel TSS promoter other --samplesLabel H3K27me3

plotHeatmap -m $project_folder/results/fragileNucs/H3K27me3/H3K27me3_FNgroups_refPoint.mat.gz --outFileName $project_folder/results/fragileNucs/H3K27me3/H3K27me3_FNgroups_heatmap_refPoint.pdf --dpi 400 --regionsLabel TSS promoter other --samplesLabel H3K27me3

#H3K9me3
computeMatrix reference-point -S /Users/sara/data/D.Melanogaster_genome/Public_datasets/modEncode_313_H3K9me3/signal_data_files/H3K9me3:Cell-Line=S2-DRSC#Developmental-Stage=Late-Embryonic-stage#Tissue=Embryo-derived-cell-line:ChIP-chip:Rep-1::Dmel_r5.32:modENCODE_313:repset.3506194.smoothedM.bigwig -R $project_folder/results/DANPOS/ChIP-H4/nucMACC_characterization/fragileNucs_annotated_80bp-TSS.bed $project_folder/results/DANPOS/ChIP-H4/nucMACC_characterization/fragileNucs_annotated_promotersOnly.bed $project_folder/results/DANPOS/ChIP-H4/nucMACC_characterization/fragileNucs_annotated_not-TSS.bed -a 1500 -b 1500 -o $project_folder/results/fragileNucs/H3K9me3/H3K9me3_FNgroups_refPoint.mat.gz -p 20 --missingDataAsZero --smartLabels --referencePoint center --samplesLabel H3K9me3

plotProfile -m $project_folder/results/fragileNucs/H3K9me3/H3K9me3_FNgroups_refPoint.mat.gz --outFileName $project_folder/results/fragileNucs/H3K9me3/H3K9me3_FNgroups_refPoint.pdf --plotHeight 20 --plotWidth 30 --dpi 400 --regionsLabel TSS promoter other --samplesLabel H3K9me3

plotHeatmap -m $project_folder/results/fragileNucs/H3K9me3/H3K9me3_FNgroups_refPoint.mat.gz --outFileName $project_folder/results/fragileNucs/H3K9me3/H3K9me3_FNgroups_heatmap_refPoint.pdf --dpi 400 --regionsLabel TSS promoter other --samplesLabel H3K9me3

#HP1b
computeMatrix reference-point -S /Users/sara/data/D.Melanogaster_genome/Public_datasets/modEncode_941_HP1b/signal_data_files/HP1b:Cell-Line=S2-DRSC#Developmental-Stage=Late-Embryonic-stage#Tissue=Embryo-derived-cell-line:ChIP-chip:Rep-1::Dmel_r5.32:modENCODE_941:repset.5014217.smoothedM.bigwig -R $project_folder/results/DANPOS/ChIP-H4/nucMACC_characterization/fragileNucs_annotated_80bp-TSS.bed $project_folder/results/DANPOS/ChIP-H4/nucMACC_characterization/fragileNucs_annotated_promotersOnly.bed $project_folder/results/DANPOS/ChIP-H4/nucMACC_characterization/fragileNucs_annotated_not-TSS.bed -a 1500 -b 1500 -o $project_folder/results/fragileNucs/HP1b/HP1b_FNgroups_refPoint.mat.gz -p 20 --missingDataAsZero --smartLabels --referencePoint center --samplesLabel HP1b

plotProfile -m $project_folder/results/fragileNucs/HP1b/HP1b_FNgroups_refPoint.mat.gz --outFileName $project_folder/results/fragileNucs/HP1b/HP1b_FNgroups_refPoint.pdf --plotHeight 20 --plotWidth 30 --dpi 400 --regionsLabel TSS promoter other --samplesLabel HP1b

plotHeatmap -m $project_folder/results/fragileNucs/HP1b/HP1b_FNgroups_refPoint.mat.gz --outFileName $project_folder/results/fragileNucs/HP1b/HP1b_FNgroups_heatmap_refPoint.pdf --dpi 400 --regionsLabel TSS promoter other --samplesLabel HP1b

############# per nuc group ####################
computeMatrix reference-point -S /Users/sara/data/D.Melanogaster_genome/Public_datasets/modEncode_305_H3K4me3/signal_data_files/H3K4me3:Developmental-Stage=Late-Embryonic-stage#Tissue=Embryo-derived-cell-line:ChIP-chip:Rep-1::Dmel_r5.32:modENCODE_305:repset.3211336.smoothedM.bigwig /Users/sara/data/D.Melanogaster_genome/Public_datasets/modEncode_296_H3K27ac/signal_data_files/H3K27ac_S2_modEncode296.smooth.bigwig /Users/sara/data/D.Melanogaster_genome/Public_datasets/modEncode_298_H3K27me3/signal_data_files/H3K27me3:Cell-Line=S2-DRSC#Developmental-Stage=Late-Embryonic-stage#Tissue=Embryo-derived-cell-line:ChIP-chip:Rep-1::Dmel_r5.32:modENCODE_298:repset.4621697.smoothedM.bigwig /Users/sara/data/D.Melanogaster_genome/Public_datasets/modEncode_313_H3K9me3/signal_data_files/H3K9me3:Cell-Line=S2-DRSC#Developmental-Stage=Late-Embryonic-stage#Tissue=Embryo-derived-cell-line:ChIP-chip:Rep-1::Dmel_r5.32:modENCODE_313:repset.3506194.smoothedM.bigwig /Users/sara/data/D.Melanogaster_genome/Public_datasets/modEncode_941_HP1b/signal_data_files/HP1b:Cell-Line=S2-DRSC#Developmental-Stage=Late-Embryonic-stage#Tissue=Embryo-derived-cell-line:ChIP-chip:Rep-1::Dmel_r5.32:modENCODE_941:repset.5014217.smoothedM.bigwig -R $project_folder/results/DANPOS/ChIP-H4/nucMACC_characterization/fragileNucs_annotated_not-TSS.bed -a 1500 -b 1500 -o $project_folder/results/fragileNucs/otherFN_refPoint.mat.gz -p 20 --missingDataAsZero --smartLabels --referencePoint center

plotProfile -m $project_folder/results/fragileNucs/otherFN_refPoint.mat.gz --outFileName $project_folder/results/fragileNucs/otherFN_refPoint.pdf --plotHeight 20 --plotWidth 30 --dpi 400 --samplesLabel H3K4me3 H3K27ac H3K27me3 H3K9me3 HP1b --regionsLabel "other FN"

plotHeatmap -m $project_folder/results/fragileNucs/otherFN_refPoint.mat.gz -o $project_folder/results/fragileNucs/otherFN_heatmap_refPoint.pdf --dpi 400 --samplesLabel H3K4me3 H3K27ac H3K27me3 H3K9me3 HP1b --regionsLabel "other FN"

#TSS FN
computeMatrix reference-point -S /Users/sara/data/D.Melanogaster_genome/Public_datasets/modEncode_305_H3K4me3/signal_data_files/H3K4me3:Developmental-Stage=Late-Embryonic-stage#Tissue=Embryo-derived-cell-line:ChIP-chip:Rep-1::Dmel_r5.32:modENCODE_305:repset.3211336.smoothedM.bigwig /Users/sara/data/D.Melanogaster_genome/Public_datasets/modEncode_296_H3K27ac/signal_data_files/H3K27ac_S2_modEncode296.smooth.bigwig /Users/sara/data/D.Melanogaster_genome/Public_datasets/modEncode_298_H3K27me3/signal_data_files/H3K27me3:Cell-Line=S2-DRSC#Developmental-Stage=Late-Embryonic-stage#Tissue=Embryo-derived-cell-line:ChIP-chip:Rep-1::Dmel_r5.32:modENCODE_298:repset.4621697.smoothedM.bigwig /Users/sara/data/D.Melanogaster_genome/Public_datasets/modEncode_313_H3K9me3/signal_data_files/H3K9me3:Cell-Line=S2-DRSC#Developmental-Stage=Late-Embryonic-stage#Tissue=Embryo-derived-cell-line:ChIP-chip:Rep-1::Dmel_r5.32:modENCODE_313:repset.3506194.smoothedM.bigwig /Users/sara/data/D.Melanogaster_genome/Public_datasets/modEncode_941_HP1b/signal_data_files/HP1b:Cell-Line=S2-DRSC#Developmental-Stage=Late-Embryonic-stage#Tissue=Embryo-derived-cell-line:ChIP-chip:Rep-1::Dmel_r5.32:modENCODE_941:repset.5014217.smoothedM.bigwig -R $project_folder/results/DANPOS/ChIP-H4/nucMACC_characterization/fragileNucs_annotated_80bp-TSS.bed -a 1500 -b 1500 -o $project_folder/results/fragileNucs/TSS-FN_refPoint.mat.gz -p 20 --missingDataAsZero --smartLabels --referencePoint center

plotProfile -m $project_folder/results/fragileNucs/TSS-FN_refPoint.mat.gz --outFileName $project_folder/results/fragileNucs/TSS-FN_refPoint.pdf --plotHeight 20 --plotWidth 30 --dpi 400 --samplesLabel H3K4me3 H3K27ac H3K27me3 H3K9me3 HP1b --regionsLabel "TSS FN"

plotHeatmap -m $project_folder/results/fragileNucs/TSS-FN_refPoint.mat.gz -o $project_folder/results/fragileNucs/TSS-FN_heatmap_refPoint.pdf --dpi 400 --samplesLabel H3K4me3 H3K27ac H3K27me3 H3K9me3 HP1b --regionsLabel "TSS FN"

#Promoter TSS
computeMatrix reference-point -S /Users/sara/data/D.Melanogaster_genome/Public_datasets/modEncode_305_H3K4me3/signal_data_files/H3K4me3:Developmental-Stage=Late-Embryonic-stage#Tissue=Embryo-derived-cell-line:ChIP-chip:Rep-1::Dmel_r5.32:modENCODE_305:repset.3211336.smoothedM.bigwig /Users/sara/data/D.Melanogaster_genome/Public_datasets/modEncode_296_H3K27ac/signal_data_files/H3K27ac_S2_modEncode296.smooth.bigwig /Users/sara/data/D.Melanogaster_genome/Public_datasets/modEncode_298_H3K27me3/signal_data_files/H3K27me3:Cell-Line=S2-DRSC#Developmental-Stage=Late-Embryonic-stage#Tissue=Embryo-derived-cell-line:ChIP-chip:Rep-1::Dmel_r5.32:modENCODE_298:repset.4621697.smoothedM.bigwig /Users/sara/data/D.Melanogaster_genome/Public_datasets/modEncode_313_H3K9me3/signal_data_files/H3K9me3:Cell-Line=S2-DRSC#Developmental-Stage=Late-Embryonic-stage#Tissue=Embryo-derived-cell-line:ChIP-chip:Rep-1::Dmel_r5.32:modENCODE_313:repset.3506194.smoothedM.bigwig /Users/sara/data/D.Melanogaster_genome/Public_datasets/modEncode_941_HP1b/signal_data_files/HP1b:Cell-Line=S2-DRSC#Developmental-Stage=Late-Embryonic-stage#Tissue=Embryo-derived-cell-line:ChIP-chip:Rep-1::Dmel_r5.32:modENCODE_941:repset.5014217.smoothedM.bigwig -R $project_folder/results/DANPOS/ChIP-H4/nucMACC_characterization/fragileNucs_annotated_promotersOnly.bed -a 1500 -b 1500 -o $project_folder/results/fragileNucs/Promoter-FN_refPoint.mat.gz -p 20 --missingDataAsZero --smartLabels --referencePoint center

plotProfile -m $project_folder/results/fragileNucs/Promoter-FN_refPoint.mat.gz --outFileName $project_folder/results/fragileNucs/Promoter-FN_refPoint.pdf --plotHeight 20 --plotWidth 30 --dpi 400 --samplesLabel H3K4me3 H3K27ac H3K27me3 H3K9me3 HP1b --regionsLabel "Promoter FN"

plotHeatmap -m $project_folder/results/fragileNucs/Promoter-FN_refPoint.mat.gz -o $project_folder/results/fragileNucs/Promoter-FN_heatmap_refPoint.pdf --dpi 400 --samplesLabel H3K4me3 H3K27ac H3K27me3 H3K9me3 HP1b --regionsLabel "Promoter FN"

#CTCF
computeMatrix reference-point -S /Users/sara/data/D.Melanogaster_genome/Public_datasets/modEncode_2639_CTCF/signal_data_files/CTCF:Cell-Line=S2-DRSC#Developmental-Stage=Late-Embryonic-stage#Tissue=Embryo-derived-cell-line:ChIP-seq:Rep-1:input:Dmel_r5.32:modENCODE_2639.bigwig -R $project_folder/results/DANPOS/ChIP-H4/nucMACC_characterization/fragileNucs_annotated_80bp-TSS.bed $project_folder/results/DANPOS/ChIP-H4/nucMACC_characterization/fragileNucs_annotated_promotersOnly.bed $project_folder/results/DANPOS/ChIP-H4/nucMACC_characterization/fragileNucs_annotated_not-TSS.bed -a 1500 -b 1500 -o $project_folder/results/fragileNucs/CTCF/CTCF_FNgroups_refPoint.mat.gz -p 20 --missingDataAsZero --smartLabels --referencePoint center --samplesLabel CTCF

plotProfile -m $project_folder/results/fragileNucs/CTCF/CTCF_FNgroups_refPoint.mat.gz --outFileName $project_folder/results/fragileNucs/CTCF/CTCF_FNgroups_refPoint.pdf --plotHeight 20 --plotWidth 30 --dpi 400 --regionsLabel TSS promoter other --samplesLabel CTCF

plotHeatmap -m $project_folder/results/fragileNucs/CTCF/CTCF_FNgroups_refPoint.mat.gz --outFileName $project_folder/results/fragileNucs/CTCF/CTCF_FNgroups_heatmap_refPoint.pdf --dpi 400 --regionsLabel TSS promoter other --samplesLabel CTCF

#RNApol-II
computeMatrix reference-point -S /Users/sara/data/D.Melanogaster_genome/Public_datasets/modEncode_329_RNApol-II/signal_data_files/RNA-polymerase-II:Cell-Line=S2-DRSC#Developmental-Stage=Late-Embryonic-stage#Tissue=Embryo-derived-cell-line:ChIP-chip:Rep-1::Dmel_r5.32:modENCODE_329:repset.4162648.smoothedM.bigwig -R $FN $exp_TSS $not_exp_TSS -a 1500 -b 1500 -o $project_folder/results/fragileNucs/RNApol-II/RNApol-II_FNgroups_refPoint.mat.gz -p 20 --missingDataAsZero --smartLabels --referencePoint center --samplesLabel RNApol-II

plotProfile -m $project_folder/results/fragileNucs/RNApol-II/RNApol-II_FNgroups_refPoint.mat.gz --outFileName $project_folder/results/fragileNucs/RNApol-II/RNApol-II_FNgroups_refPoint.pdf --plotHeight 20 --plotWidth 30 --dpi 400 --regionsLabel TSS promoter other --samplesLabel RNApol-II

plotHeatmap -m $project_folder/results/fragileNucs/RNApol-II/RNApol-II_FNgroups_refPoint.mat.gz --outFileName $project_folder/results/fragileNucs/RNApol-II/RNApol-II_FNgroups_heatmap_refPoint.pdf --dpi 400 --regionsLabel TSS promoter other --samplesLabel RNApol-II

#RNApol-II on promoters
computeMatrix reference-point -S /Users/sara/data/D.Melanogaster_genome/Public_datasets/modEncode_329_RNApol-II/signal_data_files/RNA-polymerase-II:Cell-Line=S2-DRSC#Developmental-Stage=Late-Embryonic-stage#Tissue=Embryo-derived-cell-line:ChIP-chip:Rep-1::Dmel_r5.32:modENCODE_329:repset.4162648.smoothedM.bigwig -R /Volumes/Sara_RAID_20TB/Sara_data_2018/MNase_project/GSM3274626_GRO-seq/FN_promoters.bed /Volumes/Sara_RAID_20TB/Sara_data_2018/MNase_project/GSM3274626_GRO-seq/other_promoters.bed -a 1500 -b 1500 -o $project_folder/results/fragileNucs/RNApol-II/RNApol-II_promoters_refPoint.mat.gz -p 20 --missingDataAsZero --smartLabels --referencePoint center --samplesLabel RNApol-II

plotProfile -m $project_folder/results/fragileNucs/RNApol-II/RNApol-II_promoters_refPoint.mat.gz --outFileName $project_folder/results/fragileNucs/RNApol-II/RNApol-II_promoters_refPoint.pdf --plotHeight 20 --plotWidth 30 --dpi 400 --regionsLabel "FN promoters" "other promoters" --samplesLabel RNApol-II

plotHeatmap -m $project_folder/results/fragileNucs/RNApol-II/RNApol-II_promoters_refPoint.mat.gz --outFileName $project_folder/results/fragileNucs/RNApol-II/RNApol-II_promoters_heatmap_refPoint.pdf --dpi 400 --regionsLabel "FN promoters" "other promoters" --samplesLabel RNApol-II


#H3K4me1
computeMatrix reference-point -S /Users/sara/data/D.Melanogaster_genome/Public_datasets/modEncode_304_H3K4me1/signal_data_files/H3K4me1:Cell-Line=S2-DRSC#Developmental-Stage=Late-Embryonic-stage#Tissue=Embryo-derived-cell-line:ChIP-chip:Rep-1::Dmel_r5.32:modENCODE_304:repset.3211335.smoothedM.bigwig -R $project_folder/results/DANPOS/ChIP-H4/nucMACC_characterization/fragileNucs_annotated_80bp-TSS.bed $project_folder/results/DANPOS/ChIP-H4/nucMACC_characterization/fragileNucs_annotated_promotersOnly.bed $project_folder/results/DANPOS/ChIP-H4/nucMACC_characterization/fragileNucs_annotated_not-TSS.bed -a 1500 -b 1500 -o $project_folder/results/fragileNucs/H3K4me1/H3K4me1_FNgroups_refPoint.mat.gz -p 20 --missingDataAsZero --smartLabels --referencePoint center --samplesLabel H3K4me1

plotProfile -m $project_folder/results/fragileNucs/H3K4me1/H3K4me1_FNgroups_refPoint.mat.gz --outFileName $project_folder/results/fragileNucs/H3K4me1/H3K4me1_FNgroups_refPoint.pdf --plotHeight 20 --plotWidth 30 --dpi 400 --regionsLabel TSS promoter other --samplesLabel H3K4me1

plotHeatmap -m $project_folder/results/fragileNucs/H3K4me1/H3K4me1_FNgroups_refPoint.mat.gz --outFileName $project_folder/results/fragileNucs/H3K4me1/H3K4me1_FNgroups_heatmap_refPoint.pdf --dpi 400 --regionsLabel TSS promoter other --samplesLabel H3K4me1


#########################################################
######### MONO-NUCS: Hypo/Hyper accessible nucs #########
#########################################################

#1. 5 histone modifications
computeMatrix reference-point -S /Users/sara/data/D.Melanogaster_genome/Public_datasets/modEncode_304_H3K4me1/signal_data_files/H3K4me1:Cell-Line=S2-DRSC#Developmental-Stage=Late-Embryonic-stage#Tissue=Embryo-derived-cell-line:ChIP-chip:Rep-1::Dmel_r5.32:modENCODE_304:repset.3211335.smoothedM.bigwig /Users/sara/data/D.Melanogaster_genome/Public_datasets/modEncode_296_H3K27ac/signal_data_files/H3K27ac_S2_modEncode296.smooth.bigwig /Users/sara/data/D.Melanogaster_genome/Public_datasets/modEncode_3761_H3K4me3/signal_data_files/H3K4me3:Cell-Line=S2-DRSC#Developmental-Stage=Late-Embryonic-stage#Tissue=Embryo-derived-cell-line:ChIP-chip:Rep-1::Dmel_r5.32:modENCODE_3761:repset.15664090.smoothedM.bigwig /Users/sara/data/D.Melanogaster_genome/Public_datasets/modEncode_313_H3K9me3/signal_data_files/H3K9me3:Cell-Line=S2-DRSC#Developmental-Stage=Late-Embryonic-stage#Tissue=Embryo-derived-cell-line:ChIP-chip:Rep-1::Dmel_r5.32:modENCODE_313:repset.3506194.smoothedM.bigwig /Users/sara/data/D.Melanogaster_genome/Public_datasets/modEncode_298_H3K27me3/signal_data_files/H3K27me3:Cell-Line=S2-DRSC#Developmental-Stage=Late-Embryonic-stage#Tissue=Embryo-derived-cell-line:ChIP-chip:Rep-1::Dmel_r5.32:modENCODE_298:repset.4621697.smoothedM.bigwig -R $project_folder/results/DANPOS/ChIP-H4/nucMACC_characterization/highMACC_monoNucs.bed $project_folder/results/DANPOS/ChIP-H4/nucMACC_characterization/lowMACC_monoNucs.bed -a 1500 -b 1500 -o $project_folder/results/Hypo-Hyper_AccesNucs/5_histoneMod/5_histoneMod_NucGroups_refPoint.mat.gz -p 20 --missingDataAsZero --smartLabels --referencePoint center --samplesLabel H3K4me1 H3K27ac  H3K4me3 H3K9me3 H3K27me3

#5 histone modifications
plotHeatmap -m $project_folder/results/Hypo-Hyper_AccesNucs/5_histoneMod/5_histoneMod_NucGroups_refPoint.mat.gz --outFileName $project_folder/results/Hypo-Hyper_AccesNucs/5_histoneMod/5_histoneMod_NucGroups_refPoint.pdf --dpi 400 --regionsLabel Hyper Hypo --samplesLabel H3K4me1 H3K27ac H3K4me3 H3K9me3 H3K27me3

#H3K4me1 only
plotHeatmap -m $project_folder/results/Hypo-Hyper_AccesNucs/H3K4me1/H3K4me1_NucGroups_refPoint.mat.gz --outFileName $project_folder/results/Hypo-Hyper_AccesNucs/H3K4me1/H3K4me1_NucGroups_refPoint.pdf --dpi 400 --regionsLabel Hyper Hypo --samplesLabel H3K4me1

#HP1b
computeMatrix reference-point -S /Users/sara/data/D.Melanogaster_genome/Public_datasets/modEncode_941_HP1b/signal_data_files/HP1b:Cell-Line=S2-DRSC#Developmental-Stage=Late-Embryonic-stage#Tissue=Embryo-derived-cell-line:ChIP-chip:Rep-1::Dmel_r5.32:modENCODE_941:repset.5014217.smoothedM.bigwig -R $project_folder/results/DANPOS/ChIP-H4/nucMACC_characterization/fragileNucs_annotated_80bp-TSS.bed $project_folder/results/DANPOS/ChIP-H4/nucMACC_characterization/fragileNucs_annotated_promotersOnly.bed $project_folder/results/DANPOS/ChIP-H4/nucMACC_characterization/fragileNucs_annotated_not-TSS.bed -a 1500 -b 1500 -o $project_folder/results/fragileNucs/HP1b/HP1b_FNgroups_refPoint.mat.gz -p 20 --missingDataAsZero --smartLabels --referencePoint center --samplesLabel HP1b

plotHeatmap -m $project_folder/results/fragileNucs/HP1b/HP1b_FNgroups_refPoint.mat.gz --outFileName $project_folder/results/fragileNucs/HP1b/HP1b_FNgroups_heatmap_refPoint.pdf --dpi 400 --regionsLabel TSS promoter other --samplesLabel HP1b


#H3K4me2
computeMatrix reference-point -S /Users/sara/data/D.Melanogaster_genome/Public_datasets/H3K4me2/GSM2756642_s2_h3k4me2_input.bw /Users/sara/data/D.Melanogaster_genome/Public_datasets/H3K4me2/GSM2756643_s2_h3k4me2.bw /Users/sara/data/D.Melanogaster_genome/Public_datasets/H3K4me2/GSM2756644_s2_h3k4me2_rep2.bw -R $project_folder/results/DANPOS/ChIP-H4/nucMACC_characterization/fragileNucs_annotated_80bp-TSS.bed $project_folder/results/DANPOS/ChIP-H4/nucMACC_characterization/fragileNucs_annotated_promotersOnly.bed $project_folder/results/DANPOS/ChIP-H4/nucMACC_characterization/fragileNucs_annotated_not-TSS.bed -a 1500 -b 1500 -o $project_folder/results/fragileNucs/H3K4me2/H3K4me2_FNgroups_refPoint.mat.gz -p 20 --missingDataAsZero --referencePoint center --samplesLabel "Input" "H3K4me2 rep1" "H3K4me2 rep2"

plotHeatmap -m $project_folder/results/fragileNucs/H3K4me2/H3K4me2_FNgroups_refPoint.mat.gz --outFileName $project_folder/results/fragileNucs/H3K4me2/H3K4me2_FNgroups_heatmap_refPoint.pdf --dpi 400 --regionsLabel TSS promoter other --samplesLabel "Input" "H3K4me2 rep1" "H3K4me2 rep2"


 
