#!/bin/sh

bedtools intersect -a path/to/RUN/10_sub-nucMACC/unStable_subNucs.bed -b path/to/project_foldert/public_data/GSM3274626_GRO-seq/GRO-seq_minusStrand_peaks.bed -wa -u > path/to/RUN/10_sub-nucMACC/unStable_subNucs_overlap-M1BP-minus.bed

bedtools intersect -a path/to/RUN/10_sub-nucMACC/unStable_subNucs.bed -b path/to/project_foldert/public_data/GSM3274626_GRO-seq/GRO-seq_plusStrand_peaks.bed -wa -u > path/to/RUN/10_sub-nucMACC/unStable_subNucs_overlap-M1BP-plus.bed


#Combine replicates (take max score)
macs2 cmbreps -i path/to/project_foldert/public_data/GSM3274626_GRO-seq/GSM3274626_SRX026303.flybase.minus.bedGraph path/to/project_foldert/public_data/GSM3274626_GRO-seq/GSM3274631_SRX026391.flybase.minus.bedGraph -m max -o path/to/project_foldert/public_data/GSM3274626_GRO-seq/GRO-seq_minusStrand.bedGraph

macs2 cmbreps -i path/to/project_foldert/public_data/GSM3274626_GRO-seq/GSM3274626_SRX026303.flybase.plus.bedGraph path/to/project_foldert/public_data/GSM3274626_GRO-seq/GSM3274631_SRX026391.flybase.plus.bedGraph -m max -o path/to/project_foldert/public_data/GSM3274626_GRO-seq/GRO-seq_plusStrand.bedGraph

#Call peaks on merged files
macs2 bdgpeakcall -i path/to/project_foldert/public_data/GSM3274626_GRO-seq/GRO-seq_plusStrand.bedGraph --outdir path/to/project_foldert/public_data/GSM3274626_GRO-seq/ -o GRO-seq_plusStrand_peaks.bed

macs2 bdgpeakcall -i path/to/project_foldert/public_data/GSM3274626_GRO-seq/GRO-seq_minusStrand.bedGraph --outdir path/to/project_foldert/public_data/GSM3274626_GRO-seq/ -o GRO-seq_minusStrand_peaks.bed


######### M1BP ########
#Call peaks on merged files
macs2 bdgpeakcall -i /Users/sara/data/D.Melanogaster_genome/Public_datasets/M1BP_ChIP-seq/GSM1208162_M1BP.bedGraph --outdir /Users/sara/data/D.Melanogaster_genome/Public_datasets/M1BP_ChIP-seq/ -o M1BP_peaks.bed

#Intersect with FN
bedtools intersect -a path/to/RUN/10_sub-nucMACC/unStable_subNucs.bed -b /Users/sara/data/D.Melanogaster_genome/Public_datasets/M1BP_ChIP-seq/M1BP_peaks.bed -wa -u > path/to/RUN/10_sub-nucMACC/unStable_subNucs_overlap-M1BP.bed

#Intersect with expressed genes containg FNs
bedtools intersect -a path/to/RUN/10_sub-nucMACC/Expressed_genes_TSS_with_FN.bed -b /Users/sara/data/D.Melanogaster_genome/Public_datasets/M1BP_ChIP-seq/M1BP_peaks.bed -wa -u > path/to/RUN/10_sub-nucMACC/expressed_unStable_subNucs_overlap-M1BP.bed

bedtools intersect -a path/to/RUN/10_sub-nucMACC/Expressed_genes_TSS_with_FN.bed -b /Users/sara/data/D.Melanogaster_genome/Public_datasets/M1BP_ChIP-seq/M1BP_peaks.bed -wa -u -f 0.2 > path/to/RUN/10_sub-nucMACC/expressed_unStable_subNucs_overlap-M1BP_stringent.bed

#Not overlapping
bedtools intersect -a path/to/RUN/10_sub-nucMACC/Expressed_genes_TSS_with_FN.bed -b /Users/sara/data/D.Melanogaster_genome/Public_datasets/M1BP_ChIP-seq/M1BP_peaks.bed -v > path/to/RUN/10_sub-nucMACC/expressed_unStable_subNucs_NotOverlap-M1BP.bed


#Intersect with expressed genes without FNs
bedtools intersect -a path/to/RUN/10_sub-nucMACC/Expressed_genes_TSS_withOut_FN.bed -b /Users/sara/data/D.Melanogaster_genome/Public_datasets/M1BP_ChIP-seq/M1BP_peaks.bed -wa -u > path/to/RUN/10_sub-nucMACC/expressed_without-unStable_subNucs_overlap-M1BP.bed

bedtools intersect -a path/to/RUN/10_sub-nucMACC/Expressed_genes_TSS_withOut_FN.bed -b /Users/sara/data/D.Melanogaster_genome/Public_datasets/M1BP_ChIP-seq/M1BP_peaks.bed -wa -u -f 0.2 > path/to/RUN/10_sub-nucMACC/expressed_without-unStable_subNucs_overlap-M1BP_stringent.bed

#Not overlapping
bedtools intersect -a path/to/RUN/10_sub-nucMACC/Expressed_genes_TSS_withOut_FN.bed -b /Users/sara/data/D.Melanogaster_genome/Public_datasets/M1BP_ChIP-seq/M1BP_peaks.bed -v > path/to/RUN/10_sub-nucMACC/expressed_without-unStable_subNucs_NotOverlap-M1BP.bed


######## final files
wc -l path/to/RUN/10_sub-nucMACC/expressed_unStable_subNucs_NotOverlap-M1BP.bed path/to/RUN/10_sub-nucMACC/expressed_unStable_subNucs_overlap-M1BP.bed path/to/RUN/10_sub-nucMACC/expressed_without-unStable_subNucs_NotOverlap-M1BP.bed path/to/RUN/10_sub-nucMACC/expressed_without-unStable_subNucs_overlap-M1BP.bed