#!/bin/sh

#Run as: ./00_MNase_script.sh

#The nucMACC pipeline requires raw paired-end MNase-seq data (.fastq format) with at least two MNase titrations per sample
#Project related variables need to be set by the user, before running the script
#Replicates should be pooled before running the script or processed seperately

#OUTPUT COLORS
GREEN="\033[1;32m"
PURPLE="\033[1;35m"
BLUE="\033[1;34m"
RESET="\033[1;0m"

set -e

echo "$PURPLE--------- MNase-seq analysis ---------$RESET"

echo "$PURPLE--------- Set project related variables ---------$RESET"

fasta=/path/to/genome.fa
annotation=/path/to/genes.gtf
project_folder=/path/to/project_folder
samples="$project_folder/raw_data/*_1.fastq.gz"
genome=/path/to/Bowtie2Index/genome
blacklisted=/path/to/dm3-blacklist+Chromosomes.bed
chromSizes=/path/to/dm3.genome.sizes
DANPOS_genes=/path/to/DM3_RefSeq_DANPOS_genePredictions.csv
DANPOS_expressed=/path/to/DM3_RefSeq_genePredictions_DANPOS_ExpressedGenes.csv
DANPOS_silent=/path/to/DM3_RefSeq_genePredictions_DANPOS_SilentGenes2.csv
DANPOS=/path/to/danpos_v3.1/danpos.py
Trimmomatic=/Users/sara/tools/Trimmomatic-0.38
genomeSize="162367812"
genome_igv="dm3"

echo  "$PURPLE The project folder: $RESET" $project_folder
echo  "$PURPLE Annotation file: $RESET" $annotation
echo  "$PURPLE Genome: $RESET" $genome
echo  "$PURPLE Genome fasta file: $RESET" $fasta
echo  "$PURPLE Samples: $RESET" $(basename ${samples/%_1.fastq.gz/})

cd $project_folder

echo  "$PURPLE--------- START OF ANALYSIS ---------$RESET"

mkdir -p raw_data processed_data/trimmed processed_data/aligned logs scripts results/DANPOS results/nucMACC results/fastQC

echo  "$GREEN--------- Quality Control ---------$RESET"

fastqc -t 10 -o $project_folder/results/fastQC $project_folder/raw_data/*.fastq.gz
multiqc $project_folder/results/fastQC/ -o $project_folder/results/fastQC/

echo  "$BLUE--------- Adapter trimming and read pairing ---------$RESET"

cd $project_folder/raw_data/

for file in *_1.fastq.gz; do echo "Processing file: " $file;  java -jar $Trimmomatic/trimmomatic-0.38.jar PE -phred33 -threads 20 $file ${file/%_1.fastq.gz/_2.fastq.gz} ${file/%_1.fastq.gz/_paired_1.fastq.gz} ${file/%_1.fastq.gz/_unpaired_1.fastq.gz} ${file/%_1.fastq.gz/_paired_2.fastq.gz} ${file/%_1.fastq.gz/_unpaired_2.fastq.gz} ILLUMINACLIP:$Trimmomatic/adapters/NGS_contaminants.fa:2:30:10 SLIDINGWINDOW:10:25 TRAILING:3 LEADING:3 2>> $project_folder/logs/Trimmomatic_log.txt; done

mv $project_folder/raw_data/*paired* $project_folder/processed_data/trimmed/

echo  "$BLUE--------- Alignment to the refrence genome ---------$RESET"

cd $project_folder/processed_data/trimmed/

for file in *_paired_1.fastq.gz; do echo "Processing sample: " $file  " + " ${file/%_1.fastq.gz/_2.fastq.gz}; bowtie2 --threads 16 --local -t --very-sensitive-local --no-discordant -x $genome -1 $file -2 ${file/%_1.fastq.gz/_2.fastq.gz} 2>> $project_folder/logs/Bowtie2_alignment_stats.txt | samtools view -bS -f 2 -@ 8 - > $project_folder/processed_data/aligned/${file/%_1.fastq.gz/_PE_aligned.bam}; done

for file in $project_folder/processed_data/aligned/*_aligned.bam; do echo "Processing sample: " $file; samtools sort -@ 20 $file > ${file/.bam/.sorted.bam} ; done

for file in $project_folder/processed_data/aligned/*.sorted.bam; do echo "Processing sample: " $file; samtools index -b $file ;done

echo "$BLUE----- Mapping quality control -----$RESET"

cd $project_folder

mkdir -p $project_folder/results/qualimap/raw_data/

#Make a file with information on sample name and file name
cd $project_folder/processed_data/aligned/

gfind -name \*.sorted.bam -printf "%f\t%f\n" > $project_folder/processed_data/aligned/MNase_samples.txt

qualimap multi-bamqc --java-mem-size=16G -d $project_folder/processed_data/aligned/MNase_samples.txt -outdir $project_folder/results/qualimap/raw_data/ -outformat PDF -r

echo "$BLUE----- Filter low mapping quality reads -----$RESET"

cd $project_folder/processed_data/aligned

for file in *.sorted.bam; do echo $file; alignmentSieve -b $file -o ${file/.sorted.bam/_qFilt.bam} --minMappingQuality 30 --filterMetrics ${file/.sorted.bam/_qFiltLog.txt} -p 18; done

echo "$BLUE----- Check the quality again with qualiMap -----$RESET"

mkdir $project_folder/results/qualimap/qFilt

gfind -name \*_qFilt.bam -printf "%f\t%f\n" > $project_folder/processed_data/aligned/MNase_samples_qFilt.txt

qualimap multi-bamqc --java-mem-size=16G -d $project_folder/processed_data/aligned/MNase_samples_qFilt.txt -outdir $project_folder/results/qualimap/qFilt -outformat PDF -r

echo "$BLUE----- Filter only monoNucs -----$RESET"

#Based on the insert size plot (qualimap) --> the mononuc sizes are XXX - XXX nt (normally 140-200 nt)

cd  $project_folder/processed_data/aligned/pooled/

for file in *_paired_PE_aligned.sorted.bam; do echo "Processing file: " $file;  alignmentSieve -b $file -o ${file/.sorted.bam/_len140-200.bam} --minMappingQuality 30 --filterMetrics ${file/.sorted.bam/_sizeFiltLog.txt} -p 18  --minFragmentLength 140 --maxFragmentLength 200 --blackListFileName $blacklisted; done

for file in *_aligned_len140-200.bam; do samtools sort -@ 20 $file > ${file/.bam/.sorted.bam} ; done

for file in *_aligned_len140-200.sorted.bam; do samtools index -b $file ;done

mkdir -p $project_folder/processed_data/aligned/pooled/len140-200/

mv $project_folder/processed_data/aligned/pooled/*len140-200* $project_folder/processed_data/aligned/pooled/len140-200/

echo "$BLUE----- Filter only subNucs -----$RESET"

#Based on the insert size plot (qualimap) --> the subnuc sizes are XXX - XXX nt (normally 50-139 nt)

for file in *_paired_PE_aligned.sorted.bam; do echo "Processing file: " $file;  alignmentSieve -b $file -o ${file/.sorted.bam/_len50-139.bam} --minMappingQuality 30 --filterMetrics ${file/.sorted.bam/_sizeFiltLog.txt} -p 18  --minFragmentLength 50 --maxFragmentLength 139 --blackListFileName $blacklist; done

for file in *_aligned_len50-139.bam; do samtools sort -@ 20 $file > ${file/.bam/.sorted.bam} ; done

for file in *_aligned_len50-139.sorted.bam; do samtools index -b $file ;done

mkdir -p $project_folder/processed_data/aligned/pooled/len50-139/

mv $project_folder/processed_data/aligned/pooled/*len50-139* $project_folder/processed_data/aligned/pooled/len50-139/


echo "$BLUE----- Pooled all mono-nucleosome titrations -----$RESET"

samtools merge -@ 20 $project_folder/processed_data/aligned/pooled/len140-200/MNase_ChIP-H4_pooled_paired_PE_aligned_len140-200.bam $project_folder/processed_data/aligned/pooled/len140-200/MNase_ChIP-H4_1.5U_paired_PE_aligned_len140-200.sorted.bam $project_folder/processed_data/aligned/pooled/len140-200/MNase_ChIP-H4_6.25U_paired_PE_aligned_len140-200.sorted.bam $project_folder/processed_data/aligned/pooled/len140-200/MNase_ChIP-H4_25U_paired_PE_aligned_len140-200.sorted.bam $project_folder/processed_data/aligned/pooled/len140-200/MNase_ChIP-H4_100U_paired_PE_aligned_len140-200.sorted.bam

samtools index -b $project_folder/processed_data/aligned/pooled/len140-200/*pooled*

echo "$BLUE----- Call mono-nucleosomes -----$RESET"

cd $project_folder/processed_data/aligned/pooled/len140-200/

for file in *_paired_PE_aligned_len140-200.sorted.bam; do mkdir $project_folder/results/DANPOS/len140-200/${file/_paired_PE_aligned_len140-200.sorted.bam/}; done

for file in *_paired_PE_aligned_len140-200.sorted.bam; do echo $file; python3 $DANPOS dpos $file -m 1 --extend 70 -c $genomeSize -u 0 -z 1 -a 1 -e 1 -o $project_folder/results/DANPOS/len140-200/${file/_paired_PE_aligned_len140-200.sorted.bam/} > $project_folder/results/DANPOS/len140-200/${file/_paired_PE_aligned_len140-200.sorted.bam/}/DANPOS_stats.txt ; done

python3 $DANPOS dpos *pooled* -m 1 --extend 70 -c $genomeSize -u 0 -z 1 -a 1 -e 1 -o $project_folder/results/DANPOS/len140-200/ >$project_folder/results/DANPOS/len140-200/DANPOS_stats.txt

echo "$BLUE----- Call sub-nucleosomes -----$RESET"

mkdir -p $project_folder/results/DANPOS/len50-139

cd $project_folder/processed_data/aligned/pooled/

for file in *_paired_PE_aligned_len50-139.sorted.bam; do mkdir $project_folder/results/DANPOS/len50-139/${file/_paired_PE_aligned_len50-139.sorted.bam/}; done

for file in *_paired_PE_aligned_len50-139.sorted.bam; do echo $file; python3 $DANPOS dpos $file -m 1 --extend 70 -c $genomeSize -u 0 -z 70 -e 1 -o $project_folder/results/DANPOS/len50-139/${file/_paired_PE_aligned_len50-139.sorted.bam/}/ > $project_folder/results/DANPOS/len50-139/${file/_paired_PE_aligned_len50-139.sorted.bam/}/DANPOS_stats.txt ; done

echo "$BLUE----- mono-nuc profile over TSS of all genes (DANPOS) -----$RESET"

mkdir $project_folder/results/DANPOS/len140-200/profile/

cd $project_folder/results/DANPOS/len140-200/profile/

wig_files=$(ls $project_folder/results/DANPOS/len140-200/*U/pooled/*.wig | paste -s -d, -)

python3 $DANPOS profile --genefile_paths $DANPOS_expressed,$DANPOS_silent --genomic_sites TSS $wig_files

echo "$BLUE----- mono-nuc profile over TSS of all genes (DeepTools) -----$RESET"
wig_files=$(ls /Users/sara/data/MACC_project/nextflow_out/H4-ChIP_pre/05_MONO-NUCS_PROFILE/*bw | paste -s -)

computeMatrix reference-point -S $wig_files -R ${DANPOS_expressed/.csv/.bed} ${DANPOS_silent/.csv/.bed} -o monoNuc_matrix.gz --referencePoint TSS -b 1000 -a 1000 -p 20 --smartLabels

plotProfile --dpi 400 --matrixFile monoNuc_matrix.gz -o monoNuc_profile.pdf --plotType lines --perGroup --regionsLabel "Expressed genes" "Silent genes" --samplesLabel "1.5U" "100U" "25U" "6.25U" --legendLocation upper-right --plotHeigh 10 --plotWidth 10 --yAxisLabel "Relative distance" --plotTitle "Mono-nucleosomes TSS profile" --colors "#00d455" "#ff00cc" "#8800aa" "#00ccff"

#00d455 green
#00ccff blue
#8800aa purple
#ff00cc pink

echo "$BLUE----- sub-nuc profile over TSS of all genes (DANPOS) -----$RESET"

mkdir $project_folder/results/DANPOS/len50-139/profile/

cd $project_folder/results/DANPOS/len50-139/profile/

wig_files=$(ls $project_folder/results/DANPOS/len50-139/*U/pooled/*.wig | paste -s -d, -)

python3 $DANPOS profile --genefile_paths $DANPOS_expressed,$DANPOS_silent --genomic_sites TSS $wig_files

echo "$BLUE----- sub-nuc profile over TSS of all genes (DeepTools) -----$RESET"
wig_files=$(ls /Users/sara/data/MACC_project/nextflow_out/H4-ChIP_pre/06_SUB-NUCS_PROFILE/*bw | paste -s -)

computeMatrix reference-point -S $wig_files -R ${DANPOS_expressed/.csv/.bed} ${DANPOS_silent/.csv/.bed} -o subNuc_matrix.gz --referencePoint TSS -b 1000 -a 1000 -p 20 --smartLabels

plotProfile --dpi 400 --matrixFile subNuc_matrix.gz -o subNuc_profile.pdf --plotType lines --perGroup --regionsLabel ExpressedGenes SilentGenes --samplesLabel "1.5U" "100U" "25U" "6.25U" --legendLocation upper-right --plotHeigh 10 --plotWidth 10 --yAxisLabel "Relative distance" --plotTitle "Sub-nucleosome TSS profile" --colors "#00d455" "#ff00cc" "#8800aa" "#00ccff"

#00d455 green
#00ccff blue
#8800aa purple
#ff00cc pink

echo "$BLUE----- wig to bigWig conversion -----$RESET"

for file in $project_folder/results/DANPOS/len*/pooled/*.wig; do echo $file; wigToBigWig -clip $file $chromSizes ${file/%.wig/.bigwig}; done

for file in $project_folder/results/DANPOS/len*/pooled/*.bigwig; do echo ${file/_paired_PE_aligned_len50-139.sorted.Fnor.smooth.bw/.bigwig}; mv $file ${file/_paired_PE_aligned_len50-139.sorted.Fnor.smooth.bw/_sub-nucleosomes.bigwig};done

#for file in $project_folder/results/DANPOS/len*/pooled/*.wig; do echo $file; igvtools toTDF $file ${file/%.wig/.tdf} $genome_igv ; done

echo "$BLUE----- Mono-nucs: Count reads per nucleosome -----$RESET"

#1. xls to SAF file
echo -e "GeneID\tChr\tStart\tEnd\tStrand\tFuzzinessScore" > MNase_len140-200_nucPositions.saf

awk '{FS=OFS="\t"} NR > 1  {printf("nuc%06d%s\n",NR,"\t"$1"\t"$2"\t"$3"\t"".""\t"$5)}' $project_folder/results/DANPOS/len140-200/pooled/*.positions.xls >> $project_folder/results/DANPOS/len140-200/pooled/MNase_len140-200_nucPositions.saf

awk '{FS=OFS="\t"} NR > 1 {printf($1"\t"$2"\t"$3"\t""nuc%06d%s\n",NR,"\t"$5"\t"".")}' $project_folder/results/DANPOS/len140-200/pooled/*.positions.xls >> $project_folder/results/DANPOS/len140-200/pooled/MNase_len140-200_nucPositions.bed
 
gsed 's/-//g' -i $project_folder/results/DANPOS/len140-200/pooled/MNase_len140-200_nucPositions.saf > $project_folder/results/DANPOS/len140-200/pooled/MNase_len140-200_nucPositions.saf
 
#2. read count featureCount

featureCounts -F SAF -a $project_folder/results/DANPOS/len140-200/pooled/MNase_len140-200_nucPositions.saf -o $project_folder/results/DANPOS/len140-200/pooled/MNase_len140-200_readCounts.csv -T 20 -p -B -C --largestOverlap $project_folder/results/DANPOS/len140-200/processed_data/aligned/*1min*.sorted.bam  $project_folder/results/DANPOS/len140-200/processed_data/aligned/*3min*.sorted.bam $project_folder/results/DANPOS/len140-200/processed_data/aligned/*9min*.sorted.bam $project_folder/results/DANPOS/len140-200/processed_data/aligned/*27min*.sorted.bam

#3. add GC% to read count

awk '{FS=OFS="\t"} NR > 2 {print $2,$3,$4,$1,$6,$5,$7,$8,$9,$10}' $project_folder/results/DANPOS/len140-200/pooled/MNase_len140-200_readCounts.csv | bedtools nuc -fi $fasta -bed - | cut -f 1-10,12 > $project_folder/results/DANPOS/len140-200/pooled/MNase_len140-200_readCounts_wGC.csv

echo "$BLUE----- Sub-nucs: Count reads per nucleosome -----$RESET"

#1. xls to SAF file

echo -e "GeneID\tChr\tStart\tEnd\tStrand\tFuzzinessScore" > $project_folder/results/DANPOS/len50-139/pooled/MNase_len50-139_nucPositions.saf

awk '{FS=OFS="\t"} NR > 1 {printf("nuc%06d%s\n",NR,"\t"$1"\t"$2"\t"$3"\t"".""\t"$5)}' $project_folder/results/DANPOS/len50-139/pooled/*.positions.xls >> $project_folder/results/DANPOS/len50-139/pooled/MNase_len50-139_nucPositions.saf

 awk '{FS=OFS="\t"} NR > 1 {printf($1"\t"$2"\t"$3"\t""nuc%06d%s\n",NR,"\t"$5"\t"".")}' $project_folder/results/DANPOS/len50-139/pooled/*.positions.xls >> $project_folder/results/DANPOS/len50-139/pooled/MNase_len50-139_nucPositions.bed

#2. read count featureCount

featureCounts -F SAF -a $project_folder/results/DANPOS/len50-139/pooled/MNase_len50-139_nucPositions.saf -o $project_folder/results/DANPOS/len50-139/pooled/MNase_len50-139_readCounts.csv -T 20 -p -B -C --largestOverlap $project_folder/processed_data/aligned/len50-139/*1min*.sorted.bam  $project_folder/processed_data/aligned/len50-139/*3min*.sorted.bam $project_folder/processed_data/aligned/len50-139/*9min*.sorted.bam $project_folder/processed_data/aligned/len50-139/*27min*.sorted.bam

#3. add GC% to read count
awk '{FS=OFS="\t"} NR > 2 {print $2,$3,$4,$1,$6,$5,$7,$8,$9,$10}' $project_folder/results/DANPOS/len50-139/pooled/MNase_len50-139_readCounts.csv | bedtools nuc -fi $fasta -bed - | cut -f 1-10,12 > $project_folder/results/DANPOS/len50-139/pooled/MNase_len50-139_readCounts_wGC.csv


echo "$BLUE----- Nucleosome accessability score calculation: monoNucs -----$RESET"

filepath=$(pwd)

Rscript $project_folder/scripts/01_nuc-accessability-score.R $project_folder/results/DANPOS/len140-200/pooled/MNase_len140-200_readCounts_wGC.csv $project_folder/results/DANPOS/len140-200/pooled/

echo "$BLUE----- Nucleosome stability score calculation: subNucs -----$RESET"

Rscript $project_folder/scripts/02_nuc-stability-score.R $project_folder/results/DANPOS/len50-139/pooled/MNase_len50-139_readCounts_wGC.csv $project_folder/results/DANPOS/len50-139/pooled/

echo "$BLUE----- Score cut-off -----$RESET"

for file in $project_folder/results/DANPOS/len*/pooled/nucMACC/nucMACC_scores.csv; do echo $file; awk '{FS=OFS="\t"} NR > 1 {print $1,$2,$3,$13,$4}' $file > ${file/%.csv/.bedgraph}; done

for file in $project_folder/results/DANPOS/len*/pooled/nucMACC/nucMACC_scores.csv; do echo $file; awk '{FS=OFS="\t"} NR >1 {print $1,$2,$3,$4,$13,"."}' $file > ${file/%.csv/.bed}; done

#OPTIONAL
#Rscript -e "rmarkdown::render('$project_folder/scripts/03_specialNucs-plots_H3.Rmd',params=list(args = myarg))"

echo  "$BLUE--------- DONE ---------$RESET"
