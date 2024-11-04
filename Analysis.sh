#!/bin/bash

Raw_reads="/home/santhoshhegde/Bioinfo_Task/Raw_Reads"
Ref_genome="/home/santhoshhegde/Bioinfo_Task/Ref_Genome/GRCm39.primary_assembly.genome.fa"
Ref_gtf="/home/santhoshhegde/Bioinfo_Task/Ref_Genome/gencode.vM35.basic.annotation.gtf"

mkdir Out_Analysis

out_dir="/home/santhoshhegde/Bioinfo_Task/Out_Analysis"

##################################### Assessing data quality using FastQC  #################################################################
mkdir $out_dir/QC

/home/gnomeadmin/software/FastQC/fastqc -o $out_dir/"QC" $Raw_reads/*.fastq.gz -t 32

###########################################  Quality and Adapter Trimming using fastp  ###########################################################

mkdir $out_dir/Trimmed
while read p;
do
fastp -i $Raw_reads/$p"_R1.fastq.gz" -I $Raw_reads/$p"_R2.fastq.gz" -o $out_dir/"Trimmed"/$p"_R1_trimmed.fastq.gz" -O $out_dir/"Trimmed"/$p"_R2_trimmed.fastq.gz" -h $out_dir/$p".html" -p -R $out_dir/$p"_Report.txt" -w 32 -p
done <Sample_names.txt

######################################  Indexing Reference Genome  ###############################################################################
/home/santhoshhegde/Transcriptome_Demo/Refseq/Softwares/hisat2-2.2.1-Linux_x86_64/hisat2-2.2.1/hisat2-build $Ref_genome Mouse


############################################ Read Mapping to Reference genome using hisat2 #######################################################

mkdir $out_dir/Mapping
while read p;
do
/home/santhoshhegde/Transcriptome_Demo/Refseq/Softwares/hisat2-2.2.1-Linux_x86_64/hisat2-2.2.1/hisat2 -x Ref_Genome/Mouse -1 $out_dir/"Trimmed"/$p"_R1_trimmed.fastq.gz" -2 $out_dir/"Trimmed"/$p"_R2_trimmed.fastq.gz" -S $out_dir/Mapping/$p".sam" --summary-file $out_dir/Mapping/$p"_mapping_stats.txt" -p 32
done <Sample_names.txt

###################################### Read Mapping to Reference genome using hisat2 (with removing multi mapped reads) ############################
mkdir $out_dir/Mapping/Filtered_multimap_reads
while read p;
do
/home/santhoshhegde/Transcriptome_Demo/Refseq/Softwares/hisat2-2.2.1-Linux_x86_64/hisat2-2.2.1/hisat2 -x Ref_Genome/Mouse -1 $out_dir/Trimmed/$p"_R1_trimmed.fastq.gz" -2 $out_dir/Trimmed/$p"_R2_trimmed.fastq.gz" -S $out_dir/Mapping/Filtered_multimap_reads/$p".sam" --summary-file $out_dir/Mapping/Filtered_multimap_reads/$p"_mapping_stats.txt" -k 1 -p 46
done <Sample_names.txt


################################### Processing Alignments using Samtools (Version: 0.1.19-44428cd) #################################################
mkdir $out_dir/Bam_files
while read p;
do
samtools view -Sb -@ 32 -o $out_dir/Bam_files/$p".bam" $out_dir/Mapping/Filtered_multimap_reads/$p".sam"
samtools sort -@ 32 $out_dir/Bam_files/$p".bam" $out_dir/Bam_files/$p"_sorted"
samtools index $out_dir/Bam_files/$p"_sorted.bam"
done <Sample_names.txt

################################### Read Counts using Featurecounts #################################################
featureCounts -T 3 -p -g gene_name -a Ref_Genome/gencode.vM35.basic.annotation.gtf -o Counts Bam_files/*.bam

################################### Differential Expression Analysis using DESeq2 #################################################
mv Counts_gene_names Count_matrix.txt
echo -e "Sample_Name\ttissue\ttime\nHeart_ZT0_1\tHeart\t0_hrs\nHeart_ZT0_2\tHeart\t0_hrs\nHeart_ZT12_1\tHeart\t12_hrs\nHeart_ZT12_2\tHeart\t12_hrs\nLiver_ZT0_1\tLiver\t0_hrs\nLiver_ZT0_2\tLiver\t0_hrs\nLiver_ZT12_1\tLiver\t12_hrs\nLiver_ZT12_2\tLiver\t12_hrs" > Sample_Info.txt

Rscript DESEq2.r


