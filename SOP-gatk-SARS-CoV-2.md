#### Author  

AMBARISH KUMAR

er.ambarish@gmail.com

ambari73_sit@jnu.ac.in

----------------------------------------------------------------------------------------------------------------------------------------

This is standard operating protocol for genomic variant detection using GTK3. This will be effective and usefull for getting SARS-CoV-2 genome variants.

----------------------------------------------------------------------------------------------------------------------------------------

It uses paired-end illumina SARS-CoV-2 RNASEQ reads and reference genome.

Tool description

- Bowtie2

- SAMTools

- PICARD

- GATK3


----------------------------------------------------------------------------------------------------------------------------------------

Purpose is to generate CWL script for SOP for genomic variant detection using GATK3. SARS-CoV-2 virus genome will be a case study for it.

It uses illumina RNASEQ reads and genome sequence. 

--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Standard operating protocol for variant calling using GATK
Software installation and environment setup 
      - Access .jar file for each GATK and PICARD tool from their respective folder.
      - Export absolute path of bowtie executables using bash variable $PATH.
## Mapping to the reference genome
Mutated read data are mapped against reference genome using aligner - bowtie2.
Mapping steps can further be split into indexing and alignment steps.

### Indexing 
Indexing is compressing and assigning genomic position to genome sequence for fast and efficient traversal of genomic bases. Indexing is performed for GATK using the bowtie command bowtie2-build.
Command-line usage:

bowtie2-build reference_genome index_base_name

Where

bowtie2-build - Bowtie command to build index of reference genome.

reference_genome - fasta format reference genome.

index_base_name - Base name of generated index file.

Alignment 

Alignment is matching read bases with that of reference genome sequence. It is performed for GATK using bowtie aligner as data preparation step. 

Command-line usage:

bowtie2 -q -x index -1 left_read -2 right_read -S alignment_file

Where

bowtie2 - Bowtie command for alignment.

q - bowtie option for fastq input file format.

index - prefix for bowtie index.

1 - option for left-end reads. 

2 - option for right-end reads.

S - option for .sam format alignment file.

left_read - reads generated from forward strand.

right_read - reads generated from reverse strand.

alignment_file - alignment file in .sam format.
      
--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

###### Data pre-processing using PICARD tools

--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
 
## Add read groups, sort, mark duplicates, and create index

In this step, we merge all read groups with different  sequencing lane reads and read group header @RG and assign a single read group header.

Command-line usage:

java -jar  AddOrReplaceReadGroups.jar -I input_file -O output_file_name -RGID rgid_id -RGLB rglb_string -RGPL rgpl_string -RGSM rgsm _string -RGPU rgpu_string

Where

AddOrReplaceReadGroup.jar - Picard tool to merge all read groups and assign single read group ID.

input_file - Alignment file in .sam format.

output_file - Alignment file in .sam format with merged read group ID.

rglb_string   -   Read Group Library Required. 

rgpl_string   -   Read Group platform (e.g. illumina, solid) Required

rgsm_string   -   Read Group sample name Required

rgpu_string   -   Read Group platform unit Required

## Sort and Index

This step sorts and indexes the reads, and converts between SAM to BAM format. The indexed bam file is further coordinate sorted using the PICARD tool SortSam.jar

Command-line usage:

java -jar SortSam.jar I=input_file O=output_file CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT SORT_ORDER=coordinate

Where

Input_file -  Sorted bam file

Output_file  - The sorted BAM or SAM output file

## Mark the duplicate reads

Duplicate sequenced reads are marked and removed using PICARD tool MarkDuplicates.jar

Command-line usage:

java -jar MarkDuplicates.jar I=input_file O=output_file CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT M=output_metrics 

Where

input_file - sorted bam input file

output_file - markdup bam file

output_metrices - File to write duplication metrics to Required

## Split'N'Trim and reassign mapping qualities

Hardclip the hanging part of reads into intronic region and assign reads mapping quality of 60 replacing 255  which is understood and interpreted by GATK

Command-line usage:

java -jar GenomeAnalysisTK.jar -T SplitNCigarReads -R reference_sequence -I input_file -o out_put -rf ReassignOneMappingQuality -RMQF rmqf_value -RMQT rmsqt_value -U ALLOW_N_CIGAR_READS

Where

reference_sequence - reference sequence in fasta format.

input_file - markedup .bam file.

outputfile - split .bam output file.

rf - allows read filter and mapping quality. 

rmqf_value - RMQF value 255

rmqt_value - RMQT value 60

## Variant calling

HaplotypeCaller - calls all plausible haplotypes and detect variants. 

Command-line usage:

java -jar GenomeAnalysisTK.jar -T HaplotypeCaller -R reference_sequence  -I input_file -o output_file -dontUseSoftClippedBases -stand_call_conf stand_call_conf_value -stand_emit_conf stand_emit_conf_value

Where

input_file   -   split bam file.

output_file  -  vcf file as output.

stand_call_conf_value   -   20.0 confidence interval value

stand_emit_conf_value   -   20.0 confidence interval value

stand_call_conf   - The minimum phred-scaled confidence threshold at which variants        should be called

stand_call_emit   -   The minimum phred-scaled confidence threshold at which variants should be emitted (and filtered with LowQual if less than the calling threshold)

## Variant filtering

Filter out the low quality variants.

Command-line usage

java -jar GenomeAnalysisTK.jar -T VariantFiltration -R reference_genome -V  input file -window window_size -cluster cluster_size -filterName filter_name -filter "FS > 30.0" -filterName filter_name -filter "QD < 2.0" -o out_put

Where

reference_genome   -  reference sequence in fasta file format.

input_file  -   .vcf file i.e variant set to used as an input.

window_size   -   window size in which to evaluate clustered snp’s.

cluster_size   -   the no of snp’s which makeup a cluster.

filter_name   -   name of the filter.

out_put  -   filtered variants or VCF’s.
--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

### Command-line implementation 

#### Input dataset

sars-cov-2.fasta - sars-cov-2 reference genome.
sars-cov-2-reads_1.fq - left-end RNASEQ reads.
sars-cov-2-reads_2.fq - right-end RNASEQ reads.

#### Preparation of alignment files

export PATH=$PATH:/absolute path to bowtie executables/

bowtie2-build sars-cov-2.fasta sars-cov-2

bowtie2 -q -x sars-cov-2 -1 sars-cov-2-reads_1.fq -2 sars-cov-2-reads_2.fq -S sars-cov-2.sam

#### Preparation of index and dictionary files

samtools faidx sars-cov-2.fasta

--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

###### Data pre-processing using PICARD tools

--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

java -jar CreateSequenceDictionary.jar R=sars-cov-2.fasta O=sars-cov-2.dict



java -jar AddOrReplaceReadGroups.jar I=sars-cov-2-mutant.sam O=sars-cov-2-mutantsorted.bam RGID=1 RGLB=445_LIB RGPL=illumina RGSM=RNA RGPU=illumina

java -jar  SortSam.jar I=sars-cov-2-mutantsorted.bam O=sars-cov-2-mutantsort2.bam CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT SO=coordinate

java -jar  MarkDuplicates.jar I=sars-cov-2-mutantsort2.bam O=sars-cov-2-mutantmarkdup.bam CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT M=output.metrics 

--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

###### GATK command-line execution

--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

java -jar GenomeAnalysisTK.jar -T SplitNCigarReads -R sars-cov-2.fasta -I sars-cov-2-mutantmarkdup.bam -o sars-cov-2-mutantsplit.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS

java -jar GenomeAnalysisTK.jar -T HaplotypeCaller -R sars-cov-2.fasta -I sars-cov-2-mutantsplit.bam -o sars-cov-2-mutant.vcf -dontUseSoftClippedBases -stand_call_conf 20.0 -stand_emit_conf 20.0

java -jar GenomeAnalysisTK.jar -T VariantFiltration -R sars-cov-2.fasta -V  sars-cov-2-mutant.vcf -window 35 -cluster 3 -filterName FS -filter "FS > 30.0" -filterName QD -filter "QD < 2.0" -o sars-cov-2-mutantfilter.vcf 

java -jar GenomeAnalysisTK.jar -T SelectVariants -R sars-cov-2.fasta -V sars-cov-2-mutantfilter.vcf -o sars-cov-2-indel.vcf -selectType indel

java -jar GenomeAnalysisTK.jar -T SelectVariants -R sars-cov-2.fasta -V sars-cov-2-mutantfilter.vcf -o sars-cov-2-snp.vcf -selectType snp

--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

###### Final output file

--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
sars-cov-2-snp.vcf - file containing filtered SNP.

sars-cov-2-indel.vcf - file containing filtered INDEL.

--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
