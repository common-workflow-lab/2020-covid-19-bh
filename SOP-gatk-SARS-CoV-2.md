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

bowtie2-build <reference genome> <index base name>

Where

bowtie2-build - Bowtie command to build index of reference genome.

<reference genome> - fasta format reference genome.

<index base name> - Base name of generated index file.

Alignment 

Alignment is matching read bases with that of reference genome sequence. It is performed for GATK using bowtie aligner as data preparation step. 

Command-line usage:

bowtie2 -q -x <index> -1 <left read> -2 <right read> -S <alignment file>

Where

bowtie2 - Bowtie command for alignment.

q - bowtie option for fastq input file format.

<index> - prefix for bowtie index.

1 - option for left-end reads. 

2 - option for right-end reads.

S - option for .sam format alignment file.

<left read> - reads generated from forward strand.

<right read> - reads generated from reverse strand.

<alignment file> - alignment file in .sam format.
 
## Add read groups, sort, mark duplicates, and create index

In this step, we merge all read groups with different  sequencing lane reads and read group header @RG and assign a single read group header.

Command-line usage:

java -jar  AddOrReplaceReadGroups.jar - I <input file> - O <output file name>  - RGID <rgid id> - RGLB <rglb string>  - RGPL <rgpl string> - RGSM <rgsm string> - RGPU <rgpu string>

Where

AddOrReplaceReadGroup.jar - Picard tool to merge all read groups and assign single read group ID.

<input file> - Alignment file in .sam format.

<output file> - Alignment file in .sam format with merged read group ID.

<rglb string>   -   Read Group Library Required. 

<rgpl string>   -   Read Group platform (e.g. illumina, solid) Required

<rgsm string>   -   Read Group sample name Required

<rgpu string>   -   Read Group platform unit Required

## Sort and Index
This step sorts and indexes the reads, and converts between SAM to BAM format. The indexed bam file is further coordinate sorted using the PICARD tool SortSam.jar

Command-line usage:

java -jar SortSam.jar I=<input file> O=<output file> CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT SORT_ORDER=coordinate

Where

<Input file> -  Sorted bam file

<Output file>  - The sorted BAM or SAM output file

## Mark the duplicate reads

Duplicate sequenced reads are marked and removed using PICARD tool MarkDuplicates.jar

Command-line usage:

java -jar MarkDuplicates.jar I=<input file> O=<output file> CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT M=<output metrics> 

Where

<input file> - sorted bam input file

<output file> - markdup bam file

<output metrices> - File to write duplication metrics to Required

## Split'N'Trim and reassign mapping qualities

Hardclip the hanging part of reads into intronic region and assign reads mapping quality of 60 replacing 255  which is understood and interpreted by GATK

Command-line usage:

java -jar GenomeAnalysisTK.jar -T SplitNCigarReads -R <reference sequence> -I <inputfile> -o <output> -rf ReassignOneMappingQuality -RMQF <rmqf value> -RMQT <rmsqt value> -U ALLOW_N_CIGAR_READS

Where

<reference sequence> - reference sequence in fasta format.

<input file> - markedup .bam file.

<outputfile> - split .bam output file.

<rf> - allows read filter and mapping quality. 

<rmqf value> - RMQF value 255

<rmqt value> - RMQT value 60

## Variant calling

HaplotypeCaller - calls all plausible haplotypes and detect variants. 

Command-line usage:

java -jar GenomeAnalysisTK.jar -T HaplotypeCaller -R <reference sequence>  -I <inputfile> -o <outputfile> -dontUseSoftClippedBases -stand_call_conf <stand call conf value> -stand_emit_conf <stand emit conf  value>

Where

<inputfile>   -   split bam file.

<outputfile>  -  vcf file as output.

<stand call conf value>   -   20.0 confidence interval value

<stand emit conf value>   -   20.0 confidence interval value

-stand_call_conf   - The minimum phred-scaled confidence threshold at which variants        should be called

-stand_call_emit   -   The minimum phred-scaled confidence threshold at which variants should be emitted (and filtered with LowQual if less than the calling threshold)

## Variant filtering

Filter out the low quality variants.

Command-line usage

java -jar GenomeAnalysisTK.jar -T VariantFiltration -R <reference genome> -V  <input file> -window <window size> -cluster <cluster size> -filterName <filter name> -filter "FS > 30.0" -filterName <filter name> -filter "QD < 2.0" -o <output>

Where

<reference genome>   -  reference sequence in fasta file format.

<input file>  -   .vcf file i.e variant set to used as an input.

<window size>   -   window size in which to evaluate clustered snp’s.

<cluster size>   -   the no of snp’s which makeup a cluster.

<filter name>   -   name of the filter.

<output>  -   filtered variants or VCF’s.
