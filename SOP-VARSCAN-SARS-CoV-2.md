## Author

AMBARISH KUMAR

er.ambarish@gmail.com

ambari73_sit@jnu.ac.in

SOP for genomic variant discovery using VARSCAN would be beneficial in detecting SARS2 genomic variations.

It utilises paired-end illunina RNASEQ reads and reference genome for SARS2 virus.

## Mapping to reference genome

Mutated read data are mapped against reference genome using bowtie.
Mapping steps can further be split into indexing and alignment steps.

#### Indexing 

Command-line usage:

bowtie2-build reference_genome index_base_name

where

bowtie2-build - Bowtie command to build index of reference genome.

reference_genome - fasta format reference genome.

index_base_name - Base name of generated index file.

#### Alignment 

Command-line usage:

bowtie2 -q -x index -1 left_read -2 right_read -S alignment_file

where

bowtie2 - Bowtie command for alignment.

index - prefix for bowtie index.

left_read - reads generated from forward strand.

right_read - reads generated from reverse strand.

alignment_file - alignment file in .sam format.

#### Conversion of .sam file format to .bam file format

Samtools view utility is used to convert alignment file format from .sam format to .bam format.

Command-line usage:

samtools view -bS input_file > output_file

where

input_file - Alignment file in .sam format.

output_file - Alignment file in .bam format.

#### Sorting of alignment file

Samtools sort utility sorts bam format alignment file as per coordinate order.

Command-line usage:

samtools sort input_file base_name

where

input_file - alignment file in .bam format.

base_name - Base name for sorted .bam format alignment file.

#### Indexing of sorted file

Samtools index utility prepare index of sorted file.

Command-line usage:

samtools index input_file

where

input_file - sorted alignment file in .bam format.

#### Pileup file creation from sorted bam file

Samtools utility mpileup generates pileup file for each genomic base position. 

Command-line usage:

samtools mpileup -B -f reference_genome sorted_bam_file > pileup_file

where

reference_genome - Fasta format reference genome.

sorted_bam_file - Coordinate sorted alignment file in .bam format.

pileup_file - pileup file for each genomic position.
