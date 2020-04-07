### Author

AMBARISH KUMAR

er.ambarish@gmail.com

ambari73_sit@jnu.ac.in



Standard Operating Protocol using SAMTools serves the purpose of genome variants discovery i.e SNPs and INDELS. It will be useful in analysing genomic variations into SARS-CoV-2 genome.



It utilizes paired-end illumina RNASEQ datasets and reference genome of SARS-CoV-2 virus.



### Indexing 

Command-line usage:

bowtie2-build <reference genome> <index base name>

where

bowtie2-build - Bowtie command to build index of reference genome.

<reference genome> - fasta format reference genome.

<index base name> - Base name of generated index file.

Alignment 

Command-line usage:

bowtie2 -q -x <index> -1 <left read> -2 <right read> -S <alignment file>

where

bowtie2 - Bowtie command for alignment.

<index> - prefix for bowtie index.

<left read> - reads generated from forward strand.

<right read> - reads generated from reverse strand.

<alignment file> - alignment file in .sam format.

### Conversion of .sam file format to .bam file format

Samtools view utility is used to convert alignment file format from .sam format to .bam format.

Command-line usage:

samtools view -bS <input file> > <output file>

where

<input file> - Alignment file in .sam format.

<output file> - Alignment file in .bam format.

### Sorting of alignment file

Samtools sort utility sorts bam format alignment file as per coordinate order.

Command-line usage:

samtools sort <input file> <base name>

where

<input file> - alignment file in .bam format.

<base name> - Base name for sorted .bam format alignment file.

### Indexing of sorted file

Samtools index utility prepare index of sorted file.

Command-line usage:

samtools index <input file>

where

<input file> - sorted alignment file in .bam format.

### Pileup file creation from sorted bam file

Samtools utility mpileup generates pileup file for each genomic base position. 

Command-line usage:

samtools mpileup -B -f <reference genome> <sorted bam file> > <pileup file>

where

<reference genome> - Fasta format reference genome.

<sorted bam file> - Coordinate sorted alignment file in .bam format.

<pileup file> - pileup file for each genomic position.

### Variants calling

Bcftools view utility calls variants from samtools generated pileup file.

Command-line usage:

bcftools view -vcg <input file> > <vcf output>

where

<input file> - samtools generated pileup file.

<vcf output> - .vcf format file containing raw output.


###### Command line implementation.


#### Input dataset

SARS-CoV-2.fasta - SARS-CoV-2 reference genome.
SARS-CoV-2-reads_1.fq - left-end reads.
SARS-CoV-2-reads_2.fq - right-end reads.
 
#### Preparation of alignment file

bowtie2-build SARS-CoV-2.fasta SARS-CoV-2

bowtie2 -q -x SARS-CoV-2 -1 SARS-CoV-2-reads_1.fq -2 SARS-CoV-2-reads_2.fq -S SARS-CoV-2-mutant.sam

#### Pre-processing of alignment-file

samtools view -bS SARS-CoV-2-mutant.sam > SARS-CoV-2-mutant.bam

samtools sort SARS-CoV-2-mutant.bam SARS-CoV-2-mutantsorted

samtools index SARS-CoV-2-mutantsorted.bam

#### Pileup-file generation

samtools mpileup -uf SARS-CoV-2.fasta SARS-CoV-2-mutantsorted.bam | bcftools view -vcg - > SARS-CoV-2-mutantraw.vcf


Output file - SARS-CoV-2-mutantraw.vcf is combined file containing SNPs and INDELs.


