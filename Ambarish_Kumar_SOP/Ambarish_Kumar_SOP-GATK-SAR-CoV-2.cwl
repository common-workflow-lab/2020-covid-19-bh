#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: Workflow

doc: |
  Author: AMBARISH KUMAR er.ambarish@gmail.com & ambari73_sit@jnu.ac.in
  This is a proposed standard operating procedure for genomic variant detection using GATK4.
  It is hoped to be effective and useful for getting SARS-CoV-2 genome variants.
  
  It uses Illumina RNASEQ reads and genome sequence.

inputs:
  sars_cov_2_reference_genome:
    type: File
    format: edam:format_1929  # FASTA

  rnaseq_left_reads:
    type: File
    format: edam:format_1930  # FASTQ

  rnaseq_right_reads:
    type: File
    format: edam:format_1930  # FASTQ

steps:
  index_reference_genome_with_bowtie2:
    run: ../tools/bowtie2/bowtie2_build.cwl
    in:
      reference_in: sars_cov_2_reference_genome
      bt2_index_base:
        valueFrom: "sars-cov-2"
    out: [ indices ]

  align_rnaseq_reads_to_genome:
    run: ../tools/bowtie2/bowtie2_align.cwl
    in:
      indices_file: index_reference_genome_with_bowtie2/indices
      filelist: rnaseq_left_reads
      filelist_mates: rnaseq_right_reads
      output_filename:
        valueFrom: sars-cov-2.sam
    out: [ output ]

  index_reference_genome_with_samtools:
    run: ../tools/samtools/samtools_faidx.cwl
    in:
      sequences: sars_cov_2_reference_genome
    out: [sequences_with_index]

  create_sequence_dictionary:
    run: ../tools/picard/picard_CreateSequenceDictionary.cwl
    in:
      REFERENCE: index_reference_genome_with_samtools/sequences_with_index
    out: [ sequences_with_dictionary ]

  update_read_group:
    run: ../tools/picard/picard_AddOrReplaceReadGroups.cwl
    in:
      INPUT: align_rnaseq_reads_to_genome/output
      OUTPUT:
        valueFrom: sars-cov-2-newreadgroups.bam
      RGID:
        valueFrom: "1"
      RGLB:
        valueFrom: 445_LIB
      RGPL:
        valueFrom: illumina
      RGSM:
        valueFrom: RNA
      RGPU:
        valueFrom: illumina
      SORT_ORDER:
        valueFrom: coordinate
    out: [ sequences_with_new_read_group ]
 
  mark_duplicates:
    run: ../tools/picard/picard_markdup.cwl
    in:
      bam_sorted: update_read_group/sequences_with_new_read_group
    out: [ bam_duprem ]

  split_N_cigar_reads:
    run: ../tools/GATK/GATK-SplitNCigarReads.cwl
    in:
      reference: create_sequence_dictionary/sequences_with_dictionary
      reads: mark_duplicates/bam_duprem
      output_filename:
        valueFrom: sars-cov-2-mutantsplit.bam
      read_filter:
        valueFrom: ReassignOneMappingQuality
    out: [ output ]
 
outputs: []
