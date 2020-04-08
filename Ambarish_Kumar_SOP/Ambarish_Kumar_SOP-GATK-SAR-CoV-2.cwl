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
  index_reference_genome:
    run: ../tools/bowtie2/bowtie2_build.cwl
    in:
      reference_in: sars_cov_2_reference_genome
      bt2_index_base:
        valueFrom: "sars-cov-2"
    out: [ indices ]

  align_rnaseq_reads_to_genome:
    run: ../tools/bowtie2/bowtie2_align.cwl
    in:
      indices_folder: index_reference_genome/indices
      filelist: rnaseq_left_reads
      filelist_mates: rnaseq_right_reads
      output_filename:
        valueFrom: sars-cov-2.sam
    out: [ output ]


outputs: []
