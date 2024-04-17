#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

hints:
  DockerRequirement:
    dockerPull: quay.io/biocontainers/mafft:7.458--h516909a_0
  ResourceRequirement:
    coresMin: 8
    ramMin: 40000

baseCommand: mafft

inputs:
  save_memory:
    type: boolean?
    doc: for long genomic sequences
    inputBinding:
      prefix: --memsave

  no_save_memory:
    type: boolean?
    inputBinding:
      prefix: --nomemsave

  anysymbol:
    type: boolean
    doc: for unusual characters
    inputBinding:
      prefix: --anysymbol

  sequences:
    label: Sequences to align
    format: edam:format_1929  # FASTA
    type: File
    inputBinding:
      position: 1

arguments:
  - --auto
  - prefix: --thread
    valueFrom: $(runtime.cores)

stdout: $(inputs.sequences.nameroot).alignment.fasta

outputs:
  alignment:
    type: File
    format: edam:format_1929  # FASTA
    streamable: true
    outputBinding:
      glob: $(inputs.sequences.nameroot).alignment.fasta

$namespaces:
  edam: http://edamontology.org/
# $schemas: [ "http://edamontology.org/EDAM_1.18.owl" ]
