#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

hints:
  DockerRequirement:
    dockerPull: quay.io/biocontainers/mafft:7.458--h516909a_0

baseCommand: mafft

inputs:
  anysymbol:
    type: boolean
    doc: for unusual characters
    inputBinding:
      prefix: --anysymbol

  sequences:
    label: Sequences to align
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
   outputBinding:
     glob: $(inputs.sequences.nameroot).alignment.fasta

$namespaces:
  edam: http://edamontology.org/
$schemas:
  - http://edamontology.org/EDAM_1.18.owl

