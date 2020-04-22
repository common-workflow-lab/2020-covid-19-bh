#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: Workflow

inputs:
  sequences:
    type: File
    format: edam:format_1929  # FASTA

steps:
  align:
    run: mafft.cwl
    in:
      sequences: sequences
      anysymbol:
        default: true
    out: [ alignment ]

  build_consensus_tree:
    run: iqtree.cwl
    in:
      alignments: align/alignment
      ultrafast_bootstrap_replicates:
        default: 1000
      optimize_ufboot:
        default: true
      substitution_model:
        default: HKY
    out: [ consensus_tree ]

outputs:
  consensus_tree:
    type: File
    outputSource: build_consensus_tree/consensus_tree

$namespaces:
  edam: http://edamontology.org/
$schemas:
  - http://edamontology.org/EDAM_1.18.owl
