#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: Workflow

inputs:
  sequences:
    type: File
    format: edam:format_1929  # FASTA
  ultrasfast_max_iterations:
    type: int?
    default: 100000

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
      ultrasfast_max_iterations: ultrasfast_max_iterations
      optimize_ufboot:
        default: true
      substitution_model:
        default: HKY
    out: [ result_tree, log ]

  fail_for_warning:
    run:
      class: CommandLineTool
      hints:
        DockerRequirement:
          dockerPull: quay.io/biocontainers/iqtree:1.6.9--he860b03_1
      requirements:
        InlineJavascriptRequirement: {}
      inputs:
        log: File
      baseCommand: grep
      arguments:
        - "WARNING: bootstrap analysis did not converge. You should rerun with higher number of iterations (-nm option)"
      successCodes: [ 1 ]
      permanentFailCodes: [ 0 ]
      outputs:
        warning: stdout
    in:
      log: build_consensus_tree/log
    out: [ warning ]

outputs:
  result_tree:
    type: File
    outputSource: build_consensus_tree/result_tree
  warning:
    type: File
    outputSource: fail_for_warning/warning

$namespaces:
  edam: http://edamontology.org/
$schemas: [ "http://edamontology.org/EDAM_1.18.owl" ]
