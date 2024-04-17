#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

s:license: http://www.apache.org/licenses/LICENSE-2.0

requirements:
  InitialWorkDirRequirement:
    listing:
      - $(inputs.alignments)
hints:
  ResourceRequirement:
    coresMin: 8
    ramMin: 30720
  DockerRequirement:
    dockerPull: quay.io/biocontainers/iqtree:1.6.9--he860b03_1
  SoftwareRequirement:
    packages:
      iqtree:
        version: ["1.6.9"]
        specs: 
          - https://bio.tools/iq-tree
          - https://identifiers.org/RRID/RRID:SCR_017254

baseCommand: iqtree

arguments: 
 - prefix: -nt
   valueFrom: AUTO
 - prefix: -ntmax
   valueFrom: $(runtime.cores)
 - prefix: -seed
   valueFrom: "1"

# TODO: what about these options?
# http://www.phylo.org/index.php/rest/iqtree_xsede.html
# vparam.specify_runtype_=2
# vparam.specify_runtype_ - (Excl) - Specify the nrun type; 1 (Model Selection); 2 (Tree Inference); no default.
# vparam.specify_numparts_=1
# vparam.specify_numparts_ - (Integer) - How many partitions does your data set have.

inputs:
  alignments:
    type: File
    format: edam:format_1929  # FASTA
    inputBinding:
      prefix: -s

  optimize_ufboot:
    type: boolean?
    doc: Optimize UFBoot trees by NNI on bootstrap alignment
    inputBinding:
      prefix: -bnni

  ultrasfast_max_iterations:
    type: int?
    label: "Maximum number of iterations (default: 1000)"
    inputBinding:
      prefix: -nm

  ultrafast_bootstrap_replicates:
    type: int?
    doc: |
      Specify the number of replicates to this option to choose the
      Ultrafast bootstrapping appoach ("-bb")
    inputBinding:
      prefix: -bb

  single_branch_test_replicates:
    type: int?
    doc: |
      Specify the number of replicates to this option to choose the
      SH-like approximate likelihood ratio test (SH-aLRT)

      Specifying 0 replicates chooses the
      Parametric aLRT test (Anisimova and Gascuel 2006)
    label: number of replicates
    inputBinding:
      prefix: -alrt

  substitution_model:
    type: string?
    inputBinding:
      prefix: -m
    doc: |
      .                 DNA: HKY (default), JC, F81, K2P, K3P, K81uf, TN/TrN, TNef,
      .                      TIM, TIMef, TVM, TVMef, SYM, GTR, or 6-digit model
      .                      specification (e.g., 010010 = HKY)
      .             Protein: LG (default), Poisson, cpREV, mtREV, Dayhoff, mtMAM,
      .                      JTT, WAG, mtART, mtZOA, VT, rtREV, DCMut, PMB, HIVb,
      .                      HIVw, JTTDCMut, FLU, Blosum62, GTR20, mtMet, mtVer, mtInv
      .     Protein mixture: C10,...,C60, EX2, EX3, EHO, UL2, UL3, EX_EHO, LG4M, LG4X
      .              Binary: JC2 (default), GTR2
      .     Empirical codon: KOSI07, SCHN05
      .   Mechanistic codon: GY (default), MG, MGK, GY0K, GY1KTS, GY1KTV, GY2K,
      .                      MG1KTS, MG1KTV, MG2K
      .Semi-empirical codon: XX_YY where XX is empirical and YY is mechanistic model
      .      Morphology/SNP: MK (default), ORDERED, GTR
      .      Lie Markov DNA: One of the following, optionally prefixed by RY, WS or MK:
      .                      1.1,  2.2b, 3.3a, 3.3b,  3.3c,
      .                      3.4,  4.4a, 4.4b, 4.5a,  4.5b,
      .                      5.6a, 5.6b, 5.7a, 5.7b,  5.7c,
      .                      5.11a,5.11b,5.11c,5.16,  6.6,
      .                      6.7a, 6.7b, 6.8a, 6.8b,  6.17a,
      .                      6.17b,8.8,  8.10a,8.10b, 8.16,
      .                      8.17, 8.18, 9.20a,9.20b,10.12,
      .                      10.34,12.12
      .      Non-reversible: STRSYM (strand symmetric model, synonymous with WS6.6)
      .      Non-reversible: UNREST (most general unrestricted model, functionally equivalent to 12.12)
      .      Models can have parameters appended in brackets.
      .          e.g. '-mRY3.4{0.2,-0.3}+I' specifies parameters for
      .          RY3.4 model but leaves proportion of invariant sites
      .          unspecified. '-mRY3.4{0.2,-0.3}+I{0.5} gives both.
      .          When this is done, the given parameters will be taken
      .          as fixed (default) or as start point for optimization
      .          (if -optfromgiven option supplied)

      .       Otherwise: Name of file containing user-model parameters
      .                  (rate parameters and state frequencies)

outputs: 
  result_tree:
    label: Maximum-likelihood tree
    type: File
    outputBinding:
      glob: $(inputs.alignments.basename).treefile
  distances:
    label: Likelihood distances
    type: File?
    outputBinding:
      glob: $(inputs.alignments.basename).mldist
  model_information:
    type: File?
    outputBinding:
       glob: $(inputs.alignments.basename).model.gz
  split_supports:
    type: File?
    label: Ultrafast bootstrap approximation split support values
    outputBinding:
       glob: $(inputs.alignments.basename).splits.nex
  ufboot_trees:
    type: File?
    label: Ultrafast bootstrap approximation UFBoot trees
    outputBinding:
       glob: $(inputs.alignments.basename).ufboot
  consensus_tree:
    type: File?
    label: Ultrafast bootstrap approximation consensus tree
    outputBinding:
       glob: $(inputs.alignments.basename).contree
  report:
    label: IQ-TREE Report
    type: File?
    outputBinding:
       glob: $(inputs.alignments.basename).iqtree
  alignment_with_unique_sequences:
    type: File?
    label: alignment with unique sequences
    outputBinding:
      glob: $(inputs.alignments.basename).uniqueseq.phy
  log:
    type: File?
    outputBinding:
      glob: $(inputs.alignments.basename).log

$namespaces:
  s: http://schema.org/
  edam: http://edamontology.org/

# $schemas:
# - "http://schema.org/version/latest/schema.rdf"
# - "http://edamontology.org/EDAM_1.18.owl"
