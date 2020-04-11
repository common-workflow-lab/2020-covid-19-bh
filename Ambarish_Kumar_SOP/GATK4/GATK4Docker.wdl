workflow variantcall {

  
  File refIndex
  File refDict
  File referenceGenome
  String gatk_docker
  String gatk_path
  String name
 
  call Alignment{
    input:
         ReferenceGenome=referenceGenome,
         sampleName=name,
         index=name
}
  call AddOrReplaceReadGroups{
    input:
          inputSAM=Alignment.rawSAM,
          sampleName=name,
          docker=gatk_docker,
          gatk_path=gatk_path
          
 }
  call SortSam {
    input: 
          inputBAM=AddOrReplaceReadGroups.rawBAM,
          sampleName=name,
          docker=gatk_docker,
          gatk_path=gatk_path

  }
  call ReferenceSeqIndex{
   input:
         ReferenceGenome=referenceGenome,
         sampleName=name
}
  call ReferenceSeqDictionary {
    input:
           docker=gatk_docker,
           gatk_path=gatk_path,
           ReferenceGenome=referenceGenome,
           sampleName=name
} 
  call MarkDuplicates{
   input:
        inputBAM=SortSam.rawBAM,
        docker=gatk_docker,
        gatk_path=gatk_path,
        sampleName=name 
 }
  call SplitNCigarReads{
   input:
      inputBAM=MarkDuplicates.rawBAM,
      RefIndex=refIndex,
      RefDict=refDict,
      docker=gatk_docker,
      gatk_path=gatk_path,
      ReferenceGenome=referenceGenome,
      sampleName=name
}
 
  call HaplotypeCaller{
   input:
       inputBAM=SplitNCigarReads.rawBAM,
       RefIndex=refIndex,
       RefDict=refDict,
       docker=gatk_docker,
       gatk_path=gatk_path,
       ReferenceGenome=referenceGenome,
       sampleName=name
       
 }
           
 call VariantFilteration{
  input:
      mutantVCF=HaplotypeCaller.rawVCF,
       RefIndex=refIndex,
       RefDict=refDict,
       docker=gatk_docker,
       gatk_path=gatk_path,
       ReferenceGenome=referenceGenome,
       sampleName=name
      
      
 }
 call SelectSNPs {
  input:
   mutantVCF=VariantFilteration.rawVCF,
   docker=gatk_docker,
   gatk_path=gatk_path,
   ReferenceGenome=referenceGenome,
   RefIndex=refIndex,
   RefDict=refDict,
   sampleName=name
 }
 call SelectINDELs {
  input:
  mutantVCF=VariantFilteration.rawVCF,
  docker=gatk_docker,
  gatk_path=gatk_path,
  ReferenceGenome=referenceGenome,
  RefIndex=refIndex,
  RefDict=refDict,
  sampleName=name
}
}
task Alignment {
 
  File leftFastq
  File rightFastq
  File ReferenceGenome
  String sampleName
  String index

  command {
          export PATH=$PATH:/home/cbl/20/
          bowtie2-build ${ReferenceGenome} ${index} 
          
          bowtie2 -q -x ${index} -1 ${leftFastq} -2 ${rightFastq} -S ${sampleName}.sam           
}


output {
    
    File rawSAM = "${sampleName}.sam"
  }
}

task ReferenceSeqIndex{
 File ReferenceGenome
 String sampleName
 
 command {
   samtools faidx ${ReferenceGenome} > ${sampleName}.fasta.fai
 }



output {
  File refIndex = "${sampleName}.fasta.fai"  
 }
}

task ReferenceSeqDictionary{
 String docker
 String gatk_path
 File ReferenceGenome
 String sampleName
 command {
    java -jar ${gatk_path} CreateSequenceDictionary -R ${ReferenceGenome} -O ${sampleName}.fasta.dict 


 }

runtime
{
  docker:docker
}

output {
 
 File refDict = "${sampleName}.fasta.dict"  
 } 
}

task AddOrReplaceReadGroups {
  String docker
  String gatk_path
  File inputSAM 
  String sampleName

  command {
   java -jar ${gatk_path} AddOrReplaceReadGroups -I ${inputSAM} -O ${sampleName}.bam -RGID 1 -RGLB 445_LIB -RGPL illumina -RGSM RNA -RGPU illumina
      
  }

runtime
{
  docker:docker
}

  output {
    File rawBAM = "${sampleName}.bam"
  }
}

task SortSam {
 String docker
 String gatk_path
 File inputBAM
 String sampleName
 
 command {
     java -jar ${gatk_path} SortSam -I ${inputBAM} -O ${sampleName}.bam -CREATE_INDEX true -VALIDATION_STRINGENCY LENIENT -SO coordinate}

runtime
{
  docker:docker
}

 
 output {
    File rawBAM = "${sampleName}.bam"
  }
}

task MarkDuplicates {
 String docker
 String gatk_path
 File inputBAM
 String sampleName

 command {
   java -jar ${gatk_path} MarkDuplicates -I ${inputBAM} -O ${sampleName}.markdup.bam -CREATE_INDEX true -VALIDATION_STRINGENCY LENIENT -M output.metrics
}

runtime
{
  docker:docker
}

 output {
  File rawBAM = "${sampleName}.markdup.bam"
 }
} 

task SplitNCigarReads {
 String docker
 String gatk_path
 File inputBAM
 File ReferenceGenome
 File RefIndex
 File RefDict
 String sampleName

 command {
  java -jar ${gatk_path} SplitNCigarReads -R ${ReferenceGenome} -I ${inputBAM} -O ${sampleName}.split.bam  
}

runtime
{
  docker:docker
}

 output {
  File rawBAM = "${sampleName}.split.bam"
}
}

task HaplotypeCaller {
 String docker
 String gatk_path
 File inputBAM
 File ReferenceGenome
 File RefIndex
 File RefDict
 String sampleName

 command {
                samtools index ${inputBAM} > ${sampleName}.split.bam.bai
                      
             java -jar ${gatk_path} HaplotypeCaller -R ${ReferenceGenome} -I ${inputBAM} -O ${sampleName}.mutant.vcf 
}

runtime
{
  docker:docker
}

 output {
  File rawVCF = "${sampleName}.mutant.vcf"
 }
}

task VariantFilteration {
  String docker
  String gatk_path
  File mutantVCF
  File ReferenceGenome
  File RefIndex
  File RefDict
  String sampleName
  

 command {
         java -jar ${gatk_path} IndexFeatureFile -F ${mutantVCF}
      java -jar ${gatk_path} VariantFiltration -R ${ReferenceGenome} -V ${mutantVCF} -window 35 -cluster 3 -O ${sampleName}.mutantfilter.vcf
}

runtime
{
  docker:docker
}

 output {
 
 File rawVCF = "${sampleName}.mutantfilter.vcf"
 }
 
}

task SelectSNPs {
  String docker
  String gatk_path
  File mutantVCF
  File ReferenceGenome
  File RefIndex
  File RefDict
  String sampleName

 command {
   java -jar ${gatk_path} SelectVariants -R ${ReferenceGenome} -V ${mutantVCF} -O ${sampleName}.mutantsnp.vcf -select-type-to-include SNP
}

runtime
{
  docker:docker
}

 output {
 File rawVCF = "${sampleName}.mutantsnp.vcf"
}
}

task SelectINDELs {
  String docker
  String gatk_path
  File mutantVCF
  File ReferenceGenome
  File RefIndex
  File RefDict
  String sampleName

 command {
  java -jar ${gatk_path} SelectVariants -R ${ReferenceGenome} -V ${mutantVCF} -O ${sampleName}.mutantindel.vcf -select-type-to-include INDEL
}

runtime
{
  docker:docker
}

 output {
 File rawVCF = "${sampleName}.mutantindel.vcf"
}
}

