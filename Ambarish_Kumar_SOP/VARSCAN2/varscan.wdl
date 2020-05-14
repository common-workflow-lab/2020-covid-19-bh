workflow varscan {

  String varscanPath
  File referenceGenome
  String name
 
  call Alignment {
    input:
         ReferenceGenome=referenceGenome,
         sampleName=name,
         index=name
}

 call samtoolsView {
    input:
         inputSAM=Alignment.rawSAM,
         sampleName=name
}

 call samtoolsSort {
    input:
         inputBAM=samtoolsView.rawBAM,
         sampleName=name
    
}

 call samtoolsIndex {
    input:
        inputBAM=samtoolsSort.rawBAM,
        sampleName=name
}

 call samtoolsMpileup {
   input:
        ReferenceGenome=referenceGenome,
        inputBAM=samtoolsSort.rawBAM,
        sampleName=name
}

 call mpileup2snp {
   input:
        inputMpileup=samtoolsMpileup.rawMpileup,
        sampleName=name,
        varscanPath=varscanPath
} 
        }

task Alignment {
 
  File leftFastq
  File rightFastq
  File ReferenceGenome
  String sampleName
  String index

  command {
          export PATH=$PATH:/home/ngsap2/Downloads/bowtie2-2.4.1
          bowtie2-build ${ReferenceGenome} ${index} 
          
          bowtie2 -q -x ${index} -1 ${leftFastq} -2 ${rightFastq} -S ${sampleName}.sam           
}



output {
    
    File rawSAM = "${sampleName}.sam"
  }
}

task samtoolsView {
 
 File inputSAM
 String sampleName
 command {
            samtools view -bS ${inputSAM} > ${sampleName}.bam
}

output {
 
 File rawBAM = "${sampleName}.bam"
}

}
task samtoolsSort {

 File inputBAM
 String sampleName

command {

  samtools sort ${inputBAM} ${sampleName}.sorted

}

output {
 
 File rawBAM = "${sampleName}.sorted.bam"

}

}

task samtoolsIndex {

 File inputBAM
 String sampleName

command {

  samtools index ${inputBAM} > ${sampleName}.bam.bai
  
 }

output {
 
 File rawBAM = "${sampleName}.bam.bai"

}
}

task samtoolsMpileup {
 
 File inputBAM
 File ReferenceGenome
 String sampleName 

command {

samtools index ${inputBAM}

samtools mpileup -B -f ${ReferenceGenome} ${inputBAM} > ${sampleName}.mpileup 

}

output {
 
File rawMpileup = "${sampleName}.mpileup" 

}

}
task mpileup2snp
{

File inputMpileup
String sampleName
String varscanPath


command {

java -jar ${varscanPath} mpileup2snp ${inputMpileup} > ${sampleName}.vcf

}

output{

File snpVCF = "${sampleName}.vcf"

}
 
}


