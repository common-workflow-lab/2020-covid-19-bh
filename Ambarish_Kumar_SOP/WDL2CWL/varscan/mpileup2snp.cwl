#!/usr/bin/env cwl-runner
# This tool description was generated automatically by wdl2cwl ver. 0.2

{
    "baseCommand": [],
    "outputs": [
        {
            "outputBinding": {
                "glob": "$(inputs.sampleName).vcf"
            },
            "type": "File",
            "id": "snpVCF"
        }
    ],
    "inputs": [
        {
            "type": "File",
            "id": "inputMpileup"
        },
        {
            "type": "string",
            "id": "sampleName"
        },
        {
            "type": "string",
            "id": "varscanPath"
        }
    ],
    "requirements": [
        {
            "class": "ShellCommandRequirement"
        },
        {
            "class": "InlineJavascriptRequirement"
        }
    ],
    "arguments": [
        {
            "shellQuote": false,
            "valueFrom": "java -jar $(inputs.varscanPath) mpileup2snp $(inputs.inputMpileup.path) > $(inputs.sampleName).vcf"
        }
    ],
    "cwlVersion": "v1.0",
    "class": "CommandLineTool",
    "id": "mpileup2snp"
}