# rnaSeqQC

QC metrics for RNASeq data

## Overview

## Dependencies

* [samtools 1.9](https://github.com/samtools/samtools)
* [picard 2.21.2](https://broadinstitute.github.io/picard/command-line-overview.html)
* [production-tools-python 1.1.2](https://bitbucket.oicr.on.ca/projects/GSI/repos/production-tools-python/)
* [bwa 0.7.17](https://github.com/lh3/bwa/releases/download/v0.7.17/bwa-0.7.17.tar.bz2)
* [bam-qc-metrics 0.2.5](https://github.com/oicr-gsi/bam-qc-metrics.git)


## Usage

### Cromwell
```
java -jar cromwell.jar run rnaSeqQC.wdl --inputs inputs.json
```

### Inputs

#### Required workflow parameters:
Parameter|Value|Description
---|---|---
`bwaMem.refFasta`|String|Path to human genome FASTA reference
`bwaMem.modules`|String|required environment modules
`picard.refFlat`|String|Path to Picard flatfile reference
`picard.refFasta`|String|Path to human genome FASTA reference
`picard.modules`|String|Required environment modules


#### Optional workflow parameters:
Parameter|Value|Default|Description
---|---|---|---
`bamFile`|File?|None|Input BAM file on which to compute QC metrics
`inputFastqs`|Array[Pair[Pair[File,File],String]]?|None|Array of pairs of fastq files together with RG information strings
`outputFileNamePrefix`|String|"rnaSeqQC"|Prefix for output files
`strandSpecificity`|String|"NONE"|Indicates if we have strand-specific data, could be NONE or empty, default: NONE


#### Optional task parameters:
Parameter|Value|Default|Description
---|---|---|---
`star.indexBam_timeout`|Int|48|hours before task timeout
`star.indexBam_modules`|String|"picard/2.19.2"|modules for running indexing job
`star.indexBam_jobMemory`|Int|12|Memory allocated indexing job
`star.runStar_timeout`|Int|72|hours before task timeout
`star.runStar_jobMemory`|Int|64|Memory allocated for this job
`star.runStar_threads`|Int|6|Requested CPU threads
`star.runStar_peOvMMp`|Float|0.1|maximum proportion of mismatched bases in the overlap area
`star.runStar_peOvNbasesMin`|Int|12|minimum number of overlap bases to trigger mates merging and realignment
`star.runStar_chimOutJunForm`|Int?|None|flag to add metadata to chimeric junction output for functionality with starFusion - 1 for metadata, 0 for no metadata
`star.runStar_chimNonchimScoDMin`|Int|10|to trigger chimeric detection, the drop in the best non-chimeric alignment score with respect to the read length has to be greater than this value
`star.runStar_chimMulmapNmax`|Int|20|maximum number of chimeric multi-alignments
`star.runStar_chimScoJunNonGTAG`|Int|-4|penalty for a non-GTAG chimeric junction
`star.runStar_chimMulmapScoRan`|Int|3|the score range for multi-mapping chimeras below the best chimeric score
`star.runStar_alignIntMax`|Int|100000|maximum intron size
`star.runStar_alignMatGapMax`|Int|100000|maximum gap between two mates
`star.runStar_alignSJDBOvMin`|Int|10|minimum overhang for annotated spliced alignments
`star.runStar_chimJunOvMin`|Int|12|minimum overhang for a chimeric junction
`star.runStar_chimSegmin`|Int|12|minimum length of chimeric segment length
`star.runStar_multiMax`|Int|-1|multiMax parameter for STAR
`star.runStar_saSparsed`|Int|2|saSparsed parameter for STAR
`star.runStar_uniqMAPQ`|Int|255|Score for unique mappers
`star.runStar_modules`|String|"star/2.7.3a hg38-star-index100/2.7.3a"|modules for running STAR
`star.runStar_addParam`|String?|None|Additional STAR parameters
`star.runStar_genereadSuffix`|String|"ReadsPerGene.out"|ReadsPerGene file suffix
`star.runStar_chimericjunctionSuffix`|String|"Chimeric.out"|Suffix for chimeric junction file
`star.runStar_transcriptomeSuffix`|String|"Aligned.toTranscriptome.out"|Suffix for transcriptome-aligned file
`star.runStar_starSuffix`|String|"Aligned.sortedByCoord.out"|Suffix for sorted file
`star.runStar_genomeIndexDir`|String|"$HG38_STAR_INDEX100_ROOT/"|Directory with STAR index files
`bamqc.bamqcSuffix`|String|"bamqc.json"|Suffix for output file
`bamqc.modules`|String|"bam-qc-metrics/0.2.5"|required environment modules
`bamqc.jobMemory`|Int|16|Memory allocated for this job
`bamqc.threads`|Int|4|Requested CPU thread
`bamqc.timeout`|Int|4|hours before task timeout
`bwaMem.contamSuffix`|String|"contaminationBwaFlagstat.txt"|Suffix for output file
`bwaMem.threads`|Int|4|Requested CPU threads
`bwaMem.jobMemory`|Int|16|Memory allocated for this job
`bwaMem.timeout`|Int|4|hours before task timeout
`countUniqueReads.uniqueReadsSuffix`|String|"uniqueReads.txt"|Suffix for output file
`countUniqueReads.modules`|String|"samtools/1.9"|required environment modules
`countUniqueReads.jobMemory`|Int|16|Memory allocated for this job
`countUniqueReads.threads`|Int|4|Requested CPU threads
`countUniqueReads.timeout`|Int|4|hours before task timeout
`picard.picardMem`|Int|6000|Memory to run picard JAR, in MB
`picard.picardSuffix`|String|"picardCollectRNASeqMetrics.txt"|Suffix for output file
`picard.jobMemory`|Int|64|Memory allocated for this job
`picard.threads`|Int|4|Requested CPU threads
`picard.timeout`|Int|4|hours before task timeout
`collate.collatedSuffix`|String|"collatedMetrics.json"|Suffix for output file
`collate.modules`|String|"production-tools-python/2"|required environment modules
`collate.jobMemory`|Int|16|Memory allocated for this job
`collate.threads`|Int|4|Requested CPU threads
`collate.timeout`|Int|4|hours before task timeout


### Outputs

Output | Type | Description
---|---|---
`result`|File|JSON file of collated rnaSeqQC output


## Commands
 This section lists command(s) run by rnaSeqQC workflow
 
 * Running rnaSeqQC
 
 Run QC analysis on RNAseq data
 
 ```
     write_fast_metrics.py
     -b INPUT_BAM 
     -o RESULTS_NAME
 
 ```
 Alternative to sorting, collate used here to realign with bwa mem and run samtools flagstat:
 
 ```
     samtools collate 
     -O 
     --reference REFERENCE_FASTA
     INPUT_BAM
     | samtools fastq - 
     | bwa mem -M -t 8 -p 
     BWA_INDEXES_ROOTNAME
     - 
     | 
     samtools flagstat - > RESULTS_NAME
 
 ```
 
 Collate Results:
 
 ```
     COLLATION_SCRIPT
     --bamqc BAM_QC_RESULTS
     --contam CONTAMINATION_RESULTS
     --picard PICARD_RESULTS 
     --unique-reads UNIQUE_READS
     --out RESULT_NAME
     STRANDEDNESS_OPTION
 ```
 
 Count the number of primary alignments:
 
 ```
     samtools view -F 256 INPUT_BAM | wc -l > RESULT_FILE
 
 ```
 
 Collect RNASEQ metrics with picard:
 
 ```
     java -Xmx~[PICARD_MEMORY]M 
     -jar picard.jar CollectRnaSeqMetrics 
     I=INPUT_BAM 
     O=RESULT_NAME 
     STRAND_SPECIFICITY=STRAND_SPECIFICITY 
     REF_FLAT=REF_FLAT 
     REFERENCE_SEQUENCE=REF_FASTA 
     VALIDATION_STRINGENCY=SILENT
 
 ```
 ## Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .

_Generated with generate-markdown-readme (https://github.com/oicr-gsi/gsi-wdl-tools/)_
