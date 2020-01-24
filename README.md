# RNASeqQC

QC metrics for RNASeq data

## Overview

## Dependencies

* [samtools 1.9](https://github.com/samtools/samtools)
* [picard 2.21.2](https://broadinstitute.github.io/picard/command-line-overview.html)
* [production-tools-python 1.0.1](https://bitbucket.oicr.on.ca/projects/GSI/repos/production-tools-python/)
* [bwa 0.7.17](https://github.com/lh3/bwa/releases/download/v0.7.17/bwa-0.7.17.tar.bz2)
* [bam-qc-metrics 0.2.3](https://github.com/oicr-gsi/bam-qc-metrics.git)


## Usage

### Cromwell
```
java -jar cromwell.jar run RNASeqQC.wdl --inputs inputs.json
```

### Inputs

#### Required workflow parameters:
Parameter|Value|Description
---|---|---
`bamFile`|File|Input BAM file on which to compute QC metrics
`bwaRef`|File|Ribosomal reference file in FASTA format, for alignment by BWA
`refFlat`|File|Reference flat file for Picard CollectRNASeqMetrics
`humanGenomeRef`|File|Human genome FASTA reference for Picard CollectRNASeqMetrics
`humanGenomeRefIndex`|File|Human genome FASTA reference index for Picard CollectRNASeqMetrics
`picardJarDir`|String|Directory containing the Picard JAR


#### Optional workflow parameters:
Parameter|Value|Default|Description
---|---|---|---
`outputFileNamePrefix`|String|"RNASeqQC"|Prefix for output files


#### Optional task parameters:
Parameter|Value|Default|Description
---|---|---|---
`bamqc.bamqcSuffix`|String|"bamqc.json"|Suffix for output file
`bamqc.modules`|String|"bam-qc-metrics/0.2.3"|required environment modules
`bamqc.jobMemory`|Int|16|Memory allocated for this job
`bamqc.threads`|Int|4|Requested CPU threads
`bamqc.timeout`|Int|4|hours before task timeout
`bamToFastq.suffixAll`|String|"all.fastq"|Suffix for FASTQ file of all reads
`bamToFastq.suffixR1`|String|"R1.fastq"|Suffix for FASTQ file of read 1
`bamToFastq.suffixR2`|String|"R2.fastq"|Suffix for FASTQ file of read 2
`bamToFastq.modules`|String|"samtools/1.9"|required environment modules
`bamToFastq.jobMemory`|Int|16|Memory allocated for this job
`bamToFastq.threads`|Int|4|Requested CPU threads
`bamToFastq.timeout`|Int|4|hours before task timeout
`bwaMem.contamSuffix`|String|"contaminationBwaFlagstat.txt"|Suffix for output file
`bwaMem.modules`|String|"samtools/1.9 bwa/0.7.17"|required environment modules
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
`picard.strandSpecificity`|String|"NONE"|String to denote strand specificity for Picard
`picard.modules`|String|"picard/2.21.2"|required environment modules
`picard.jobMemory`|Int|64|Memory allocated for this job
`picard.threads`|Int|4|Requested CPU threads
`picard.timeout`|Int|4|hours before task timeout
`collate.collatedSuffix`|String|"collatedMetrics.json"|Suffix for output file
`collate.modules`|String|"production-tools-python/0"|required environment modules
`collate.jobMemory`|Int|16|Memory allocated for this job
`collate.threads`|Int|4|Requested CPU threads
`collate.timeout`|Int|4|hours before task timeout


### Outputs

Output | Type | Description
---|---|---
`result`|File|JSON file of collated RNASeqQC output


## Niassa + Cromwell

This WDL workflow is wrapped in a Niassa workflow (https://github.com/oicr-gsi/pipedev/tree/master/pipedev-niassa-cromwell-workflow) so that it can used with the Niassa metadata tracking system (https://github.com/oicr-gsi/niassa).

* Building
```
mvn clean install
```

* Testing
```
mvn clean verify \
-Djava_opts="-Xmx1g -XX:+UseG1GC -XX:+UseStringDeduplication" \
-DrunTestThreads=2 \
-DskipITs=false \
-DskipRunITs=false \
-DworkingDirectory=/path/to/tmp/ \
-DschedulingHost=niassa_oozie_host \
-DwebserviceUrl=http://niassa-url:8080 \
-DwebserviceUser=niassa_user \
-DwebservicePassword=niassa_user_password \
-Dcromwell-host=http://cromwell-url:8000
```

## Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .

_Generated with wdl_doc_gen (https://github.com/oicr-gsi/wdl_doc_gen/)_
