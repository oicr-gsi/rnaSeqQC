version 1.0

workflow RNASeqQC {

    input {
	File bamFile
	File bwaRef
	File collationScript
	File refFlat
	String picardJarDir
	String outputFileNamePrefix = "RNASeqQC"
    }
    
    call bamqc {
	input:
	bamFile = bamFile,
	outputFileNamePrefix = outputFileNamePrefix
    }

    call bamToFastq {
	input:
	bamFile = bamFile,
	outputFileNamePrefix = outputFileNamePrefix
    }

    call bwaMem {
	input:
	fastqR1=bamToFastq.fastqR1,
	fastqR2=bamToFastq.fastqR2,
	bwaRef=bwaRef,
	outputFileNamePrefix = outputFileNamePrefix
    }

    call countUniqueReads {
	input:
	bamFile = bamFile,
	outputFileNamePrefix = outputFileNamePrefix
    }
    
    call picard {
	input:
	bamFile = bamFile,
	refFlat = refFlat,
	picardJarDir = picardJarDir,
	outputFileNamePrefix = outputFileNamePrefix
    }

    call collate {
	input:
	collationScript = collationScript,
	bamqc = bamqc.result,
	contam = bwaMem.result,
	picard = picard.result,
	uniqueReads = countUniqueReads.result,
	outputFileNamePrefix = outputFileNamePrefix
    }
    
    output {
	File result = collate.collatedResults
    }

    meta {
	author: "Iain Bancarz"
	email: "ibancarz@oicr.on.ca"
	description: "QC metrics for RNASeq data"
	dependencies: [
	{
	    name: "samtools/1.9",
	    url: "https://github.com/samtools/samtools"
	},
	{
	    name: "picard/2.21.2",
	    url: "https://broadinstitute.github.io/picard/command-line-overview.html"
	},
	{
	    name: "production-tools-python/1.0.1",
	    url: "https://bitbucket.oicr.on.ca/projects/GSI/repos/production-tools-python/"
	},
	{
	    name: "bwa/0.7.17",
	    url: "https://github.com/lh3/bwa/releases/download/v0.7.17/bwa-0.7.17.tar.bz2"
	},
	{
	    name: "bam-qc-metrics/0.2.3",
	    url: "https://github.com/oicr-gsi/bam-qc-metrics.git"
	},
	]
    }
}


task bamqc {

    input {
	File bamFile
	String outputFileNamePrefix
	String bamqcSuffix = "bamqc.json"
    }

    parameter_meta {
	bamFile: "Input BAM file of aligned RNASeqQC data"
	outputFileNamePrefix: "Prefix for output file"
	bamqcSuffix: "Suffix for output file"
    }

    String resultName = "~{outputFileNamePrefix}.~{bamqcSuffix}"
    
    command <<<
	write_fast_metrics.py \
	-b ~{bamFile} \
	-o ~{resultName} \
    >>>

    output {
	File result = "~{resultName}"
    }

    meta {
	output_meta: {
            result: "JSON file containing BAMQC metrics",
	}
    }
}

task bamToFastq {

    input {
	File bamFile
	String outputFileNamePrefix
	String suffixAll = "all.fastq"
	String suffixR1 = "R1.fastq"
	String suffixR2 = "R2.fastq"
    }

    parameter_meta {
	bamFile: "Input BAM file of aligned RNASeqQC data"
	outputFileNamePrefix: "Prefix for output files"
	suffixAll: "Suffix for FASTQ file of all reads"
	suffixR1: "Suffix for FASTQ file of read 1"
	suffixR2: "Suffix for FASTQ file of read 2"
    }

    String allFastq = "~{outputFileNamePrefix}.~{suffixAll}"
    String R1 = "~{outputFileNamePrefix}.~{suffixR1}"
    String R2 = "~{outputFileNamePrefix}.~{suffixR2}"
    
    command <<<
	samtools bam2fq ~{bamFile} > ~{allFastq} && \
	cat ~{allFastq} | grep '^@.*/1$' -A 3 --no-group-separator > ~{R1} && \
	cat ~{allFastq} | grep '^@.*/2$' -A 3 --no-group-separator > ~{R2}
    >>>

    output {
	File fastqR1 = "~{R1}"
	File fastqR2 = "~{R2}"
    }

    meta {
	output_meta: {
            fastqR1: "FASTQ file for read 1",
	    fastqR2: "FASTQ file for read 2"
	}
    }
}

task bwaMem {

    input {
	File fastqR1
	File fastqR2
	File bwaRef
	String outputFileNamePrefix
	String contamSuffix = "contaminationBwaFlagstat.txt"
    }

    parameter_meta {
	fastqR1: "FASTQ file for read 1"
	fastqR2: "FASTQ file for read 2"
	bwaRef: "Ribosomal reference file for alignment by BWA"
	outputFileNamePrefix: "Prefix for output file"
	contamSuffix: "Suffix for output file"
    }

    String resultName = "~{outputFileNamePrefix}.~{contamSuffix}"

    command <<<
	bwa mem \
	-M \
	-t 8 \
	~{bwaRef} \
	~{fastqR1} \
	~{fastqR2} \
	| \
	samtools view -S -b - \
	| \
	samtools flagstat - > ~{resultName}
    >>>

    output {
	File result = "~{resultName}"
    }

    meta {
	output_meta: {
            result: "Text file with results of running 'samtools flagstat' on BWA output"
	}
    }
}

task collate {

    input {
	File collationScript
	File bamqc
	File contam
	File picard
	File uniqueReads
	String outputFileNamePrefix
	String collatedSuffix = "collatedMetrics.json"
    }

    parameter_meta {
	collationScript: "Python script to collate intermediate results"
	bamqc: "JSON output from bamqc task"
	contam: "Text output from ribosomal contamination check by bwaMem task"
	picard: "Text output from picard task"
	uniqueReads: "Text output from uniqueReads task"
	outputFileNamePrefix: "Prefix for output file"
	collatedSuffix: "Suffix for output file"
    }

    String resultName = "~{outputFileNamePrefix}.~{collatedSuffix}"
    
    command <<<
	~{collationScript} \
	--bamqc ~{bamqc} \
	--contam ~{contam} \
	--picard ~{picard} \
	--unique-reads ~{uniqueReads} \
	--out ~{resultName}
    >>>
    
    output {
	File collatedResults="~{resultName}"
    }

    meta {
	output_meta: {
            result: "JSON file of collated RNASeqQC output"
	}
    }
}

task countUniqueReads {

    input {
	File bamFile
	String outputFileNamePrefix
	String uniqueReadsSuffix = "uniqueReads.txt"
    }

    parameter_meta {
	bamFile: "Input BAM file of aligned RNASeqQC data"
	outputFileNamePrefix: "Prefix for output file"
	uniqueReadsSuffix: "Suffix for output file"
    }

    String resultName = "~{outputFileNamePrefix}.~{uniqueReadsSuffix}"

    command <<<
	samtools view -F 256 ~{bamFile} \
	| wc -l \
	> ~{resultName}
    >>>
    
    output {
	File result = "~{resultName}"
    }

    meta {
	output_meta: {
            result: "Text file with unique read count"
	}
    }
}

task picard {

    input {
	File bamFile
	String picardJarDir
	File refFlat
	String outputFileNamePrefix
	Int picardMem=6000
	String picardSuffix = "picardCollectRNASeqMetrics.txt"
	String strandSpecificity="NONE"
    }

    parameter_meta {
	bamFile: "Input BAM file of aligned RNASeqQC data"
	picardJarDir: "Path of directory containing the Picard JAR"
	refFlat: "Flat reference file required by Picard"
	outputFileNamePrefix: "Prefix for output file"
	picardMem: "Memory to run picard JAR, in MB"
	picardSuffix: "Suffix for output file"
	strandSpecificity: "String to denote strand specificity for Picard"
    }

    String resultName = "~{outputFileNamePrefix}.~{picardSuffix}"
    
    command <<<
	java -Xmx~{picardMem}M \
	-jar ~{picardJarDir}picard.jar CollectRnaSeqMetrics \
	I=~{bamFile} \
	O=~{resultName} \
	STRAND_SPECIFICITY=~{strandSpecificity} \
	REF_FLAT=~{refFlat}
    >>>

    output {
	File result = "~{resultName}"
    }

    meta {
	output_meta: {
            result: "Text file with Picard output"
	}
    }
}
