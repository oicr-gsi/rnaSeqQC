version 1.0

workflow rnaSeqQC {

    input {
	File bamFile
	String outputFileNamePrefix = "rnaSeqQC"
	String strandSpecificity = "NONE"
    }

    parameter_meta {
	bamFile: "Input BAM file on which to compute QC metrics"
	outputFileNamePrefix: "Prefix for output files"
    }
    
    call bamqc {
	input:
	bamFile = bamFile,
	outputFileNamePrefix = outputFileNamePrefix
    }

    call bwaMem {
	input:
	bamFile = bamFile,
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
	outputFileNamePrefix = outputFileNamePrefix,
	strandSpecificity = strandSpecificity
    }

    call collate {
	input:
	bamqc = bamqc.result,
	contam = bwaMem.result,
	picard = picard.result,
	uniqueReads = countUniqueReads.result,
	outputFileNamePrefix = outputFileNamePrefix,
	strandSpecificity = strandSpecificity
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
	    name: "production-tools-python/1.1.2",
	    url: "https://bitbucket.oicr.on.ca/projects/GSI/repos/production-tools-python/"
	},
	{
	    name: "bwa/0.7.17",
	    url: "https://github.com/lh3/bwa/releases/download/v0.7.17/bwa-0.7.17.tar.bz2"
	},
	{
	    name: "bam-qc-metrics/0.2.5",
	    url: "https://github.com/oicr-gsi/bam-qc-metrics.git"
	}
	]
    }
}


task bamqc {

    input {
	File bamFile
	String outputFileNamePrefix
	String bamqcSuffix = "bamqc.json"
	String modules = "bam-qc-metrics/0.2.5"
	Int jobMemory = 16
	Int threads = 4
	Int timeout = 4
    }

    parameter_meta {
	bamFile: "Input BAM file of aligned rnaSeqQC data"
	outputFileNamePrefix: "Prefix for output file"
	bamqcSuffix: "Suffix for output file"
	modules: "required environment modules"
	jobMemory: "Memory allocated for this job"
	threads: "Requested CPU threads"
	timeout: "hours before task timeout"
    }

    String resultName = "~{outputFileNamePrefix}.~{bamqcSuffix}"
    
    command <<<
	write_fast_metrics.py \
	-b ~{bamFile} \
	-o ~{resultName} \
    >>>

    runtime {
	modules: "~{modules}"
	memory:  "~{jobMemory} GB"
	cpu:     "~{threads}"
	timeout: "~{timeout}"
    }

    output {
	File result = "~{resultName}"
    }

    meta {
	output_meta: {
            result: "JSON file containing BAMQC metrics"
	}
    }
}

task bwaMem {

    input {
	File bamFile
	String outputFileNamePrefix
	String refFasta
	String modules
	String contamSuffix = "contaminationBwaFlagstat.txt"
	Int threads = 4
	Int jobMemory = 16
	Int timeout = 4
    }

    parameter_meta {
	bamFile: "Input BAM file of aligned rnaSeqQC data"
	outputFileNamePrefix: "Prefix for output file"
	refFasta: "Path to human genome FASTA reference"
	modules: "required environment modules"
	contamSuffix: "Suffix for output file"
	jobMemory: "Memory allocated for this job"
	threads: "Requested CPU threads"
	timeout: "hours before task timeout"
    }

    String resultName = "~{outputFileNamePrefix}.~{contamSuffix}"
    String rrnaRefName = "human_all_rRNA.fasta"

    # $RNASEQQC_RIBOSOME_GRCH38_BWA_INDEX_ROOT in module rnaseqqc-ribosome-grch38-bwa-index

    command <<<
	set -e
	set -o pipefail
	samtools collate \
	-O \
	--reference ~{refFasta} \
	~{bamFile} \
	| samtools fastq - \
	| bwa mem \
	-M \
	-t 8 \
	-p \
	$RNASEQQC_RIBOSOME_GRCH38_BWA_INDEX_ROOT/~{rrnaRefName} \
	- \
	| \
	samtools flagstat - > ~{resultName}
    >>>

    output {
	File result = "~{resultName}"
    }

    runtime {
	modules: "~{modules}"
	memory:  "~{jobMemory} GB"
	cpu:     "~{threads}"
	timeout: "~{timeout}"
    }

    meta {
	output_meta: {
            result: "Text file with results of running 'samtools flagstat' on BWA output"
	}
    }
}

task collate {

    input {
	File bamqc
	File contam
	File picard
	File uniqueReads
	String outputFileNamePrefix
	String strandSpecificity
	String collatedSuffix = "collatedMetrics.json"
	String modules = "production-tools-python/2"
	Int jobMemory = 16
	Int threads = 4
	Int timeout = 4
    }

    parameter_meta {
	bamqc: "JSON output from bamqc task"
	contam: "Text output from ribosomal contamination check by bwaMem task"
	picard: "Text output from picard task"
	uniqueReads: "Text output from uniqueReads task"
	outputFileNamePrefix: "Prefix for output file"
	collatedSuffix: "Suffix for output file"
	modules: "required environment modules"
	jobMemory: "Memory allocated for this job"
	threads: "Requested CPU threads"
	timeout: "hours before task timeout"
    }

    String resultName = "~{outputFileNamePrefix}.~{collatedSuffix}"
    String collationScript = "rnaseqqc-collate"
    String strandOption = if strandSpecificity=="NONE" then "--no-strand-specificity " else ""

    command <<<
	~{collationScript} \
	--bamqc ~{bamqc} \
	--contam ~{contam} \
	--picard ~{picard} \
	--unique-reads ~{uniqueReads} \
	--out ~{resultName} \
	~{strandOption}
    >>>

    runtime {
	modules: "~{modules}"
	memory:  "~{jobMemory} GB"
	cpu:     "~{threads}"
	timeout: "~{timeout}"
    }

    output {
	File collatedResults="~{resultName}"
    }

    meta {
	output_meta: {
            collatedResults: "JSON file of collated rnaSeqQC output"
	}
    }
}

task countUniqueReads {

    input {
	File bamFile
	String outputFileNamePrefix
	String uniqueReadsSuffix = "uniqueReads.txt"
	String modules = "samtools/1.9"
	Int jobMemory = 16
	Int threads = 4
	Int timeout = 4
    }

    parameter_meta {
	bamFile: "Input BAM file of aligned rnaSeqQC data"
	outputFileNamePrefix: "Prefix for output file"
	uniqueReadsSuffix: "Suffix for output file"
	modules: "required environment modules"
	jobMemory: "Memory allocated for this job"
	threads: "Requested CPU threads"
	timeout: "hours before task timeout"
    }

    String resultName = "~{outputFileNamePrefix}.~{uniqueReadsSuffix}"

    command <<<
	samtools view -F 256 ~{bamFile} \
	| wc -l \
	> ~{resultName}
    >>>

    runtime {
	modules: "~{modules}"
	memory:  "~{jobMemory} GB"
	cpu:     "~{threads}"
	timeout: "~{timeout}"
    }
    
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
	String outputFileNamePrefix
	String refFlat
	String refFasta
	String modules
	Int picardMem=6000
	String picardSuffix = "picardCollectRNASeqMetrics.txt"
	String strandSpecificity="NONE"
	Int jobMemory = 64
	Int threads = 4
	Int timeout = 4
    }

    parameter_meta {
	bamFile: "Input BAM file of aligned rnaSeqQC data"
	outputFileNamePrefix: "Prefix for output file"
	refFlat: "Path to Picard flatfile reference"
	refFasta: "Path to human genome FASTA reference"
	modules: "Required environment modules"
	picardMem: "Memory to run picard JAR, in MB"
	picardSuffix: "Suffix for output file"
	strandSpecificity: "String to denote strand specificity for Picard"
	jobMemory: "Memory allocated for this job"
	threads: "Requested CPU threads"
	timeout: "hours before task timeout"
    }

    String resultName = "~{outputFileNamePrefix}.~{picardSuffix}"

    # Environment variables from modulefiles:
    # $PICARD_ROOT <- picard
    # $HG38_ROOT <- hg38
    # $HG38_REFFLAT_ROOT <- hg38-refflat

    # VALIDATION_STRINGENCY=SILENT prevents BAM parsing errors with the given REFERENCE_SEQUENCE

    # check if HG19_ROOT or HG38_ROOT variable is set by environment module
    command <<<
	java -Xmx~{picardMem}M \
	-jar $PICARD_ROOT/picard.jar CollectRnaSeqMetrics \
	I=~{bamFile} \
	O=~{resultName} \
	STRAND_SPECIFICITY=~{strandSpecificity} \
	REF_FLAT=~{refFlat} \
	REFERENCE_SEQUENCE=~{refFasta} \
	VALIDATION_STRINGENCY=SILENT
    >>>

    runtime {
	modules: "~{modules}"
	memory:  "~{jobMemory} GB"
	cpu:     "~{threads}"
	timeout: "~{timeout}"
    }

    output {
	File result = "~{resultName}"
    }

    meta {
	output_meta: {
            result: "Text file with Picard output"
	}
    }
}
