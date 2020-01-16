version 1.0

workflow RNASeqQC {

    input {
	File bamFile
	File bwaRef
	File collationScript
	File refFlat
	String picardJarDir
    }
    
    call bamqc {
	input:
	bamFile = bamFile
    }

    call bamToFastq {
	input:
	bamFile = bamFile
    }

    call bwaMem {
	input:
	fastqR1=bamToFastq.fastqR1,
	fastqR2=bamToFastq.fastqR2,
	bwaRef=bwaRef
    }

    call countUniqueReads {
	input:
	bamFile = bamFile
    }
    
    call picard {
	input:
	bamFile = bamFile,
	refFlat = refFlat,
	picardJarDir = picardJarDir
    }

    call collate {
	input:
	collationScript = collationScript,
	bamqc = bamqc.result,
	contam = bwaMem.result,
	picard = picard.result,
	uniqueReads = countUniqueReads.result
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
    }

    parameter_meta {
	bamFile: "Input BAM file of aligned RNASeqQC data"
    }

    String resultName = "bamqc.json"
    
    command <<<
	write_fast_metrics.py \
	-b ~{bamFile} \
	-o ~{resultName} \
    >>>

    output {
	File result = "~{resultName}"
    }
}

task bamToFastq {

    input {
	File bamFile
    }

    parameter_meta {
	bamFile: "Input BAM file of aligned RNASeqQC data"
    }

    String allFastq = "all.fastq"
    String R1 = "R1.fastq"
    String R2 = "R2.fastq"
    
    command <<<
	samtools bam2fq ~{bamFile} > ~{allFastq} && \
	cat ~{allFastq} | grep '^@.*/1$' -A 3 --no-group-separator > ~{R1} && \
	cat ~{allFastq} | grep '^@.*/2$' -A 3 --no-group-separator > ~{R2}
    >>>

    output {
	File fastqR1 = "~{R1}"
	File fastqR2 = "~{R2}"
    }
}

task bwaMem {

    input {
	File fastqR1
	File fastqR2
	File bwaRef
    }

    parameter_meta {
	fastqR1: "FASTQ file for read 1"
	fastqR2: "FASTQ file for read 2"
	bwaRef: "Ribosomal reference file for alignment by BWA"
    }

    String resultName = "contamination_summary.txt"

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
}

task collate {

    input {
	File collationScript
	File bamqc
	File contam
	File picard
	File uniqueReads
    }

    parameter_meta {
	collationScript: "Python script to collate intermediate results"
	bamqc: "JSON output from bamqc task"
	contam: "Text output from ribosomal contamination check by bwaMem task"
	picard: "Text output from picard task"
	uniqueReads: "Text output from uniqueReads task"
    }

    String resultName = "RNASeqQC.json"
    
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
}

task countUniqueReads {

    input {
	File bamFile
    }

    parameter_meta {
	bamFile: "Input BAM file of aligned RNASeqQC data"
    }

    String resultName = "unique_reads.txt"

    command <<<
	samtools view -F 256 ~{bamFile} \
	| wc -l \
	> ~{resultName}
    >>>
    
    output {
	File result = "~{resultName}"
    }
}

task picard {

    input {
	File bamFile
	String picardJarDir
	File refFlat
	Int picardMem=6000
	String strandSpecificity="NONE"
    }

    parameter_meta {
	bamFile: "Input BAM file of aligned RNASeqQC data"
	picardJarDir: "Path of directory containing the Picard JAR"
	refFlat: "Flat reference file required by Picard"
	picardMem: "Memory to run picard JAR, in MB"
	strandSpecificity: "String to denote strand specificity for Picard"
    }

    String resultName = "CollectRNASeqMetrics.txt"
    
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
}
