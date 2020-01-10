version 1.0

workflow RNASeqQC {

    input {
	File bamFile
	File bwaRef
	File collationScript
	String picardJarDir
	File refFlat
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
	String resultsJSON = collate.collatedJSON
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

    File resultFile = "bamqc.json"
    
    command <<<
	write_fast_metrics.py \
	-b ~{bamFile} \
	-o ~{resultFile} \
    >>>

    output {
	File result = resultFile
    }
}

task bamToFastq {

    input {
	File bamFile
    }

    File allFastq = "all.fastq"
    File R1 = "R1.fastq"
    File R2 = "R2.fastq"
    
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

    File resultFile = "contamination_summary.txt"

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
	samtools flagstat - > ~{resultFile}
    >>>

    output {
	File result = resultFile
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

    File resultFile = "RNASeqQC.json"
    
    command <<<
	~{collationScript} \
	--bamqc ~{bamqc} \
	--contam ~{contam} \
	--picard ~{picard} \
	--unique-reads ~{uniqueReads} \
	--out ~{resultFile}
    >>>
    
    output {
	File collatedJSON=resultFile
    }
}

task countUniqueReads {

    input {
	File bamFile
    }

    File resultFile = "unique_reads.txt"

    command <<<
	samtools view -F 256 ~{bamFile} \
	| wc -l \
	> ~{resultFile}
    >>>
    
    output {
	File result = resultFile
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

    File resultFile = "CollectRNASeqMetrics.txt"
    
    command <<<
	java -Xmx~{picardMem}M \
	-jar ~{picardJarDir}picard.jar CollectRnaSeqMetrics \
	I=~{bamFile} \
	O=~{resultFile} \
	STRAND_SPECIFICITY=~{strandSpecificity} \
	REF_FLAT=~{refFlat}
    >>>

    output {
	File result = resultFile
    }
}

