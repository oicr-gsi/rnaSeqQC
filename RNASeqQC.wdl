version 1.0

workflow RNASeqQC {

    input {
	File bamPath
	File bwaRef
	File collationScript
	String picardJarDir
	File refFlat
    }
    
    call bamqc {
	input:
	bamPath = bamPath
    }

    call bamToFastq {
	input:
	bamPath = bamPath
    }

    call bwaMem {
	input:
	fastqR1=bamToFastq.fastqR1,
	fastqR2=bamToFastq.fastqR2,
	bwaRef=bwaRef
    }

    call countUniqueReads {
	input:
	bamPath = bamPath
    }
    
    call picard {
	input:
	bamPath = bamPath,
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
	File resultsJSON = collate.collatedJSON
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
	File bamPath
    }

    File resultName = "bamqc.json"
    
    command <<<
	write_fast_metrics.py \
	-b ~{bamPath} \
	-o ~{resultName} \
    >>>

    output {
	File result = "~{resultName}"
    }
}

task bamToFastq {

    input {
	File bamPath
    }

    File allFastq = "all.fastq"
    File R1 = "R1.fastq"
    File R2 = "R2.fastq"
    
    command <<<
	samtools bam2fq ~{bamPath} > ~{allFastq} && \
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

    File resultName = "contamination_summary.txt"

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

    File resultName = "RNASeqQC.json"
    
    command <<<
	~{collationScript} \
	--bamqc ~{bamqc} \
	--contam ~{contam} \
	--picard ~{picard} \
	--unique-reads ~{uniqueReads} \
	--out ~{resultName}
    >>>
    
    output {
	File collatedJSON="~{resultName}"
    }
}

task countUniqueReads {

    input {
	File bamPath
    }

    File resultName = "unique_reads.txt"

    command <<<
	samtools view -F 256 ~{bamPath} \
	| wc -l \
	> ~{resultName}
    >>>
    
    output {
	File result = "~{resultName}"
    }
}

task picard {

    input {
	File bamPath
	String picardJarDir
	File refFlat
	Int picardMem=6000
	File strandSpecificity="NONE"
    }

    File resultName = "CollectRNASeqMetrics.txt"
    
    command <<<
	java -Xmx~{picardMem}M \
	-jar ~{picardJarDir}picard.jar CollectRnaSeqMetrics \
	I=~{bamPath} \
	O=~{resultName} \
	STRAND_SPECIFICITY=~{strandSpecificity} \
	REF_FLAT=~{refFlat}
    >>>

    output {
	File result = "~{resultName}"
    }
}

