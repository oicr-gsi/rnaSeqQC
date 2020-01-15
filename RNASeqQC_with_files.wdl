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
    
}


task bamqc {

    input {
	File bamFile
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
