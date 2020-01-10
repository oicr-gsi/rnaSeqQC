version 1.0

workflow RNASeqQCTestSimple {

    input {
	File bamFile
    }
    
    call countUniqueReads {
	input:
	bamFile = bamFile
    }

    output {
	File result = countUniqueReads.result
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

