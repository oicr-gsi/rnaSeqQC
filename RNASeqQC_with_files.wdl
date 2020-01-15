version 1.0

workflow RNASeqQC {

    input {
	File bamFile
    }
    
    call bamqc {
	input:
	bamFile = bamFile
    }

    output {
	File result = bamqc.result
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
