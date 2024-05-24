
process porechop {
    input:
        path input_fastq
    
    output:
        path "${params.experiment_name}_trim.fastq.gz", emit: fastq
        path "porechop.log", emit: logs

    script:
        """
        eval "\$(micromamba shell hook --shell=bash)"
        
        # Run the script
        porechop_abi -abi -i $input_fastq -o chopped.fastq.gz > porechop.log

        # trim T's at the beginning of the reads
        cutadapt -j 0 -g ^T{30} -e 0.2 --poly-a -m ${params.min_len} -M ${params.max_len} -o ${params.experiment_name}_trim.fastq.gz chopped.fastq.gz
        """
}