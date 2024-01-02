
process porechop {
    input:
        path input_fastq
    
    output:
        path "trimmed.fastq.gz", emit: fastq
        path "porechop.log", emit: logs

    script:
        """
        # Run the script
        porechop_abi -abi -i $input_fastq -o chopped.fastq.gz > porechop.log

        # trim T's at the beginning of the reads
        cutadapt -j 0 -g ^T{30} -e 0.2 --poly-a -m ${params.min_len} -o trimmed.fastq.gz chopped.fastq.gz
        """
}

process rna2dna {
    input:
        path input_fastq
    
    output:
        path "dna.fastq.gz"

    script:
        """
        zcat ${input_fastq} | perl -pe 's/U/T/g if \$. % 4 == 2' | gzip -c > dna.fastq.gz
        """
}