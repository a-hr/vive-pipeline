
process merge_fastqs {
    input:
        path fastqs
    output:
        path "merged.fastq.gz"
    script:

    """
    first_file=\$(basename "${fastqs[0]}")
    extension=\${first_file##*.}

    if [ \$extension == 'gz' ]
    then
        cat ${fastqs} > merged.fastq.gz
    else
        cat ${fastqs} | gzip > merged.fastq.gz
    fi
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
