
process summary_secondary {
    publishDir "${params.output_dir}", mode: 'copy'

    input:
        path original_fastq
        path trimmed_fastq
        path target_bam
        path target_unmapped_fastq
        path seconday_bam
        path unknown_fastq
    
    output:
        path "summary.txt"

    script:
        """
        if [[ $original_fastq == *.gz ]]
        then
            total_reads=\$(zcat $original_fastq | wc -l)
        else
            total_reads=\$(cat $original_fastq | wc -l)
        fi
        total_reads=\$((total_reads / 4))

        trim_reads=\$(zcat $trimmed_fastq | wc -l)
        trim_reads=\$((trim_reads / 4))

        for bam in *.bam
        do
            samtools index \$bam
        done

        target_reads=\$(samtools view -c $target_bam)

        target_unmapped_reads=\$(zcat $target_unmapped_fastq | wc -l)
        target_unmapped_reads=\$((target_unmapped_reads / 4))

        secondary_reads=\$(samtools view -c $seconday_bam)

        unknown_reads=\$(zcat $unknown_fastq | wc -l)
        unknown_reads=\$((unknown_reads / 4))

        echo "total_reads: \$total_reads" > summary.txt
        echo "-------------------------" >> summary.txt
        echo "trim_reads: \$trim_reads" >> summary.txt
        echo "target_reads: \$target_reads" >> summary.txt
        echo "target_unmapped_reads: \$target_unmapped_reads" >> summary.txt
        echo "secondary_reads: \$secondary_reads" >> summary.txt
        echo "unknown_reads: \$unknown_reads" >> summary.txt
        """
}

process summary_base {
        publishDir "${params.output_dir}", mode: 'copy'

    input:
        path original_fastq
        path trimmed_fastq
        path target_bam
        path target_unmapped_fastq
    
    output:
        path "summary.txt"

    script:
        """
        if [[ $original_fastq == *.gz ]]
        then
            total_reads=\$(zcat $original_fastq | wc -l)
        else
            total_reads=\$(cat $original_fastq | wc -l)
        fi
        total_reads=\$((total_reads / 4))

        trim_reads=\$(zcat $trimmed_fastq | wc -l)
        trim_reads=\$((trim_reads / 4))

        samtools index $target_bam
        target_reads=\$(samtools view -c $target_bam)

        target_unmapped_reads=\$(zcat $target_unmapped_fastq | wc -l)
        target_unmapped_reads=\$((target_unmapped_reads / 4))

        echo "total_reads: \$total_reads" > summary.txt
        echo "-------------------------" >> summary.txt
        echo "trim_reads: \$trim_reads" >> summary.txt
        echo "target_reads: \$target_reads" >> summary.txt
        echo "target_unmapped_reads: \$target_unmapped_reads" >> summary.txt
        """
}