
process count_mapping_reads {
    publishDir "${params.output_dir}", mode: 'copy'

    input:
        path original_fastq
        path trimmed_fastq
        path target_bam
        path human_bam
        path unknown_bam
    
    output:
        path "reads_table.txt"

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

        target_reads=\$(samtools view -c -F 260 $target_bam)
        human_reads=\$(samtools view -c $human_bam)
        unknown_reads=\$(samtools view -c $unknown_bam)

        echo "total_reads: \$total_reads" > reads_table.txt
        echo "-------------------------" >> reads_table.txt
        echo "trim_reads: \$trim_reads" >> reads_table.txt
        echo "target_reads: \$target_reads" >> reads_table.txt
        echo "human_reads: \$human_reads" >> reads_table.txt
        echo "unknown_reads: \$unknown_reads" >> reads_table.txt
        """
}