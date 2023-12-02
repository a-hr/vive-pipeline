
process count_mapping_reads {
    publishDir "${params.output_dir}"

    input:
        path original_fastq
        path target_bam
        path human_bam
        path unknown_bam
    
    output:
        path "reads_table.txt"

    script:
        """
        total_reads=$(zcat $original_fastq | wc -l)
        total_reads=$((total_reads / 4))

        samtools index -M *.bam

        target_reads=$(samtools view -c $target_bam)
        human_reads=$(samtools view -c $human_bam)
        unknown_reads=$(samtools view -c $unknown_bam)

        echo "total_reads: $total_reads" > reads_table.txt
        echo "-------------------------" >> reads_table.txt
        echo "target_reads: $target_reads" >> reads_table.txt
        echo "human_reads: $human_reads" >> reads_table.txt
        echo "unknown_reads: $unknown_reads" >> reads_table.txt
        """
}