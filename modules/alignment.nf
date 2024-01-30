
process align {
    publishDir "${params.output_dir}/bams",
        enabled: params.get_bams,
        mode: 'copy'

    input:
        path reference
        path fastq
        val out_name

    output:
        path "${out_name}_mapped.bam", emit: "mapped_bam"
        path "${out_name}_unmapped.fastq.gz", emit: "unmapped_fastq"
        path "${out_name}_mapped.bam.bai", emit: "mapped_bai"

    script:
        """
        # Aligning fastq
        minimap2 -2 -t 8 -ax splice -c --MD --secondary=no $reference $fastq \\
            | samtools view -b - \\
            | samtools sort -o aligned.bam -

        samtools index aligned.bam

        # Separating mapped and unmapped reads
        samtools view -b -F 2820 \\
            -U unmapped.bam \\
            -o ${out_name}_mapped.bam aligned.bam

        # Indexing
        samtools index ${out_name}_mapped.bam

        # Converting unmapped reads to fastq
        samtools bam2fq unmapped.bam | gzip > ${out_name}_unmapped.fastq.gz
        """
}

process bam2fastq {
    input:
        path bam
    output:
        path "*.fastq.gz"

    script:
        """
        filename=\$(basename $bam .bam)
        samtools fastq -n $bam | gzip > \$filename.fastq.gz
        """
}