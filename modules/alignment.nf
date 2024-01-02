
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
        path "${out_name}_unmapped.bam", emit: "unmapped_bam"
        path "${out_name}_mapped.bam.bai", emit: "mapped_bai"
        path "${out_name}_unmapped.bam.bai", emit: "unmapped_bai"

    script:
        """
        # Aligning fastq
        minimap2 -t 8 -ax map-ont -c --MD $reference $fastq \\
            | samtools view -b - \\
            | samtools sort -o aligned.bam -

        samtools index aligned.bam

        # Separating mapped and unmapped reads
        samtools view -b -F 4 -o ${out_name}_mapped.bam aligned.bam
        samtools view -b -f 4 -o ${out_name}_unmapped.bam aligned.bam

        # Indexing
        for bam in *.bam
        do
            samtools index \$bam
        done
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