
process align {
    publishDir "${params.output_dir}",
         saveAs: { filename -> "${out_name}_${filename}" }, 
         enabled: params.get_bams

    input:
        path reference
        path fastq
        val out_name

    output:
        path "mapped.bam", emit: "mapped_bam"
        path "unmapped.bam", emit: "unmapped_bam"
        path "mapped.bam.bai", emit: "mapped_bai"
        path "unmapped.bam.bai", emit: "unmapped_bai"

    script:
        """
        # Aligning fastq
        minimap2 -t 8 -ax map-ont -c --MD $reference $fastq \\
            | samtools view -b - \\
            | samtools sort -o aligned.bam -

        samtools index aligned.bam

        # Separating mapped and unmapped reads
        samtools view -b -F 4 -o mapped.bam aligned.bam
        samtools view -b -f 4 -o unmapped.bam aligned.bam

        # Indexing
        samtools index -M *.bam
        """
}

process bam2fastq {
    input:
        path bam
    output:
        path "sample.fastq.gz"

    script:
        """
        samtools fastq -n $bam | gzip > sample.fastq.gz
        """
}