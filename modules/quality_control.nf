process fastqc {
    input:
        path fastq
    output:
        path "fastqc_output/"

    """
    mkdir fastqc_output
    fastqc $fastq \\
        -o fastqc_output/ \\
        -t ${task.cpus}
    """
}

process multiqc {
    publishDir "${params.output_dir}"

    input:
        path porechop_logs
        path fastqc_logs
    output:
        path "*.html"
    
    """
    multiqc . -n ${params.run_name}_multiqc_report.html
    """
}
