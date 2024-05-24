
process seqAnalysis {
    publishDir "${params.output_dir}/splice_sites", mode: 'copy'
    
    input:
        path reference_fa
    output:
        path "*.tsv", emit: "tables"
        path "*.png", emit: "plots"
    script:
    """
    eval "\$(micromamba shell hook --shell=bash)"

    mkdir -p tmp/
    
    # extract the reference sequence
    seq=`grep -v '^>' ${reference_fa} | tr -d '\\n'`

    # get splice score tables
    get_splice_sites.R \\
        ${params.score5} \\
        ${params.score3} \\
        tmp/ \\
        \$seq \\
        ${params.me2x5} \\
        ${params.splice5seqs} \\
        ${params.splice_models}
    
    # generate plot
    plot_ss.R 5ss_scores.tsv 3ss_scores.tsv
    """
}