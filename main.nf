#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { merge_fastqs; rna2dna} from './modules/utils'
include { porechop } from './modules/porechop'
include { fastqc; multiqc } from './modules/quality_control'
include { bam2fastq; align as target_align; align as secondary_align} from './modules/alignment'
include { summary_secondary; summary_base } from './modules/final_stats'
include { seqAnalysis } from './modules/alternative_splicing'

// Define the input parameters
Channel
    .fromPath(params.input_fastqs + "/*.fastq*", checkIfExists: true)
    .set { input_fastq }

Channel
    .fromPath(params.plasmid_ref_fa, checkIfExists: true)
    .set { plasmid_ref_fa }

if (params.assess_secondary) {
    Channel
        .fromPath(params.secondary_ref_fa, checkIfExists: true)
        .set { secondary_ref_fa }
}

// Define the workflow
workflow {
    // Merge the input fastq files
    input_fastq | collect \
        | merge_fastqs \
        | set { merged_fastq }

    // Convert to DNA if RNA
    dna_fastq = params.is_rna ? rna2dna(merged_fastq) : merged_fastq

    // Run porechop to remove adapter sequences
    porechop(dna_fastq)
    clean_fastq = porechop.out.fastq

    // Run fastqc
    fastqc(clean_fastq)

    // Align to the plasmid reference
    target_align(plasmid_ref_fa, clean_fastq, "target")

    // Assess contamination
    if (params.assess_secondary) {
        // Retrieve unmapped reads from the target alignment
        unmapped_fastq = target_align.out.unmapped_fastq

        // Align unmapped reads to the human reference
        secondary_align(secondary_ref_fa, unmapped_fastq, "human")

        summary_secondary(
            dna_fastq,
            clean_fastq,
            target_align.out.mapped_bam,
            target_align.out.unmapped_fastq,
            secondary_align.out.mapped_bam,
            secondary_align.out.unmapped_fastq
        )
    } else {
        summary_base(
            dna_fastq,
            clean_fastq,
            target_align.out.mapped_bam,
            target_align.out.unmapped_fastq,
        )
    }

    // Alternative splicing analysis
    seqAnalysis(plasmid_ref_fa)
    
    // Run multiqc
    multiqc(porechop.out.logs, fastqc.out)
}
