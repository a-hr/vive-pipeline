#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { porechop; rna2dna } from './modules/porechop'
include { fastqc; multiqc } from './modules/quality_control'
include { bam2fastq; align as target_align; align as human_align} from './modules/alignment'
include { count_mapping_reads } from './modules/final_stats'

// Define the input parameters
Channel
    .fromPath(params.input_fastq)
    .set { input_fastq }

Channel
    .fromPath(params.plasmid_ref_fa)
    .set { plasmid_ref_fa }

if (params.assess_contamination) {
    Channel
        .fromPath(params.human_ref_fa)
        .set { human_ref_fa }
}

// Define the workflow
workflow {
    // Run porechop to remove adapter sequences
    input_fastq = params.is_rna ? rna2dna(input_fastq) : input_fastq
    porechop(input_fastq)
    clean_fastq = porechop.out.fastq

    // Run fastqc
    fastqc(clean_fastq)

    // Align to the plasmid reference
    target_align(plasmid_ref_fa, clean_fastq, "target")

    // Assess contamination
    if (params.assess_contamination) {
        // Retrieve unmapped reads from the target alignment
        unmapped_fastq = target_align.out.unmapped_fastq

        // Align unmapped reads to the human reference
        human_align(human_ref_fa, unmapped_fastq, "human")

        count_mapping_reads(
            input_fastq,
            clean_fastq,
            target_align.out.mapped_bam,
            target_align.out.unmapped_fastq,
            human_align.out.mapped_bam,
            human_align.out.unmapped_fastq
        )
    }

    // Alternative splicing analysis
    // TODO
    
    // Run multiqc
    // logs = Channel.of([porechop.out.logs, fastqc.out])
    multiqc(porechop.out.logs, fastqc.out)
}
