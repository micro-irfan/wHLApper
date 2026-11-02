#! /usr/bin/env nextflow

nextflow.enable.dsl=2

/*
 * ---------------------------------------------------------------------------
 *  minimap2: realign reads to the HLA decoy reference
 * ---------------------------------------------------------------------------
 */

// Reads -> FASTQ (keeping tags) -> minimap2 map-ont -> sorted BAM (adapted from minimap2_alignment)
process MINIMAP2_ALIGN {
    tag "${meta.id}"
    label 'hlala'

    input:
        tuple val(meta), path(xam), path(xamIdx)
        path index          // .mmi or FASTA of the HLA decoy reference
        path reference      // reference the input reads are aligned to (for CRAM decoding)
        path reference_fai

    output:
        tuple val(meta), path("${meta.id}.realigned.bam"), path("${meta.id}.realigned.bam.bai")

    script:
    def cram_input    = xam.name.endsWith('.cram') ? "--reference ${reference}" : ""
    def sort_threads  = Math.max(1, Math.min((task.cpus / 4) as int, 4))
    def map_threads   = Math.max(1, task.cpus - sort_threads)
    """
    samtools fastq -@ 1 -T '*' ${cram_input} ${xam} \\
        | minimap2 -t ${map_threads} -ayL -x map-ont ${index} - \\
        | samtools sort -@ ${sort_threads} -o ${meta.id}.realigned.bam -

    samtools index -@ ${task.cpus} ${meta.id}.realigned.bam
    """
}
