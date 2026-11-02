#! /usr/bin/env nextflow

nextflow.enable.dsl=2

/*
 * ---------------------------------------------------------------------------
 *  HLA*LA typing and per-sample QC
 * ---------------------------------------------------------------------------
 */

process HLALA {
    tag "${meta.id}"
    label 'hlala'
    publishDir "${params.outdir}/${meta.id}", mode: 'copy', saveAs: { fn -> fn.endsWith('_hla') ? 'hla' : null }

    input:
        tuple val(meta), path(xam), path(xamIdx)
        path graph_dir
        path reference
        path reference_fai

    output:
        tuple val(meta), path("${meta.id}_hla/R1_bestguess_G.txt"), emit: bestguess
        path("${meta.id}_hla"), emit: hla_dir

    script:
    // HLA-LA (BamTools) expects <file>.bai / <file>.crai next to the alignment
    def idx_ext = xam.name.endsWith('.cram') ? 'crai' : 'bai'
    """
    [ -e ${xam}.${idx_ext} ] || ln -s ${xamIdx} ${xam}.${idx_ext}
    mkdir -p hlala_work

    HLA-LA.pl --BAM ${xam} \\
              --customGraphDir ${graph_dir} \\
              --sampleID ${meta.id} \\
              --maxThreads ${task.cpus} \\
              --workingDir \$PWD/hlala_work \\
              --longReads ont2d \\
              --samtools_T ${reference}

    mv hlala_work/${meta.id}/hla ${meta.id}_hla
    """
}

process HLALA_QC {
    tag "${meta.id}"
    label 'python'
    publishDir "${params.outdir}/${meta.id}", mode: 'copy'

    input:
        tuple val(meta), path(bestguess)

    output:
        path("${meta.id}.hla_qc.tsv")

    script:
    """
    hlala_qc.py -s ${meta.id} -r ${bestguess} --min-cov ${params.min_coverage} --min-q1 ${params.min_q1}
    """
}
