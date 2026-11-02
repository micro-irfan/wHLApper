#! /usr/bin/env nextflow

nextflow.enable.dsl=2

/*
 * ---------------------------------------------------------------------------
 *  Cohort summary matrices
 * ---------------------------------------------------------------------------
 */

process COMPILE {
    label 'python'
    publishDir "${params.outdir}/summary", mode: 'copy'

    input:
        path(qc_tables)
        path(truth)

    output:
        path("*.tsv", includeInputs: false)

    script:
    def truth_arg = truth.name != 'NO_TRUTH' ? "--truth ${truth}" : ''
    """
    compile_hla.py --qc ${qc_tables} ${truth_arg} --outdir .
    """
}
