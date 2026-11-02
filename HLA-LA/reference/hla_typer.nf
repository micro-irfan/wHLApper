#! /usr/bin/env nextflow

nextflow.enable.dsl=2
version = "Version 1.0"

include { samtools_index } from './common.nf'
include { filter_reads } from './hla_typing.nf'
include { minimap2_alignment } from './common.nf'
include { hla_typing_core } from './hla_typing.nf'
include { process_hla_la } from './hla_typing.nf'

workflow hlaTyper {
    take:
        bamIdx
        refFasta
        refIdx

    emit:
        results

    main:
        refFasHLA = file("${params.s3_LongRead_references}" + '/GRCh38_full_analysis_set_plus_decoy_hla.fa')
        refIdxHLA = file("${params.s3_LongRead_references}" + '/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai')
        refMMIHLA = file("${params.s3_LongRead_references}" + '/GRCh38_full_analysis_set_plus_decoy_hla.mmi')
        prgGraphDir = Channel.fromPath("${params.s3_LongRead_references}" + "/graph/PRG_MHC_GRCh38_withIMGT/PRG_MHC_GRCh38_withIMGT/")

        tmpBam = filter_reads(bamIdx, refFasta, refIdx)
        realignedCram = minimap2_alignment(tmpBam, refMMIHLA, refFasHLA, refIdxHLA)
        realignedCramIdx = samtools_index(realignedCram)

        results = hla_typing_core(realignedCramIdx.combine(prgGraphDir), refFasHLA, refIdxHLA)
        summary = process_hla_la(results)
    
}

// NXF_VER=22.10.0 nextflow run main.nf --threads 8 --input_list test.txt -c nextflow.config --ecr_repo 682905754632.dkr.ecr.ap-southeast-1.amazonaws.com