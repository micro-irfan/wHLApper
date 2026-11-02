#! /usr/bin/env nextflow

nextflow.enable.dsl=2

process filter_reads {
    tag "LongRead_${run_id}_${sample_id}_filter_reads"
    cpus 8
    memory 16.GB
    container "${params.ecr_repo}/nala/minimap2:2.22-samtools"
    
    input: 
        tuple val(run_id), val(sample_id), path(xam), path(xamIdx)
        path refFasta
        path refIdx
    
    output:
        tuple val(run_id), val(sample_id), path("tmp.merge.bam")

    script:
    def cram_input = xam.name.endsWith('.cram') ? "-C -T ${refFasta}" : ""
    """
    samtools view -@ ${task.cpus} -o tmp_mapped.bam ${cram_input} ${xam} chr6
    samtools view -@ ${task.cpus} -f 4 -o tmp_unmapped.bam ${cram_input} ${xam}
    samtools merge -@ ${task.cpus} -o tmp.merge.bam tmp_mapped.bam tmp_unmapped.bam 
    """
}


process process_hla_la {
    tag "LongRead_${run_id}_${sample_id}_process_hla_la"
    publishDir "${params.outdir}/${sample_id}/hla-la/", mode:"copy"
    container "${params.ecr_repo}/nala/longreadpy:${params.longreadpyver}"
    cpus 1

    input:
        tuple val(run_id), val(sample_id), path(hla_results)

    output:
        tuple val(run_id), val(sample_id), path("${sample_id}.hla_typing_results.txt")

    script:
    """
    python /opt/process_hla_results.py -s ${sample_id} -r ${hla_results}
    """
}


// Create a process for HLA-LA
process hla_typing_core {
    tag "LongRead_${run_id}_${sample_id}_hla_typing"
    cpus params.threads
    memory 63.GB
    container "682905754632.dkr.ecr.ap-southeast-1.amazonaws.com/nala/hla-la:v1.0.3"
    publishDir "${params.outdir}/${sample_id}/hla-la/", mode:"copy"

    input:
        tuple val(run_id), val(sample_id), path(cram), path(cramFileIndex), path(prgGraphDir)
        path(refFasta)
        path(refFastaIndex)

    output:
        tuple val(run_id), val(sample_id), path("${sample_id}_results/${sample_id}/hla/R1_bestguess_G.txt")

    script:
    """
    mkdir -p ${sample_id}_results
    
    HLA-LA.pl --BAM ${cram} \
              --customGraphDir ${prgGraphDir} \
              --sampleID ${sample_id} \
              --maxThreads ${task.cpus} \
              --workingDir ${sample_id}_results \
              --longReads ont2d \
              --samtools_T ${refFasta}       
    """
}

