#! /usr/bin/env nextflow

nextflow.enable.dsl=2
version = "Version 1.0"

process rename_bam {
    tag "LongRead_${run_id}_${sample_id}_rename_bam"
    
    input:
        tuple val(run_id), val(sample_id), path(xam)

    output:
        tuple val(run_id), val(sample_id), path("*_1hr.cram")

    shell:
    """
    mv ${xam} ${sample_id}_1hr.cram
    """

}

process samtools_merge {
    tag "LongRead_${run_id}_${sample_id}_merge_bam"
    cpus 8
    memory 4.GB
    container "${params.ecr_repo}/nala/minimap2:2.22-samtools"
    publishDir "${params.outdir}/merge/", mode:"copy"
    
    input:
        tuple val(run_id), val(sample_id), path(xam1), path(xai1), path(xam2), path(xai2)
        file ref 
        file refIdx
    
    output: 
        tuple val(run_id), val(sample_id), path("*.merged.cram"), path("*.merged.cram.crai")
    
    script:

    def cram_params = xam1.name.endsWith('.cram') ? "--reference ${ref}" : ""

    """ 
    samtools merge -@ ${task.cpus} -c ${cram_params} -O cram tmp.cram ${xam1} ${xam2}
    samtools addreplacerg -r 'ID:${sample_id}' -r 'SM:${sample_id}' -o ${sample_id}.merged.cram tmp.cram
    samtools index -@ ${task.cpus} ${sample_id}.merged.cram
    """
}

workflow merge_bam {
    take:
        alignedTAS
        alignedWGS
        ref
        refIdx

    emit:
        merge_out

    main:
        merge_out = samtools_merge(alignedTAS.combine(alignedWGS, by:[0,1]), ref, refIdx)

}

process cram_to_bam {
    tag "LongRead_${run_id}_${sample_id}_cram_to_bam"
    cpus 8
    memory 4.GB
    container "${params.ecr_repo}/nala/minimap2:2.22-samtools"
    
    input:
        tuple val(run_id), val(sample_id), path(cram), path(crai)
        file ref
        file refIdx

    
    output:
        tuple val(run_id), val(sample_id), path("${sample_id}.bam"), path("${sample_id}.bam.bai")
    
    script:
    """
    samtools view -@ ${task.cpus} -b -T ${ref} ${cram} -o ${sample_id}.bam
    samtools index -@ ${task.cpus} ${sample_id}.bam
    """
} 


process samtools_index {
    tag "LongRead_${run_id}_${sample_id}_samtools_index"
    cpus 4
    memory 2.GB
    container "${params.ecr_repo}/nala/minimap2:2.22-samtools"
    
    input:
        tuple val(run_id), val(sample_id), path(xam)
    
    output:
        tuple val(run_id), val(sample_id), path(xam), path("${xam}.*ai")
    
    script:
    """
    samtools index -@ ${task.cpus} $xam
    """
} 

process minimap2_alignment {
    tag "LongRead_${run_id}_${sample_id}_minimap2_alignment"
    cpus params.threads
    memory {
        refMMI.size() > 1e9 ? "31 GB" : "11 GB"
    }
    container "${params.ecr_repo}/nala/minimap2:2.22-samtools"
    publishDir "${params.outdir}/alignment/", mode:"copy"
    
    input:
        tuple val(run_id), val(sample_id), path(bam)
        path refMMI
        path refFasta
        path refIdx
        
    output:
        tuple val(run_id), val(sample_id), path("*.cram"), path("*cram.crai")

    script:
        int sorting_threads = Math.min((task.cpus / 3) as int, 3)
        int mapping_threads = params.threads - sorting_threads

    """
    samtools fastq -@ 1 -T '*' ${bam} \
        | minimap2 -t ${mapping_threads} -ayL -x map-ont ${refMMI} - \
        | samtools sort -@ ${sorting_threads - 1} - \
        | samtools view -T ${refFasta} -C -o ${sample_id}.cram -

    samtools index ${sample_id}.cram
    """
}

workflow alignment {
    take:
        refMMI
        refFa
        refFaIdx
        ubam

    emit:
        bam

    main:
        bam = minimap2_alignment(ubam, refMMI, refFa, refFaIdx)
}

