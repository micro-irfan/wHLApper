#! /usr/bin/env nextflow

nextflow.enable.dsl=2

/*
 * ---------------------------------------------------------------------------
 *  samtools: merge, index and HLA region extraction
 * ---------------------------------------------------------------------------
 */

process MERGE_BAMS {
    tag "${meta.id}"
    label 'hlala'

    input:
        tuple val(meta), path(xams, stageAs: 'input/*')
        path reference
        path reference_fai

    output:
        tuple val(meta), path("${meta.id}.merged.bam"), path("${meta.id}.merged.bam.bai")

    script:
    """
    samtools merge -@ ${task.cpus} --reference ${reference} -o ${meta.id}.merged.bam input/*
    samtools index -@ ${task.cpus} ${meta.id}.merged.bam
    """
}

process SAMTOOLS_INDEX {
    tag "${meta.id}"
    label 'hlala'

    input:
        tuple val(meta), path(xam)

    output:
        tuple val(meta), path(xam), path("${xam}.{bai,crai}")

    script:
    """
    samtools index -@ ${task.cpus} ${xam}
    """
}

// Keep chr6, chr6 alt and HLA-* contigs plus unmapped reads (adapted from filter_reads)
process EXTRACT_HLA_REGION {
    tag "${meta.id}"
    label 'hlala'

    input:
        tuple val(meta), path(xam), path(xamIdx)
        path reference
        path reference_fai

    output:
        tuple val(meta), path("${meta.id}.hla.bam"), path("${meta.id}.hla.bam.bai")

    script:
    def cram_input = xam.name.endsWith('.cram') ? "-T ${reference}" : ""
    """
    contigs=\$(samtools view -H ${cram_input} ${xam} \\
        | awk -F'\\t' '\$1=="@SQ" { for (i = 2; i <= NF; i++) if (\$i ~ /^SN:/) { sn = substr(\$i, 4);
              if (sn == "chr6" || sn == "6" || sn ~ /^chr6_.*_alt\$/ || sn ~ /^HLA-/) print sn } }' \\
        | tr '\\n' ' ')
    if [ -z "\$contigs" ]; then
        echo "No chr6/HLA contigs found in ${xam} header" >&2
        exit 1
    fi

    samtools view -@ ${task.cpus} -b -F 4 ${cram_input} -o tmp_mapped.bam ${xam} \$contigs
    samtools view -@ ${task.cpus} -b -f 4 ${cram_input} -o tmp_unmapped.bam ${xam}
    samtools merge -@ ${task.cpus} -o ${meta.id}.hla.bam tmp_mapped.bam tmp_unmapped.bam
    samtools index -@ ${task.cpus} ${meta.id}.hla.bam
    rm -f tmp_mapped.bam tmp_unmapped.bam
    """
}
