#!/usr/bin/env nextflow

/*
 * wHLApper - Nextflow wrapper for HLA*LA (ONT)
 *
 * Input : one aligned BAM/CRAM, or a samplesheet (CSV/TSV):
 *           sample_id, bam (file or folder of BAMs to merge), [HLA-A_1, HLA-A_2, HLA-B_1, HLA-B_2]
 * Output: per-sample HLA*LA calls + QC, and cohort matrices in ${outdir}/summary
 */

nextflow.enable.dsl = 2

def helpMessage() {
    log.info """
    wHLApper - HLA*LA wrapper (ONT)

    Usage:
      nextflow run main.nf -profile docker --input samplesheet.csv --reference ref.fa \\
          --hla_reference GRCh38_full_analysis_set_plus_decoy_hla.fa --graph_dir /path/PRG_MHC_GRCh38_withIMGT

    Required:
      --input          Aligned .bam/.cram, or samplesheet .csv/.tsv
      --reference      Reference FASTA the BAMs were aligned to (.fai alongside)
      --graph_dir      HLA*LA graph folder (e.g. PRG_MHC_GRCh38_withIMGT)
      --hla_reference  HLA decoy FASTA reads are realigned to (.fai alongside),
                       e.g. GRCh38_full_analysis_set_plus_decoy_hla.fa

    Optional:
      --hla_reference_mmi  Prebuilt minimap2 index of --hla_reference (otherwise built on the fly)
      --sample_id      Sample name in single-BAM mode (default: BAM basename)
      --outdir         Output folder [${params.outdir}]
      --min_coverage   AverageCoverage must be > this [${params.min_coverage}]
      --min_q1         Q1 must be > this [${params.min_q1}]
      --skip_extract   Pass the whole BAM to HLA*LA instead of chr6/HLA contigs + unmapped
      --threads        Threads for HLA*LA [${params.threads}]

    Profiles: docker | singularity | conda  (optionally combined with: slurm)
    """.stripIndent()
}

/*
 * ---------------------------------------------------------------------------
 *  Modules
 * ---------------------------------------------------------------------------
 */

include { MERGE_BAMS; SAMTOOLS_INDEX; EXTRACT_HLA_REGION } from './module/samtools'
include { MINIMAP2_ALIGN                                 } from './module/minimap2'
include { HLALA; HLALA_QC                                } from './module/hlala'
include { COMPILE                                        } from './module/compile'

/*
 * ---------------------------------------------------------------------------
 *  Input parsing
 * ---------------------------------------------------------------------------
 */

def findIndex(f) {
    def s = f.toString()
    def candidates = s.endsWith('.cram')
        ? ["${s}.crai", s.replaceAll(/\.cram$/, '.crai')]
        : ["${s}.bai", s.replaceAll(/\.bam$/, '.bai')]
    return candidates.collect { file(it) }.find { it.exists() }
}

def parseInput(String input) {
    def inFile = file(input, checkIfExists: true)
    def name = inFile.name.toLowerCase()

    // Single BAM/CRAM
    if (name ==~ /.*\.(bam|cram)$/) {
        def id = params.sample_id ?: inFile.name.replaceAll(/\.(bam|cram)$/, '')
        return [[id: id, files: [inFile], index: findIndex(inFile), type: 'file', truth: [:]]]
    }

    if (!(name ==~ /.*\.(csv|tsv|txt)$/)) {
        error("--input must be a .bam/.cram file or a .csv/.tsv samplesheet: ${input}")
    }

    def sep = name.endsWith('.csv') ? ',' : '\t'
    def seen = [] as Set
    return inFile.splitCsv(header: true, sep: sep).collect { raw ->
        def row = raw.collectEntries { k, v -> [(k.toString().trim().toLowerCase()), v?.toString()?.trim()] }
        def id = row['sample_id'] ?: row['sample'] ?: row['id']
        def p  = row['bam'] ?: row['bam_path'] ?: row['path'] ?: row['bam_dir']
        if (!id || !p) error("Samplesheet row is missing sample_id or bam: ${raw}")
        if (!seen.add(id)) error("Duplicate sample_id in samplesheet: ${id}")

        def truth = [:]
        row.each { k, v ->
            def m = (k =~ /^hla-?([a-z0-9]+)_([12])$/)
            if (m.matches() && v) truth["HLA-${m.group(1).toUpperCase()}_${m.group(2)}".toString()] = v
        }

        def path = file(p, checkIfExists: true)
        if (path.isDirectory()) {
            def xams = path.listFiles().findAll { it.name ==~ /.*\.(bam|cram)$/ }.sort { it.name }
            if (!xams) error("No .bam/.cram files found in folder for ${id}: ${path}")
            if (xams.size() == 1) {
                return [id: id, files: xams, index: findIndex(xams[0]), type: 'file', truth: truth]
            }
            return [id: id, files: xams, index: null, type: 'dir', truth: truth]
        }
        return [id: id, files: [path], index: findIndex(path), type: 'file', truth: truth]
    }
}

/*
 * ---------------------------------------------------------------------------
 *  Workflow
 * ---------------------------------------------------------------------------
 */

workflow {
    main:
    // Top-level statements aren't allowed by the strict syntax parser, so --help
    // lives inside the workflow
    if (params.help) {
        helpMessage()
        exit 0
    }

    if (!params.input)       error("Missing --input (BAM/CRAM or samplesheet)")
    if (!params.reference)     error("Missing --reference (FASTA)")
    if (!params.graph_dir)     error("Missing --graph_dir (HLA*LA graph folder)")
    if (!params.hla_reference) error("Missing --hla_reference (HLA decoy FASTA, e.g. GRCh38_full_analysis_set_plus_decoy_hla.fa)")

    def reference     = file(params.reference, checkIfExists: true)
    def reference_fai = file("${params.reference}.fai", checkIfExists: true)
    def hla_reference = file(params.hla_reference, checkIfExists: true)
    def hla_fai       = file("${params.hla_reference}.fai", checkIfExists: true)
    def hla_index     = params.hla_reference_mmi ? file(params.hla_reference_mmi, checkIfExists: true) : hla_reference
    def graph_dir     = file(params.graph_dir, checkIfExists: true)

    def samples = parseInput(params.input)
    log.info "wHLApper: ${samples.size()} sample(s) -> ${params.outdir}"

    // Truth table (only when any sample has expected HLA types)
    def truthKeys = samples.collectMany { it.truth.keySet() as List }.unique()
    def truthLines = [(['sample_id'] + truthKeys).join('\t')] +
        samples.findAll { it.truth }.collect { s -> ([s.id] + truthKeys.collect { s.truth[it] ?: '' }).join('\t') }
    def ch_truth = truthKeys
        ? channel.of(truthLines.join('\n') + '\n').collectFile(name: 'truth.tsv')
        : channel.value(file("${projectDir}/assets/NO_TRUTH"))

    // Merge folders, index where needed
    def ch_samples = channel.fromList(samples).branch { s ->
        merge:     s.type == 'dir'
        indexed:   s.index != null
        unindexed: true
    }

    MERGE_BAMS(ch_samples.merge.map { s -> [[id: s.id], s.files] }, reference, reference_fai)
    SAMTOOLS_INDEX(ch_samples.unindexed.map { s -> [[id: s.id], s.files[0]] })

    def ch_aligned = MERGE_BAMS.out
        .mix(SAMTOOLS_INDEX.out)
        .mix(ch_samples.indexed.map { s -> [[id: s.id], s.files[0], s.index] })

    def ch_extracted = params.skip_extract
        ? ch_aligned
        : EXTRACT_HLA_REGION(ch_aligned, reference, reference_fai)

    // Realign to the HLA decoy reference expected by the HLA*LA graph
    MINIMAP2_ALIGN(ch_extracted, hla_index, reference, reference_fai)

    HLALA(MINIMAP2_ALIGN.out, graph_dir, hla_reference, hla_fai)
    HLALA_QC(HLALA.out.bestguess)
    COMPILE(HLALA_QC.out.collect(), ch_truth)

    onComplete:
    log.info "wHLApper finished: ${workflow.success ? 'OK' : 'FAILED'} - results in ${params.outdir}"
}
