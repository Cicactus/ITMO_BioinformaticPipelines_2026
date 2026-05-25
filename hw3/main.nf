nextflow.enable.dsl = 2

// INPUTS
params.sra_id   = "DRR030302"
params.reads    = null        // --reads 'path/*_{1,2}.fastq' for local files
params.ref      = null        // --ref 'path/to/ref.fa' for provided ref
params.outdir   = "results"

include { DOWNLOAD_SRA; TRIMMING; ASSEMBLY; INDEX_REF; MAPPING; PLOT_COVERAGE } from './modules/processes.nf'
include { FASTQC as FASTQC_RAW } from './modules/processes.nf'
include { FASTQC as FASTQC_TRIM } from './modules/processes.nf'
include { BCFTOOLS_MPILEUP } from './modules/bcftools/mpileup/main'

// WORKFLOW
workflow {
    // Data
    read_ch = params.reads ? 
        Channel.fromFilePairs(params.reads, checkIfExists: true) : 
        DOWNLOAD_SRA(params.sra_id)

    // Quality Control
    FASTQC_RAW(read_ch, "raw")

    // Trimming
    TRIMMING(read_ch)
    FASTQC_TRIM(TRIMMING.out, "trimmed")

    // Assembly
    if (params.ref) {
        ref_ch = Channel.fromPath(params.ref, checkIfExists: true).collect()
    } else {
        ASSEMBLY(TRIMMING.out)
        ref_ch = ASSEMBLY.out.map { it[1] }.collect()
    }

    // Downstream Analysis
    MAPPING(TRIMMING.out, ref_ch)
    
    // Variant Calling
    mpileup_feed_ch = MAPPING.out.map { sample_id, bam, bai ->
        return [ [id: sample_id], [bam, bai] ]
    }

    ref_feed_ch = ref_ch.map { fasta ->
        return [ [id: fasta.baseName], fasta ]
    }

    BCFTOOLS_MPILEUP(mpileup_feed_ch, ref_feed_ch, false)

    // Plotting the
    PLOT_COVERAGE(MAPPING.out)
}