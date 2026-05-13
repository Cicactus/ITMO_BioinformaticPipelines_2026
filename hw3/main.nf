nextflow.enable.dsl = 2

// INPUTS
params.sra_id   = "DRR030302"
params.reads    = null        // --reads 'path/*_{1,2}.fastq' for local files
params.ref      = null        // --ref 'path/to/ref.fa' for provided ref
params.outdir   = "results"

include { DOWNLOAD_SRA; TRIMMING; ASSEMBLY; INDEX_REF; MAPPING; PLOT_COVERAGE; VARIANT_CALLING } from './modules/processes.nf'
include { FASTQC as FASTQC_RAW } from './modules/processes.nf'
include { FASTQC as FASTQC_TRIM } from './modules/processes.nf'

// WORKFLOW
workflow {
    read_ch = params.reads ? Channel.fromFilePairs(params.reads, checkIfExists: true) : DOWNLOAD_SRA(params.sra_id)

    // Raw QC
    FASTQC_RAW(read_ch, "qc_raw")

    // Trimming
    TRIMMING(read_ch)

    // Trimmed QC
    FASTQC_TRIM(TRIMMING.out, "qc_trimmed") 

    // Assembly 
    if (params.ref) {
        ref_file = file(params.ref)
    } else {
        ASSEMBLY(TRIMMING.out)
        ref_file = ASSEMBLY.out.map { it[1] }
    }

    ref_ch = Channel.fromPath(ref_file).first()
    index_ch = INDEX_REF(ref_ch).first()

    // Mapping
    MAPPING(TRIMMING.out, ref_ch)

    VARIANT_CALLING(MAPPING.out, ref_ch)

    // Plotting
    PLOT_COVERAGE(MAPPING.out)
}