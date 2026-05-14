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
    read_ch = (params.reads ? Channel.fromFilePairs(params.reads, checkIfExists: true) : DOWNLOAD_SRA(params.sra_id)).first()
    
    // Raw QC
    FASTQC_RAW(read_ch, "raw")

    // Trimming
    TRIMMING(read_ch)

    // Trimmed QC
    FASTQC_TRIM(TRIMMING.out, "trimmed")

    // Assembly 
    def ref_input_ch
    if (params.ref) {
        ref_input_ch = Channel.fromPath(params.ref, checkIfExists: true)
    } else {
        ASSEMBLY(TRIMMING.out)
        ref_input_ch = ASSEMBLY.out.map { it[1] }
    }
    ref_ch = ref_input_ch.first()
    index_ch = INDEX_REF(ref_ch)

    // Mapping
    MAPPING(TRIMMING.out, ref_ch)

    VARIANT_CALLING(MAPPING.out, ref_ch)

    // Plotting
    PLOT_COVERAGE(MAPPING.out)
}