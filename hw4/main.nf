nextflow.enable.dsl = 2

// =========================================================================================
// PARAMETERS & IMPORTS
// =========================================================================================
params.input  = "data/sra_samplesheet.csv" // Points directly to our new SRA sheet
params.ref    = null                       
params.outdir = "results"

include { 
    DOWNLOAD_SRA; 
    TRIMMING; 
    ASSEMBLY; 
    MAPPING; 
    PLOT_COVERAGE;
    FILTER_VARIANTS 
} from './modules/processes.nf'

include { FASTQC as FASTQC_RAW } from './modules/processes.nf'
include { FASTQC as FASTQC_TRIM } from './modules/processes.nf'
include { BCFTOOLS_MPILEUP } from './modules/bcftools/mpileup/main'

// Named closures
def parseSraManifest(row) {
    def meta = [
        id    : row.id?.trim(),
        repeat: row.repeat?.trim() ?: '1',
        group : row.group?.trim() ?: 'unspecified_cohort'
    ]
    return tuple(row.sra_id?.trim(), meta)
}

// WORKFLOW
workflow {

    if (params.input) {
        //  Read all samples into one manifest channel
        ch_manifest = Channel.fromPath(params.input, checkIfExists: true)
            .splitCsv(header: true)
            .map { row -> parseSraManifest(row) } // Updated to call the function explicitly

        // Extract raw IDs to trigger background SRA internet fetches
        ch_sra_ids = ch_manifest.map { sra_id, meta -> sra_id }
        ch_fetched_reads = DOWNLOAD_SRA(ch_sra_ids) // Emits: [sra_id, [fastq_files]]

        // Join internet files back to your structural tracking metadata
        ch_staged_inputs = ch_fetched_reads
            .join(ch_manifest) // Matches on matching sra_id string
            .map { sra_id, reads, meta -> tuple(meta.id, reads) }

        // Retain metadata lookup cross-reference key
        ch_metadata_lookup = ch_manifest.map { sra_id, meta -> tuple(meta.id, meta) }

    } else {
        error "Please specify a valid SRA tracking sheet via --input"
    }

    // Upstream Processing
    FASTQC_RAW(ch_staged_inputs, "raw")
    TRIMMING(ch_staged_inputs)
    FASTQC_TRIM(TRIMMING.out, "trimmed")

    ch_mapping_input = Channel.empty()
    ref_feed_ch      = Channel.empty()

    if (params.ref) {
        ref_ch = Channel.fromPath(params.ref, checkIfExists: true).collect()
        ch_mapping_input = TRIMMING.out.combine(ref_ch)
        
        // Structure for BCFTOOLS: [ [id: baseName], fasta, [] ]
        ref_feed_ch = ref_ch.map { fasta -> tuple([id: fasta.baseName], fasta, []) }
    } else {
        ASSEMBLY(TRIMMING.out) // Emits: [sample_id, scaffolds.fasta]
        ch_mapping_input = TRIMMING.out.join(ASSEMBLY.out)
        
        ref_feed_ch = ASSEMBLY.out.map { sample_id, fasta -> tuple([id: 'denovo_ref'], fasta, []) }.first()
    }

    // Pass the correctly scoped channel into MAPPING
    MAPPING(ch_mapping_input)

    // Split / isolate your samples by group using .subMap()
    ch_mapped_with_meta = MAPPING.out
        .join(ch_metadata_lookup) // Matches on sample_id
        .map { sample_id, bam, bai, meta -> 
            return tuple(meta.subMap(['group']), bam, bai)
        }

    // Accumulate parallelized cohorts dynamically using groupTuple
    ch_grouped_for_analysis = ch_mapped_with_meta
        .groupTuple(by: 0)
        .map { grouping_key, bams, bais ->
            def meta = [id: grouping_key.group, group: grouping_key.group]
            return tuple(meta, [bams, bais].flatten(), [], [])
        }

    // Multi-sample Variant Calling & Converge to one channel
    BCFTOOLS_MPILEUP(ch_grouped_for_analysis, ref_feed_ch, false)

    // Run final variant analysis utilizing your newly added module
    FILTER_VARIANTS(BCFTOOLS_MPILEUP.out.vcf)

    PLOT_COVERAGE(MAPPING.out)
}