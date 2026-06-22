nextflow.enable.dsl = 2

// INPUTS
params.sra_id   = "DRR030302"
params.reads    = null        // --reads 'path/*_{1,2}.fastq' for local files
params.ref      = null        // --ref 'path/to/ref.fa' for provided ref
params.outdir   = "results"

// PROCESSES
process DOWNLOAD_SRA {
    tag "$sra_id"
    publishDir "${params.outdir}/raw_data", mode: 'copy'

    input:
    val sra_id

    output:
    tuple val(sra_id), path("${sra_id}_{1,2}.fastq")

    script:
    """
    PREFIX=\$(echo ${sra_id} | cut -c 1-6)
    
    URL1="https://ftp.sra.ebi.ac.uk/vol1/fastq/\$PREFIX/${sra_id}/${sra_id}_1.fastq.gz"
    URL2="https://ftp.sra.ebi.ac.uk/vol1/fastq/\$PREFIX/${sra_id}/${sra_id}_2.fastq.gz"

    echo "--- Downloading from: \$URL1"

    wget --no-check-certificate \$URL1 -O ${sra_id}_1.fastq.gz
    wget --no-check-certificate \$URL2 -O ${sra_id}_2.fastq.gz

    gunzip ${sra_id}_1.fastq.gz
    gunzip ${sra_id}_2.fastq.gz
    """
}

process FASTQC_RAW {
    tag "$sample_id"
    publishDir "${params.outdir}/qc", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    path "*.html"

    script:
    """
    fastqc ${reads}
    """
}

process TRIMMING {
    tag "$sample_id"
    publishDir "${params.outdir}/trimmed", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_R{1,2}_trim.fq")

    script:
    """
    fastp -i ${reads[0]} -I ${reads[1]} \
          -o ${sample_id}_R1_trim.fq -O ${sample_id}_R2_trim.fq \
          --html ${sample_id}_fastp.html
    """
}

process FASTQC_TRIM {
    tag "$sample_id"
    publishDir "${params.outdir}/qc_trim", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    path "*.html"

    script:
    """
    fastqc ${reads}
    """
}

process ASSEMBLY {
    tag "$sample_id"
    publishDir "${params.outdir}/assembly", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("scaffolds.fasta")

    script:
    """
    spades.py -1 ${reads[0]} -2 ${reads[1]} -o . -t ${task.cpus} --only-assembler
    """
}

process MAPPING {
    tag "$sample_id"
    publishDir "${params.outdir}/mapping", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)
    path ref

    output:
    tuple val(sample_id), path("aligned.bam"), path("aligned.bam.bai")

    script:
    """
    bwa index ${ref}
    bwa mem -t ${task.cpus} ${ref} ${reads[0]} ${reads[1]} | \
    samtools view -Sb - | samtools sort -o aligned.bam
    samtools index aligned.bam
    """
}

process PLOT_COVERAGE {
    tag "$sample_id"
    publishDir "${params.outdir}/plots", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    path "${sample_id}_coverage.png"

    script:
    """
    #!/usr/bin/env python3
    import pandas as pd
    import matplotlib.pyplot as plt
    import subprocess
    import os

    with open("coverage.txt", "w") as f:
        subprocess.run(["samtools", "depth", "${bam}"], stdout=f, check=True)

    if os.path.exists('coverage.txt') and os.path.getsize('coverage.txt') > 0:
        df = pd.read_csv('coverage.txt', sep='\t', names=['Ref', 'Pos', 'Depth'])

        longest_scaffold = df['Ref'].value_counts().idxmax()
        df_filtered = df[df['Ref'] == longest_scaffold].copy()

        plt.figure(figsize=(12, 5))
        plt.fill_between(df_filtered['Pos'], df_filtered['Depth'], color="skyblue", alpha=0.4, label='Depth')
        plt.plot(df_filtered['Pos'], df_filtered['Depth'], color="Slateblue", lw=1)
        
        avg_cov = df_filtered['Depth'].mean()
        plt.axhline(y=avg_cov, color='red', linestyle='--', label=f'Avg Depth: {avg_cov:.1f}x')
        
        plt.title(f"Coverage Plot: {longest_scaffold} (${sample_id})")
        plt.xlabel("Position on Scaffold")
        plt.ylabel("Read Depth")
        plt.legend()
        plt.grid(axis='y', alpha=0.3)
        
        plt.savefig("${sample_id}_coverage.png")
    else:
        open("${sample_id}_coverage.png", 'w').close()
        print("!!! No coverage data found to plot.")
    """
}

// WORKFLOW
workflow {
    if (params.reads) {
        read_ch = Channel.fromFilePairs(params.reads, checkIfExists: true)
    } else {
        read_ch = DOWNLOAD_SRA(params.sra_id)
    }

    // Raw QC
    FASTQC_RAW(read_ch)

    // Trimming
    TRIMMING(read_ch)

    // Trimmed QC
    FASTQC_TRIM(TRIMMING.out) 

    // Assembly 
    if (params.ref) {
        ref_ch = Channel.fromPath(params.ref)
    } else {
        ASSEMBLY(TRIMMING.out)
        ref_ch = ASSEMBLY.out.map { it[1] } 
    }

    // Mapping
    MAPPING(TRIMMING.out, ref_ch)

    // Plotting
    PLOT_COVERAGE(MAPPING.out)
}