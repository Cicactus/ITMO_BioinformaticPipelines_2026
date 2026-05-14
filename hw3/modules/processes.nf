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

process FASTQC {
    tag "$sample_id"
    publishDir "${params.outdir}/qc/${subdir}", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)
    val subdir

    output:
    path "*.html"

    script:
    """
    fastqc ${reads}
    """
}

process TRIMMING {
    tag "$sample_id"
    
    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_R{1,2}_trim.fq")

    script:
    def (r1, r2) = reads
    """
    fastp -i ${r1} -I ${r2} \
          -o ${sample_id}_R1_trim.fq -O ${sample_id}_R2_trim.fq \
          --html ${sample_id}_fastp.html
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

process INDEX_REF {
    tag "$ref"
    
    input: 
    path ref
    
    output: 
    path "${ref}*"
    
    script:
    """
    bwa index ${ref}
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

process VARIANT_CALLING {
    tag "$sample_id"
    publishDir "${params.outdir}/variants", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)
    path ref

    output:
    path "${sample_id}.vcf"

    script:
    """
    bcftools mpileup -f ${ref} ${bam} | bcftools call -mv -Ob -o ${sample_id}.vcf
    """
}