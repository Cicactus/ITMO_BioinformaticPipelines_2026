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
    # FIX: Let Nextflow safely inject the string, and let tail extract the final digit
    LAST_DIGIT=\$(echo -n "${sra_id}" | tail -c 1)

    # Default to the standard 9-character flat layout
    URL1="https://ftp.sra.ebi.ac.uk/vol1/fastq/\$PREFIX/${sra_id}/${sra_id}_1.fastq.gz"
    URL2="https://ftp.sra.ebi.ac.uk/vol1/fastq/\$PREFIX/${sra_id}/${sra_id}_2.fastq.gz"

    # Use --spider to check if the flat URL exists. If it fails (404), fall back to the padded folder layout
    if ! wget --spider --no-check-certificate \$URL1 &>/dev/null; then
        URL1="https://ftp.sra.ebi.ac.uk/vol1/fastq/\$PREFIX/00\$LAST_DIGIT/${sra_id}/${sra_id}_1.fastq.gz"
        URL2="https://ftp.sra.ebi.ac.uk/vol1/fastq/\$PREFIX/00\$LAST_DIGIT/${sra_id}/${sra_id}_2.fastq.gz"
    fi

    # Perform actual downloads
    wget --no-check-certificate \$URL1 -O ${sra_id}_1.fastq.gz
    wget --no-check-certificate \$URL2 -O ${sra_id}_2.fastq.gz
    
    gunzip ${sra_id}_1.fastq.gz
    gunzip ${sra_id}_2.fastq.gz
    """

    stub:
    """
    touch ${sra_id}_1.fastq
    touch ${sra_id}_2.fastq
    echo "STUB: Mocked fastq pairs for SRA run ${sra_id}"
    """
}

process FASTQC {
    tag "$sample_id"
    publishDir { "${params.outdir}/qc/${subdir}" }, mode: 'copy'

    input:
    tuple val(sample_id), path(reads)
    val subdir

    output:
    path "*.html"

    script:
    """
    fastqc ${reads}
    """

    stub:
    """
    touch ${sample_id}_fastqc.html
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

    stub:
    """
    touch ${sample_id}_R1_trim.fq
    touch ${sample_id}_R2_trim.fq
    """
}

process ASSEMBLY {
    tag "$sample_id"
    
    input:
    tuple val(sample_id), path(reads)
    
    output:
    tuple val(sample_id), path("${sample_id}_scaffolds.fasta")

    script:
    """
    spades.py -1 ${reads[0]} -2 ${reads[1]} -o . -t 4 --only-assembler
    
    # FIX: Change '//' to '#' so Bash treats this as a comment
    # If real spades runs, we would rename the default output:
    mv scaffolds.fasta ${sample_id}_scaffolds.fasta
    """

    stub:
    """
    touch ${sample_id}_scaffolds.fasta
    """
}

process MAPPING {
    tag "$sample_id"

    input:
    // Receives the unified 3-element tuple: [sample_id, [reads1, reads2], scaffolds]
    tuple val(sample_id), path(reads), path(scaffolds)

    output:
    tuple val(sample_id), path("${sample_id}.bam"), path("${sample_id}.bam.bai")

    script:
    """
    # Index only this specific reference fasta
    bwa index ${scaffolds}

    # Map reads cleanly against its matched reference
    bwa mem -t 4 ${scaffolds} ${reads[0]} ${reads[1]} | \
    samtools view -Sb - | \
    samtools sort -o ${sample_id}.bam

    samtools index ${sample_id}.bam
    """

    stub:
    """
    touch ${sample_id}.bam
    touch ${sample_id}.bam.bai
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

    stub:
    """
    touch ${sample_id}_coverage.png
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

process FILTER_VARIANTS {
    tag "$meta.group"
    label 'process_low'
    publishDir "${params.outdir}/filtered_variants", mode: 'copy'

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.filtered.vcf.gz"), emit: vcf

    script:
    """
    bcftools filter \
        -e "QUAL < 20 || DP < 10" \
        -O z \
        -o ${meta.group}.filtered.vcf.gz \
        ${vcf}
    """

    stub:
    """
    touch ${meta.group}.filtered.vcf.gz
    echo "STUB EXECUTION: Filtered variants completed for group channel: ${meta.group}"
    """
}