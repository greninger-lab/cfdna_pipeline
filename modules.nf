// Adapter trim with Trimmomatic
process trim_reads { 
    // container "quay.io/biocontainers/fastp:0.20.1--h2e03b76_1"
    container "srynobio/fastp:latest"

	// Retry on fail at most three times 
    errorStrategy 'retry'
    maxRetries 1

    input:
        tuple file(R1), file(R2), val(reference) // from input_read_ch
        file ADAPTERS
    output: 
        tuple env(base), val(reference), file("*R1_37bp.fastq.gz"), file("*R2_37bp.fastq.gz"), file("*summary.csv") // into Trim_out_ch
    
    publishDir "${params.OUTDIR}37bp_trimmed", mode: 'copy', pattern: '*_37bp.fastq.gz'

    script:
    """
    #!/bin/bash

    base=\$(echo ${R1} | cut -f1 -d'_')
    fastp -w ${task.cpus} -i \$base*R1*.fastq.gz -I \$base*R2*.fastq.gz -o \$base'_R1_37bp.fastq.gz' -O \$base'_R2_37bp.fastq.gz' -b 37 -B 37 > fastp_out.txt 2>&1

    trimmed_reads=\$(grep "reads passed filter:" fastp_out.txt | cut -f2 -d":" | cut -f2 -d" ")

    echo "sample_name,raw_reads,trimmed_reads,human_reads,dedup_human_reads,virus_reads,dedup_virus_reads" > \$base'_summary.csv'
    num_r1_untrimmed=\$(gunzip -c ${R1} | wc -l)
    num_r2_untrimmed=\$(gunzip -c ${R2} | wc -l)
    num_untrimmed=\$((\$((num_r1_untrimmed + num_r2_untrimmed))/4))
    printf "\$base,\$num_untrimmed,\$trimmed_reads," >> \$base'_summary.csv'

    """
}

// Use bowtie2 to filter against human and virus
process align_samples {
    container "quay.io/biocontainers/bowtie2:2.4.4--py37h13ad519_0"

    // Retry on fail at most three times 
    errorStrategy 'retry'
    maxRetries 1

    input:
        tuple val(base), val(reference), file("${base}_R1_37bp.fastq.gz"), file("${base}_R2_37bp.fastq.gz"), file("${base}_summary.csv")// from Trim_out_ch
        tuple file(GRCH37_1), file(GRCH37_2), file(GRCH37_3), file(GRCH37_4), file(GRCH37_5), file(GRCH37_6)
        tuple file(AC_000008_1), file(AC_000008_2), file(AC_000008_3), file(AC_000008_4), file(AC_000008_5), file(AC_000008_6)
        tuple file(AF157706_1), file(AF157706_2), file(AF157706_3), file(AF157706_4), file(AF157706_5), file(AF157706_6)
        tuple file(MF681662_1), file(MF681662_2), file(MF681662_3), file(MF681662_4), file(MF681662_5), file(MF681662_6)
        tuple file(NC_007605_1), file(NC_007605_2), file(NC_007605_3), file(NC_007605_4), file(NC_007605_5), file(NC_007605_6)
        tuple file(MK883610_1), file(MK883610_2), file(MK883610_3), file(MK883610_4), file(MK883610_5), file(MK883610_6)
        tuple file(NC_001348_1), file(NC_001348_2), file(NC_001348_3), file(NC_001348_4), file(NC_001348_5), file(NC_001348_6)
        tuple file(NC_010956_1), file(NC_010956_2), file(NC_010956_3), file(NC_010956_4), file(NC_010956_5), file(NC_010956_6)
        tuple file(NC_001798_1), file(NC_001798_2), file(NC_001798_3), file(NC_001798_4), file(NC_001798_5), file(NC_001798_6)
        tuple file(NC_006273_1), file(NC_006273_2), file(NC_006273_3), file(NC_006273_4), file(NC_006273_5), file(NC_006273_6)
    output:
        tuple val(base), val(reference), file("${base}_human.sam"), file("${base}_${reference}.sam"), file("${base}_summary.csv")

    publishDir "${params.OUTDIR}37bp_human", mode: 'copy', pattern: '*matched_rRNA*.fastq.gz'

    script:
    """
    #!/bin/bash

    echo "Aligning to ${reference}..."

    bowtie2 -p ${task.cpus} -x GRCh37_index -1 ${base}_R1_37bp.fastq.gz -2 ${base}_R2_37bp.fastq.gz --local --no-unal -p 10 -S ${base}_human.sam 
    bowtie2 -p ${task.cpus} -x ${reference}_index -1 ${base}_R1_37bp.fastq.gz -2 ${base}_R2_37bp.fastq.gz --local --no-unal -p 6 -S ${base}_${reference}.sam
    """
}

// Convert sams to bams, process bams, and extract relevant info
process bam_processing { 
    container "quay.io/vpeddu/lava_image:latest"

    // Retry on fail at most three times 
    errorStrategy 'retry'
    maxRetries 3

    input:
        tuple val(base), val(reference), file("${base}_human.sam"), file("${base}_${reference}.sam"), file("${base}_summary.csv")
    output:
        tuple val(base), val(reference), file ("*")

    publishDir "${params.OUTDIR}37bp_bams", mode: 'copy', pattern: "*${reference}*.bam"
    publishDir "${params.OUTDIR}37bp_human", mode: 'copy', pattern: "*human*.bam"
    publishDir "${params.OUTDIR}summary_files", mode: 'copy', pattern: '*final_summary.csv'
    publishDir "${params.OUTDIR}fragment_lengths", mode: 'copy', pattern: '*dedup_inserts.txt'

    script:
    """
    #!/bin/bash

    /usr/local/miniconda/bin/samtools view -@ 6 -Sb ${base}_${reference}.sam | /usr/local/miniconda/bin/samtools sort -o ${base}_${reference}_sorted.bam 
    /usr/local/miniconda/bin/samtools view -@ 6 -Sb ${base}_human.sam | /usr/local/miniconda/bin/samtools sort -o ${base}_human_sorted.bam

    java -jar /usr/bin/picard.jar MarkDuplicates INPUT=${base}_${reference}_sorted.bam OUTPUT=${base}_${reference}_deduplicated.bam VERBOSITY=ERROR REMOVE_DUPLICATES=TRUE VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=TRUE METRICS_FILE=${base}_deduplicated_metrics.txt
    java -jar /usr/bin/picard.jar MarkDuplicates INPUT=${base}_human_sorted.bam OUTPUT=${base}_human_deduplicated.bam VERBOSITY=ERROR REMOVE_DUPLICATES=TRUE VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=TRUE METRICS_FILE=${base}_human_deduplicated_metrics.txt

    human_bam_reads=\$(/usr/local/miniconda/bin/samtools flagstat ${base}_human_sorted.bam | grep "mapped (" | awk '{print \$1}')
    human_dedup_bam_reads=\$(/usr/local/miniconda/bin/samtools flagstat ${base}_human_deduplicated.bam | grep "mapped (" | awk '{print \$1}')
    bam_reads=\$(/usr/local/miniconda/bin/samtools flagstat ${base}_${reference}_sorted.bam | grep "mapped (" | awk '{print \$1}')
    dedup_bam_reads=\$(/usr/local/miniconda/bin/samtools flagstat ${base}_${reference}_deduplicated.bam | grep "mapped (" | awk '{print \$1}')

    /usr/local/miniconda/bin/samtools view -f66 ${base}_${reference}_deduplicated.bam | cut -f9 | awk '{print sqrt(\$0^2)}' > ${base}_${reference}_dedup_inserts.txt

    cp ${base}_summary.csv ${base}_final_summary.csv
    printf "\$human_bam_reads,\$human_dedup_bam_reads,\$bam_reads,\$dedup_bam_reads" >> ${base}_final_summary.csv

    """
} 