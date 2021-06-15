#!/usr/bin/env nextflow
/*
========================================================================================
                                cfDNA analysis
========================================================================================\
 #### Homepage / Documentation
https://github.com/greninger-lab/cfdna
----------------------------------------------------------------------------------------
*/

// Using the Nextflow DSL-2 to account for the logic flow of this workflow
nextflow.preview.dsl=2

def helpMessage() {
    log.info"""
    cfDNA analysis pipeline

    Usage: 
    An example command for running the pipeline is as follows:
    nextflow run greninger-lab/cfdna \\
        --METADATA      Three-column CSV in the format: path to R1, path to R2, reference. [REQUIRED]
        --OUTDIR        Output directory. [REQUIRED]
        
    """.stripIndent()
}

////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
/*                                                    */
/*          SET UP CONFIGURATION VARIABLES            */
/*                                                    */
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

// Show help message
params.help = false
if (params.help){
    helpMessage()
    exit 0
}

METADATA_FILE = file(params.METADATA)
params.OUTDIR= false

// if OUTDIR not set
if (params.OUTDIR == false) {
    println( "Must provide an output directory with --OUTDIR") 
    exit(1)
}
// Make sure OUTDIR ends with trailing slash
if (!params.OUTDIR.endsWith("/")){
   params.OUTDIR = "${params.OUTDIR}/"
}

/////////////////////
// Reference files //
/////////////////////
ADAPTERS = file("s3://clomp-reference-data/tool_specific_data/Tpallidum_WGS/All_adapters.fa")

AC_000008 = file("s3://clomp-reference-data/tool_specific_data/cfdna_index/AC_000008.fasta")
AC_000008_1 = file("s3://clomp-reference-data/tool_specific_data/cfdna_index/AC000008_index.1.bt2")
AC_000008_2 = file("s3://clomp-reference-data/tool_specific_data/cfdna_index/AC000008_index.2.bt2")
AC_000008_3 = file("s3://clomp-reference-data/tool_specific_data/cfdna_index/AC000008_index.3.bt2")
AC_000008_4 = file("s3://clomp-reference-data/tool_specific_data/cfdna_index/AC000008_index.4.bt2")
AC_000008_5 = file("s3://clomp-reference-data/tool_specific_data/cfdna_index/AC000008_index.rev.1.bt2")
AC_000008_6 = file("s3://clomp-reference-data/tool_specific_data/cfdna_index/AC000008_index.rev.2.bt2")
AF157706 = file("s3://clomp-reference-data/tool_specific_data/cfdna_index/AF157706.fasta")
AF157706_1 = file("s3://clomp-reference-data/tool_specific_data/cfdna_index/AF157706_index.1.bt2")
AF157706_2 = file("s3://clomp-reference-data/tool_specific_data/cfdna_index/AF157706_index.2.bt2")
AF157706_3 = file("s3://clomp-reference-data/tool_specific_data/cfdna_index/AF157706_index.3.bt2")
AF157706_4 = file("s3://clomp-reference-data/tool_specific_data/cfdna_index/AF157706_index.4.bt2")
AF157706_5 = file("s3://clomp-reference-data/tool_specific_data/cfdna_index/AF157706_index.rev.1.bt2")
AF157706_6 = file("s3://clomp-reference-data/tool_specific_data/cfdna_index/AF157706_index.rev.2.bt2")
GRCH37 = file("s3://clomp-reference-data/tool_specific_data/cfdna_index/GRCh37_genomic.fna.gz")
GRCH37_1 = file("s3://clomp-reference-data/tool_specific_data/cfdna_index/GRCh37_index.1.bt2")
GRCH37_2 = file("s3://clomp-reference-data/tool_specific_data/cfdna_index/GRCh37_index.2.bt2")
GRCH37_3 = file("s3://clomp-reference-data/tool_specific_data/cfdna_index/GRCh37_index.3.bt2")
GRCH37_4 = file("s3://clomp-reference-data/tool_specific_data/cfdna_index/GRCh37_index.4.bt2")
GRCH37_5 = file("s3://clomp-reference-data/tool_specific_data/cfdna_index/GRCh37_index.rev.1.bt2")
GRCH37_6 = file("s3://clomp-reference-data/tool_specific_data/cfdna_index/GRCh37_index.rev.2.bt2")
MF681662 = file("s3://clomp-reference-data/tool_specific_data/cfdna_index/MF681662.fasta")
MF681662_1 = file("s3://clomp-reference-data/tool_specific_data/cfdna_index/MF681662_index.1.bt2")
MF681662_2 = file("s3://clomp-reference-data/tool_specific_data/cfdna_index/MF681662_index.2.bt2")
MF681662_3 = file("s3://clomp-reference-data/tool_specific_data/cfdna_index/MF681662_index.3.bt2")
MF681662_4 = file("s3://clomp-reference-data/tool_specific_data/cfdna_index/MF681662_index.4.bt2")
MF681662_5 = file("s3://clomp-reference-data/tool_specific_data/cfdna_index/MF681662_index.rev.1.bt2")
MF681662_6 = file("s3://clomp-reference-data/tool_specific_data/cfdna_index/MF681662_index.rev.2.bt2")
MK883610 = file("s3://clomp-reference-data/tool_specific_data/cfdna_index/MK883610.fasta")
MK883610_1 = file("s3://clomp-reference-data/tool_specific_data/cfdna_index/MK883610_index.1.bt2")
MK883610_2 = file("s3://clomp-reference-data/tool_specific_data/cfdna_index/MK883610_index.2.bt2")
MK883610_3 = file("s3://clomp-reference-data/tool_specific_data/cfdna_index/MK883610_index.3.bt2")
MK883610_4 = file("s3://clomp-reference-data/tool_specific_data/cfdna_index/MK883610_index.4.bt2")
MK883610_5 = file("s3://clomp-reference-data/tool_specific_data/cfdna_index/MK883610_index.rev.1.bt2")
MK883610_6 = file("s3://clomp-reference-data/tool_specific_data/cfdna_index/MK883610_index.rev.2.bt2")
NC_001348 = file("s3://clomp-reference-data/tool_specific_data/cfdna_index/NC_001348.fasta")
NC_001348_1 = file("s3://clomp-reference-data/tool_specific_data/cfdna_index/NC001348_index.1.bt2")
NC_001348_2 = file("s3://clomp-reference-data/tool_specific_data/cfdna_index/NC001348_index.2.bt2")
NC_001348_3 = file("s3://clomp-reference-data/tool_specific_data/cfdna_index/NC001348_index.3.bt2")
NC_001348_4 = file("s3://clomp-reference-data/tool_specific_data/cfdna_index/NC001348_index.4.bt2")
NC_001348_5 = file("s3://clomp-reference-data/tool_specific_data/cfdna_index/NC001348_index.rev.1.bt2")
NC_001348_6 = file("s3://clomp-reference-data/tool_specific_data/cfdna_index/NC001348_index.rev.2.bt2")
NC_010956 = file("s3://clomp-reference-data/tool_specific_data/cfdna_index/NC_010956.fasta")
NC_010956_1 = file("s3://clomp-reference-data/tool_specific_data/cfdna_index/NC_010956_index.1.bt2")
NC_010956_2 = file("s3://clomp-reference-data/tool_specific_data/cfdna_index/NC_010956_index.2.bt2")
NC_010956_3 = file("s3://clomp-reference-data/tool_specific_data/cfdna_index/NC_010956_index.3.bt2")
NC_010956_4 = file("s3://clomp-reference-data/tool_specific_data/cfdna_index/NC_010956_index.4.bt2")
NC_010956_5 = file("s3://clomp-reference-data/tool_specific_data/cfdna_index/NC_010956_index.rev.1.bt2")
NC_010956_6 = file("s3://clomp-reference-data/tool_specific_data/cfdna_index/NC_010956_index.rev.2.bt2")
NC_001798 = file("s3://clomp-reference-data/tool_specific_data/cfdna_index/NC_001798.fasta")
NC_001798_1 = file("s3://clomp-reference-data/tool_specific_data/cfdna_index/NC001798_index.1.bt2")
NC_001798_2 = file("s3://clomp-reference-data/tool_specific_data/cfdna_index/NC001798_index.2.bt2")
NC_001798_3 = file("s3://clomp-reference-data/tool_specific_data/cfdna_index/NC001798_index.3.bt2")
NC_001798_4 = file("s3://clomp-reference-data/tool_specific_data/cfdna_index/NC001798_index.4.bt2")
NC_001798_5 = file("s3://clomp-reference-data/tool_specific_data/cfdna_index/NC001798_index.rev.1.bt2")
NC_001798_6 = file("s3://clomp-reference-data/tool_specific_data/cfdna_index/NC001798_index.rev.2.bt2")
NC_006273 = file("s3://clomp-reference-data/tool_specific_data/cfdna_index/NC_006273.fasta")
NC_006273_1 = file("s3://clomp-reference-data/tool_specific_data/cfdna_index/NC006273_index.1.bt2")
NC_006273_2 = file("s3://clomp-reference-data/tool_specific_data/cfdna_index/NC006273_index.2.bt2")
NC_006273_3 = file("s3://clomp-reference-data/tool_specific_data/cfdna_index/NC006273_index.3.bt2")
NC_006273_4 = file("s3://clomp-reference-data/tool_specific_data/cfdna_index/NC006273_index.4.bt2")
NC_006273_5 = file("s3://clomp-reference-data/tool_specific_data/cfdna_index/NC006273_index.rev.1.bt2")
NC_006273_6 = file("s3://clomp-reference-data/tool_specific_data/cfdna_index/NC006273_index.rev.2.bt2")
NC_007605 = file("s3://clomp-reference-data/tool_specific_data/cfdna_index/NC_007605.fasta")
NC_007605_1 = file("s3://clomp-reference-data/tool_specific_data/cfdna_index/NC007605_index.1.bt2")
NC_007605_2 = file("s3://clomp-reference-data/tool_specific_data/cfdna_index/NC007605_index.2.bt2")
NC_007605_3 = file("s3://clomp-reference-data/tool_specific_data/cfdna_index/NC007605_index.3.bt2")
NC_007605_4 = file("s3://clomp-reference-data/tool_specific_data/cfdna_index/NC007605_index.4.bt2")
NC_007605_5 = file("s3://clomp-reference-data/tool_specific_data/cfdna_index/NC007605_index.rev.1.bt2")
NC_007605_6 = file("s3://clomp-reference-data/tool_specific_data/cfdna_index/NC007605_index.rev.2.bt2")


input_read_ch = Channel
    .fromPath(METADATA_FILE)
    .splitCsv(header:true)
    .map{ row-> tuple(file(row.Sample_R1), file(row.Sample_R2), (row.Ref)) }

// Throws exception if paths in METADATA are not valid
Channel
    .fromPath(METADATA_FILE)
    .splitCsv(header:true)
    .map{row-> (file(row.Sample_R1, checkIfExists:true))}.ifEmpty{error "Check metadata file - R1 file does not exist."}
    //.map{row-> (file(row.Sample).isEmpty())}
    //.filter{ it == false}.subscribe{println it}

//
// Import processes
// 

include trim_reads from './modules'
include align_samples from './modules'
 include bam_processing from './modules'


////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
/*                                                    */
/*                 RUN THE WORKFLOW                   */
/*                                                    */
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

workflow {
    trim_reads (
        input_read_ch,
        ADAPTERS
    )
    align_samples (
        trim_reads.out[0],
        tuple(GRCH37_1, GRCH37_2, GRCH37_3, GRCH37_4, GRCH37_5, GRCH37_6),
        tuple(AC_000008_1, AC_000008_2, AC_000008_3, AC_000008_4, AC_000008_5, AC_000008_6),
        tuple(AF157706_1, AF157706_2, AF157706_3, AF157706_4, AF157706_5, AF157706_6),
        tuple(MF681662_1, MF681662_2, MF681662_3, MF681662_4, MF681662_5, MF681662_6),
        tuple(NC_007605_1, NC_007605_2, NC_007605_3, NC_007605_4, NC_007605_5, NC_007605_6),
        tuple(MK883610_1, MK883610_2, MK883610_3, MK883610_4, MK883610_5, MK883610_6),
        tuple(NC_001348_1, NC_001348_2, NC_001348_3, NC_001348_4, NC_001348_5, NC_001348_6),
        tuple(NC_010956_1, NC_010956_2, NC_010956_3, NC_010956_4, NC_010956_5, NC_010956_6),
        tuple(NC_001798_1, NC_001798_2, NC_001798_3, NC_001798_4, NC_001798_5, NC_001798_6),
        tuple(NC_006273_1, NC_006273_2, NC_006273_3, NC_006273_4, NC_006273_5, NC_006273_6)
    )
    bam_processing (
        align_samples.out[0]
    )
}