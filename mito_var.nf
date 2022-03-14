#!/usr/bin/env nextflow

params.ss = 'samplesheet.csv'
params.ref = "$baseDir/dm6.fa"
params.fai = "$baseDir/dm6.fa.fai"
params.dict = "$baseDir/dm6.dict"
params.outdir = "$baseDir/outdir"

// Check input path parameters to see if they exist
checkPathParamList = [
    params.ss, params.ref, params.fai, params.dict
]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Read in samplesheet and get bams
Channel
    .fromPath(params.ss, checkIfExists: true)
    .splitCsv(header:true)
    .map{row -> file(row.bam)}
    .set{ bam_ch }

//Set up reference channels
Channel
    .fromPath(params.ref)
    .first() // Converts to a value channel to avoid consuming the reference
    .set{ref_ch}

Channel
    .fromPath(params.fai)
    .first() // Converts to a value channel to avoid consuming the reference
    .set{fai_ch}

Channel
    .fromPath(params.dict)
    .first() // Converts to a value channel to avoid consuming the reference
    .set{dict_ch}


process mutect2_mito {
    publishDir params.outdir, mode:'copy'

    input:
    path bam from bam_ch
    path ref from ref_ch
    path fai from fai_ch
    path dict from dict_ch

    output:
    path "${bam}.vcf.gz" into vcf_ch
    path "${bam}.vcf.gz.tbi" into vcf_index_ch
    path "${bam}.vcf.gz.stats" into vcf_stats_ch

    script:
    """
    gatk Mutect2 -R $ref -L chrM --mitochondria-mode -I $bam -O "${bam}.vcf.gz"
    """
}

process mutect2_filter {
    publishDir params.outdir, mode:'copy'

    input:
    path vcf from vcf_ch
    path vcf_index from vcf_index_ch
    path vcf_stats from vcf_stats_ch
    path ref from ref_ch
    path fai from fai_ch
    path dict from dict_ch

    output:
    path "${vcf}.filtered.vcf.gz" into vcf_filt_ch

    script:
    """
    gatk FilterMutectCalls -R $ref -V $vcf --mitochondria-mode -O "${vcf}.filtered.vcf.gz"
    """
}
