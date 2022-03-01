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

// Read in samplesheet and get bams (two channels required)
Channel
    .fromPath(params.ss, checkIfExists: true)
    .splitCsv(header:true)
    .map{row -> file(row.bam)}
    .set{ bam_ch1 }

Channel
    .fromPath(params.ss, checkIfExists: true)
    .splitCsv(header:true)
    .map{row -> file(row.bam)}
    .set{ bam_ch2 }

// Set up channel of bam indexes for pileup
Channel
    .fromPath(params.ss, checkIfExists: true)
    .splitCsv(header:true)
    .map{row -> file(row.bam + '.bai', checkIfExists: true)}
    .set{ bai_ch }

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

// Run mpileup to call all possible variant sites which GATK can then be forced to consider
// This allows checking of evidence for every site and allows sites which were otherwise missed through
process mpileup {
    publishDir params.outdir, mode:'copy'

    input:
    path bam from bam_ch1
    path bai from bai_ch
    path ref from ref_ch
    path fai from fai_ch

    output:
    path "${bam}.mpileup.vcf.gz" into mpileup_vcf_ch
    path "${bam}.mpileup.vcf.gz.tbi" into mpileup_vcf_tbi_ch

    script:
    """
    bcftools mpileup -r chrM -O z -f $ref -o "${bam}.mpileup.vcf.gz" $bam
    tabix -p vcf "${bam}.mpileup.vcf.gz"
    """

}

process mutect2_mito {
    publishDir params.outdir, mode:'copy'

    input:
    path mpileup_vcf from mpileup_vcf_ch
    path mpileup_vcf_tbi from mpileup_vcf_tbi_ch
    path bam from bam_ch2
    path ref from ref_ch
    path fai from fai_ch
    path dict from dict_ch

    output:
    path "${bam}.mutect2.vcf.gz" into vcf_ch
    path "${bam}.mutect2.vcf.gz.tbi" into vcf_index_ch
    path "${bam}.mutect2.vcf.gz.stats" into vcf_stats_ch

    script:
    """
    #gatk Mutect2 -R hg38.fa -L chrM -I MAid261_AK136.md.bam -O MAid261_AK136.md.bam.normalmode.vcf --alleles MAid261_AK136.md.bam.bcftools.vcf.gz
    #gatk Mutect2 -R $ref -L chrM --mitochondria-mode -I $bam -O "${bam}.vcf.gz"
    # Add --mitochondria-mode back in ?????
    gatk Mutect2 -R $ref -L chrM -I $bam -O "${bam}.mutect2.vcf.gz" --alleles $mpileup_vcf
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
