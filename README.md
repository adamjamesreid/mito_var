# mito_var

## Map reads with nf-core/sarek
`slurm_sub.py -j sarek_nvs024_ak nextflow run nf-core/sarek -r 2.7.1 --input NVS024_AK_samplesheet.tsv -c /mnt/home3/ajr236/gurdon.config --genome dm6 --save_bam_mapped --outdir sarek_nvs024_ak`

## Run mito_var pipeline (runs mutect2 in mitochondrial mode )
`slurm_sub.py -j nvs024_ak_mitovar nextflow run mito_var.nf --ss NVS024_AK_bams_samplesheet.csv --ref /mnt/beegfs/home1/reid/ajr236/projects/hansong/ref/dm6.fa --fai /mnt/beegfs/home1/reid/ajr236/projects/hansong/ref/dm6.fa.fai --dict /mnt/
beegfs/home1/reid/ajr236/projects/hansong/ref/dm6.dict --outdir mito_var_outdir  -with-singularity /mnt/home3/nextflow/nfcore-sarek-2.7.1.img -c /mnt/home3/nextflow/gurdon.config`

## Further filter variants based on location in sequence (not in the control region), PASS, and allele frequency (<0.9), then count them
`python count_variants.py NVS024_AK_vcfs_samplesheet.csv > counts.txt`

## Counts were then collected across sequencing experiments (by hand into one file) and combined with mitochondrial depth calculations (calculated as below and combined into one file by hand) to provide a data file for analysis in R
### Calculate mito coverage
`for x in *.md.bam; do mosdepth -c chrM ${x}_chrM $x; done`

`python ~/code_development/mito_var/parse_mosdepth.py *summary.txt  > chrM_depth.txt`

`combine_tables.py <all_counts> <all coverage>`

### R analysis
GBP0004_analysis.R
