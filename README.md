# mito_var

## Map reads with nf-core/sarek
`slurm_sub.py -j sarek_nvs024_ak nextflow run nf-core/sarek -r 2.7.1 --input NVS024_AK_samplesheet.tsv -c /mnt/home3/ajr236/gurdon.config --genome dm6 --save_bam_mapped --outdir sarek_nvs024_ak`

## Run mito_var pipeline (runs mutect2 in mitochondrial mode )
slurm_sub.py -j nvs024_ak_mitovar nextflow run mito_var.nf --ss NVS024_AK_bams_samplesheet.csv --ref /mnt/beegfs/home1/reid/ajr236/projects/hansong/ref/dm6.fa --fai /mnt/beegfs/home1/reid/ajr236/projects/hansong/ref/dm6.fa.fai --dict /mnt/
beegfs/home1/reid/ajr236/projects/hansong/ref/dm6.dict --outdir mito_var_outdir  -with-singularity /mnt/home3/nextflow/nfcore-sarek-2.7.1.img -c /mnt/home3/nextflow/gurdon.config
