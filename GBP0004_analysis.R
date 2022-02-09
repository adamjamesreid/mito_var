# Analyse GBP0004 results
setwd('/mnt/home1/reid/ajr236/projects/hansong/GBP0004/analysis')

x<-read.table("counts_coverage_combined.tsv", header=TRUE)

# Set up ordering for genotypes
x$genotype <- factor(x$genotype , levels=c("WT", "RECdelMTS", "RECKO"))
x$age <- factor(x$age , levels=c("Young", "Old"))

# chrM coverage histogram
hist(x$chrM_depth, breaks=20)

# chrM by coverage boxplot - NVS036 has much higher coverage, 
#n.b. it also has weird GC content, lower whole genome coverage and lots of duplicates
boxplot(chrM_depth~batch, data=x)

# coverage is higher in young samples (could downsample everything)
boxplot(chrM_depth~age, data=x)
median(subset(x, x$age=="Young")$chrM_depth)
median(subset(x, x$age=="Old")$chrM_depth)
mean(subset(x, x$age=="Young")$chrM_depth)
mean(subset(x, x$age=="Old")$chrM_depth)

# WT looks a bit higher in coverage
boxplot(chrM_depth~genotype, data=x)
boxplot(chrM_depth~genotype+age, data=x)

# Total variants (log) by age and genotype
boxplot(log(variants)~genotype, data=x)
boxplot(log(variants)~age, data=x)
boxplot(log(variants)~age+genotype, data=x)

# Variants types (log) by age and genotype
boxplot(log(snps)~genotype, data=x)
boxplot(log(indels)~genotype, data=x)

# Variant details for NVS024
boxplot(variants~genotype+age, data=subset(x, x$batch=="NVS024"), cex.axis=0.7)
boxplot(snps~genotype+age, data=subset(x, x$batch=="NVS024"), cex.axis=0.7)
boxplot(indels~genotype+age, data=subset(x, x$batch=="NVS024"), cex.axis=0.7)

# Variant details for NVS036
boxplot(variants~genotype+age, data=subset(x, x$batch=="NVS036"), cex.axis=0.7)
boxplot(snps~genotype+age, data=subset(x, x$batch=="NVS036"), cex.axis=0.7)
boxplot(indels~genotype+age, data=subset(x, x$batch=="NVS036"), cex.axis=0.7)

# Variant details for NVS040
boxplot(variants~genotype+age, data=subset(x, x$batch=="NVS040"), cex.axis=0.7)
boxplot(snps~genotype+age, data=subset(x, x$batch=="NVS040"), cex.axis=0.7)
boxplot(indels~genotype+age, data=subset(x, x$batch=="NVS040"), cex.axis=0.7)

# Variants vs read depth (coverage > 2500 looks to result in larger numbers of SNPs)
plot(x$chrM_depth, x$variants)
plot(x$chrM_depth, x$variants, xlim=c(0,2500))
