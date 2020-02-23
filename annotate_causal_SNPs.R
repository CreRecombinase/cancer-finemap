
source('GWAS_funs_libs.R')

### Annotate SNPs with PIP > 0.1
### This script is gross, but it works. It basically just calls all the functions I made in GWAS_funs_libs.R

# load SuSiE results and normalized summary statistics
load('BRCA/BRCA_susie_L1.Robj')
sumstats <- vroom::vroom('~/GERMLINE/SUMMARY_STATISTICS/susie/BRCA_Sumstats_SuSiE_Ready.txt.gz', delim = '\t', col_names = T)
merged.sumstats <- merge_susie_sumstats(susie_res, sumstats)

# Which annotation did each SNP come from?
top.snps <- merged.sumstats[merged.sumstats$susie_pip > 0.1 & merged.sumstats$CS == 1, ]
annot.df <- suppressMessages(vroom::vroom('~/GERMLINE/SUMMARY_STATISTICS/torus/BRCA_Torus_annotations.txt.gz', col_names = T))
annot.df <- annot.df[annot.df$snp %in% top.snps$snp, ]
colnames(annot.df) <- gsub(colnames(annot.df), pattern = '.bed_d', replacement = '') # clean up annotation names
annot.df$merged <- apply(annot.df, MARGIN = 1, FUN = function(x){return(paste0(names(x)[which(x==1)], collapse = ';'))})
top.snps <- inner_join(top.snps, annot.df[ ,c("snp","merged")], by = "snp")

# Which gene is each SNP inside of?
genes <- suppressMessages(vroom::vroom('~/GERMLINE/refGenome/GENES_COORDs_b37.txt', col_names = T))
top.snps <- add.nearby.gene(top.snps, genes, dist=0, gene_type='inside_gene')

# Which cancer driver is near each SNP?
drivers <- readLines('~/GERMLINE/ANNOTATIONS/cancer_driver_genes/AllDrivers_FDR_0.10.txt')
top.snps <- add.nearby.gene(top.snps, genes[genes$gene %in% drivers, ], dist = 500000, gene_type='nearby_driver')

# Which immune gene is near each SNP?
immune.genes <- suppressMessages(vroom::vroom('~/GERMLINE/refGenome/immune_genes_cords.txt', col_names = T))
top.snps <- add.nearby.gene(top.snps, immune.genes, dist = 500000, gene_type='nearby_immune')

# Which promoter is each SNP in?
tss.df <- as_tibble(data.frame(chr=genes$chr, start=genes$start-2000, end=genes$start+1000, gene=genes$gene))
top.snps <- add.nearby.gene(top.snps, tss.df, dist = 0, gene_type='in_promoter')

# Which SNPs are eQTLs in breast and whole blood?
gtex <- suppressMessages(vroom::vroom('~/GERMLINE/GTEx_Analysis_v7_eQTL/Breast_Mammary_var_genes.txt.gz', col_names = T, delim = '\t'))
top.snps <- add.gtex.annotation(top.snps, gtex, tissue_type = 'breast_gtex_eGenes')
gtex <- suppressMessages(vroom::vroom('~/GERMLINE/GTEx_Analysis_v7_eQTL/Whole_Blood_var_genes.txt.gz', col_names = T, delim = '\t'))
top.snps <- add.gtex.annotation(top.snps, gtex, tissue_type = 'whole_blood_gtex_eGenes')


# What trait is found for each SNP in GWAS Catalog?
gwas <- suppressMessages(vroom::vroom('~/GERMLINE/GWAS_catalog/gwas_catalog_v1.0-associations_e98_r2020-02-08.tsv',col_names = T))
top.snps <- add.gwas_catalog(top.snps, gwas)

# What is the LocusZoom for each SNP?
top.snps <- add.locuszoom.link(top.snps)

# What region is each SNP in?
top.snps <- add.region(top.snps)
top.snps$chrPos <- paste0('chr',top.snps$chr,':',top.snps$pos)

# select relevant columns
top.snps.kable <- top.snps[ , c('region','snp','chrPos','inside_gene','in_promoter','nearby_driver','nearby_immune','breast_gtex_eGenes','whole_blood_gtex_eGenes','merged','gwas_catalog','zscore','susie_pip','locuszoom')]
colnames(top.snps.kable) <- c('Region','SNP','Position','Gene Inside','Promoter Inside','Nearby Drivers (500kb)','Nearby Immune Gene (500kb)','Breast eGenes','Whole-Blood eGenes','Annotation','GWAS Catalog','ZScore','PIP','locuszoom')

# do some formatting..
top.snps.kable <- top.snps.kable %>% replace(. == 'NA' | . == "", NA)
top.snps.kable <- top.snps.kable[order(top.snps.kable$Region), ]
top.snps.kable$PIP <- round(top.snps.kable$PIP, 2)
top.snps.kable$ZScore <- round(top.snps.kable$ZScore, 2)

# write results
vroom::vroom_write(top.snps.kable, path = '~/GERMLINE/cancer-finemap/BRCA/BRCA_finemapped_L1_pip0.1_annotated.txt', delim = '\t', col_names = T)


