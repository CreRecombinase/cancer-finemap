
source('GWAS_funs_libs.R')

# proccessing BRCA GWAS 
## Normalization
brca.gwas <- vroom::vroom('../SUMMARY_STATISTICS/brca_raw_onco2_sumstats.txt.gz', delim = '\t', col_names = T)
bigsnp.1kg <- snp_attach(rdsfile = '~/GERMLINE/1000G/EUR_variable_1kg.rds')
snplist <- vroom::vroom('~/GERMLINE/1000G/1kg_snp_list.txt.gz', delim = '\n', col_names = F)
# re-name columns
# compute z-scores
# filter INFO < 0.9
# keep SNPs with rsIDs in 1000G
cleaned.brca.gwas <- clean_sumstats(sumstats = brca.gwas,
                                    snps.to.keep = snplist,
                                    cols.to.keep = c('chr','position_b37','bcac_onco2_beta','bcac_onco2_se','a0','a1','rsIDs','bcac_onco2_P1df_Wald','bcac_onco2_r2'))

# match to reference panel
cleaned.brca.gwas <- merge.bigsnp.gwas(cleaned.brca.gwas, bigsnp.1kg)


# assign each snp to an LD block
cleaned.brca.gwas <- assign.locus.snp(cleaned.sumstats = cleaned.brca.gwas, ldBed = '../ANNOTATIONS/Euro_LD_Chunks.bed')

#brca.cleaned.pruned <- prune.regions(cleaned.brca.gwas, bigsnp.1kg)

# save for Torus use
vroom::vroom_write(cleaned.brca.gwas[,c('snp','locus','zscore')], 
                   delim = '\t', 
                   col_names = F, 
                   path = '../SUMMARY_STATISTICS/BRCA_Torus_sumstats.txt.gz')

# Add annotations
annotations <- c('Cancer_Drivers.bed','BRCA_OCRs.bed','Conserved_LindbladToh.bed','BRCA_ac.bed','Immune_OCRs.bed','Immune_h3k27ac.bed')
annotations <- paste0('~/GERMLINE/ANNOTATIONS/BRCA_run/',annotations)

cleaned.brca.gwas_annots <- annotator(cleaned.brca.gwas, annotations = annotations)

vroom::vroom_write(cleaned.brca.gwas_annots[,-c(1:6,8:13)], 
                   delim = '\t', 
                   col_names = T, 
                   path = '../SUMMARY_STATISTICS/BRCA_Torus_annotations.txt.gz')

# Torus (run them in parallel)
system('~/torus_src/torus -d ../SUMMARY_STATISTICS/BRCA_Torus_sumstats.txt.gz -annot ../SUMMARY_STATISTICS/BRCA_Torus_annotations.txt.gz --load_zval -dump_prior brca_torus_prior > BRCA/torus_enrich.txt')
system('~/torus_src/torus -d ../SUMMARY_STATISTICS/BRCA_Torus_sumstats.txt.gz -annot ../SUMMARY_STATISTICS/BRCA_Torus_annotations.txt.gz --load_zval -qtl > BRCA/loci_qval.txt')

system('sed -i"" -e "s/.bed.1//g" BRCA/torus_enrich.txt')
system('cat brca_torus_prior/* > BRCA/BRCA_SNP_PIP.txt')
system('rm -r brca_torus_prior')

# add torus PIP to summary stats
pip <- as_tibble(data.table::fread('BRCA/BRCA_SNP_PIP.txt', sep=' ', header=F))
colnames(pip) <- c('snp','torus_pip')
cleaned.brca.gwas <- inner_join(cleaned.brca.gwas, pip, by='snp')

# keep loci at 10% FDR
chunk.fdr <- read.delim('BRCA/loci_qval.txt', sep='',header=F, stringsAsFactors = F)
chunks <- chunk.fdr$V2[chunk.fdr$V3 < 0.1]

# choose chunks that pass FDR
cleaned.brca.gwas <- cleaned.brca.gwas[cleaned.brca.gwas$locus %in% chunks, ]

# ready for fine mapping
vroom::vroom_write(cleaned.brca.gwas, 
                   delim = '\t', 
                   col_names = T, 
                   path = '../SUMMARY_STATISTICS/BRCA_Sumstats_SuSiE_Ready.txt.gz')


# run susie on computing cluster

