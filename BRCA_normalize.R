
source('GWAS_funs_libs.R')

# proccessing BRCA GWAS 
## Normalization
brca.gwas <- vroom::vroom('../SUMMARY_STATISTICS/brca_raw_onco2_sumstats.txt.gz', delim = '\t', col_names = T)
snplist <- vroom::vroom('/project2/xinhe/1kg/bigsnpr/1kg_snp_list.txt.gz', delim = '\n', col_names = F)

cleaned.brca.gwas <- clean_sumstats(sumstats = brca.gwas,
                                    SNPs = snplist, 
                                    cols.to.keep = c('chr','position_b37','bcac_onco2_beta','bcac_onco2_se','a0','a1','rsIDs','bcac_onco2_P1df_Wald'))

cleaned.annot.brca.gwas <- assign.locus.snp(cleaned.sumstats = cleaned.brca.gwas,
                                            ldBed = '../ANNOTATIONS/Euro_LD_Chunks.bed')

vroom::vroom_write(cleaned.annot.brca.gwas[,c('snp','locus','zscore')], 
                   delim = '\t', 
                   col_names = F, 
                   path = '../SUMMARY_STATISTICS/BRCA_Torus_sumstats.txt.gz')

# Add annotations
annotations <- c('Cancer_Drivers.bed','BRCA_OCRs.bed','Conserved_LindbladToh.bed','BRCA_ac.bed','Immune_OCRs.bed','Immune_h3k27ac.bed')
annotations <- paste0('../ANNOTATIONS/BRCA_run/',annotations)

brca.gwas.annots <- annotator(cleaned.annot.brca.gwas, annotations = annotations)

vroom::vroom_write(brca.gwas.annots[,-c(1:6,8,9,10)], 
                   delim = '\t', 
                   col_names = T, 
                   path = '../SUMMARY_STATISTICS/BRCA_Torus_annotations.txt.gz')

# Torus
system('~/torus_src/torus -d ../SUMMARY_STATISTICS/BRCA_Torus_sumstats.txt.gz -annot ../SUMMARY_STATISTICS/BRCA_Torus_annotations.txt.gz --load_zval -dump_prior brca_torus_prior > BRCA/torus_enrich.txt')
system('~/torus_src/torus -d ../SUMMARY_STATISTICS/BRCA_Torus_sumstats.txt.gz -annot ../SUMMARY_STATISTICS/BRCA_Torus_annotations.txt.gz --load_zval -qtl > BRCA/loci_qval.txt')

system('sed -i"" -e "s/.bed.1//g" BRCA/torus_run.txt')
system('cat brca_torus_prior/* > BRCA/BRCA_SNP_PIP.txt')
system('rm -r brca_torus_prior')

# add torus PIP to summary stats
pip <- as_tibble(data.table::fread('BRCA/BRCA_SNP_PIP.txt', sep=' ', header=F))
colnames(pip) <- c('snp','torus_pip')
cleaned.annot.brca.gwas <- inner_join(cleaned.annot.brca.gwas, pip, by='snp')

# keep loci at 10% FDR
chunk.fdr <- read.delim('BRCA/loci_qval.txt', sep='',header=F, stringsAsFactors = F)
chunks <- chunk.fdr$V2[chunk.fdr$V3 < 0.1]

cleaned.annot.brca.gwas <- cleaned.annot.brca.gwas[cleaned.annot.brca.gwas$locus %in% chunks, ]

# add genotype index
bigsnp.1kg <- snp_attach(rdsfile = 'bigsnpr/EUR2.1kg.rds')
brca.gwas.1kg.merged <- merge.bigsnp.gwas(cleaned.annot.brca.gwas, bigsnp.1kg)

vroom::vroom_write(brca.gwas.1kg.merged, 
                   delim = '\t', 
                   col_names = T, 
                   path = 'summary_statistics/BRCA_Sumstats_SuSiE_Ready.txt.gz')



