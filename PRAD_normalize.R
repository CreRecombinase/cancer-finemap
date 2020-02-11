
source('GWAS_funs_libs.R')

# proccessing PRAD GWAS 

prad.gwas <- vroom::vroom('../SUMMARY_STATISTICS/prad_raw_sumstats.txt.gz', delim = '\t', col_names = T)
snplist <- vroom::vroom('/project2/xinhe/1kg/bigsnpr/1kg_snp_list.txt.gz', col_names = F)

prad.gwas$chr <- as.numeric(sapply(strsplit(prad.gwas$varname, split='_'), function(x){x[1]}))
cleaned.prad.gwas <- clean_sumstats(sumstats = prad.gwas,
                                    SNPs = snplist, 
                                    cols.to.keep = c('chr','Position','OR','SE','Effect','Baseline','SNP','Pnorm'))

cleaned.annot.prad.gwas <- assign.locus.snp(cleaned.prad.gwas,
                                            ldBed = '../ANNOTATIONS/Euro_LD_Chunks.bed')

vroom::vroom_write(cleaned.annot.prad.gwas[,c('snp','locus','zscore')], 
                   delim = '\t', 
                   col_names = F, 
                   path = '../SUMMARY_STATISTICS/PRAD_Torus_sumstats.txt.gz')


#annotations <- list.files(path = 'bed_annotations/BRCA_run/', pattern = '*.bed', full.names = T)
annotations <- c('PRAD_h3k27ac.bed','Cancer_Drivers.bed','PRAD_OCRs.bed','Conserved_LindbladToh.bed','Immune_OCRs.bed','Immune_h3k27ac.bed')
annotations <- paste0('../ANNOTATIONS/PRAD_run/',annotations)

prad.gwas.annots <- annotator(cleaned.annot.prad.gwas, annotations = annotations)

vroom::vroom_write(prad.gwas.annots[,-c(1:6,8,9,10)], 
                   delim = '\t', 
                   col_names = T, 
                   path = '../SUMMARY_STATISTICS/PRAD_Torus_annotations.txt.gz')

system('/project2/xinhe/software/dap/torus_src/torus -d ../SUMMARY_STATISTICS/PRAD_Torus_sumstats.txt.gz -annot ../SUMMARY_STATISTICS/PRAD_Torus_annotations.txt.gz --load_zval -dump_prior prad_torus_prior > PRAD/torus_enrich.txt')
system('/project2/xinhe/software/dap/torus_src/torus -d ../SUMMARY_STATISTICS/PRAD_Torus_sumstats.txt.gz -annot ../SUMMARY_STATISTICS/PRAD_Torus_annotations.txt.gz --load_zval -qtl > PRAD/loci_qval.txt')

system('sed -i "s/.bed.1//g" PRAD/torus_run.txt')
system('cat prad_torus_prior/* > PRAD/PRAD_SNP_PIP.txt')
system('rm -r prad_torus_prior')
# prepare for susie

# add torus PIP to summary stats
pip <- as_tibble(data.table::fread('PRAD/PRAD_SNP_PIP.txt', sep=' ', header=F))
colnames(pip) <- c('snp','torus_pip')
cleaned.annot.prad.gwas <- inner_join(cleaned.annot.prad.gwas, pip, by='snp')

# keep loci that pass 10% FDR
chunk.fdr <- read.delim('PRAD/loci_qval.txt',sep='',header=F, stringsAsFactors = F)
chunks <- chunk.fdr$V2[chunk.fdr$V3 < 0.1]

cleaned.annot.prad.gwas <- cleaned.annot.prad.gwas[cleaned.annot.prad.gwas$locus %in% chunks, ]

# reference panel for computing LD

bigsnp.1kg <- snp_attach(rdsfile = '/project2/xinhe/1kg/bigsnpr/EUR_variable_1kg.rds')
prad.gwas.1kg.merged <- merge.bigsnp.gwas(cleaned.annot.prad.gwas, bigsnp.1kg)

vroom::vroom_write(prad.gwas.1kg.merged, 
                   delim = '\t', 
                   col_names = T, 
                   path = '../SUMMARY_STATISTICS/PRAD_Sumstats_SuSiE_Ready.txt.gz')

