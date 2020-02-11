
source('GWAS_funs_libs.R')

# proccessing BRCA GWAS 

mel.gwas <- as_tibble(data.table::fread('../SUMMARY_STATISTICS/melanoma_raw_sumstats.txt.gz', sep = '\t', header=T))
snplist <- vroom::vroom('/project2/xinhe/1kg/bigsnpr/1kg_snp_list.txt.gz', delim = '\n', col_names = F)

mel.gwas$beta <- log(mel.gwas$OR)
cleaned.mel.gwas <- clean_sumstats(sumstats = mel.gwas,
                                    SNPs = snplist, 
                                    cols.to.keep = c('CHR','BP','beta','SE','A1','A2','SNPID_UKB','P'))

cleaned.annot.mel.gwas <- assign.locus.snp(cleaned.sumstats = cleaned.mel.gwas,
                                            ldBed = '../ANNOTATIONS/Euro_LD_Chunks.bed')

vroom::vroom_write(cleaned.annot.mel.gwas[,c('snp','locus','zscore')], 
                   delim = '\t', 
                   col_names = F, 
                   path = '../SUMMARY_STATISTICS/Melanoma_Torus_sumstats.txt.gz')


#annotations <- list.files(path = 'bed_annotations/BRCA_run/', pattern = '*.bed', full.names = T)
annotations <- c('Driver_OCRs.bed','Cancer_Drivers.bed','MEL_OCRs.bed','Conserved_LindbladToh.bed','Immune_OCRs.bed','Immune_h3k27ac.bed')
annotations <- paste0('../ANNOTATIONS/MEL_run/',annotations)

mel.gwas.annots <- annotator(cleaned.annot.mel.gwas, annotations = annotations)

vroom::vroom_write(mel.gwas.annots[,-c(1:6,8,9,10)], 
                   delim = '\t', 
                   col_names = T, 
                   path = '../SUMMARY_STATISTICS/Melanoma_Torus_annotations.txt.gz')

system('/project2/xinhe/software/dap/torus_src/torus -d ../SUMMARY_STATISTICS/Melanoma_Torus_sumstats.txt.gz -annot ../SUMMARY_STATISTICS/Melanoma_Torus_annotations.txt.gz --load_zval -dump_prior mel_torus_prior > melanoma/torus_enrich.txt')
system('/project2/xinhe/software/dap/torus_src/torus -d ../SUMMARY_STATISTICS/Melanoma_Torus_sumstats.txt.gz -annot ../SUMMARY_STATISTICS/Melanoma_Torus_annotations.txt.gz --load_zval -qtl > melanoma/loci_qval.txt')

system('sed -i "s/.bed.1//g" melanoma/torus_run.txt')
system('cat mel_torus_prior/* > melanoma/Melanoma_SNP_PIP.txt')
system('rm -r mel_torus_prior')

# add torus PIP to summary stats
pip <- as_tibble(data.table::fread('melanoma/Melanoma_SNP_PIP.txt', sep=' ', header=F))
colnames(pip) <- c('snp','torus_pip')
cleaned.annot.mel.gwas <- inner_join(cleaned.annot.mel.gwas, pip, by='snp')

# keep loci at 10% FDR
chunk.fdr <- read.delim('melanoma/loci_qval.txt', sep='',header=F, stringsAsFactors = F)
chunks <- chunk.fdr$V2[chunk.fdr$V3 < 0.1]

cleaned.annot.mel.gwas <- cleaned.annot.mel.gwas[cleaned.annot.mel.gwas$locus %in% chunks, ]

# reference panel for computing LD

bigsnp.1kg <- snp_attach(rdsfile = '/project2/xinhe/alan/Cancer/torus/1kg/EUR.1kg.rds')
mel.gwas.1kg.merged <- merge.bigsnp.gwas(cleaned.annot.mel.gwas, bigsnp.1kg)

vroom::vroom_write(mel.gwas.1kg.merged, 
                   delim = '\t', 
                   col_names = T, 
                   path = '../SUMMARY_STATISTICS/Melanoma_Sumstats_SuSiE_Ready.txt.gz')

