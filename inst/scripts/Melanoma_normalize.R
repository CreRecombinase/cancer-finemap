
source("GWAS_funs_libs.R")

# proccessing MELANOMA GWAS
sumstats <- "~/GERMLINE/SUMMARY_STATISTICS/melanoma_raw_sumstats.txt.gz"
snplist <- "~/GERMLINE/PLINK/1kg_snp_list.txt.gz"
plink <- "~/GERMLINE/PLINK/EUR_variable_1kg.rds"

mel.gwas <- as_tibble(data.table::fread(sumstats, sep = "\t", header = T))
snplist <- vroom::vroom(snplist, delim = "\n", col_names = F)

mel.gwas$beta <- log(mel.gwas$OR)
cleaned.mel.gwas <- clean_sumstats(
  sumstats = mel.gwas,
  SNPs = snplist,
  cols.to.keep = c("CHR", "BP", "beta", "SE", "A1", "A2", "SNPID_UKB", "P")
)

cleaned.annot.mel.gwas <- assign.locus.snp(
  cleaned.sumstats = cleaned.mel.gwas,
  ldBed = "~/GERMLINE/ANNOTATIONS/Euro_LD_Chunks.bed"
)

vroom::vroom_write(cleaned.annot.mel.gwas[, c("snp", "locus", "zscore")],
  delim = "\t",
  col_names = F,
  path = "~/GERMLINE/SUMMARY_STATISTICS/ldsc_ready/Melanoma_Torus_sumstats.txt.gz"
)


annotations <- c("Driver_OCRs.bed", "Cancer_Drivers.bed", "MEL_OCRs.bed", "Conserved_LindbladToh.bed", "Immune_h3k27ac.bed")
annotations <- paste0("~/GERMLINE/ANNOTATIONS/MEL_run/", annotations)

mel.gwas.annots <- annotator(cleaned.annot.mel.gwas, annotations = annotations)

vroom::vroom_write(mel.gwas.annots[, -c(1:6, 8, 9, 10)],
  delim = "\t",
  col_names = T,
  path = "~/GERMLINE/SUMMARY_STATISTICS/torus_ready/Melanoma_Torus_annotations.txt.gz"
)

system("~/torus_src/torus -d ../SUMMARY_STATISTICS/torus_ready/Melanoma_Torus_sumstats.txt.gz -annot ../SUMMARY_STATISTICS/torus_ready/Melanoma_Torus_annotations.txt.gz --load_zval -dump_prior mel_torus_prior > melanoma/torus_enrich.txt")
system("~/torus_src/torus -d ../SUMMARY_STATISTICS/torus_ready/Melanoma_Torus_sumstats.txt.gz -annot ../SUMMARY_STATISTICS/torus_ready/Melanoma_Torus_annotations.txt.gz --load_zval -qtl > melanoma/loci_qval.txt")

system('sed -i"" -e "s/.bed.1//g" melanoma/torus_enrich.txt')
system("cat mel_torus_prior/* > melanoma/Melanoma_SNP_PIP.txt")
system("rm -r mel_torus_prior")

# add torus PIP to summary stats
pip <- as_tibble(data.table::fread("melanoma/Melanoma_SNP_PIP.txt", sep = " ", header = F))
colnames(pip) <- c("snp", "torus_pip")
cleaned.annot.mel.gwas <- inner_join(cleaned.annot.mel.gwas, pip, by = "snp")

# keep loci at 10% FDR
chunk.fdr <- read.delim("melanoma/loci_qval.txt", sep = "", header = F, stringsAsFactors = F)
chunks <- chunk.fdr$V2[chunk.fdr$V3 < 0.1]

cleaned.annot.mel.gwas <- cleaned.annot.mel.gwas[cleaned.annot.mel.gwas$locus %in% chunks, ]

# reference panel for computing LD

bigsnp.1kg <- snp_attach(rdsfile = plink)
mel.gwas.1kg.merged <- merge.bigsnp.gwas(cleaned.annot.mel.gwas, bigsnp.1kg)

vroom::vroom_write(mel.gwas.1kg.merged,
  delim = "\t",
  col_names = T,
  path = "../SUMMARY_STATISTICS/susie_ready/Melanoma_Sumstats_SuSiE_Ready.txt.gz"
)

# prune regions

mel.gwas.1kg.merged.pruned <- prune.regions(sumstats = mel.gwas.1kg.merged, ref_panel = bigsnp.1kg)

vroom::vroom_write(mel.gwas.1kg.merged.pruned,
  delim = "\t",
  col_names = T,
  path = "../SUMMARY_STATISTICS/susie_ready/Melanoma_Sumstats_SuSiE_Ready_Pruned.txt.gz"
)

# run susie
susie_res <- list()
for (i in 1:length(chunks)) {
  print(paste0(i, " out of ", length(chunks)))
  susie_res[[as.character(chunks[i])]] <- run.susie(mel.gwas.1kg.merged.pruned, bigsnp.1kg, chunks[i], L = 1, prior = T)
}

save(susie_res, file = "melanoma/melanoma_susie_L1.Robj")

susie_res_unif <- list()
for (i in 1:length(chunks)) {
  print(paste0(i, " out of ", length(chunks)))
  susie_res_unif[[as.character(chunks[i])]] <- run.susie(mel.gwas.1kg.merged.pruned, bigsnp.1kg, chunks[i], L = 1, prior = F)
}

save(susie_res_unif, file = "melanoma/melanoma_susie_L1_UNIFORM.Robj")
