
suppressMessages(library(tidyverse))
suppressMessages(library(bigsnpr))
suppressMessages(library(susieR))
suppressMessages(library(ggplot2))
suppressMessages(library(rtracklayer))
suppressMessages(library(biomaRt))
suppressMessages(library(GenomicRanges))

# Cleans summary statistics and normalizes them for downstream analysis
clean_sumstats <- function(sumstats, SNPs, cols.to.keep){
  
  stopifnot(!is.null(sumstats))
  stopifnot(!is.null(SNPs))
  stopifnot(length(cols.to.keep) == 8)
  
  chr <- cols.to.keep[1]
  pos <- cols.to.keep[2]
  beta <- cols.to.keep[3]
  se <- cols.to.keep[4]
  a0 <- cols.to.keep[5]
  a1 <- cols.to.keep[6]
  rs <- cols.to.keep[7]
  pval <- cols.to.keep[8]
  colnames(SNPs) <- rs
  
  # keep SNPs in 1kg
  sumstats <- inner_join(sumstats, SNPs, by=rs)
  # Extract relevant columns
  clean.sumstats <- sumstats[ ,c(chr, pos, beta, se, a0, a1, rs, pval)]
  colnames(clean.sumstats) <- c('chr','pos','beta','se','a0','a1','snp', 'pval')
  
  clean.sumstats$beta <- as.numeric(clean.sumstats$beta)
  clean.sumstats$se <- as.numeric(clean.sumstats$se)
  
  # Compute Zscores
  zscore <- clean.sumstats$beta/clean.sumstats$se
  clean.sumstats <- clean.sumstats[!is.na(zscore),]
  zscore <- zscore[!is.na(zscore)]
  clean.sumstats['zscore'] <- zscore
  
  return(clean.sumstats)
}

# assigns ref and alternative allele to dataframe of chromosome/positions and allele1/allele2
assign.ref.alt <- function(df, refGenome, chr, pos, a0, a1){
  
  subdf <- df[,c(chr, pos, a0, a1)]
  colnames(subdf) <- c('CHR','POS','ALLELE1','ALLELE2')
  
  dfRanges <- GRanges(seqnames = subdf$CHR, ranges = IRanges(start = subdf$POS, end = subdf$POS))
  subdf <- subdf %>% 
    group_by(CHR) %>% 
    mutate(REF = unlist(strsplit(as.character(refGenome[[CHR[1]]][POS]), split="")))
  
  a0 <- subdf$ALLELE1
  a1 <- subdf$ALLELE2
  a0.is.alt <- which(a0 != subdf$REF)
  a1.is.alt <- which(a1 != subdf$REF)
  
  subdf[a0.is.alt, 'ALT'] <- a0[a0.is.alt]
  subdf[a1.is.alt, 'ALT'] <- a1[a1.is.alt]
  
  df['ref'] <- subdf['REF']
  df['alt'] <- subdf['ALT']
  
  return(df)
}


# Assigns each SNP to one ld-block
assign.locus.snp <- function(cleaned.sumstats, ldBed){
  
  ld <- vroom::vroom(ldBed, col_names = F)
  ldRanges <- GRanges(seqnames = ld$X1, ranges = IRanges(start = ld$X2, end = ld$X3, names = ld$X4))
  ldRanges <- plyranges::mutate(ldRanges, locus=names(ldRanges))
  
  snpRanges <- GRanges(seqnames = cleaned.sumstats$chr, 
                       ranges   = IRanges(start = cleaned.sumstats$pos, 
                                          end   = cleaned.sumstats$pos,
                                          names = cleaned.sumstats$snp))
  
  snpRanges <- plyranges::mutate(snpRanges, snp=names(snpRanges))
  
  snp.ld.overlap <- plyranges::join_overlap_inner(snpRanges, ldRanges)
  snp.ld.block <- as_tibble(snp.ld.overlap@elementMetadata)
  snp.ld.block <- snp.ld.block[!duplicated(snp.ld.block$snp), ] # some SNPs are in multiple ld-blocks due to edge of ld blocks
  cleaned.annot.sumstats <- inner_join(cleaned.sumstats, snp.ld.block, 'snp')
  
  return(cleaned.annot.sumstats)
}


# Each annotation gets assigned SNPs based on overlap
annotator <- function(gwas, annotations){
  
  snpRanges <- GRanges(seqnames = gwas$chr, 
                       ranges = IRanges(start = gwas$pos, end = gwas$pos))
  snpRanges <- plyranges::mutate(snpRanges, snp=gwas$snp)
  
  for(f in annotations){
    
    name <- paste0(basename(f),'_d')
    curr <- import(f, format='bed')
    subdf <- subsetByOverlaps(snpRanges, curr)
    snpsIn <- unique(subdf$snp)
    
    gwas <- mutate(gwas, !!name := ifelse(snp %in% snpsIn,1,0))
  }
  return(gwas)
}

merge.bigsnp.gwas <- function(gwas, bigSNP){
  
  map <- bigSNP$map
  snp_info <- map[,c('chromosome','physical.pos','allele1','allele2')]
  colnames(snp_info) <- c('chr','pos','a0','a1')
  
  matched.gwas <- as_tibble(bigsnpr::snp_match(gwas, 
                                               snp_info, 
                                               strand_flip = T)) %>% 
    dplyr::rename(og_index = `_NUM_ID_.ss`) %>% 
    dplyr::rename(bigSNP_index = `_NUM_ID_`) %>%
    mutate(zscore = beta/se)
  
  return(matched.gwas)
  
}

############## SUSIE 
run.susie <- function(sumstats, ref_panel, ldchunk, L, prior){
  
  sub.sumstats <- sumstats[sumstats$locus == ldchunk, ]
  if(nrow(sub.sumstats) > 1){
    X <- ref_panel$genotypes[ ,sub.sumstats$bigSNP_index]
    X <- scale(X, center = T, scale = T)
    zhat <- sub.sumstats$zscore
    R <- cov2cor((crossprod(X) + tcrossprod(zhat))/nrow(X))
    if(prior){
      res <- susie_rss(z = zhat, prior_weights = sub.sumstats$torus_pip, R = R, L = L, check_z = FALSE)
    }
    else{
      res <- susie_rss(z = zhat, R = R, L = L)
    }
    return(res)
  }
}

# merges susie results with original summary statistics data frame
# ASSUMES L = 1! ONLY ONE CREDIBLE SET PER LOCUS! 
merge_susie_sumstats <- function(susie_results, sumstats){
  
  sumstats$susie_pip <- 0
  sumstats$CS <- 0
  loci <- names(susie_results)
  
  for(l in loci){
    n.snps <- length(susie_results[[l]]$pip)
    sumstats[sumstats$locus == as.numeric(l), "susie_pip"] <- susie_results[[l]]$pip
    
    snps.in.cs <- rep(0, n.snps)
    if(!is.null(susie_results[[l]]$sets$cs)){
      snps.in.cs[unlist(susie_results[[l]]$sets$cs$L1)] <- 1
    }
    sumstats[sumstats$locus == as.numeric(l), "CS"] <- snps.in.cs
  }
  return(sumstats)
}


add.gtex.annotation <- function(sumstats, gtex){
  
  stopifnot('chrom' %in% colnames(gtex))
  stopifnot('position' %in% colnames(gtex))
  stopifnot('Name' %in% colnames(gtex))
  stopifnot('Tissue' %in% colnames(gtex))
  
  # load GTEx v7 significant snp-gene pairs
  gtex$id <- paste0(gtex$chrom, '_', gtex$position)
  gtex <- gtex[, c('id', 'Name','Tissue')]
  sumstats$id <- paste0(sumstats$chr, '_', sumstats$pos)
  
  # first left join, to get annotations of fine-mapped snps
  annot.sumstats <- left_join(x = sumstats, y = gtex, by='id')
  # collapse annotation into one
  annots <- annot.sumstats %>% group_by(id) %>% summarise(eGenes = paste(unique(Name), collapse = '\n'), Tissues=paste(unique(Tissue), collapse = '\n'))
  # join again to get original df with new annotation
  annot.sumstats <- inner_join(sumstats, annots, by = 'id') %>% dplyr::select(-id)
  
  return(annot.sumstats)
  
}


add.nearby.gene <- function(sumstats, gene.ranges){
  
  snp.ranges <- GRanges(seqnames = sumstats$chr, 
                               ranges = IRanges(start = sumstats$pos, 
                                                end   = sumstats$pos))
  snp.ranges <- plyranges::mutate(snp.ranges, snp = sumstats$snp)
  
  nearby.genes <- plyranges::join_overlap_inner(gene.ranges, snp.ranges, maxgap=10000) %>%
    as_tibble() %>%
    group_by(snp) %>%
    summarise(nearby_genes = paste(name, collapse='\n'))
  
  annot.sumstats <- left_join(sumstats, nearby.genes, by='snp')
  return(annot.sumstats)
}

add.consequence <- function(sumstats){
  
  stopifnot('snp' %in% colnames(sumstats))
  
  snp.list <- sumstats$snp
  snp.mart <- useMart(biomart = "ENSEMBL_MART_SNP",
                      dataset = "hsapiens_snp")
  
  result <- getBM(attributes = c("refsnp_id","consequence_type_tv"),
                  filters    = "snp_filter", 
                  values     = snp.list, 
                  mart       = snp.mart) %>% as_tibble %>% group_by(refsnp_id) %>% summarise(consequence = paste0(consequence_type_tv, collapse = '\n'))
  colnames(result) <- c('snp', 'consequence')
  annot.sumstats <- left_join(sumstats, result, by='snp')
  return(annot.sumstats)
}



