
## suppressMessages(library(tidyverse))
## suppressMessages(library(bigsnpr))
## suppressMessages(library(susieR))
## suppressMessages(library(ggplot2))
## suppressMessages(library(GenomicRanges))

##' @title
##'
##'
##' @param seqname
##' @param start
##' @param end
##' @return
##' @export
make_ranges <- function(seqname, start, end) {
  return(GenomicRanges::GRanges(seqnames = seqname, ranges = IRanges::IRanges(start = start, end = end)))
}

##' @title
##'
##'
##' @param locus
##' @param ld_chunk_file
##' @return
##' @export
locus.lookup <- function(locus, ld_chunk_file = "~/GERMLINE/ANNOTATIONS/Euro_LD_Chunks.bed") {
  euro_ld <- read.delim(ld_chunk_file, sep = "\t", stringsAsFactors = F, header = F)
  currLocus <- euro_ld[euro_ld$V4 == locus, ]
  chrom <- currLocus$V1
  start <- currLocus$V2
  end <- currLocus$V3
  return(list(locus = locus, chrom = chrom, start = start, end = end))
}

# Cleans summary statistics and normalizes them for downstream analysis
##' @title
##'
##'
##' @param sumstats
##' @param cols.to.keep
##' @return
##' @export
clean_sumstats <- function(sumstats, cols.to.keep) {
  stopifnot(!is.null(sumstats))
  stopifnot(length(cols.to.keep) == 9)

  chr <- cols.to.keep[1]
  pos <- cols.to.keep[2]
  beta <- cols.to.keep[3]
  se <- cols.to.keep[4]
  a0 <- cols.to.keep[5]
  a1 <- cols.to.keep[6]
  rs <- cols.to.keep[7]
  pval <- cols.to.keep[8]
  info <- cols.to.keep[9]

  # keep SNPs in 1kg
  # sumstats <- dplyr::inner_join(sumstats, snps.to.keep, by=rs)
  # Extract relevant columns
  clean.sumstats <- sumstats[, c(chr, pos, beta, se, a0, a1, rs, pval, info)]
  colnames(clean.sumstats) <- c("chr", "pos", "beta", "se", "a0", "a1", "snp", "pval", "info")

  # filter imputation quality
  clean.sumstats <- clean.sumstats[!is.na(clean.sumstats$info), ]
  clean.sumstats <- clean.sumstats[clean.sumstats$info > 0.9, ]

  # keep SNPs with rs ids
  clean.sumstats <- clean.sumstats[startsWith(clean.sumstats$snp, "rs"), ]

  clean.sumstats$beta <- as.numeric(clean.sumstats$beta)
  clean.sumstats$se <- as.numeric(clean.sumstats$se)

  # Compute Zscores
  zscore <- clean.sumstats$beta / clean.sumstats$se
  clean.sumstats <- clean.sumstats[!is.na(zscore), ]
  zscore <- zscore[!is.na(zscore)]
  clean.sumstats["zscore"] <- zscore

  return(clean.sumstats)
}

# assigns ref and alternative allele to dataframe of chromosome/positions and allele1/allele2
##' @title
##'
##'
##' @param df
##' @param refGenome
##' @param chr
##' @param pos
##' @param a0
##' @param a1
##' @return
##' @export
assign.ref.alt <- function(df, refGenome, chr, pos, a0, a1) {
  subdf <- df[, c(chr, pos, a0, a1)]
  colnames(subdf) <- c("CHR", "POS", "ALLELE1", "ALLELE2")

  dfRanges <- GenomicRanges::GRanges(seqnames = subdf$CHR, ranges = IRanges::IRanges(start = subdf$POS, end = subdf$POS))
  subdf <- subdf %>%
    dplyr::group_by(CHR) %>%
    dplyr::mutate(REF = unlist(strsplit(as.character(refGenome[[CHR[1]]][POS]), split = "")))

  a0 <- subdf$ALLELE1
  a1 <- subdf$ALLELE2
  a0.is.alt <- which(a0 != subdf$REF)
  a1.is.alt <- which(a1 != subdf$REF)

  subdf[a0.is.alt, "ALT"] <- a0[a0.is.alt]
  subdf[a1.is.alt, "ALT"] <- a1[a1.is.alt]

  df["ref"] <- subdf["REF"]
  df["alt"] <- subdf["ALT"]

  return(df)
}


# Assigns each SNP to one ld-block
##' @title
##'
##'
##' @param cleaned.sumstats
##' @param ldBed
##' @return
##' @export
assign.locus.snp <- function(cleaned.sumstats, ldBed) {
  ld <- vroom::vroom(ldBed, col_names = F)
  ldRanges <- plyranges::make_ranges(ld$X1, ld$X2, ld$X3)
  ldRanges <- plyranges::mutate(ldRanges, locus = ld$X4)

  snpRanges <- GenomicRanges::GRanges(
    seqnames = cleaned.sumstats$chr,
    ranges = IRanges::IRanges(
      start = cleaned.sumstats$pos,
      end = cleaned.sumstats$pos,
      names = cleaned.sumstats$snp
    )
  )

  snpRanges <- plyranges::mutate(snpRanges, snp = names(snpRanges))

  snp.ld.overlap <- plyranges::join_overlap_inner(snpRanges, ldRanges)
  snp.ld.block <- tibble::as_tibble(snp.ld.overlap@elementMetadata)
  snp.ld.block <- snp.ld.block[!duplicated(snp.ld.block$snp), ] # some SNPs are in multiple ld-blocks due to edge of ld blocks
  cleaned.annot.sumstats <- dplyr::inner_join(cleaned.sumstats, snp.ld.block, "snp")

  return(cleaned.annot.sumstats)
}


# Each annotation gets assigned SNPs based on overlap
##' @title
##'
##'
##' @param gwas
##' @param annotations
##' @return
##' @export
annotator <- function(gwas, annotations) {
  snpRanges <- plyranges::make_ranges(gwas$chr, gwas$pos, gwas$pos)
  snpRanges <- plyranges::mutate(snpRanges, snp = gwas$snp)

  for (f in annotations) {
    name <- paste0(basename(f), "_d")
    curr <- rtracklayer::import(f, format = "bed")
    subdf <- subsetByOverlaps(snpRanges, curr)
    snpsIn <- unique(subdf$snp)

    gwas <- mutate(gwas, !!name := ifelse(snp %in% snpsIn, 1, 0))
  }
  return(gwas)
}
##' @title
##'
##'
##' @param gwas
##' @param bigSNP
##' @return
##' @export
merge.bigsnp.gwas <- function(gwas, bigSNP) {
  map <- bigSNP$map
  snp_info <- map[, c("chromosome", "physical.pos", "allele1", "allele2")]
  colnames(snp_info) <- c("chr", "pos", "a0", "a1")

  matched.gwas <- tibble::as_tibble(bigsnpr::snp_match(gwas,
    snp_info,
    strand_flip = T,
    match.min.prop = 1
  )) %>%
    rename(og_index = `_NUM_ID_.ss`) %>%
    rename(bigSNP_index = `_NUM_ID_`) %>%
    mutate(zscore = beta / se)

  return(matched.gwas)
}
##' @title
##'
##'
##' @param sumstats
##' @param ref_panel
##' @return
##' @export
prune.regions <- function(sumstats, ref_panel) {
  df_list <- list()
  r2_thresh <- 0.1
  for (l in unique(sumstats$locus)) {
    sub.sumstats <- sumstats[sumstats$locus == l, ]

    topSnpPos <- which.max(sub.sumstats$susie_pip)
    topSnp <- sub.sumstats$bigSNP_index[topSnpPos]
    topSnpG <- ref_panel$genotypes[, topSnp]
    G <- ref_panel$genotypes[, sub.sumstats$bigSNP_index]

    r2 <- as.vector(cor(topSnpG, G, method = "pearson"))
    subRegions <- sub.sumstats[r2 > r2_thresh, ]
    subRegions$r2 <- r2[r2 > r2_thresh]
    df_list[[l]] <- subRegions
  }
  return(Reduce(rbind, df_list))
}

# SUSIE related functions
##' @title
##'
##'
##' @param sumstats
##' @param ref_panel
##' @param ldchunk
##' @param L
##' @param prior
##' @return
##' @export
run.susie <- function(sumstats, ref_panel, ldchunk, L, prior) {
  sub.sumstats <- sumstats[sumstats$locus == ldchunk, ]
  if (nrow(sub.sumstats) > 1) {
    X <- ref_panel$genotypes[, sub.sumstats$bigSNP_index]
    X <- scale(X, center = T, scale = T)
    zhat <- sub.sumstats$zscore
    R <- cov2cor((crossprod(X) + tcrossprod(zhat)) / nrow(X))
    if (prior) {
      res <- susieR::susie_rss(z = zhat, prior_weights = sub.sumstats$torus_pip, R = R, L = L, check_z = FALSE)
    }
    else {
      res <- susieR::susie_rss(z = zhat, R = R, L = L)
    }
    return(res)
  }
}

# merges susie results with original summary statistics data frame
# ASSUMES L = 1! ONLY ONE CREDIBLE SET PER LOCUS!
##' @title
##'
##'
##' @param susie_results
##' @param sumstats
##' @return
##' @export
merge_susie_sumstats <- function(susie_results, sumstats) {
  sumstats$susie_pip <- 0
  sumstats$CS <- 0
  loci <- names(susie_results)

  for (l in loci) {
    n.snps <- length(susie_results[[l]]$pip)
    sumstats[sumstats$locus == as.numeric(l), "susie_pip"] <- susie_results[[l]]$pip

    snps.in.cs <- rep(0, n.snps)
    if (!is.null(susie_results[[l]]$sets$cs)) {
      snps.in.cs[unlist(susie_results[[l]]$sets$cs$L1)] <- 1
    }
    sumstats[sumstats$locus == as.numeric(l), "CS"] <- snps.in.cs
  }
  return(sumstats)
}

# Annotations for causal SNPs (apply these after fine-mapping!)
##' @title
##'
##'
##' @param sumstats
##' @param gtex
##' @param tissue_type
##' @return
##' @export
add.gtex.annotation <- function(sumstats, gtex, tissue_type = "") {
  stopifnot("chr" %in% colnames(gtex))
  stopifnot("pos" %in% colnames(gtex))
  stopifnot("name" %in% colnames(gtex))

  # load GTEx v7 significant snp-gene pairs
  gtex$var_id <- paste0(gtex$chr, "_", gtex$pos)
  gtex <- gtex[, c("var_id", "name")]
  sumstats$var_id <- paste0(sumstats$chr, "_", sumstats$pos)

  # first left join, to get annotations of fine-mapped snps
  annot.sumstats <- dplyr::left_join(x = sumstats, y = gtex, by = "var_id")
  # collapse annotation into one
  annots <- annot.sumstats %>%
    dplyr::group_by(var_id) %>%
    dplyr::summarise(!!tissue_type := paste(unique(name), collapse = ";"))
  # join again to get original df with new annotation
  annot.sumstats <- dplyr::inner_join(sumstats, annots, by = "var_id") %>%
    dplyr::select(-var_id)

  return(annot.sumstats)
}

##' @title
##'
##'
##' @param sumstats
##' @param gene.df
##' @param dist
##' @param gene_type
##' @return
##' @export
add.nearby.gene <- function(sumstats, gene.df, dist = 10000, gene_type) {
  gene.ranges <- plyranges::make_ranges(gene.df$chr, gene.df$start, gene.df$end)
  gene.ranges <- plyranges::mutate(gene.ranges, name = gene.df$gene)

  snp.ranges <- GenomicRanges::GRanges(
    seqnames = sumstats$chr,
    ranges = IRanges::IRanges(
      start = sumstats$pos,
      end = sumstats$pos
    )
  )
  snp.ranges <- plyranges::mutate(snp.ranges, snp = sumstats$snp)

  nearby.genes <- plyranges::join_overlap_inner(gene.ranges, snp.ranges, maxgap = dist) %>%
    tibble::as_tibble() %>%
    dplyr::group_by(snp) %>%
    dplyr::summarise(!!gene_type := paste(name, collapse = ";"))

  annot.sumstats <- left_join(sumstats, nearby.genes, by = "snp")
  return(annot.sumstats)
}
##' @title
##'
##'
##' @param sumstats
##' @return
##' @export
add.consequence <- function(sumstats) {
  stopifnot("snp" %in% colnames(sumstats))

  snp.list <- sumstats$snp
  snp.mart <- biomaRt::useMart(
    biomart = "ENSEMBL_MART_SNP",
    dataset = "hsapiens_snp"
  )

  result <- suppressMessages(biomaRt::getBM(
    attributes = c("refsnp_id", "consequence_type_tv"),
    filters = "snp_filter",
    values = snp.list,
    mart = snp.mart
  ) %>% tibble::as_tibble %>% dplyr::group_by(refsnp_id) %>% dplyr::summarise(consequence = paste0(consequence_type_tv, collapse = ";")))

  colnames(result) <- c("snp", "consequence")
  annot.sumstats <- left_join(sumstats, result, by = "snp")
  return(annot.sumstats)
}

##' @title
##'
##'
##' @param sumstats
##' @param gwas
##' @return
##' @export
add.gwas_catalog <- function(sumstats, gwas) {
  gwas <- gwas[, c("SNPS", "DISEASE/TRAIT")]
  colnames(gwas) <- c("snp", "other_gwas_trait")

  annot.sumstats <- left_join(sumstats, gwas, by = "snp") %>%
    dplyr::group_by(snp) %>%
    dplyr::summarise(gwas_catalog = paste0(unique(other_gwas_trait), collapse = ";"))

  sumstats <- left_join(sumstats, annot.sumstats, by = "snp")
  return(sumstats)
}


# BRCA locus zoom link
##' @title
##'
##'
##' @param sumstats
##' @return
##' @export
add.locuszoom.link <- function(sumstats) {
  links <- paste0("https://my.locuszoom.org/gwas/844282/region/?chrom=", sumstats$chr, "&start=", sumstats$pos - 1000, "&end=", sumstats$pos + 1000)
  sumstats$locuszoom <- links
  return(sumstats)
}
##' @title
##'
##'
##' @param sumstats
##' @return
##' @export
add.region <- function(sumstats) {
  loci <- sumstats$locus
  L <- length(loci)
  region <- rep(0, L)
  for (l in 1:L) {
    res <- locus.lookup(loci[l])
    region[l] <- paste0("chr", res$chr, ":", round(res$start / 10^6, 1), "-", round(res$end / 10^6, 1), "M")
  }
  sumstats$region <- region
  return(sumstats)
}

# Hi-C related functions

# HiC data obtained from
# Jung I, Schmitt A, Diao Y, Lee AJ et al. A compendium of promoter-centered long-range chromatin interactions in the human genome. Nat Genet 2019
# does not contain genes for the promoters
# This function assigns genes to the promoters
##' @title
##'
##'
##' @param hic
##' @param gene_cords
##' @return
##' @export
promoters_to_genes <- function(hic, gene_cords) {
  stopifnot("chr_prom" %in% colnames(hic))
  stopifnot("start_prom" %in% colnames(hic))
  stopifnot("end_prom" %in% colnames(hic))

  stopifnot("chr" %in% colnames(gene_cords))
  stopifnot("start" %in% colnames(gene_cords))
  stopifnot("end" %in% colnames(gene_cords))

  # Maps Genes to Promoters
  # GRanges object for promoters
  promoter_df <- hic[, c("chr_prom", "start_prom", "end_prom")] %>% distinct() # keep unique promoters
  promoter.ranges <- plyranges::make_ranges(promoter_df$chr_prom, promoter_df$start_prom, promoter_df$end_prom)
  promoter.ranges <- plyranges::mutate(promoter.ranges,
    promoter_id = paste0(promoter_df$chr_prom, "_", promoter_df$start_prom, "_", promoter_df$end_prom)
  )

  # GRanges object for gene TSS
  ## a genes promoter should be 2kb upstream and 1kb downstream of gene TSS
  gene.ranges <- plyranges::make_ranges(gene_cords$chr, gene_cords$start - 2000, gene_cords$start + 1000)
  gene.ranges <- plyranges::mutate(gene.ranges, gene = gene_cords$gene)

  # overlap to get promoter - gene dictionary
  promoter.genes.link <- plyranges::join_overlap_inner(gene.ranges, promoter.ranges)
  promoter.genes.link <- tibble::as_tibble(promoter.genes.link) %>% dplyr::select(c(gene, promoter_id))

  # add gene to original hic data
  hic$promoter_id <- paste0(hic$chr_prom, "_", hic$start_prom, "_", hic$end_prom)
  hic_with_gene <- dplyr::inner_join(x = hic, y = promoter.genes.link, by = "promoter_id")

  return(hic_with_gene)
}
##' @title
##'
##'
##' @param sumstats
##' @param hic_with_genes
##' @return
##' @export
compute_gene_pip <- function(sumstats, hic_with_genes) {
  stopifnot("chr_enh" %in% colnames(hic_with_genes))
  stopifnot("start_enh" %in% colnames(hic_with_genes))
  stopifnot("end_enh" %in% colnames(hic_with_genes))
  # Map SNPs to Genes via Enhancers
  enhancer.ranges <- plyranges::make_ranges(hic_with_genes$chr_enh, hic_with_genes$start_enh, hic_with_genes$end_enh)
  enhancer.ranges <- plyranges::mutate(enhancer.ranges, gene = hic_with_genes$gene)
  enhancer.ranges <- plyranges::mutate(enhancer.ranges,
    enhancer_id = paste0(hic_with_genes$chr_enh, "_", hic_with_genes$start_enh, "_", hic_with_genes$end_enh)
  )

  snp.ranges <- plyranges::make_ranges(sumstats$chr, sumstats$pos, sumstats$pos)
  snp.ranges <- plyranges::mutate(snp.ranges, snp = sumstats$snp, susie_pip = sumstats$susie_pip)

  hic.snp.intersect <- plyranges::join_overlap_intersect(snp.ranges, enhancer.ranges)
  gene_pip <- tibble::as_tibble(hic.snp.intersect) %>%
    dplyr::group_by(gene) %>%
    dplyr::summarise(PIP = sum(susie_pip))

  return(gene_pip)
}
