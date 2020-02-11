library(ggplot2)
library(plotly)
library(patchwork)

locus.lookup <- function(locus){
  euro_ld <- read.delim('~/GERMLINE/ANNOTATIONS/Euro_LD_Chunks.bed',sep = '\t', stringsAsFactors = F, header = F)
  currLocus <- euro_ld[euro_ld$V4 == locus, ]
  chrom <- currLocus$V1
  start <- currLocus$V2
  end <- currLocus$V3
  return(list(locus=locus, chrom=chrom, start=start,end=end))
}

numeric.track <- function(df, pos_col, track_col, locus_col, locus_to_plot, log.track=FALSE, track.name='x'){
  
  stopifnot(track_col %in% colnames(df))
  stopifnot(locus_col %in% colnames(df))
  stopifnot(pos_col %in% colnames(df))
  stopifnot(is.numeric(unlist(df[track_col])))
  
  locusInfo <- locus.lookup(locus_to_plot)
  locus <- locusInfo$locus
  chrom <- locusInfo$chrom
  start <- locusInfo$start
  end <- locusInfo$end
  
  sub.df <- df[df[locus_col] == locus_to_plot, ]
  pos <- unlist(sub.df[pos_col])
  currTrack <- unlist(sub.df[track_col])
  if(log.track){
    currTrack <- log10(currTrack)
  }
  
  df <- data.frame(pos = pos, currTrack = currTrack)
  
  p1 <- ggplot(df, aes(x = pos, y=currTrack)) + 
    geom_point(na.rm=TRUE) + 
    theme_bw() + 
    theme(legend.position = 'none', 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x=element_blank()) + 
    ylab(track.name) + xlab('') +
    xlim(c(start, end))
  
  return(p1)
}

annotation.track <- function(bed_annotations, locus){
  
  locusInfo <- locus.lookup(locus)
  
  chrom <- locusInfo$chrom
  start <- locusInfo$start
  end <- locusInfo$end
  
  y <- 0
  track.df <- data.frame(matrix(NA, nrow=0, ncol=6))
  colnames(track.df) <- c("X1","X2","X3","Y1","Y2","annot")
  for(annots in bed_annotations){
    name <- sub(basename(annots), pattern = '.bed', replacement = '')
    currAnnot <- vroom::vroom(annots, delim = '\t', col_names = F)
    currAnnotChrom <- currAnnot[currAnnot$X1==chrom, ]
    
    currAnnotChrom$Y1 <- y
    currAnnotChrom$Y2 <- y+1
    currAnnotChrom$annot <- name
    
    track.df <- rbind(track.df, currAnnotChrom)
    y <- y + 1
  }
  
  p <- ggplot() + 
    geom_rect(data=track.df, mapping=aes(xmin=X2, xmax=X3, ymin=Y1, ymax=Y2, fill=annot), alpha=1, na.rm=TRUE) + 
    theme_void() +
    theme(plot.title = element_text(size=10), legend.position='bottom') +
    xlim(c(start, end))
  
  return(p)
}

gene.track <- function(locus){
  
  gene_cords <- vroom::vroom('~/GERMLINE/refGenome/GENES_COORDs_b37.txt', delim = '\t', col_names = F)
  locusInfo <- locus.lookup(locus)
  
  chrom <- locusInfo$chrom
  start <- locusInfo$start
  end <- locusInfo$end
  
  currChrom <- gene_cords[gene_cords$X1==chrom & gene_cords$X2>start & gene_cords$X3<end, ]
  currChrom <- currChrom[order(currChrom$X2), ] # sort by start
  
  inChain <- c(1)
  outChain <- c()
  currEnd <- currChrom$X3[1]
  
  for(i in 2:nrow(currChrom)){
    if(currChrom$X2[i] > currEnd){
      inChain <- c(inChain, i)
      currEnd <- currChrom$X3[i]
    }else{
      outChain <- c(outChain, i)
    }
  }
  
  currChrom[inChain, "Y1"] <- 0
  currChrom[inChain, "Y2"] <- 0.2
  
  currChrom[outChain, "Y1"] <- 0.8
  currChrom[outChain, "Y2"] <- 1
  
  p <- ggplot() +
       geom_rect(data=currChrom, mapping=aes(xmin=X2, xmax=X3, ymin=Y1, ymax=Y2), color = "blue", fill="lightblue", alpha=1, na.rm=TRUE) +
       #geom_text(data=currChrom, aes(x=X2+(X3-X2)/2, y=Y2+0.1, label=X4), size=3, na.rm = T) + 
       ggrepel::geom_text_repel(data=currChrom, mapping = aes(x=X2+(X3-X2)/2, y=Y1+(Y2-Y1)/2, label=X4), na.rm = T,segment.colour = "black", min.segment.length = 0) + 
       theme_void() +
       theme(plot.title = element_text(size=10), legend.position='bottom') +
       xlim(c(start, end))
  
  return(p)
  
}

locusPlotter <- function(df, locus, annotations){
  
  sub.df <- df[df$locus == locus, ]
  pval <- numeric.track(sub.df, pos_col = 'pos', track_col = 'pval', locus_col = 'locus', locus_to_plot = locus, log.track = TRUE, track.name=expression(log[10]~pval)) + ggtitle(locus)
  PIP <- numeric.track(sub.df, pos_col = 'pos', track_col = 'susie_pip', locus_col = 'locus', locus_to_plot = locus, log.track = FALSE, track.name='PIP')
  
  genes <- suppressMessages(gene.track(locus))
  annot.track <- suppressMessages(suppressMessages(annotation.track(annotations, locus)))
  
  topPipPos <- as.numeric(sub.df[which.max(sub.df$susie_pip), 'pos'])

  Z <- wrap_plots(pval, PIP, annot.track, genes, nrow=4, heights = c(2,2,4,2))
  
  Z <- Z & geom_vline(xintercept = topPipPos, color='black', linetype="dashed")
  
  Z_zoomed <- suppressMessages(Z & xlim(c(topPipPos-100000, topPipPos+100000)) & theme(legend.position = 'none'))
  
  final <- wrap_plots(Z, Z_zoomed, ncol = 2, widths = c(6,4))
  
  return(final)
}


