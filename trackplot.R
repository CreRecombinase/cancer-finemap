library(ggplot2)
library(plotly)
library(patchwork)
library(wesanderson)

locus.lookup <- function(locus){
  euro_ld <- read.delim('~/GERMLINE/ANNOTATIONS/Euro_LD_Chunks.bed',sep = '\t', stringsAsFactors = F, header = F)
  currLocus <- euro_ld[euro_ld$V4 == locus, ]
  chrom <- currLocus$V1
  start <- currLocus$V2
  end <- currLocus$V3
  return(list(locus=locus, chrom=chrom, start=start,end=end))
}

pval.track <- function(df, locus_to_plot){
  
  locusInfo <- locus.lookup(locus_to_plot)
  locus <- locusInfo$locus
  chrom <- locusInfo$chrom
  start <- locusInfo$start
  end <- locusInfo$end

  pos <- df$pos
  currTrack <- log10(df$pval)
  r2val <- df$r2
  
  df <- data.frame(pos = pos, currTrack = currTrack, r2=r2val)
  p1 <- ggplot(df, aes(x = pos, y=currTrack, color=r2)) + 
    geom_point(na.rm=TRUE) + 
    theme_bw() + 
    theme(legend.position = 'right', 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x=element_blank()) + 
    scale_color_gradientn(colors=wes_palette(n=5, name="Zissou1")) +
    xlab('') + 
    ylab(expression(log[10]~pval)) +
    ggtitle(locus) + 
    xlim(c(start, end))
  
  return(p1)
}

pip.track <- function(df, locus_to_plot){
  
  locusInfo <- locus.lookup(locus_to_plot)
  locus <- locusInfo$locus
  chrom <- locusInfo$chrom
  start <- locusInfo$start
  end <- locusInfo$end
  
  pos <- df$pos
  currTrack <- df$susie_pip
  inCS <- factor(df$CS)
  
  df <- data.frame(pos = pos, currTrack = currTrack, Credible_Set=inCS)
  p1 <- ggplot(df, aes(x = pos, y=currTrack, color=Credible_Set)) + 
    geom_point(na.rm=TRUE) + 
    theme_bw() + 
    theme(legend.position = 'none', 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x=element_blank()) + 
    scale_color_manual(values=c('black','green')) +
    ylab('PIP') + 
    xlab('') +
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
  colnames(track.df) <- c("X1","X2","X3","Y1","Y2","Annotations")
  for(annots in bed_annotations){
    name <- sub(basename(annots), pattern = '.bed', replacement = '')
    currAnnot <- vroom::vroom(annots, delim = '\t', col_names = F)
    currAnnotChrom <- currAnnot[currAnnot$X1==chrom, ]
    
    currAnnotChrom$Y1 <- y
    currAnnotChrom$Y2 <- y+1
    currAnnotChrom$Annotations <- name
    
    track.df <- rbind(track.df, currAnnotChrom)
    y <- y + 1
  }
  
  p <- ggplot() + 
    geom_rect(data=track.df, mapping=aes(xmin=X2, xmax=X3, ymin=Y1, ymax=Y2, fill=Annotations), alpha=1, na.rm=TRUE) + 
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

  pval <- pval.track(df = sub.df, locus_to_plot = locus)
  PIP <- pip.track(df = sub.df, locus_to_plot = locus)
  #genes <- suppressMessages(gene.track(locus))
  annot.track <- suppressMessages(suppressMessages(annotation.track(annotations, locus)))
  
  topPipPos <- as.numeric(sub.df[which.max(sub.df$susie_pip), 'pos'])
  
  print(sub.df[which.max(sub.df$susie_pip), c('chr','pos')])
  
  Z <- wrap_plots(pval, PIP, annot.track, nrow=3, heights = c(3,3,4))
  Z <- suppressMessages(Z & geom_vline(xintercept = topPipPos, color='black', linetype="dashed") & xlim(c(topPipPos-50000, topPipPos+50000)))
  
  #Z_zoomed <- suppressMessages(Z & xlim(c(topPipPos-100000, topPipPos+100000)) & theme(legend.position = 'none'))
  #final <- wrap_plots(Z, Z_zoomed, ncol = 2, widths = c(6,4))
  return(Z)
}

