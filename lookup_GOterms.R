
# simple script demonstrating GO IDs -> gene entrez ID -> gene symbol
library(annotate)
library(org.Hs.eg.db)
library(GO.db)

go.terms<- Term(GOTERM) 
goID <- GOID(GOTERM[go.terms == 'immune response to tumor cell' | go.terms == 'immune system process'])

gene.entrez <- get(goID, org.Hs.egGO2ALLEGS)

immune.genes <- unname(getSYMBOL(gene.entrez, data = 'org.Hs.eg'))
immune.genes <- enframe(x = immune.genes, value = 'gene', name = NULL)

b37.genes <- vroom::vroom('~/GERMLINE/refGenome/GENES_COORDs_b37.txt', col_names = T)

immune.genes.cords <- inner_join(b37.genes, immune.genes, by='gene') %>% distinct()
vroom::vroom_write(immune.genes.cords, path = '~/GERMLINE/refGenome/immune_genes_cords.txt', delim = '\t', col_names = T)
