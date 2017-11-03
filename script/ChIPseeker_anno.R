#! /usr/bin/env Rscript

args <- commandArgs(T)
sqlit <- args[1]
gene <- args[2]
name <- args[3]

# sqlit <- '../table//tair10_txdb.sqlite'
# name <- '../../peaks/Zrf/Zrf_2_peaks'

library(GenomicFeatures)
library(ChIPseeker)
txdb <- loadDb(sqlit)

peak <- readPeakFile(paste0(name, '.bed'))
peakAnno <- annotatePeak(peak, TxDb=txdb)
anno <- as.data.frame(peakAnno)
colnames(anno)[6:7] = c('peak', 'score')

genes <- read.table(gene)
colnames(genes) <- c('transcriptId', 'geneSymbol')
anno <- merge(anno, genes)
anno <- anno[, c('seqnames','start','end','width','strand','peak','score','annotation','geneChr','geneStart','geneEnd','geneLength','geneStrand','geneId','transcriptId','geneSymbol','distanceToTSS')]

write.table(anno, paste0(name, '_anno.xls'), sep='\t', quote=F, row.names=F)

pdf(paste0(name, '_annoPie.pdf'), width=7, height=5)
plotAnnoPie(peakAnno)
dev.off()

