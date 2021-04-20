#!/usr/bin/env Rscript

# take arguments from command line
args = commandArgs(trailingOnly=TRUE)

if (length(args)<=1) {
    stop("<peak file> <gmt data folder>.n", call.=FALSE)
}

# load libraries
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)

#library(TxDb.Hsapiens.UCSC.mm38.knownGene)
#library(org.Mm.eg.db)

library(ChIPseeker)
library(patchwork)
library(clusterProfiler)
library(plotly)

pfile <- args[1]
gmtfolder <- args[2]

# Annotate
peakAnno <- annotatePeak( pfile, tssRegion=c(-3000, 3000), TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene, annoDb="org.Hs.eg.db", verbose = 0)

# plot
pdf( paste0(pfile,"_chip_annotation.pdf"), width=10, height=8 )

p1 <- plotAnnoBar(peakAnno) + ggtitle("Distribution of Genomics Features")

p2 <- plotDistToTSS(peakAnno, title="Distribution of peaks relative to TSS")

p1/p2

peakAnno.genes <- as.data.frame(peakAnno)[,c("geneId","ENSEMBL","SYMBOL","GENENAME","distanceToTSS")]

n.10k <- dim(peakAnno.genes[abs(peakAnno.genes$distanceToTSS)<=10000,])[1]
n.5k <- dim(peakAnno.genes[abs(peakAnno.genes$distanceToTSS)<=5000,])[1]
n.1k <- dim(peakAnno.genes[abs(peakAnno.genes$distanceToTSS)<=1000,])[1]
n.100 <- dim(peakAnno.genes[abs(peakAnno.genes$distanceToTSS)<=100,])[1]

d <- data.frame( DAR=c("PSA pos DAR"), D10=c(n.10k), D5=c(n.5k), D1=c(n.1k), D0=c(n.100) )

fig <- plot_ly( d,
  x = ~DAR,
  y = ~D10,
  type = "bar",
  name = "10Kb distance"
)
fig <- fig %>% add_trace( y = ~D5, name = "5Kb distance")
fig <- fig %>% add_trace( y = ~D1, name = "1Kb distance")
fig <- fig %>% add_trace( y = ~D0, name = "100b distance")
fig <- fig %>% layout(yaxis = list(title = 'Count'), barmode = 'group')
fig


# GO/ other association studies

# load GMT files
# hall mark
hgmtfile <- paste0(gmtfolder, "/h.all.entrez.gmt")
hgene <- read.gmt(hgmtfile)

# ontology
c5gmtfile <- paste0(gmtfolder, "/c5.go.bp.entrez.gmt")
bpgene <- read.gmt(c5gmtfile)

# Human Phenotype Ontology
hpogmtfile <- paste0(gmtfolder, "/c5.hpo.entrez.gmt")
hpogene <- read.gmt(hpogmtfile)

# oncogenic signature
c6gmtfile <- paste0(gmtfolder, "/c6.all.entrez.gmt")
oncogene <- read.gmt(c6gmtfile)

# immuneSigDB
c7gmtfile <- paste0(gmtfolder, "/c7.immunesigdb.entrez.gmt")
imgene <- read.gmt(c7gmtfile)

# TF targets Signatures
c3gmtfile <- paste0(gmtfolder, "/c3.tft.entrez.gmt")
tftgene <- read.gmt(c3gmtfile)

# Curated Reactome Pathways
rpgmtfile <- paste0(gmtfolder, "/c2.cp.reactome.entrez.gmt")
rpgene <- read.gmt(rpgmtfile)

# Curated Wikipathways 
wpgmtfile <- paste0(gmtfolder, "/c2.cp.wikipathways.entrez.gmt")
wpgene <- read.gmt(wpgmtfile)

# Curated Biocarta
bcgmtfile <- paste0(gmtfolder, "/c2.cp.biocarta.entrez.gmt")
bcgene <- read.gmt(bcgmtfile)

# Curated KEGG
kegggmtfile <- paste0(gmtfolder, "/c2.cp.kegg.entrez.gmt")
kegggene <- read.gmt(kegggmtfile)

# Curated Canonical Pathways
cpgmtfile <- paste0(gmtfolder, "/c2.cp.entrez.gmt")
cpgene <- read.gmt(cpgmtfile)

# Curated All
c2gmtfile <- paste0(gmtfolder, "/c2.all.entrez.gmt")
curatedgene <- read.gmt(c2gmtfile)

# cutoff
pvalueCutoff = 0.05
pAdjustMethod = "BH"
minGSSize = 5
maxGSSize = 500
qvalueCutoff = 0.2

# hallmark
hallmark <- enricher(as.data.frame(peakAnno)$geneId, TERM2GENE = hgene, pvalueCutoff=pvalueCutoff, pAdjustMethod = pAdjustMethod, minGSSize = minGSSize, maxGSSize = maxGSSize, qvalueCutoff = qvalueCutoff)
if ( dim(hallmark)[1] > 0 ) {
    dotplot(hallmark, title = "MSigDB Hallmark")
}

# The GO geneset Biology Process signatures
bp <- enricher(as.data.frame(peakAnno)$geneId, TERM2GENE = bpgene, pvalueCutoff=pvalueCutoff, pAdjustMethod = pAdjustMethod, minGSSize = minGSSize, maxGSSize = maxGSSize, qvalueCutoff = qvalueCutoff)
if ( dim(bp)[1] > 0 ) {
    dotplot(bp, title = "MSigDB GO BP")
}

# The Human Phenotype Ontology signatures
hpo <- enricher(as.data.frame(peakAnno)$geneId, TERM2GENE = hpogene, pvalueCutoff=pvalueCutoff, pAdjustMethod = pAdjustMethod, minGSSize = minGSSize, maxGSSize = maxGSSize, qvalueCutoff = qvalueCutoff)
if ( dim(hpo)[1] > 0 ) {
    dotplot(hpo, title = "MSigDB Human Phenotype Ontology")
}

# Oncogenic
onco <- enricher(as.data.frame(peakAnno)$geneId, TERM2GENE = oncogene, pvalueCutoff=pvalueCutoff, pAdjustMethod = pAdjustMethod, minGSSize = minGSSize, maxGSSize = maxGSSize, qvalueCutoff = qvalueCutoff)
if ( dim(onco)[1] > 0 ) {
    dotplot(onco, title = "MSigDB Oncogenic")
}

# The ImmuneSigDB  signatures
im <- enricher(as.data.frame(peakAnno)$geneId, TERM2GENE = imgene, pvalueCutoff=pvalueCutoff, pAdjustMethod = pAdjustMethod, minGSSize = minGSSize, maxGSSize = maxGSSize, qvalueCutoff = qvalueCutoff)
if ( dim(im)[1] > 0 ) {
    dotplot(im, title = "MSigDB ImmuneSigDB")
}

# The TF targets signatures
tft <- enricher(as.data.frame(peakAnno)$geneId, TERM2GENE = tftgene, pvalueCutoff=pvalueCutoff, pAdjustMethod = pAdjustMethod, minGSSize = minGSSize, maxGSSize = maxGSSize, qvalueCutoff = qvalueCutoff)
if ( dim(tft)[1] > 0 ) {
    dotplot(tft, title = "MSigDB TF targets")
}

# The Reactome Pathways
rp <- enricher(as.data.frame(peakAnno)$geneId, TERM2GENE = rpgene, pvalueCutoff=pvalueCutoff, pAdjustMethod = pAdjustMethod, minGSSize = minGSSize, maxGSSize = maxGSSize, qvalueCutoff = qvalueCutoff)
if ( dim(rp)[1] > 0 ) {
    dotplot(rp, title = "MSigDB Reactome Pathways")
}

# The WikiPathways signatures
wp <- enricher(as.data.frame(peakAnno)$geneId, TERM2GENE = wpgene, pvalueCutoff=pvalueCutoff, pAdjustMethod = pAdjustMethod, minGSSize = minGSSize, maxGSSize = maxGSSize, qvalueCutoff = qvalueCutoff)
if ( dim(wp)[1] > 0 ) {
    rownames(wp@result) <- gsub(wp@result$ID, pattern = "%.*$", replacement = "")
    wp@result$ID <- gsub(wp@result$ID, pattern = "%.*$", replacement = "")
    wp@result$Description <- gsub(wp@result$Description, pattern = "%.*$", replacement = "")
    dotplot(wp, title = "MSigDB WikiPathways")
}

# The BioCarta Pathways
bc <- enricher(as.data.frame(peakAnno)$geneId, TERM2GENE = bcgene, pvalueCutoff=pvalueCutoff, pAdjustMethod = pAdjustMethod, minGSSize = minGSSize, maxGSSize = maxGSSize, qvalueCutoff = qvalueCutoff)
if ( dim(bc)[1] > 0 ) {
    dotplot(bc, title = "MSigDB BioCarta")
}

# The KEGG Pathways
kegg <- enricher(as.data.frame(peakAnno)$geneId, TERM2GENE = kegggene, pvalueCutoff=pvalueCutoff, pAdjustMethod = pAdjustMethod, minGSSize = minGSSize, maxGSSize = maxGSSize, qvalueCutoff = qvalueCutoff)
if ( dim(kegg)[1] > 0 ) {
    dotplot(kegg, title = "MSigDB BioCarta")
}

# Curated all
curated <- enricher(as.data.frame(peakAnno)$geneId, TERM2GENE = curatedgene, pvalueCutoff=pvalueCutoff, pAdjustMethod = pAdjustMethod, minGSSize = minGSSize, maxGSSize = maxGSSize, qvalueCutoff = qvalueCutoff)
if ( dim(curated)[1] > 0 ) {
    dotplot(curated, title = "MSigDB Curated All")
}

dev.off()
