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

library(ggplot2)
library(ChIPseeker)
library(patchwork)
library(clusterProfiler)

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

peakAnno.genes.10k <- peakAnno.genes[abs(peakAnno.genes$distanceToTSS)<=10000,]

# GO/ other association studies

# load GMT files
# hall mark
hgmtfile <- paste0(gmtfolder, "/h.all.v7.5.1.entrez.gmt")
hgene <- read.gmt(hgmtfile)

# ontology
c5gmtfile <- paste0(gmtfolder, "/c5.go.bp.v7.5.1.entrez.gmt")
bpgene <- read.gmt(c5gmtfile)

# Human Phenotype Ontology
hpogmtfile <- paste0(gmtfolder, "/c5.hpo.v7.5.1.entrez.gmt")
hpogene <- read.gmt(hpogmtfile)

# oncogenic signature
c6gmtfile <- paste0(gmtfolder, "/c6.all.v7.5.1.entrez.gmt")
oncogene <- read.gmt(c6gmtfile)

# immuneSigDB
c7gmtfile <- paste0(gmtfolder, "/c7.immunesigdb.v7.5.1.entrez.gmt")
imgene <- read.gmt(c7gmtfile)

# TF targets Signatures
c3gmtfile <- paste0(gmtfolder, "/c3.tft.v7.5.1.entrez.gmt")
tftgene <- read.gmt(c3gmtfile)

# Curated Reactome Pathways
rpgmtfile <- paste0(gmtfolder, "/c2.cp.reactome.v7.5.1.entrez.gmt")
rpgene <- read.gmt(rpgmtfile)

# Curated Wikipathways 
wpgmtfile <- paste0(gmtfolder, "/c2.cp.wikipathways.v7.5.1.entrez.gmt")
wpgene <- read.gmt(wpgmtfile)

# Curated Biocarta
bcgmtfile <- paste0(gmtfolder, "/c2.cp.biocarta.v7.5.1.entrez.gmt")
bcgene <- read.gmt(bcgmtfile)

# Curated KEGG
kegggmtfile <- paste0(gmtfolder, "/c2.cp.kegg.v7.5.1.entrez.gmt")
kegggene <- read.gmt(kegggmtfile)

# Curated Canonical Pathways
cpgmtfile <- paste0(gmtfolder, "/c2.cp.v7.5.1.entrez.gmt")
cpgene <- read.gmt(cpgmtfile)

## Curated All
#c2gmtfile <- paste0(gmtfolder, "/c2.all.entrez.gmt")
#curatedgene <- read.gmt(c2gmtfile)

# cutoff
pvalueCutoff = 0.05
pAdjustMethod = "BH"
minGSSize = 5
maxGSSize = 500
qvalueCutoff = 0.2

# hallmark
hallmark <- enricher(peakAnno.genes.10k$geneId, TERM2GENE = hgene, pvalueCutoff=pvalueCutoff, pAdjustMethod = pAdjustMethod, minGSSize = minGSSize, maxGSSize = maxGSSize, qvalueCutoff = qvalueCutoff)
if ( !is.null(hallmark) && dim(hallmark)[1] > 0 ) {
    dotplot(hallmark, title = "MSigDB Hallmark")
}

# The GO geneset Biology Process signatures
bp <- enricher(peakAnno.genes.10k$geneId, TERM2GENE = bpgene, pvalueCutoff=pvalueCutoff, pAdjustMethod = pAdjustMethod, minGSSize = minGSSize, maxGSSize = maxGSSize, qvalueCutoff = qvalueCutoff)
if ( !is.null(bp) && dim(bp)[1] > 0 ) {
    dotplot(bp, title = "MSigDB GO BP")
}

# The Human Phenotype Ontology signatures
hpo <- enricher(peakAnno.genes.10k$geneId, TERM2GENE = hpogene, pvalueCutoff=pvalueCutoff, pAdjustMethod = pAdjustMethod, minGSSize = minGSSize, maxGSSize = maxGSSize, qvalueCutoff = qvalueCutoff)
if ( !is.null(hpo) && dim(hpo)[1] > 0 ) {
    dotplot(hpo, title = "MSigDB Human Phenotype Ontology")
}

# Oncogenic
onco <- enricher(peakAnno.genes.10k$geneId, TERM2GENE = oncogene, pvalueCutoff=pvalueCutoff, pAdjustMethod = pAdjustMethod, minGSSize = minGSSize, maxGSSize = maxGSSize, qvalueCutoff = qvalueCutoff)
if ( !is.null(onco) && dim(onco)[1] > 0 ) {
    dotplot(onco, title = "MSigDB Oncogenic")
}

# The ImmuneSigDB  signatures
im <- enricher(peakAnno.genes.10k$geneId, TERM2GENE = imgene, pvalueCutoff=pvalueCutoff, pAdjustMethod = pAdjustMethod, minGSSize = minGSSize, maxGSSize = maxGSSize, qvalueCutoff = qvalueCutoff)
if ( !is.null(im) && dim(im)[1] > 0 ) {
    dotplot(im, title = "MSigDB ImmuneSigDB")
}

# The TF targets signatures
tft <- enricher(peakAnno.genes.10k$geneId, TERM2GENE = tftgene, pvalueCutoff=pvalueCutoff, pAdjustMethod = pAdjustMethod, minGSSize = minGSSize, maxGSSize = maxGSSize, qvalueCutoff = qvalueCutoff)
if ( !is.null(tft) && dim(tft)[1] > 0 ) {
    dotplot(tft, title = "MSigDB TF targets")
}

# The Reactome Pathways
rp <- enricher(peakAnno.genes.10k$geneId, TERM2GENE = rpgene, pvalueCutoff=pvalueCutoff, pAdjustMethod = pAdjustMethod, minGSSize = minGSSize, maxGSSize = maxGSSize, qvalueCutoff = qvalueCutoff)
if ( !is.null(rp) && dim(rp)[1] > 0 ) {
    dotplot(rp, title = "MSigDB Reactome Pathways")
}

# The WikiPathways signatures
wp <- enricher(peakAnno.genes.10k$geneId, TERM2GENE = wpgene, pvalueCutoff=pvalueCutoff, pAdjustMethod = pAdjustMethod, minGSSize = minGSSize, maxGSSize = maxGSSize, qvalueCutoff = qvalueCutoff)
if ( !is.null(wp) && dim(wp)[1] > 0 ) {
    rownames(wp@result) <- gsub(wp@result$ID, pattern = "%.*$", replacement = "")
    wp@result$ID <- gsub(wp@result$ID, pattern = "%.*$", replacement = "")
    wp@result$Description <- gsub(wp@result$Description, pattern = "%.*$", replacement = "")
    dotplot(wp, title = "MSigDB WikiPathways")
}

# The BioCarta Pathways
bc <- enricher(peakAnno.genes.10k$geneId, TERM2GENE = bcgene, pvalueCutoff=pvalueCutoff, pAdjustMethod = pAdjustMethod, minGSSize = minGSSize, maxGSSize = maxGSSize, qvalueCutoff = qvalueCutoff)
if ( !is.null(bc) && dim(bc)[1] > 0 ) {
    dotplot(bc, title = "MSigDB BioCarta")
}

# The KEGG Pathways
kegg <- enricher(peakAnno.genes.10k$geneId, TERM2GENE = kegggene, pvalueCutoff=pvalueCutoff, pAdjustMethod = pAdjustMethod, minGSSize = minGSSize, maxGSSize = maxGSSize, qvalueCutoff = qvalueCutoff)
if ( !is.null(kegg) && dim(kegg)[1] > 0 ) {
    dotplot(kegg, title = "MSigDB KEGG")
}

## Curated all
#curated <- enricher(peakAnno.genes.10k$geneId, TERM2GENE = curatedgene, pvalueCutoff=pvalueCutoff, pAdjustMethod = pAdjustMethod, minGSSize = minGSSize, maxGSSize = maxGSSize, qvalueCutoff = qvalueCutoff)
#if ( dim(curated)[1] > 0 ) {
#    dotplot(curated, title = "MSigDB Curated All")
#}

dev.off()
