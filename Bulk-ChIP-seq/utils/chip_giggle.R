#!/usr/bin/env Rscript

# take arguments from command line
args = commandArgs(trailingOnly=TRUE)

if (length(args)<3) {
    stop("<peak file> <giggle exec path> <giggle index folder>.n", call.=FALSE)
}

# load libraries
library(ggplot2)
library(ggrepel)

pfile <- args[1]
giggle <- args[2]
gigglepath <- args[3]

ant <- read.csv(paste0(gigglepath,"/CistromeDB.sample.annotation.txt"), sep="\t", row.names=1, stringsAsFactors = FALSE)

RunGiggle <- function(peakbed, giggle.exec, giggle.path, organism, antFile){
  outputBed = peakbed
  cmd <- paste0("sort --buffer-size 2G -k1,1 -k2,2n -k3,3n ", outputBed, " | cut -f 1,2,3 | bgzip -c > ", outputBed, ".bed.gz")
  system(cmd)
  cmd <- paste0(giggle.exec, " search -g 2700000000 -i ", giggle.path, "/giggle.", organism, " -q ", outputBed, ".bed.gz -s > ", outputBed, ".result.xls")
  system(cmd)
  resultDf <- read.table(paste0(outputBed, ".result.xls"), sep="\t", row.names=NULL, comment.char="", stringsAsFactors =  FALSE)
  resultDf <- resultDf[,-9]
  colnames(resultDf) <- c("file", "file_size", "overlaps", "odds_ratio", "fishers_two_tail", "fishers_left_tail", "fishers_right_tail", "combo_score")
  resultDf <- resultDf[resultDf$overlaps>0,]
  rownames(resultDf) <- sapply(strsplit(resultDf$file, "/"), function(x) return(gsub("_5foldPeak.bed.gz", "", rev(x)[1])))
  resultDf <- resultDf[,c("file_size", "overlaps", "combo_score", "odds_ratio", "fishers_right_tail")]
  targetDf <- merge(resultDf, antFile, by.x=0, by.y=0)
  colnames(targetDf) <- c("sample_id", "sample_peak_number", "overlap_peak_number", "giggle_score", "odds_ratio", "fisher_right_tail", "GSM_id", "species", "factor", "cell_line", "cell_type", "tissue", "disease")

  targetDf$biological_resource <- apply(targetDf, 1, function(x) return(paste0(x[7:10], collapse=";")))
  targetDf$labels <- paste0(targetDf$"cell_line",":",targetDf$"factor")
  targetDf$fisher_right_tail <- -log10(targetDf$"fisher_right_tail")
  targetDf$"odds_ratio" <- log10(targetDf$"odds_ratio")  
  
  targetDf <- targetDf[, c("sample_id", "GSM_id", "species", "factor", "biological_resource", "labels", "giggle_score", "odds_ratio", "fisher_right_tail", "sample_peak_number", "overlap_peak_number")]
  targetDf <- targetDf[order(-targetDf$giggle_score), ]
  targetDf <- targetDf[!duplicated(targetDf$factor), ]
  write.table(targetDf, paste0(outputBed, ".giggle.res.tfs.txt"), sep="\t", quote=FALSE, row.names=FALSE)
  cmd <- paste0("rm ", outputBed, ".bed.gz")
  system(cmd)
  cmd <- paste0("rm ", outputBed, ".result.xls")
  system(cmd)
  return(targetDf)
}

giggle.res <- RunGiggle( pfile, giggle, gigglepath, "GRCh38", ant)

giggle.res.top10 <- giggle.res[order(giggle.res$giggle_score,decreasing = T),][1:10,]

t <- paste0(giggle.res.top10$biological_resource,collapse = "\n")

pdf( paste0(pfile,"_giggle_cistrome.pdf"), width=10, height=8 )

p <- ggplot(data=giggle.res, aes(x=odds_ratio, y=fisher_right_tail, label=labels)) + geom_point()
p + geom_label_repel(force=20, size=2, data = giggle.res.top10) + xlab("Giggle log10 odds_ratio") + ylab("Giggle Fisher-test right-tail -log10 pvalue") + annotate(geom = 'text', label = t, x = Inf, y = -Inf, hjust = 1.2, vjust = -0.5)

dev.off()
