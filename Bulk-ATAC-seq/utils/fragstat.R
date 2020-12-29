#!/usr/bin/env Rscript

library(ggplot2)
library(optparse)

option_list = list(
  make_option(c("-i", "--c1-name"), type="character", default=NULL, 
              help="condition1 name", dest="c1name", ),
  make_option(c("-j", "--c2-name"), type="character", default=NULL, 
              help="condition2 name", dest="c2name", ),
  make_option(c("-a", "--c1-files"), type="character", default=NULL, 
              help="condition1 files, seperated by ','", dest="c1file", ),
  make_option(c("-b", "--c2-files"), type="character", default=NULL, 
              help="condition2 files, seperated by ','", dest="c2file", ),
  make_option(c("-o", "--out"), type="character", default="out.txt", 
              help="output image file name (e.g. a png file)", dest="output_fig")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if ( is.null(opt$c1name) | is.null(opt$c2name) | is.null(opt$c1file) | is.null(opt$c2file) | is.null(opt$output_fig) ){
  print_help(opt_parser)
  stop("-i, -j, -a, -b, and -o are all required!", call.=FALSE)
}

c1_files <- unlist( strsplit( opt$c1file, "," ) )
c2_files <- unlist( strsplit( opt$c2file, "," ) )
c1_name <- opt$c1name
c2_name <- opt$c2name
output_fig <- opt$output_fig

all_data <- data.frame()

for ( f in c1_files ) {
  fname <- basename( f )
  sample_name <- gsub( ".sortedByPos.rmdp.clean.unique.bam.insertl.txt", "", fname )
  group_name <- paste0("C1:",c1_name)
  tmp_d <- read.table( f, col.names = c("fraglen","count"))
  tmp_d$sample <- sample_name
  tmp_d$group <- group_name
  all_data <- rbind( all_data, tmp_d )
}

for ( f in c2_files ) {
  fname <- basename( f )
  sample_name <- gsub( ".sortedByPos.rmdp.clean.unique.bam.insertl.txt", "", fname)
  group_name <- paste0("C2:",c2_name)
  tmp_d <- read.table( f, col.names = c("fraglen","count"))
  tmp_d$sample <- sample_name
  tmp_d$group <- group_name
  all_data <- rbind( all_data, tmp_d )
}

p <- ggplot(data=all_data, aes(x=fraglen, y=count, color=sample))
p <- p + geom_line() + facet_grid(cols = vars(group))
ggsave(output_fig, p, height=4, width=10)
