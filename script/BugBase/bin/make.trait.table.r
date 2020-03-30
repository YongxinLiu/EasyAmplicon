#!/usr/bin/Rscript
#Converts the precalculated PICRUSt tables into BugBase inputs

# USAGE
# make.trait.table.r -i path_to_file(s) -o name_output
# NOTE: the directory listed as input should contain ONLY files you want to convert to a trait table
library(optparse)

option_list <- list(
  make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
              help="Print extra output [default]"),
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="path to the files you want to convert [default %default]"),
  make_option(c("-o", "--output"), type="character", default=".",
              help="output table name [default %default]")
)
opts <- parse_args(OptionParser(option_list=option_list))

#Change into the directory specified
working_dir <- opts$input
setwd(working_dir)

#Get the precalc files you want to convert
trait_files <- list.files()

#Number of otus in the set is 203,452 (GG, 97%)
#Create ouput matrix
output_trait_table <- matrix(0, 203452, length(trait_files))
output_names <- c()
for(i in 1:length(trait_files)){
  trait_name <- strsplit(trait_files[i], ".", fixed=T)[[1]][1]
  output_names <- c(output_names, trait_name)
}
colnames(output_trait_table) <- output_names

#Fill in the matrix with the percent of genes present in each otu (out of all genes in trait)
for(i in 1:length(trait_files)){
  working_table <- read.table(trait_files[i], sep='\t', header=TRUE, check.names=FALSE, comment='', row.names=1)
  if(i == 1){
    rownames(output_trait_table) <- rownames(working_table)
  } else {
    intersect_btwn <- intersect(rownames(output_trait_table),rownames(working_table))
    working_table2 <- working_table[intersect_btwn,,drop=F]
  }
  if(ncol(working_table) > 1){
      output_trait_table[,i] <- rowSums(working_table) / ncol(working_table)
  } else {
      output_trait_table[,i] <- working_table[,1]
  }
}

#Write the new table
trait_table_name <- opts$output
write.table(output_trait_table, trait_table_name, sep="\t", quote=F, 
            col.names=NA)
