#!/usr/bin/env Rscript
# Create User defined trait tables
# - Take KEGG Orthologies and make custom pathways/traits for BugBase
# - Input is a directory containing one file per trait/pathway.
#   Each trait file should contain the KO IDs in the pathway, one per 
#   line. New line character in this file should be '\n'
# - Output will be created within the same directory specified by 
#   the input
# - The pathway name will be the name of each pathway file

# USAGE
# Default
# make.user.table.r 	-i path_to_directory_of_pathway_files 

# Options
# -w 	Data is whole genome shotgun data (picked against IMG database), 
#     default is 16S 

library(optparse)

option_list <- list(
  make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
              help="Print extra output [default]"),
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="path to directory housing the files that list the KO IDs per trait 
              [default %default]"),
  make_option(c("-w", "--wgs"), action="store_true", default=FALSE,
              help="traits are for whole genome sequencing, 
              default is 16S [default %default]")
)
opts <- parse_args(OptionParser(option_list=option_list))

#Set BugBase path
my_env <- Sys.getenv(c('BUGBASE_PATH'))
if(my_env == ""){
  stop("BUGBASE_PATH not set.")
}
db_fp <- paste(my_env, "/usr", sep='/')

#Change into the directory specified
working_dir <- opts$input
setwd(working_dir)

print("loading table, this can take up to 30 minutes")

#Get the precalc table needed
if(opts$wgs){
  precalc_fp <- paste(db_fp, "WGS_KOs_precalculated.txt.gz", 
                      sep='/')
  precalc <- read.table(precalc_fp, 
                        header=T, sep='\t', quote='', row.names=1, 
                        comment='')
} else {
  precalc_fp <- paste(db_fp, "16S_KOs_precalculated.txt.bz2", 
                      sep='/')
  precalc <- read.table(precalc_fp, 
                        header=T, sep='\t', quote='', row.names=1, 
                        comment='')
}

#Get the KO ID files
trait_files <- list.files()

print("creating intermediate files")

#Make output directory
dir.create("intermediates")

#Make a new trait intermediate for each file
#These tables will be the KO count for each OTU ID
for(i in 1:length(trait_files)){
  IDs <- scan(trait_files[i], what="", sep="\n")
  file_name <- paste("intermediates", trait_files[i], sep='/')
  write.table(precalc[,colnames(precalc) %in% IDs,drop=F], 
              file=file_name, sep='\t', col.names = NA, quote=F)
}

#Make one table that combines the intermediate tables 
#This table will be the proportion of each trait covered within each OTU
#Get the precalc files you want to convert
trait_files <- list.files("intermediates/")

#Number of otus in the set is 203,451 (GG, 97%) or 4,755 (WGS, img)
#
#Create ouput matrix
print("creating custom BugBase trait table")

if(opts$wgs){
  output_trait_table <- matrix(0, 4755, length(trait_files))
  output_names <- c()
  for(i in 1:length(trait_files)){
    trait_name <- paste(trait_files[i])
    output_names <- c(output_names, trait_name)
  }
colnames(output_trait_table) <- output_names
} else {
  output_trait_table <- matrix(0, 203451, length(trait_files))
  output_names <- c()
  for(i in 1:length(trait_files)){
    trait_name <- paste(trait_files[i])
    output_names <- c(output_names, trait_name)
  }
  colnames(output_trait_table) <- output_names
}


#Fill in the matrix with the percent of genes present in each otu 
#(out of all genes in trait)
for(i in 1:length(trait_files)){
  working_table <- read.table(paste("intermediates",trait_files[i], sep='/'), 
                              sep='\t', header=TRUE, check.names=FALSE, 
                              comment='', row.names=1)
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
trait_table_name <- "Custom_BugBase_Traits.txt"
write.table(output_trait_table, trait_table_name, sep="\t", quote=F, 
            col.names=NA)
