# Loads the original otu table and 16S copy number file.
# OTU table should be in txt format.
# Uses the 16S copy number table and otu table to normalize the 16 
#    copy number
# Required: 16S table (copyNo_table_fp), otu table (otu_table_fp)
# Returns: normalized_otus

"copyNo.normalize.otu" <- function(copyNo_table, otu_table_loaded, output){
	
	#Load the 16S table and the otu table
	#16S_table is otu x copy number
	copyNo_table <- as.matrix(read.table(copyNo_table, 
										comment='', 
										sep='\t', 
										head=TRUE, 
										row=1, 
										check=F))

	#otu_table is otu x sample
	#otu_table <- as.matrix(read.table(otu_table_fp, sep='\t', head=T, row=1, 
	#  check=F, comment='', skip=1))
	otu_table <- otu_table_loaded

	#Subset 16S table and otu table to include only otus within the 
	#	otu table
	otus_keep <- intersect(rownames(otu_table),rownames(copyNo_table))
	otu_table <- otu_table[otus_keep,,drop=F]
	copyNo_table <- copyNo_table[otus_keep,,drop=F]

	#make 16S normalized otu table by dividing each count by the copy number
	#	normalized_otus is otus by samples
	normalized_otus <-  sweep(otu_table, 1, copyNo_table, "/")

	#set directory to write file
	dir <- paste(output, "normalized_otus", sep="/")
	copy_name <- paste(dir, "16s_normalized_otus.txt", sep='/')

	#write normalized otu table to file  
	write.table(normalized_otus, copy_name, sep="\t", quote=F, 
							col.names=NA) 

	return(normalized_otus=normalized_otus)
}
