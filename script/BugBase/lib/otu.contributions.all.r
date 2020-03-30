# Creates a table of which OTUs per sample are contributing to the trait 
#   predictions and also plots the OTUs as a bar chart
# Inputs: otus_contributing, otu table, taxonomy
# Returns: one plot per trait (pdf) and table of OTUs 

"otu.contributions.all.r" <- function(otus_contributing, otu_table, taxonomy,
									taxa_level=NULL, output_fp=NULL){
	#set directory
	dir <- paste(output, "otu_contributions", sep="/")
	
	contributing_name <- paste(dir, "contributing_otus.txt", sep='/')
	write.table(otus_contributing, contributing_name, sep='\t', quote=F, 
							col.names=NA)
	
	#define traits
	traits <- colnames(otus_contributing)
	
	#set taxonomy level
	if(is.null(taxa_level)){
		taxa_level <- 2
	}
	taxa_level <- as.numeric(taxa_level)
	
	#read in gg_taxonomy table
	gg_taxonomy <- read.table(taxonomy, sep="\t", row.names=1, check=F, quote='')
	
	#keep only otus in the gg taxonomy table that are in the otu table
	#these tables will be in the same otu order
	#OTUs are now rows
	otu_table <- t(otu_table)
	otus_keep <- intersect(rownames(otu_table),rownames(gg_taxonomy))
	gg_taxonomy <- gg_taxonomy[otus_keep,, drop=F]
	
	#subset the gg taxonomy to the level specified (default=2)
	names_split <- array(dim=c(length(gg_taxonomy[,1]), 7))
	otu_names <- as.character(gg_taxonomy[,1])
	for(i in 1:length(otu_names)){
		names_split[i,] <- strsplit(otu_names[i], ";", fixed=T)[[1]]
	}
	otu_names <- names_split[,taxa_level]
	for(i in 1:length(otu_names)){
	  otu_names[i] <- strsplit(otu_names[i], "__", fixed=T)[[1]][2]
	}
	names_split[,taxa_level] <- otu_names
	for(i in 1:nrow(names_split)){
		if(is.na(names_split[i,taxa_level])){
			if(taxa_level > 1){
				names_split[i, taxa_level] <- names_split[i, taxa_level -1]
			} else {
				names_split[i, taxa_level] <- "unknown"
			}
		}
	}

	#store the full otu_table
	otu_table1 <- otu_table

	#add taxonomy as the rownames in the otu table
	rownames(otu_table) <- names_split[,taxa_level]
	
	#aggregate to the same taxa
	#you must t() again, to have samples as columns
	if(ncol(otu_table) ==1){
		otu_table <- t(t(sapply(by(otu_table,rownames(otu_table),colSums),identity)))
	} else {
		otu_table <- t(sapply(by(otu_table,rownames(otu_table),colSums),identity))
	}

	#set colors for plotting and legend creation
	#get palette 1 from R ColorBrewer
	cols <- colorRampPalette(brewer.pal(9,'Set1'))
	
	#create as many colors as there are taxa, and name them with the taxa names
	cols2 <- cols(length(rownames(otu_table)))
	names(cols2) <- unique(rownames(otu_table))
	cols2 <- c(cols2,"#C0C0C0")
	names(cols2)[length(cols2)] <- "Other"

	taxa_list <- c()
	
	if(is.null(clr_trans)){
		xlabel <- "Relative Abundance"
	} else {
		xlabel <- "CLR Transformed Relative Abundance"
	}

	#create taxa summaries
	for(x in 1:length(traits)){
		trait <- traits[x]
		
		#multiple the otu_table by the boolean table for otus contributing to 
		# the trait
		positive_otu_table <- t(sweep(t(otu_table1), 2, 
									otus_contributing[,x],"*"))

		#add taxonomy as the rownames in the otu table
		rownames(positive_otu_table) <- names_split[,taxa_level]
  		
		#aggregate to the same taxa
		#you must t() again, to have samples as columns
		if(ncol(positive_otu_table) == 1){
			positive_otu_table <- t(t(sapply(by(positive_otu_table,
							rownames(positive_otu_table),colSums),identity)))
		} else {
			positive_otu_table <- t(sapply(by(positive_otu_table,
							rownames(positive_otu_table),colSums),identity))
		}
		

		#melt the otu_table and collapse by Taxa
		melted_otu_table <- melt(positive_otu_table)
		colnames(melted_otu_table) <- c("Taxa", "SampleID", "Count")
		melted_otu_table$SampleID <- as.factor(1)
		taxa_collapsed <- ddply(melted_otu_table, .(Taxa, SampleID),summarize, 
								Count = mean(Count))
		#set value for cutoff (1/10 of the highest proportion)
		max_abund <- max(taxa_collapsed$Count)
		cutoff_val <- max_abund / 10

		#call taxa that are less than cutoff of the population "Other"
		taxa_collapsed$Taxa <- as.character(taxa_collapsed$Taxa)
		taxa_collapsed$Count <- as.numeric(taxa_collapsed$Count)

		taxa_collapsed[which(taxa_collapsed$Count < cutoff_val),
								"Taxa"] <- "Other"
		#re-collapse to group the 'Others'
		taxa_collapsed <- ddply(taxa_collapsed, .(Taxa, SampleID), 
								Count = sum(Count))

		taxa_list <- c(taxa_list, unique(taxa_collapsed$Taxa))

		#make the plot
		taxa_plot <- NULL
		taxa_plot <- ggplot(taxa_collapsed, aes_string(x="SampleID",
			y="Count", fill="Taxa")) + 
			geom_bar(stat="identity", show_guide=F) + 
			labs(y = xlabel, x = "") +
			theme_classic() +
			theme(axis.text.x=element_blank(), axis.line.x = 
				element_line(colour = 'black', size=0.5, linetype='solid'), 
				axis.line.y = element_line(colour = 'black', size=0.5, 
				linetype='solid')) +
			scale_fill_manual(values=cols2)
		
		#assign pdf name
		file <- c(".pdf")
		name <- paste(trait, ".pdf", sep='')
		name <- paste(dir, name, sep="/")
		
		#make the pdf
		pdf(name, height=6,width=6)
		par(mar=c(8,4,0.5,6), oma=c(0.1,0.1,0.1,0.1), mgp=c(1.5,0.5,0))
		
		# Plot the taxa summary
		print(taxa_plot)
		dev.off()
	}

	#subset to the taxa that met the cutoff
	taxa_list <- unique(taxa_list)

	cols_keep <- cols2[which(names(cols2) %in% taxa_list)]

	#create a legend table with names, colors and coordinates
	legend_info <- matrix(0, length(cols_keep), 4)
	legend_info[,1] <- names(cols_keep)
	legend_info[,4] <- 1
	counter <- c(1:length(cols_keep)+1)
	for(i in 1:length(cols_keep)){
		legend_info[i,2] <- cols_keep[[i]]
		legend_info[i,3] <- counter[i]
	}
	colnames(legend_info) <- c("Taxa", "Color", "Y", "X")
	legend_info <- as.data.frame(legend_info)
	legend_info$X <- as.numeric(legend_info$X)
	legend_info$Y <- as.numeric(legend_info$Y)
	legend_info$Color <- as.character(legend_info$Color)
	#Add a space before name so plot legend looks nice
	legend_info$Taxa <- sub("^", " ", legend_info$Taxa )
	
	#name pdf
	taxa_legend <- c("taxa_legend.pdf")
	legend_name <- paste(dir, taxa_legend, sep="/")
	
	#make the pdf
	pdf(legend_name, height=6,width=6)
	par(mar=c(0.5,0.5,0.5,0.5), oma=c(0.1,0.1,0.1,0.1), mgp=c(1.5,0.5,0))
	
	#plot the points
	plot(legend_info$X, legend_info$Y, 
			 pch=15, 
			 cex=3, 
			 col= legend_info$Color, 
			 axes=FALSE, 
			 xlab='', 
			 ylab='')
	#add names
	text(legend_info$X, legend_info$Y, legend_info$Taxa, pos=4)
	graphics.off()
}
