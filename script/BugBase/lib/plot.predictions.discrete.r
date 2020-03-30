# Plots predictions as beeswarm plots using user-define groups and compares 
#   trait relative abundance between groups using non-parametric tests.
# Required: prediction table, mapping file, map column, output filepath
# Returns: pdfs of plots and text files of statistics

"plot.predictions.discrete" <- function(predictions, 
											map, 
											map_column,
											clr_trans,
											output_fp=NULL){

	#set directory
	dir <- paste(output, "predicted_phenotypes", sep="/")

	#cat("\n\n RESULTS: \n\n")
	#check dimensions of traits
	#cat("\nDimensions of trait table:\n")
	#cat(dim(predictions))
	#cat("Dimensions of the mapping file:\n")
	#cat(dim(map))
	
	#ensure same order of samples in map and traits
	map <- map[rownames(predictions),,drop=F]
	
	#define traits
	traits <- colnames(predictions)

	#define groups
	groups <- sort(unique(map[,map_column]))
	groups <- lapply(groups, as.character)

	#cat("\nNumber of samples in each treatment group:\n")
	#cat(table(map[,map_column]))
  
	for(x in 1:length(traits)){
		trait <- traits[x]

		#show number of samples in each body site and trait
		header <- paste("\n\n", trait, sep='')
		# cat(header)
		# cat("\nRelative Abundance with trait (mean):\n")
		# print(tapply(predictions[,trait], map[,map_column], mean))

		# cat("\nRelative Abundance with trait (median):\n")
		# print(tapply(predictions[,trait], map[,map_column], median))

		# cat("\nStandard deviation:\n")
		# print(tapply(predictions[,trait], map[,map_column], sd))

		#non-parametric tests - either two classes or multi-class
		#print to screen and file
		if(length(groups)==2){
			group.pvalue <- wilcox.test(predictions[,trait] ~ 
										map[,map_column])$p.value
		
			# cat("\np-value is:\n")
			# cat(group.pvalue,"\n")
		
			# cat("FDR-corrected p-value is:\n")
			# cat(p.adjust(group.pvalue,'fdr'), "\n")
			
			outfile <- paste(trait, "_stats.txt", sep="")

			sink(paste(dir,outfile,sep='/'))
	
			cat("\nNumber of samples in each treatment group:\n")
			print(table(map[,map_column]))

			cat("\nRelative Abundance with trait (mean):\n")
			print(tapply(predictions[,trait], map[,map_column], mean))

			cat("\nRelative Abundance with trait (median):\n")
			print(tapply(predictions[,trait], map[,map_column], median))

			cat("\nStandard deviation:\n")
			print(tapply(predictions[,trait], map[,map_column], sd))
			
			cat("\nMann-Whitney-Wilcoxon Test was performed.\n")
			cat("p-values is:\n")
			cat(group.pvalue, "\n")
			
			cat("\nFDR-corrected p-value is:\n")
			cat(p.adjust(group.pvalue,'fdr'),"\n")

			sink()
			
		} else {
			group.pvalue <- kruskal.test(predictions[,trait] ~ 
											map[,map_column])$p.value

			#run pairwise tests if > 2 categories
			pw.pvalues <- NULL
			pw.names <- NULL
			for(i in 1:(length(groups) - 1)){
				for(j in (i+1):length(groups)){
					ix.trait.i <- map[,map_column] == groups[i]
					ix.trait.j <- map[,map_column] == groups[j]
					
					pvalue <- wilcox.test(predictions[ix.trait.i,trait], 
										predictions[ix.trait.j,trait])$p.value
					
					pw.pvalues <- c(pw.pvalues, pvalue)
					test.name <- paste(groups[i], "_vs_", groups[j],sep='')
					pw.names <- c(pw.names, test.name)
				}
			}
			names(pw.pvalues) <- pw.names
			
			# cat("\nPairwise p-values are:\n")
			# cat(pw.pvalues)
			
			# cat("\nFDR-corrected pairwise p-values are:\n")
			# cat(p.adjust(pw.pvalues,'fdr'))
			
			outfile <- paste(trait, "_stats.txt", sep="")
			
			sink(paste(dir,outfile,sep='/'))
	  
			cat(header)
			cat("\nNumber of samples in each treatment group:\n")
			print(table(map[,map_column]))

			cat("\nProportion with phenotype (mean):\n")
			print(tapply(predictions[,trait], map[,map_column], mean))

			cat("\nProportion with phenotype (median):\n")
			print(tapply(predictions[,trait], map[,map_column], median))

			cat("\nStandard deviation:\n")
			print(tapply(predictions[,trait], map[,map_column], sd))
			
			cat("\nPairwise Mann-Whitney-Wilcoxon Tests were performed.\n")
			cat("Pairwise p-values are:\n")
			print(pw.pvalues,"\n", digits=5, quote=FALSE)
			
			cat("\nFDR-corrected pairwise p-values are:\n")
			print(p.adjust(pw.pvalues,'fdr'))

			cat("\nKruskal-Wallis Test was performed.\n")
			cat("\nGroup p-value is:\n")
			cat(group.pvalue, "\n")

			sink()
		}

	#get palette 1 from R ColorBrewer
	cols <- sprintf('%s97',brewer.pal(9,'Set1'))

	#assign pdf name
	file <- c(".pdf")
	name <- paste(trait, ".pdf", sep='')
	name <- paste(dir, name, sep="/")

	#save the plot as a pdf h/w 6 inches
	pdf(name, height=6,width=6)
	par(mar=c(6,4,0.5,0.5), oma=c(0.1,0.1,0.1,0.1), mgp=c(1.5,0.5,0))
	beeswarm(predictions[,trait] ~ map[,map_column], 
			corral='random',  
					cex.axis=1, 
					pch=16, 
					col=cols, 
					xlab='', 
					ylab='', 
					cex=1, 
					cex.lab=1, 
					las=2)
	bxplot(predictions[,trait] ~ map[,map_column], add=TRUE)
	if(is.null(clr_trans)){
		mtext("Relative Abundance", 2, 3)
	} else {
		mtext("CLR Transformed Relative Abundance", 2, 3)
	}
	dev.off()
	}
}
	
