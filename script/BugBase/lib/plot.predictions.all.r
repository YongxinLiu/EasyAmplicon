# Plots predictions as beeswarm plots for all samples (no groups).
# Required: prediction table, output filepath
# Returns: pdfs of plots and text files of summary statistics 
#   (no comparisons made).


"plot.predictions.all" <- function(predictions, 
									clr_trans,
									output_fp=NULL){

	#set directory
	dir <- paste(output, "predicted_phenotypes", sep="/")
	
	traits<- colnames(predictions)

	for(x in 1:length(traits)){
		trait <- traits[x]

		# #Print to screen
		# header <- paste("\n\n", trait, sep='')
		# cat(header)
		# cat("\nRelative Abundance with trait (mean):\n")
		# cat(mean(predictions[,trait]))
		# cat("\nRelative Abundance with trait (median):\n")
		# cat(median(predictions[,trait]))
		# cat("\nStandard deviation:\n")
		# cat(sd(predictions[,trait]),"\n")

		#write summary stats to output file
		outfile <- paste(trait, "_stats.txt", sep="")   
		sink(paste(dir,outfile,sep='/'))
		cat("\nRelative Abundance with trait (mean):\n")
		cat(mean(predictions[,trait]))
		cat("\nRelative Abundance with trait (median):\n")
		cat(median(predictions[,trait]))
		cat("\nStandard deviation:\n")
		cat(sd(predictions[,trait]),"\n")
		sink()

		#get palette 1 from R ColorBrewer
		cols <- sprintf('%s95',brewer.pal(9,'Set1'))

		#assign pdf name
		file <- c(".pdf")
		name <- paste(trait, ".pdf", sep='')
		name <- paste(dir, name, sep="/")

		#save the plot as a pdf h/w 6 inches
		pdf(name, height=6,width=6);
		par(mar=c(6,4,0.5,0.5), oma=c(0.1,0.1,0.1,0.1), mgp=c(1.5,0.5,0))
		beeswarm(predictions[,trait], 
					corral='random',  
					cex.axis=1, 
					pch=16, 
					col=cols, 
					xlab='', 
					ylab='', 
					cex=1, 
					cex.lab=1, 
					las=2)
		bxplot(predictions[,trait], add=TRUE)
		if(is.null(clr_trans)){
			mtext("Relative Abundance", 2, 3)
		} else {
			mtext("CLR Transformed Relative Abundance", 2, 3)
		}
		dev.off()
	}
	
}

