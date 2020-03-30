# Plots the average trait relative abundance for all samples from thresholds 
#   0 to 1.
# Inputs: prediction array, output filepath
# Returns: pdfs of plots

"plot.thresholds.all" <- function(prediction_array, 
									clr_trans, 
									output_fp=NULL){

	#set directory
	dir <- paste(output, "thresholds", sep="/")

	#get palette 1 from R ColorBrewer
	cols <- brewer.pal(9,'Set1')

	#define traits
	traits<- colnames(prediction_array)

	#define thresholds
	thresholds <- dimnames(prediction_array)[3][[1]]

	for(x in 1:length(traits)){
		trait <- traits[x]
		
		#store all the mean trait abundances        
		means <- c()
		for(i in 1:length(thresholds)){
			means <- c(means,mean(prediction_array[,trait,i]))
		}

		#assign pdf name
		file <- c(".pdf")
		name <- paste(trait, ".pdf", sep='')
		name <- paste(dir, name, sep="/")

		#set up plot parameters
		pdf(name, height=6,width=6);
		par(mar=c(6,4,0.5,6), oma=c(0.1,0.1,0.1,0.1), mgp=c(1.5,0.5,0))

		# Plot the mean relative abundance at each threshold

		plot(thresholds, means, 
				type="n", 
				cex.axis=1, 
				pch=16, 
				xlab='', 
				ylab='', 
				cex=1)
		lines(thresholds, means, col=cols[1], lwd=2)
		if(is.null(clr_trans)){
			mtext("Relative Abundance", 2, 3)
		} else {
			mtext("CLR Transformed Relative Abundance", 2, 3)
		}
		mtext("Threshold (% of category covered)", 1, 3)    
		dev.off()
	}
}
	
