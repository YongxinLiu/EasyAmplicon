# Plots the average trait relative abundance for each treatment from thresholds 
#   0 to 1.
# Inputs: prediction array, mapping file, map column, output filepath
# Returns: pdfs of plots

"plot.thresholds" <- function(prediction_array, 
								map, 
								map_column, 
								clr_trans,
								output_fp=NULL){

	#set directory
	dir <- paste(output, "thresholds", sep="/")

	#ensure same order of samples in map and traits
	map <- map[rownames(prediction_array),,drop=F]
	
	#define groups
	groups <- sort(unique(map[,map_column]))

	#get palette 1 from R ColorBrewer
	cols <- brewer.pal(9,'Set1')

	#define traits
	traits<- colnames(prediction_array)

	#define thresholds
	thresholds <- dimnames(prediction_array)[3][[1]]

	for(x in 1:length(traits)){
		trait <- traits[x]

		#assign pdf name
		file <- c(".pdf")
		name <- paste(trait, ".pdf", sep='')
		name <- paste(dir, name, sep="/")

		pdf(name, height=6,width=6);
		par(mar=c(6,4,0.5,6), oma=c(0.1,0.1,0.1,0.1), mgp=c(1.5,0.5,0))
				
		# Plot the mean relative abundance for each threshold
		if(is.null(clr_trans)){
			yvars <- as.numeric(thresholds)
			plot(thresholds, yvars, 
				type="n", 
				cex.axis=1, 
				pch=16, 
				xlab='', 
				ylab='', 
				cex=1)
			yrange <- abs(max(yvars)-min(yvars))
		} else {
			yvars <- seq(min(prediction_array[,trait,]), max(prediction_array[,trait,]), by = abs((max(prediction_array[,trait,]) - min(prediction_array[,trait,]))/(length(thresholds)-1)))
			plot(thresholds, yvars, 
				type="n", 
				cex.axis=1, 
				pch=16, 
				xlab='', 
				ylab='', 
				cex=1)
			yrange <- abs(max(yvars)-min(yvars))
		}
		for(i in 1:length(groups)){
			group <- groups[i]
			means <- c()
			ix <- map[,map_column] == group
			for(n in 1:length(thresholds)){
				means <- c(means,mean(prediction_array[ix,trait,n]))
			}
			lines(thresholds, means, col=cols[i], lwd=2)
			legend(x=1.05, y=0.8+(yrange*0.05*i), group, cex=0.80, lty=1, lwd=2, 
					col=cols[i], bty="n", xpd=TRUE)
			if(is.null(clr_trans)){
				mtext("Relative Abundance", 2, 3)
			} else {
				mtext("CLR Transformed Relative Abundance", 2, 3)
			}
			mtext("Threshold (% of category covered)", 1, 3)    
		}
		dev.off()
	}
}
	
