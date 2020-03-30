# Plots predictions as scatter plot with linear model as a line, 
#   performs a correlation test (Pearson) and reports linear model stats.
# Required: prediction table, mapping file, map column, output filepath
# Returns: pdfs of plots and text files of statistics

"plot.predictions.continuous" <- function(predictions, 
											map, 
											map_column,
											clr_trans, 
											output_fp=NULL){
	#set directory
	dir <- paste(output, "predicted_phenotypes", sep="/")

	#cat("\n\n RESULTS: \n\n")
	#check dimensions of traits
	#cat("\nDimensions of trait table:\n")
	#cat(dim(predictions),"\n")
	#cat("Dimensions of the mapping file:\n")
	#cat(dim(map),"\n")

	#make sure map column is numeric
	map[,map_column] <- map[,map_column,drop=F]
	
	#ensure same order of samples in map and traits
	map <- map[rownames(predictions),,drop=F]
	
	map[,map_column] <- as.numeric(as.character(map[,map_column]))
	
	#define traits to test
	traits<- colnames(predictions)

	for(x in 1:length(traits)){
		trait <- traits[x]

		#correlation test and linear model
		cor.tests <- cor.test(predictions[,trait], map[,map_column], method="spearman")
		lm.out <- lm(predictions[,trait] ~ map[,map_column])
	
		coeff <- summary(lm.out)$coefficients[2,]
		fstats <- summary(lm.out)$fstatistic
		p <- pf(fstats[1],fstats[2],fstats[3],lower.tail=F)
	
		header <- paste("\n\n", trait, sep='')
		# cat(header)
		# cat("\n\nLiner Model Statistics for: ")
		# cat(trait)
		# cat("\n\nCoefficients\n")
		# cat("Estimate    Std. Error    t value    Pr(>|t|)\n")
		# cat(coeff)
		# cat("\np-value\n")
		# cat(p)
		# cat("\n\nPearson's Correlation\n\n")
		# cat("Correlation Estimate:\n")
		# cat(cor.tests$estimate)
		# cat("\np-value\n")
		# cat(cor.tests$p.value)
		# cat("\n\n")
	
		outfile <- paste(trait, "_stats.txt", sep="")
	
		sink(paste(dir,outfile,sep='/'))
		cat("\n\nLiner Model Statistics for: ")
		cat(trait)
		cat("\n\nCoefficients\n")
		cat("Estimate    Std. Error    t value    Pr(>|t|)\n")
		cat(coeff)
		cat("\np-value\n")
		cat(p)
		cat("\n\nSpearman's Correlation\n\n")
		cat("Correlation Estimate:\n")
		cat(cor.tests$estimate)
		cat("\np-value\n")
		cat(cor.tests$p.value)
		cat("\n\n")
		sink()
	
		#set color color palette from RColorBrewer, 95% transparency
		cols <- sprintf('%s97',brewer.pal(9,'Set1'))
	
		#to use gradient colors
		#Pal <- colorRampPalette(c(cols[1],cols[2]))
		#map$Col <- Pal(5)[as.numeric(cut(map[,map_column],breaks = 5))]
	
		#assign pdf name
		file <- c(".pdf")
		name <- paste(trait, ".pdf", sep='')
		name <- paste(dir, name, sep="/")
		
		#now save the plot as a pdf h/w 6 inches
		pdf(name, height=6,width=6);
		par(mar=c(6,4,0.5,0.5), oma=c(0.1,0.1,0.1,0.1), mgp=c(1.5,0.5,0))
		plot(predictions[,trait] ~ map[,map_column], 
			cex.axis=1, 
			pch=16, 
			col=cols[1], 
			xlab='', 
			ylab='', 
			cex=1, 
			cex.lab=1)
		abline(lm.out, lty=2)
		if(is.null(clr_trans)){
			mtext("Relative Abundance", 2, 3)
		} else {
			mtext("CLR Transformed Relative Abundance", 2, 3)
		}
		mtext(map_column, 1, 3)
		dev.off()
	}
}


