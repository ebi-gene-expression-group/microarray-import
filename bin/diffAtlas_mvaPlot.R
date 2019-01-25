#!/usr/bin/env Rscript

# diffAtlas_mvaPlot()
#  - Create an "MvA" plot (average intensity against log2 fold change), and write it to a file.
# ARGUMENTS:
# 	- plotDataFile <- filename of differential expression stats to be used for plotting
# 	- contrastNAme <- text description of the contrast
# 	- plotFile <- filename to write the plot to
# 	- techType <- "microarray" or "rnaseq"
diffAtlas_mvaPlot <<- function(plotDataFile, contrastName, plotFile, techType) {

	# Read data
	madf <- read.delim(plotDataFile, stringsAsFactors=FALSE)

	# First filter out rows with any "NA"s in, can't use them.
	# complete.cases() function returns TRUE for rows where there is no missing
	# data (e.g. "NA"s).
	madf <- madf[complete.cases(madf),]

	if(techType == "microarray") {

		# Find minimum and maximum values for doing the ticks on the axes. Not
		# using for DESeq results for RNA-seq data because for that we'll use
		# log10 scale on x-axis.
		minLogFC <- floor(min(madf$logFC))
		maxLogFC <- ceiling(max(madf$logFC))
		maxAvgExpr <- ceiling(max(madf$avgExpr))

		# label for x-axis
		xAxisLabel = "average intensity"
	
	} else if(techType == "rnaseq") {
		
		# x-axis label different for RNA-seq data
		xAxisLabel = "average normalized count"
	}

	# load ggplot2
	library(ggplot2)

	# Plot!
	#  - Start process with ggplot() function, give it the data for x and y, tell
	#    it how to colour the points (using the "adjPval" column in the data
	#    frame). "factor(adjPval < 0.05)" gives FALSE or TRUE. The levels
	#    should be in the order: FALSE, TRUE. So later specify colour values
	#    with "scale_colour_manual" in that order.
	#  - Use geom_jitter() to reduce "overplotting" so that we can see points
	#    that overlap easier. pass it alpha value for transparency of points
	#    (lower=more transparent) and size for point size.
	mvaPlot <- ggplot(madf, aes(x=avgExpr, y=logFC, colour=factor(adjPval < 0.05))) + geom_jitter(alpha=0.8, size=1) +
		
		# Apply the colours to the points in the plot. These should be in order
		# "non-DE-colour, DE-colour". What's in the "name" argument ends up as
		# the title of the legend/key.
		scale_colour_manual(name="DE call (FDR < 0.05):", values=c("dimgrey", "red"), labels=c("non-DE", "DE")) +
		
		# Remove the default grey grid background to the plot and legend; add
		# back the axis lines.
		theme(panel.grid=element_blank(), panel.background=element_blank(), legend.key=element_blank(), axis.line=element_line()) +
		
		# Position the legend
		theme(legend.position="bottom", legend.direction="horizontal") 

	if(techType == "microarray") {
			
		# Add numbers to the axes, with minima and maxima from above.
		mvaPlot <- mvaPlot + scale_y_continuous(breaks=minLogFC:maxLogFC) + scale_x_continuous(breaks=1:maxAvgExpr)
	
	} else if(techType == "rnaseq") {
		
		# use log scale for DESeq results from RNA-seq data
		mvaPlot <- mvaPlot + scale_x_log10()
	} 
	
	# Label the axes. expression() function lets you make the "2" subscript (highly important! :)
	mvaPlot <- mvaPlot + xlab(xAxisLabel) + ylab(expression(log[2](fold~change))) +

		# Add the title (the description for this contrast). Use strwrap()
		# function to make it wrap after 70 characters.
		ggtitle(paste(strwrap(contrastName, width=70), collapse="\n"))
	
	# Save the plot to a file (type specified by filename extension). Pass
	# type="cairo" for the graphics library as the default doesn't seem to work
	# on lime. The trouble with this is that for some reason it creates a file
	# called "Rplots.pdf". Also, it doesn't make as nice an image as on my
	# MacBook (lower res).
	ggsave(mvaPlot, file=plotFile, dpi=150, type="cairo")
}


# Run with arguments if there are any, otherwise don't do anything.
args <- commandArgs(TRUE)
if(length(args) > 0) {
	do.call(diffAtlas_mvaPlot, as.list(args))
}
		

