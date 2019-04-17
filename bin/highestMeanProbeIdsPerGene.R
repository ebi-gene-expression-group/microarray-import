#!/usr/bin/env Rscript

suppressPackageStartupMessages(require(data.table))

# Get commandline arguments.
args <- commandArgs( TRUE )

if( length( args ) != 1 ) {
  stop( "\nUsage:\n\tHishestmeanProbeIdsPerGene.R <expNorm_decorated_filename>\n\n" )
}

# Get the analytics filename from the arguments.
expNormFileName <- args[ 1 ]

if( !grepl( "-normalized-expressions.tsv.decorated.tmp$", expNormFileName ) ) {
  stop( paste(   expNormFileName, "does not look like a decorated expNorm filename. Please check." ) )
}

findHighestMeanProbePerGene<-function( expNormFileName ) {
   
  fread(input=expNormFileName)->exprWAnnot
  
  exprWAnnot$mean<-rowMeans(exprWAnnot[,4:ncol(exprWAnnot)])
  
  # sort by gene id and mean, ascending
  setkeyv(exprWAnnot, c("Gene ID","mean"))
  
  # retrieve the last element for each gene (highest mean), leave only probes
  exprWAnnot[,.SD[.N],by=`Gene ID`][,DesignElementAccession]->highestMeanProbePerGene
  
  # expTTable have unique probe ids to gene name mapping with highest mean probe per gene
  exprWAnnot[DesignElementAccession %in% highestMeanProbePerGene, ] -> expTTable
  
  # remove the  Gene ID/Name and mean colum to make it undecorated 
  expTTable[,-c("Gene ID", "Gene Name", "mean")] -> expTTable
  
  ## export undecorated expNormFileName 
  write.table( expTTable, file = gsub("decorated.tmp","undecorated",expNormFileName), row.names=FALSE, quote=FALSE, sep="\t" )
}

## merge prode ids with highest mean
findHighestMeanProbePerGene( expNormFileName )

  
