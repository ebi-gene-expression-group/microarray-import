#!/usr/bin/env Rscript

args <- commandArgs( TRUE )

analyticsFilename <- args[ 1 ]

analyticsResults <- read.delim( analyticsFilename, header=TRUE )

pvalueColnames <- colnames( analyticsResults )[ grep( "p.value", colnames( analyticsResults ) ) ]

nonNAcounts <- sapply( pvalueColnames, function( x ) {
    
    length( which( !is.na( analyticsResults[[ x ]] ) ) )
} )

if( !any( nonNAcounts > 0 ) ) {
    stop( paste( "ERROR - No non-NA p-values found in", analyticsFilename ) )
}
