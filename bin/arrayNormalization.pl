#!/usr/bin/env perl
#
# This script parsed MAGE-TAB for a given experiment accession and runs
# microarray normalization by calling an R script.
#
# NB: this script is very similar to arrayQC.pl so maybe (a) they should be
# combined, or (b) we should pull out common parts to a module they both use.

use strict;
use warnings;
use 5.10.0;

# XML config parsing.
use Atlas::AtlasConfig::Reader qw( parseAtlasConfig );

# MAGE-TAB parsing.
use Atlas::Magetab4Atlas;

use File::Spec;
use File::Basename;
use File::Copy;
use File::Path;
use IPC::Cmd qw( can_run );

use Log::Log4perl;

# Log4perl config.
my $logger_config = q(
log4perl.rootlogger			         = INFO, SCREEN
log4perl.appender.SCREEN             = Log::Log4perl::Appender::Screen
log4perl.appender.SCREEN.stderr      = 0
log4perl.appender.SCREEN.layout      = Log::Log4perl::Layout::PatternLayout
log4perl.appender.SCREEN.layout.ConversionPattern = %-5p - %m%n
);

# Absolute directory path to the file storage 
my $abs_path = dirname(File::Spec->rel2abs(__FILE__));

# Initialise logger.
Log::Log4perl::init(\$logger_config);
my $logger = Log::Log4perl::get_logger;


# Experiment directory idf filename and ArrayExpress and miRBase Load directory as args
my ($atlasExperimentDir, $idfFilename, $loadDir, $miRBaseDirectory) = @ARGV;
my $exptAccession = (split '\/', $atlasExperimentDir)[-1];

unless( $atlasExperimentDir ) {
	$logger->logdie( "Please provide experiment accession directory path as an argument." );
}
unless( $idfFilename ) {
	$logger->logdie( "Please provide idfFilename as an argument." );
}
unless( $loadDir ) {
	$logger->logdie( "Please provide AE loadDir as an argument." );
}
unless( $miRBaseDirectory ) {
	$logger->logdie( "Please provide miRBase directory as an argument." );
}

# Filename of R script for normalization.
my $normalizationRscript = "$abs_path/arrayNormalization.R";
unless( can_run( $normalizationRscript ) ) {
	$logger->logdie( "Script \"$normalizationRscript\" not found. Please ensure it is in your \$PATH and you can run it.");
}

# Check that we can run R.
unless( can_run( "R" ) ) {
	$logger->logdie( "R was not found. Please ensure it is installed and you can run it." );
}

# miRBase mapped array designs -- we need to subset probes if we find one of these.
# Get an array of miRBase mapping files.
my @A_miRBaseFiles = glob( "$miRBaseDirectory/*.A-*.tsv" );

# Create a hash for easy checking.
my $H_miRBaseFileHash = {};
foreach my $miRBaseFile (@A_miRBaseFiles) {
	# Get the array design from the file name.
	(my $arrayDesign = $miRBaseFile) =~ s/.*(A-\w{4}-\d+)\.tsv/$1/;

	# Add the miRBase mapping file to the hash with the array design as key.
	$H_miRBaseFileHash->{ $arrayDesign } = $miRBaseFile;
}

# Atlas XML config file name.
my $atlasXMLconfigFile = $exptAccession . "-configuration.xml";
# Full path to XML config file.
my $atlasXMLconfigPath = File::Spec->catfile( $atlasExperimentDir, $atlasXMLconfigFile );

# Parse the config.
my $experimentConfig = parseAtlasConfig( $atlasXMLconfigPath );

# Check that the experiment type is a microarray one.
unless( $experimentConfig->get_atlas_experiment_type =~ /array/ ) {
	$logger->logdie( "This does not look like a microarray experiment. Experiment type is \"", $experimentConfig->get_atlas_experiment_type );
}

# Read the MAGE-TAB.
$logger->info( "Reading MAGE-TAB..." );
my $magetab4atlas = Atlas::Magetab4Atlas->new( "idf_filename" => $idfFilename );
$logger->info( "Read MAGE-TAB." );

# Create hash mapping assay names to raw data file names for each array design:
#
# 	$arraysToAssaysToFiles->{ <array design 1> }->{ <assay 1> } = <file 1>
# 												->{ <assay 2> } = <file 2>
#						  ->{ <array design 2> }->{ <assay 3> } = <file 3>
# 		  				  ...
# Also pass $experimentType and get back normalization mode to pass to R script.
my ($H_arraysToAssaysToFiles, $normalizationMode) = &makeArraysToAssaysToFiles( $magetab4atlas, $loadDir, $experimentConfig);

# Log how many array designs and assays were found
my $arrayCount = keys %{ $H_arraysToAssaysToFiles };
$logger->info( "Found $arrayCount array designs" );
foreach my $arrayDesign (keys %{ $H_arraysToAssaysToFiles }) {
	my $assayCount = keys %{ $H_arraysToAssaysToFiles->{ $arrayDesign }};
	$logger->info( "\t$arrayDesign: $assayCount assays." );
}
# Also what kind of normalization will be done. This is either:
# 	- "oligo" : using oligo package for Affymetrix arrays.
# 	- "agil1" : using limma pacakge for Agilent 1-colour data.
# 	- "agil2" : using limma package for Agilent 2-colour data.
$logger->info( "The normalization mode is $normalizationMode." );

# Run normalization for each array design.
foreach my $arrayDesign (keys %{ $H_arraysToAssaysToFiles }) {
	# Check if the array design has a miRBase mapping file
	# Flag
	my $miRBaseFile = 0;
	if(exists($H_miRBaseFileHash->{ $arrayDesign })) {
		$logger->info( "$arrayDesign is a microRNA array design." );
		$miRBaseFile = $H_miRBaseFileHash->{ $arrayDesign };
	}

	# Write a file to read into R
	my $tempFile = File::Spec->catfile( $ENV{ "HOME" }, "tmp", "$exptAccession"."_$arrayDesign.$$.tsv" );
	open(my $tmpFH, ">", $tempFile) or $logger->logdie( "Can't create file \"$tempFile\": $!" );

	# Write headers
	print $tmpFH "AssayName\tFilename";
	foreach my $assayName (keys %{ $H_arraysToAssaysToFiles->{ $arrayDesign }}) {
		print $tmpFH "\n$assayName\t$H_arraysToAssaysToFiles->{ $arrayDesign }->{ $assayName }";
	}
	close $tmpFH;

	# Name for normalized data file.
	my $normalizedDataFile = $exptAccession."_".$arrayDesign."-normalized-expressions.tsv.undecorated";

	$logger->info( "Running normalization in R for $exptAccession, array design $arrayDesign..." );

	# Run R script to do normalization with Bioconductor packages in R.
	# NB: Using 2>&1 means that nothing from R is printed to STDOUT. If there's an
	# error, then it's printed by the part following this line.
	my $RscriptOutput = `$normalizationRscript $tempFile $normalizationMode $normalizedDataFile $miRBaseFile 2>&1`;

	# If there was an error in R die and print the error.
	if($RscriptOutput =~ /error/i) {
		$logger->logdie( "Error encountered during normalization of $exptAccession on array $arrayDesign. Full output from R is below:
			------------------------
			$RscriptOutput
			" );
	}
	else {
		# For 2-colour data, rename files created.
		if($normalizationMode eq "agil2") {
			my $logFCfile = $exptAccession."_$arrayDesign-log-fold-changes.tsv.undecorated";
			my $aValuesFile = $normalizedDataFile.".A-values";
			my $avgIntensitiesFile = $exptAccession."_$arrayDesign-average-intensities.tsv.undecorated";

			`mv $normalizedDataFile $logFCfile`;
			`mv $aValuesFile $avgIntensitiesFile`;
		}

		$logger->info( "Normalization for array design $arrayDesign completed." );
	}

	# Delete temporary file.
	`rm $tempFile`;
}


# Subroutines

# &makeArraysToAssaysToFiles
#
# 	- Creates a hash matching each data file to each assay, for each array design.
# 	- E.g.:
# 		$H->{ <array design 1> }->{ <assay 1> } = <file 1>
# 								->{ <assay 2> } = <file 2>
# 		  ->{ <array design 2> }->{ <assay 3> } = <file 3>
# 		  ...
# Arguments:
# 	- $magetab4atlas : a Atlas::Magetab4Atlas object.
# 	- $loadDir : path to load directory containing raw data files.
# 	- $experimentConfig : Atlas::AtlasConfig::ExperimentConfig object.
sub makeArraysToAssaysToFiles {
	# Atlas::Magetab4Atlas object and path to load directory.
	my ($magetab4atlas, $loadDir, $experimentConfig ) = @_;

	my $experimentType = $experimentConfig->get_atlas_experiment_type;

	# Get all the assay names from the experiment config.
	my $allAnalytics = $experimentConfig->get_atlas_analytics;

	# Create a hash of all the assay names found in the experiment config,
	# mapped to their respective array designs.
	my $arrayDesignsToAssayNames = {};

	foreach my $analytics ( @{ $allAnalytics }) {

		my $analyticsAssays = $analytics->get_assays;

		my $arrayDesign = $analytics->get_platform;

		$arrayDesignsToAssayNames->{ $arrayDesign } = [];

		foreach my $assay ( @{ $analyticsAssays } ) {

			push @{ $arrayDesignsToAssayNames->{ $arrayDesign } }, $assay->get_name;
		}
	}

	# Ref to empty hash to fill
	my $H_arraysToAssaysToFiles = {};

	# Normalization mode.
	my $normalizationMode = 0;

	# Go through the assays...
	foreach my $assay4atlas (@{ $magetab4atlas->get_assays }) {

		# Get assay name
		my $assayName = $assay4atlas->get_name;

		# Escape any metacharacters
		my $assayNameEsc = quotemeta( $assayName );

		# Array design
		my $arrayDesign = $assay4atlas->get_array_design;

		# Check that we saw this assay in the XML config. If not, skip it.
		unless( grep { /^$assayNameEsc$/ } @{ $arrayDesignsToAssayNames->{ $arrayDesign } } ) {
			$logger->info( "Assay \"$assayName\" not found in XML config, not including in normalization." );
			next;
		}

		# Get technology from Array design file using peach API. 
		my $adfInfoUrl = "http://peach.ebi.ac.uk:8480/api/array.txt?acc=";
		my $arrayDataTech=`curl -s $adfInfoUrl$arrayDesign`;

		# Raw data filename.
		my $arrayDataFile = File::Spec->catfile( $loadDir, $assay4atlas->get_array_data_file );

		# For 1-colour array data, need to tell the R script whether this is
		# Affymetrix or Illumina or Agilent (or other -- for now we handle Affy, Agil and Illumina
		# data). Worked this out based on the Array design info stored in ADF using peach API call.
		if( $experimentType =~ /1colour/ ) {
			if( $arrayDataTech =~ /Affymetrix/i ) { $normalizationMode = "affy"; }
			elsif( $arrayDataTech =~ /Illumina/i ) { $normalizationMode = "lumi"; }
			elsif( $arrayDataTech =~ /Agilent/i ) { $normalizationMode = "agil1"; }
			else {
				$logger->logdie( "Error $arrayDataTech not found for $arrayDesign " )
			}
			} elsif( $experimentType =~ /2colour/ ) {
				$normalizationMode = "agil2";

			# Remove label name from assay name which was added by Atlas::Magetab4Atlas.
			$assayName =~ s/\.Cy\d$//;
		}

		# Add data to hash.
		$H_arraysToAssaysToFiles->{ $arrayDesign }->{ $assayName } = $arrayDataFile;
	}

	return ($H_arraysToAssaysToFiles, $normalizationMode);
}
