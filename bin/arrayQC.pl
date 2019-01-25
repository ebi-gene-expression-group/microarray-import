#!/usr/bin/env perl
#
# Script for running microarray quality assessment. This script will read
# MAGE-TAB for an experiment and work out which data files belong to which
# array design, and are associated with which factor value(s). This info is
# written to a temp file which is read by an R script to do the quality
# assessment via the arrayQualityMetrics Bioconductor package.
#
# It will read the Atlas XML config defining assay groups and contrasts and
# only run QC for assays therein. It will remove assays that fail QC, and if
# the assay group(s) they belong to no longer have enough replicates as a
# result, it will remove those assay group(s) and the contrast(s) those assay
# group(s) are used in.
# If any assays fail QC, original XML config file will be renamed to:
# 	<expt accession>-configuration.xml.beforeQC
# New XML config will be written to:
# 	<expt accession>-configuration.xml

use strict;
use warnings;

use Atlas::Magetab4Atlas;
use Atlas::AtlasConfig::Reader qw( parseAtlasConfig );
use File::Spec;
use File::Basename;
use File::Copy;
use File::Path;
use EBI::FGPT::Config qw( $CONFIG );
use Log::Log4perl;
use IPC::Cmd qw( can_run );
use Atlas::Common qw( create_atlas_site_config );

$| = 1;

my $logger_config = q(
	log4perl.rootlogger						= INFO, SCREEN
	log4perl.appender.SCREEN				= Log::Log4perl::Appender::Screen
	log4perl.appender.SCREEN.stderr			= 0
	log4j.PatternLayout.cspec.Q				= sub { return "[QC]" }
	log4perl.appender.SCREEN.layout			= Log::Log4perl::Layout::PatternLayout
	log4perl.appender.SCREEN.layout.ConversionPattern = %-5p - %Q %m%n
);

# Initialise logger.
Log::Log4perl::init( \$logger_config );
my $logger = Log::Log4perl::get_logger;

# Absolute directory path to the file storage 
my $abs_path = dirname(File::Spec->rel2abs(__FILE__));

my $atlasProdDir = $ENV{ "ATLAS_PROD" };
unless( $atlasProdDir ) {
	$logger->logdie( "ATLAS_PROD environment variable is not defined, cannot continue." );
}

my $atlasSiteConfig = create_atlas_site_config;

# Helpful message
my $usage = "Usage:
	arrayQC.pl <experiment accession>
";

# Get accession from cmdline
my $exptAccession = shift;

# If nothing was provided, print message and die.
unless($exptAccession) { die $usage; }

# Make Atlas XML config filename. This script will run from the directory
# containing the XML file.
my $atlasXMLfile = "$exptAccession-configuration.xml";

# Die if the XML file doesn't exist.
unless(-r $atlasXMLfile) {
	$logger->logdie( "Could not read file $atlasXMLfile" );
}

$logger->info( "Reading XML config from \"$atlasXMLfile\"..." );
my $experimentConfig = parseAtlasConfig( $atlasXMLfile );
$logger->info( "Successfully read XML config." );

# R script (should be in PATH).
#-------------------------------------------------- 
my $qcRscript = "$abs_path/arrayQC.R";

unless( can_run( $qcRscript ) ) {
    $logger->logdie( "$qcRscript not found. Please ensure it is installed and you can run it." );
}

# Check that R is available.
unless( can_run( "R" ) ) {
    $logger->logdie( "R not found. Please ensure it is installed and you can run it." );
}

# Path to directory with ArrayExpress/Atlas load directories underneath.
my $exptsLoadStem = File::Spec->catdir( $CONFIG->get_AE2_LOAD_DIR, "EXPERIMENT" );

# miRBase mapped array designs -- we need to subset probes if we find one of these.
# Get an array of miRBase mapping files.
my $miRBaseDirectory = File::Spec->catdir( $atlasProdDir, $atlasSiteConfig->get_mirbase_mappings_directory );
my @miRBaseFiles = glob( "$miRBaseDirectory/*.A-*.tsv" );

# Create a hash for easy checking.
my $miRBaseFileHash = {};
foreach my $miRBaseFile (@miRBaseFiles) {
	# Get the array design from the file name.
	(my $arrayDesign = $miRBaseFile) =~ s/.*(A-\w{4}-\d+)\.tsv/$1/;
	# Add the miRBase mapping file to the hash with the array design as key.
	$miRBaseFileHash->{ $arrayDesign } = $miRBaseFile;
}

# Get the pipeline (e.g. MEXP, MTAB, GEOD, ...) for this experiment.
(my $pipeline = $exptAccession) =~ s/E-(\w{4})-\d+/$1/;

# Experiment load directory and IDF filename.
my $loadDir = "$exptsLoadStem/$pipeline/$exptAccession";
my $idfFilename = "$loadDir/$exptAccession.idf.txt";

# Die if the IDF doesn't exist.
unless(-e $idfFilename) {
	die "[ERROR] Could not find IDF file $idfFilename\n";
}

# Read the MAGE-TAB.
$logger->info( "Reading MAGE-TAB from \"$idfFilename\"..." );
my $magetab4atlas = Atlas::Magetab4Atlas->new( "idf_filename" => $idfFilename );
$logger->info( "Successfully read MAGE-TAB" );

# Next need to sort the raw data files by array design and within that
# by factor value. Use a hash like:
# 	$H->{ <array design 1> }->{ <factor value(s) 1> } = [ <file 1>, <file 2>, <file 3> ]
# 	                        ->{ <factor value(s) 2> } = [ <file 4>, <file 5>, <file 6> ]
# 	  ->{ <array design 2> }->{ <factor value(s) 3> } = [ <file 7>, <file 8>, <file 9> ]
# 	  ...
# Only consider assays that are in the XML config file.
$logger->info( "Collecting factor values and raw data filenames for assays listed in XML config only..." );
my ($arraysToFactorValuesToFiles, $experimentType) = make_arrays_to_factors_to_files($magetab4atlas, $loadDir, $experimentConfig);
$logger->info( "Successfully collected factor values and raw data filenames." );


# For each array design in the experiment: Write the factor value and filename
# information into a temporary text file for each array design. Run the R QC
# script and check output for errors. Remove assays that failed QC from the
# hash representing contrasts and assay groups, as well as assay group(s) and
# contrast(s) that are no longer eligible due to insufficient replication.
# Delete full paths to load directories from HTML report(s) produced by
# arrayQualityMetrics.
# Flag to set if any failed assays are found.
my $failed;
foreach my $arrayDesign ( sort keys %{ $arraysToFactorValuesToFiles }) {
	
	$logger->info( "Running QC in R for array design \"$arrayDesign\"..." );

	# Write annotations (factor value(s), filenames, [labels]) to a temp file. This will be read by R and then deleted.
	my ($tempFile, $miRBaseFile) = write_annotations($arrayDesign, $arraysToFactorValuesToFiles, $miRBaseFileHash, $experimentType);

	# Create directory name for this experiment/array design.
	my $reportDir = $exptAccession."_$arrayDesign"."_QM";

	# Run R script.
	my $qcRscriptOutput = `$qcRscript $tempFile $experimentType $exptAccession $arrayDesign $reportDir $miRBaseFile 2>&1`;

	# Check for errors in the R output.
	if( $? ) {
		# Die here if any errors from R.
		$logger->logdie( "$exptAccession: Error during quality metrics calculation for array $arrayDesign, outout from R below.\n------------\n$qcRscriptOutput\n------------\n" );
	}
    
	# Delete the no longer needed temp file.
	`rm $tempFile`;

	$logger->info( "R process successful." );
	
	# Look for assays that failed QC and remove them from the experiment config.
	$logger->info( "Checking for assays that failed QC..." );
	($experimentConfig, $failed) = remove_rejected_assays($experimentConfig, $qcRscriptOutput, $arrayDesign);
	
	unless( $failed ) {
		$logger->info( "All assays for \"$arrayDesign\" passed QC." );
	}

	# The HTML report contains the full path to the raw data files in the load
	# directory. This looks annoying and is not necessary, so remove it. The
	# path is in index.html and arrayQualityMetrics.js.
	$logger->info( "Removing full paths to data files from QC HTML report..." );
	apply_report_fixes($reportDir, $loadDir);
	$logger->info( "Successfully removed paths from HTML report." );

	$logger->info( "Successfully finished QC for \"$arrayDesign\"" );
}

# If we are still here, rewrite XML config file without failing assays and
# contrasts without enough replicates as a result. Use assay group IDs from
# contrast (assay group pair) IDs to write back to XML.
if( $failed ) {
	
	# Rename the original config file by appending ".beforeQC" to the filename.
	$logger->info( "Renaming original XML config file to \"$atlasXMLfile.beforeQC\"" );
	`mv $atlasXMLfile $atlasXMLfile.beforeQC`;
	
	# Write the new config file.
	$logger->info( "Writing new XML config file without assays that failed QC..." );
	$experimentConfig->write_xml( "." );
	$logger->info( "Successfully written new XML config file." );

	# Because the Atlas::AtlasConfig modules write files with ".auto" on the end,
	# rename the new one so that it doesn't.
	$logger->info( "Removing \".auto\" from new XML config filename..." );
	`mv $atlasXMLfile.auto $atlasXMLfile`;
	$logger->info( "Successully finished all QC processing." );
}
# end
#####


###############
# Subroutines #
###############


# make_arrays_to_factors_to_files
#  - Creates a hash sorting out the data files by array design and then factor value.
#  - E.g.:
# 	$H->{ <array design 1> }->{ <factor value(s) 1> }->{ <assay name 1> } = <file 1>
# 	                        ->{ <factor value(s) 2> }->{ <assay name 2> } = <file 2>
# 	  ->{ <array design 2> }->{ <factor value(s) 3> }->{ <assay name 3> } = <file 3>
# 	  ...
# Arguments:
# 	- $magetab4atlas : a Atlas::Magetab4Atlas object
# 	- $loadDir : path to load directory containing raw data files.
# 	- $experimentConfig : Atlas::AtlasConfig::ExperimentConfig object.
sub make_arrays_to_factors_to_files {
	# Atlas::Magetab4Atlas object and path to load directory.
	my ($magetab4atlas, $loadDir, $experimentConfig) = @_;
	
	# Experiment type from Atlas::Magetab4Atlas will be either "one-colour array" or
	# "two-colour array". Die if it's something else.
	my $experimentType = $magetab4atlas->get_experiment_type;
	unless($experimentType eq "one-colour array" || $experimentType eq "two-colour array") {
		die "This doesn't look like a microarray experiment. Experiment type found is: $experimentType\n";
	}

	# Empty hash to store mappings between factor values, assays and
	# file names.
	my $xmlAssayNames = {};

	# Get the analytics elements from the config.
	my $allAnalytics = $experimentConfig->get_atlas_analytics;

	foreach my $analytics ( @{ $allAnalytics }) {

		my $contrasts = $analytics->get_atlas_contrasts;

		foreach my $contrast ( @{ $contrasts }) {
			
			my $testAssayGroup = $contrast->get_test_assay_group;
			my $referenceAssayGroup = $contrast->get_reference_assay_group;
			
			foreach my $assayGroup ( $testAssayGroup, $referenceAssayGroup ) {
				foreach my $assay ( @{ $assayGroup->get_assays } ) {
					$xmlAssayNames->{ $assay->get_name } = 1;
				}
			}
		}
	}
	
	# Ref to empty hash to fill.
	my $arraysToFactorValuesToFiles = {};
	# Go through the assays...
	foreach my $assay4atlas (@{ $magetab4atlas->get_assays }) {

		# Get the assay name
		my $assayName = $assay4atlas->get_name;
		# Skip this one if it's not in the XML config file
		unless(exists($xmlAssayNames->{ $assayName })) { next; }

		# Get the array design
		my $arrayDesign = $assay4atlas->get_array_design;
		# Get the factor(s) and value(s)
		my $factors = $assay4atlas->get_factors;
		# Get the raw data filename
		my $arrayDataFile = File::Spec->catfile( $loadDir, $assay4atlas->get_array_data_file );

        # Get technology from Array design file using peach API. 
        my $adfInfoUrl = $atlasSiteConfig->get_arrayexpress_adf_info_url;
        my $arrayDataTech=`curl -s $adfInfoUrl$arrayDesign`;

		# For 1-colour array data, need to tell the R script whether this is
		# Affymetrix or Agilent (or other -- for now we only handle Affy and Agil
		# data). Work this out based on the file extension for now -- another
		# possibility might be to parse the ADF, but there is no standard way
		# to record manufacturer there so that might not be ideal.

		if($experimentType eq "one-colour array") {
			if( $arrayDataTech =~ /Affymetrix/i ) { $experimentType = "affy"; }
			elsif( $arrayDataTech =~ /Illumina/i ) { $experimentType = "lumi"; }
    		elsif( $arrayDataTech =~ /Agilent/i ) { $experimentType = "agil1"; }
    		else {
    			$logger->logdie( "Error $arrayDataTech not found for $arrayDesign " )
    		}
		} elsif($experimentType eq "two-colour array") {
			$experimentType = "agil2";
		}

		# Push all factor values onto an array.
		my @factorValues = ();
		foreach my $factor ( sort keys %{ $factors } ) {
			push @factorValues, sort keys %{ $factors->{ $factor } };
		}

		# Stick the factor values together if there's more than one. If there's
		# only one this just returns the factor value by itself.
		my $factorValue = join ", ", @factorValues;
		
		$arraysToFactorValuesToFiles->{ $arrayDesign }->{ $factorValue }->{ $assayName } = $arrayDataFile;
	}
	return ($arraysToFactorValuesToFiles, $experimentType);
}


# write_annotations
# 	- Writes the raw data filenames, factor values and assay names [and labels
# 	for 2-colour] to a temporary file.
# ARGUMENTS:
# 	- $arrayDesign : ArrayExpress array design accession
# 	- $arraysToFactorValuesToFiles : reference to hash of array designs, factor values, assay names and filenames.
# 	- $miRBaseFileHash : reference to hash of array designs that are for microRNA.
# 	- $experimentType : affy, lumi, agil1, or agil2.
sub write_annotations {
	my ($arrayDesign, $arraysToFactorValuesToFiles, $miRBaseFileHash, $experimentType) = @_;

	# Check if the array design has a miRBase mapping file
	# Flag
	my $miRBaseFile = 0;
	if(exists($miRBaseFileHash->{ $arrayDesign })) {
		$miRBaseFile = $miRBaseFileHash->{ $arrayDesign };
	}
	
	# Name for temp file to write annotations to, in /tmp with process ID ($$).
	my $tempFile = "/tmp/$exptAccession"."_$arrayDesign.$$.tsv";
	# Create annotations temp file.
	open(my $tmpFH, '>', $tempFile) or die "Can't create file \"$tempFile\": $!\n";
	
	# Write factor values and corresponding raw data filename.
	# If this is two-colour data we want to include the label info in the file.
	# Create a file with headings like:
	# 	AssayName	Cy3	Cy5	FileName
	if($experimentType eq "agil2") {
		# Ref to empty hash to remember two-colour annotations and filenames.
		my $twoColourAnnotations = {};
		# Go through the factor values for this array design...
		foreach my $factorValue ( sort keys %{ $arraysToFactorValuesToFiles->{ $arrayDesign } }) {
			# Go through the assays for this factor value...
			foreach my $assayName ( sort keys %{ $arraysToFactorValuesToFiles->{ $arrayDesign }->{ $factorValue } }) {
				# Get the label (Cy3 or Cy5)
				(my $label = $assayName) =~ s/.*\.(Cy\d)$/$1/;
				# Make a version of the assay name without the label.
				(my $assayNameNoLabel = $assayName) =~ s/\.Cy\d$//;
				
				# Add the factor value to the hash for this assay for this label.
				$twoColourAnnotations->{ $assayNameNoLabel }->{ $label } = $factorValue;
				# Add the file name for this assay as well.
				$twoColourAnnotations->{ $assayNameNoLabel }->{ "filename" } = $arraysToFactorValuesToFiles->{ $arrayDesign }->{ $factorValue }->{ $assayName };
			}
		}

		# Write header.
		print $tmpFH "AssayName\tCy3\tCy5\tFileName";
		# Write the annotations.
		foreach my $assayNameNoLabel ( sort keys %{ $twoColourAnnotations }) {
			print $tmpFH "\n$assayNameNoLabel\t";
			print $tmpFH $twoColourAnnotations->{ $assayNameNoLabel }->{ "Cy3" }, "\t";
			print $tmpFH $twoColourAnnotations->{ $assayNameNoLabel }->{ "Cy5" }, "\t";
			print $tmpFH $twoColourAnnotations->{ $assayNameNoLabel }->{ "filename" };
		}
	}
	# For 1-colour arrays, don't need the Cy3/Cy5 info.
	else {
		# Write header.
		print $tmpFH "AssayName\tFactorValue\tFileName";
		# Go through the factor values for this array...
		foreach my $factorValue ( sort keys %{ $arraysToFactorValuesToFiles->{ $arrayDesign } }) {
			# Go through the assays for this factor value...
			foreach my $assayName ( sort keys %{ $arraysToFactorValuesToFiles->{ $arrayDesign }->{ $factorValue } }) {
				# Write annotations.
				print $tmpFH "\n$assayName\t$factorValue\t".$arraysToFactorValuesToFiles->{ $arrayDesign }->{ $factorValue }->{ $assayName };
			}
		}
	}
	close $tmpFH;

	return ($tempFile, $miRBaseFile);
}


# apply_report_fixes
# 	- Remove the full path to the load directory from index.html and arrayQualityMetrics.js.
# 	- Replaces the broken "hrwiter" link with the correct one in index.html.
# 	- Replaces the arrayQualityMetrics.css filename with the Atlas CSS filename in index.html.
# Arguments:
# 	- $reportDir : the directory containing the output arrayQualityMetrics.
# 	- $loadDir : the load directory path to remove.
sub apply_report_fixes {

	my ( $reportDir, $loadDir ) = @_;
	
	# We need to fix the HTML report and the JavaScript file.
	my @filesToFix = ( "$reportDir/index.html", "$reportDir/arrayQualityMetrics.js" );

	foreach my $original ( @filesToFix ) {

		# A filename for a temp file to write to, in /tmp with process ID.
		my $temp = "/tmp/$$.temp";
		
		# Open the report for reading.
		open(my $originalFH, "<", $original) or die "Can't open $original : $!\n";
		# Open the temp file for writing.
		open(my $tempFH, ">", $temp) or die "Can't open $temp : $!\n";
		
		# Go through the report line by line...
		while(defined(my $line = <$originalFH>)) {

			# Remove the load directory path if it's there.
			$line =~ s/$loadDir\/*//g;
            
            # Replace the hwriter link with the correct one if needed.
            $line =~ s/http:\/\/www\.embl\.de\/~gpau\/hwriter\/index\.html/https:\/\/cran\.r-project\.org\/web\/packages\/hwriter\/index\.html/g;

            # Replace the CSS filename with the Atlas one if neede.
            $line =~ s/arrayQualityMetrics\.css/\/gxa\/resources\/css\/quality-metrics\.css/g;

            # For the JavaScript file, we also need to fix the loop at line 87
            if( $original =~ /arrayQualityMetrics\.js$/ ) {
                if( $. == 87 ) {
                    $line = "            {\n                i++;\n            $line\n            }\n";
                }
            }
			
            # Write the line to the temp file.
			print $tempFH $line;
		}

		# Close them.
		close($originalFH);
		close($tempFH);

		# Overwrite the original report with the fixed one.
		`mv $temp $original`;
	}
}


# remove_rejected_assays
#	- Find names of assays that failed QC in the R script output and remove
#	them from the experiment config. Also remove contrasts containing them if the
#	contrast no longer has enough replicates without failed assays.
# ARGUMENTS:
# 	- $experimentConfig : Atlas::AtlasConfig::ExperimenConfig object.
# 	- $qcRscriptOutput : variable containing all output from R script (STDOUT, STDERR)
# 	- $arrayDesign : ArrayExpress array design accession
sub remove_rejected_assays {
	my ($experimentConfig, $qcRscriptOutput, $arrayDesign) = @_;
	
	# Flag to set if we see any failed assays (hence XML needs re-write).
	my $failed = 0;

	# Check for rejected assays. If there are some, check them in the contrast hash.
	if($qcRscriptOutput =~ /REJECTED ASSAYS:\t(.*)\n/) {
		# Get the string of rejected assay names separated by tabs if more than one.
		my @rejectedAssays = split "\t", $1;

		# Want to check whether assay groups containing the rejected assay(s)
		# still have enough assays without rejected ones.
		foreach my $rejected (@rejectedAssays) { 
			# Log that this assay failed.
			$logger->info( "$exptAccession: Assay \"$rejected\" failed QC and will be removed from XML config." );
			
			# Set flag
			$failed = 1;
			
			$experimentConfig->remove_assay( $rejected );
		}
	}
	return($experimentConfig, $failed);
}

