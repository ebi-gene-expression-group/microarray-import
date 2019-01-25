#!/bin/bash

# This script decorates a microarray experiment with gene name and identifier from the latest Ensembl (or miRBase - as applicable) release.
# Source script from the same (prod or test) Atlas environment as this script
scriptDir=$(cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
source ${scriptDir}/decorate_routines.sh

if [ $# -lt 1 ]; then
        echo "Usage: $0 PATH_TO_EXPERIMENT "
        echo "e.g. $0 ${ATLAS_PROD}/analysis/differential/microarray/experiments/E-MTAB-1066"
        exit 1;
fi
expPath=$1
e=`basename ${expPath}`

decorate_if_exists() {
    fileToDecorate=$1
    if [ -e "$fileToDecorate" ]; then
        decorate_microarray_file $@
        if [ $? -ne 0 ]; then
                echo "ERROR: FAILED decorate_microarray_file $@" >&2
                exit 1
        fi
    fi
}

find ${expPath} -maxdepth 1 -name "${e}_A-*-analytics.tsv.undecorated" \
    | xargs -n1 basename \
    | sed "s/${e}_//" \
    | sed "s/-analytics.tsv.undecorated//" \
    | while read -r arrayDesign ; do
    arrayDesignFile=`get_arraydesign_file ${arrayDesign}`
    if [ $? -ne 0 ]; then
        echo "ERROR: Could not find array design: $arrayDesign" >&2
        exit 1
    fi
    organism=$(get_organism_given_arraydesign_file ${arrayDesignFile} )
    if [ ! -z `echo $arrayDesignFile | grep mirbase` ]; then
        # This is a miRNA microarray experiment
        geneNameFile="${expPath}/mature.accession.tsv.aux"
        tail -n +2 ${ATLAS_PROD}/bioentity_properties/mirbase/${organism}.mature.tsv | awk -F"\t" '{print $2"\t"$1}' | sort -k 1,1 > $geneNameFile
    else
        geneNameFile=`get_geneNameFile_given_organism ${organism}`
    fi

    decorate_if_exists "${expPath}/${e}_${arrayDesign}-analytics.tsv.undecorated" "$arrayDesignFile" "$geneNameFile"
    decorate_if_exists "${expPath}/${e}_${arrayDesign}-normalized-expressions.tsv.undecorated" "$arrayDesignFile" "$geneNameFile"
    decorate_if_exists "${expPath}/${e}_${arrayDesign}-average-intensities.tsv.undecorated" "$arrayDesignFile" "$geneNameFile"
    decorate_if_exists "${expPath}/${e}_${arrayDesign}-log-fold-changes.tsv.undecorated" "$arrayDesignFile" "$geneNameFile"
    decorate_if_exists "${expPath}/${e}_${arrayDesign}-analytics.tsv.undecorated.unrounded" "$arrayDesignFile" "$geneNameFile"

done
rm -rf ${expPath}/mature.accession.tsv.aux
