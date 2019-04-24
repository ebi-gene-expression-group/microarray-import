#!/usr/bin/env bash

# Source script from the same (prod or test) Atlas environment as this script
scriptDir=$(cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
projectRoot=${scriptDir}/..

[ ! -z ${ATLASPROD_PATH+x} ] || ( echo "Env var ATLASPROD_PATH needs to be defined." && exit 1 )

# Set path (this is done at this level since this will be executed directly):
for mod in analysis exec bash_util bioentity_annotations; do
    export PATH=$ATLASPROD_PATH/$mod:$PATH
done

# source routines
source decorate_routines.sh

expPath=$1
e=$(basename $expPath)

decorate_if_exists() {
    fileToDecorate=$1
    if [ -e "$fileToDecorate" ]; then
        decorate_normalized_file $@
        if [ $? -ne 0 ]; then
                echo "ERROR: FAILED decorate_normalized_microarray_file $@" >&2
                exit 1
        fi
    fi
}

decorate_normalized_file() {
    f=$1
    arrayDesignFile=$2
    geneNameFile=$3

    echo decorate_microarray_file "$@"
    decoratedFile=`echo $f | sed 's/\.undecorated//'`
    amm -s "$projectRoot"/bioentity_annotations/decorateFile.sc \
        --geneIdFile "$arrayDesignFile" \
        --geneNameFile "$geneNameFile" \
        --source "$f" \
        | perl -e 'print scalar <>, sort <>;' \
        > $decoratedFile.decorated.tmp
    decoratedFileLength=$(wc -l "$decoratedFile.decorated.tmp" | cut -f 1 -d ' ' )
    if [ -s "$decoratedFile.decorated.tmp" ] && [ "$decoratedFileLength" -gt 1 ]; then
        return 0
    else
     return 1
    fi
}

## get normalized-expressions.tsv file
find ${expPath} -maxdepth 1 -name "${e}_A-*-normalized-expressions.tsv.undecorated" \
    | xargs -n1 basename \
    | sed "s/${e}_//" \
    | sed "s/-normalized-expressions.tsv.undecorated//" \
    | while read -r arrayDesign ; do
    arrayDesignFile=$(get_arraydesign_file ${arrayDesign})
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
            geneNameFile=$(get_geneNameFile_given_organism $organism)
    fi

    ## generate temp decorated file normalized-expressions.tsv  
    decorate_if_exists "${expPath}/${e}_${arrayDesign}-normalized-expressions.tsv.undecorated" "$arrayDesignFile" "$geneNameFile"       

    # use the temp decorated file to find highest mean of probe ids per gene
    if  [ -e "${expPath}/${e}_${arrayDesign}-normalized-expressions.tsv.decorated.tmp" ]; then
            echo "Merging probe ids with highest mean per gene"
            highestMeanProbeIdsPerGene.R "${e}_${arrayDesign}-normalized-expressions.tsv.decorated.tmp"
    else 
         echo "ERROR: ${expPath}/${e}_${arrayDesign}-normalized-expressions.tsv.decorated.tmp doesn't exist" 
         exit 1       
    fi
done    

rm -rf ${expPath}/*-normalized-expressions.tsv.decorated.tmp 
