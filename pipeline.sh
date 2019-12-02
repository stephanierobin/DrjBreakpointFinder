#! /bin/sh
red=`tput setaf 1`
green=`tput setaf 2`
yellow=`tput setaf 3`
cyan=`tput setaf 6`
bold=`tput bold`
reset=`tput sgr0`

function help {
    echo $reset
echo "====================="
echo "DrjBreakpointFinder"
echo "====================="
echo "Usage: sh pipeline.sh -r reads.fasta -g genome.fasta -i input_directory -o output_directory"
}

while getopts "r:g:i:o:h" opt; do
    case $opt in
        r)
        reads=$OPTARG
        ;;

        g)
        genome=$OPTARG
        ;;

        i)
        input=$OPTARG
        ;;

        o)
        output=$OPTARG
        ;;

        h)
        help
        exit
        ;;
    esac
done

if [[ -z "${reads}" ]]; then
    echo "${red}-r is mandatory$reset" >&2
    exit
fi
if [[ -z "${genome}" ]]; then
    echo "${red}-g is mandatory$reset" >&2
    exit
fi
if [[ -z "${input}" ]]; then
    echo "${red}-i is mandatory$reset" >&2
    exit
fi
if [[ -z "${output}" ]]; then
    echo "${red}-o is mandatory$reset" >&2
    exit
fi


## Genocluster :
source /local/env/envclustalw.sh
source /local/env/envncbi.sh
source /local/env/envr-3.5.1.sh
source /local/env/envemboss.sh

## Intermediate files or directories
blastDir=$output/blast
blastFile=$blastDir/megablast_result.m8
readLength=$input/$(basename $reads .fasta)"_rlength.tab"
breakpointDir=$output/breakpoint
drjPairs=$breakpointDir/drjPairs.tab
seqCoordPairs=$breakpointDir/seqCoordPairs.tab
vectorDir=$breakpointDir/vectors
alignDir=$breakpointDir/align
figureDir=$breakpointDir/figures
segmentationResult=$breakpointDir/segmentation.tab
multipleLog=$breakpointDir/list-multiple-exaequo.txt
pairDir=$output/drjPairs_all_segments
pairDirFinal=$output/drjPairs_merged_segments
confirmedDrjPairs=$output/drjPairs_confirmed.tab
summaryDir=$output/drjPairs_figures
alignPairDir=$output/drjPairs_alignments
vectorPairDir=$alignPairDir/vectors

mkdir $output
mkdir -p $blastDir

formatdb -i $genome -p F -o T
megablast -i $reads -d $genome -o $blastFile -m8 > $blastDir/megablast.log 2>&1
mv formatdb.log $blastDir/

echo "Megablast finished"

perl get_read_length.pl  -reads $reads -out $readLength

mkdir -p $breakpointDir


R --slave --vanilla --quiet --no-save  <<MyRScript1

source("get_fully_overlapping_triplets.R")
getFullyOverlappingTriplets(blastFile="$blastFile",readLengthFile="$readLength",coordFile="$seqCoordPairs",drjPairsFile="$drjPairs")

MyRScript1

 echo "Triplets selection finished"

 perl align_and_get_mismatch_vectors.pl -coord $seqCoordPairs -bacs $genome -reads $reads -dirFasta $breakpointDir/tmp -dirAlign $alignDir -dirVector $vectorDir

 echo "All triplet Alignments finished"


 mkdir -p $pairDir
 mkdir -p $pairDirFinal

R --slave --vanilla --quiet --no-save  <<MyRScript2

source("get_segmentation.R")
segmentASetOfSequences("$seqCoordPairs","$vectorDir","$figureDir","$segmentationResult",zoom=T,margin=20)

source("get_segmentations_by_DRJpairs.R")
getAllValidSegmentationsByDRJPair(segmentationFile="$segmentationResult",drjPairsFile="$drjPairs",pairDir="$pairDir",confirmedDrjPairFile="$confirmedDrjPairs",multipleLog="$multipleLog")

source("merge_segments_on_drj_pairs.R")
mergeAllDrjs(inputPairDir="$pairDir",outputPairDir="$pairDirFinal")

MyRScript2

echo "All segmentations finished"

mkdir -p $alignPairDir
perl align_drjPairs.pl -coord $drjPairs -bacs $genome -dirFasta $alignPairDir/tmp -dirAlign $alignPairDir/align -dirVector $vectorPairDir

echo "Alignments of DRJ pairs (similarity) finished"

mkdir -p $summaryDir

R --slave --vanilla --quiet --no-save  <<MyRScript4

source("draw_summary_scaffolds.R")
drawAllSummary(pairDir="$pairDir",outputDir="$summaryDir",readNames=FALSE,addSimi="$vectorPairDir",merge=TRUE)

MyRScript4

echo "Summary plots finished"
