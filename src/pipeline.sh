#! /bin/sh
red=`tput setaf 1`
green=`tput setaf 2`
yellow=`tput setaf 3`
cyan=`tput setaf 6`
bold=`tput bold`
reset=`tput sgr0`

function help {
    echo $reset
echo "====================================================="
echo "DrjBreakpointFinder"
echo "====================================================="
echo "Usage: ./pipeline.sh -r reads.fasta -g genome.fasta -o output_directory"
}

while getopts "r:g:o:h" opt; do
    case $opt in
        r)
        reads=$OPTARG
        ;;

        g)
        genome=$OPTARG
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
if [[ -z "${output}" ]]; then
    echo "${red}-o is mandatory$reset" >&2
    exit
fi
## Programmes
#sources="../programs"

## ma machine :
#alias clustalw='/Applications/Bioinfo/clustalw2' # ok
# blast
#PATH=$PATH:$HOME/Bin/blast-2.2.26/bin/


## Genocluster :
source /local/env/envclustalw.sh
source /local/env/envncbi.sh
#source /local/env/envR-2.14.2.sh
source /local/env/envr-3.5.1.sh
source /local/env/envemboss.sh

## Input data
#inputDir="../data/formatted-data"
#bacs="final_bacs_shortname.fasta"
#reads="HdIV_SANGER_clean_min100.fasta"
#resdir="."
resdir=$(pwd)
#resdir=".."
#genome="Hd_genome_1.0"
#scaffold_name=`basename $1 .m8`
#SCAFFOLD=$1
#scaffold_name="${SCAFFOLD%.*}"

## Intermediate files or directories
blastDir=$resdir/blast
blastFile=$resdir/megablast_result.m8
#readLength=$inputDir/$(basename $inputDir/$reads .fasta)"_rlength.tab"
readLength=$(basename $reads .fasta)"_rlength.tab"

breakpointDir=$resdir/breakpoint
#breakpointDir=$resdir/${scaffold_name}/breakpoint
drjPairs=$breakpointDir/drjPairs.tab
seqCoordPairs=$breakpointDir/seqCoordPairs.tab
vectorDir=$breakpointDir/vectors
alignDir=$breakpointDir/align
figureDir=$breakpointDir/figures
segmentationResult=$breakpointDir/segmentation.tab
multipleLog=$breakpointDir/list-multiple-exaequo.txt
pairDir=$resdir/drjPairs_all_segments
pairDirFinal=$resdir/drjPairs_merged_segments
confirmedDrjPairs=$resdir/drjPairs_confirmed.tab
summaryDir=$resdir/drjPairs_figures
alignPairDir=$resdir/drjPairs_alignments
vectorPairDir=$alignPairDir/vectors


mkdir -p $blastDir

formatdb -i $genome -p F -o T
megablast -i $reads -d $genome -o $blastFile -m8 > $blastDir/megablast.log 2>&1
mv $resdir/formatdb.log $blastDir/

echo "Megablast finished"

#perl get_read_length.pl  -reads $inputDir/$nb_reads -out $readLength

mkdir -p $breakpointDir
#echo "pwd breakpoint"
#echo $(pwd)
#echo $breakpointDir

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

source("draw_summary.R")
drawAllSummary(pairDir="$pairDir",outputDir="$summaryDir",readNames=FALSE,addSimi="$vectorPairDir",merge=TRUE)
MyRScript4

echo "Summary plots finished"
