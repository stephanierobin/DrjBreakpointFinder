#! /bin/sh

# #######################################################################################################
# DrjBreakpointFinder identifies the precise recombination breakpoints in Direct Repeat Junctions (DRJs) 
# of viral insertions in the host genome. 
#
# Authors: Stéphanie Robin, Claire Lemaitre
# #######################################################################################################

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
. /local/env/envconda.sh
conda activate  --stack /groups/bipaa/local/blast
conda activate  --stack /groups/bipaa/local/clustalw
conda activate  --stack /groups/bipaa/local/emboss
source /local/env/envr-3.5.1.sh

#src DIR
EDIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
srcDIR=$EDIR/src

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

mkdir -p $output
mkdir -p $blastDir

echo "Running Megablast..."

makeblastdb -in $genome -input_type fasta -dbtype nucl
blastn -task megablast -query $reads -db $genome -outfmt 6 -out $blastDir/megablast_result.m8 2>&1

echo "Megablast finished"

perl $srcDIR/get_read_length.pl  -reads $reads -out $readLength

mkdir -p $breakpointDir


echo "Selecting sequence triplets..."

R --slave --vanilla --quiet --no-save  <<MyRScript1

mydir = getwd()
setwd("$srcDIR")
source("get_fully_overlapping_triplets.R")
setwd(mydir)
res=getFullyOverlappingTriplets(blastFile="$blastFile",readLengthFile="$readLength",coordFile="$seqCoordPairs",drjPairsFile="$drjPairs")

MyRScript1

 echo "Triplets selection finished"

echo "Aligning triplets..."

 perl $srcDIR/align_and_get_mismatch_vectors.pl -coord $seqCoordPairs -bacs $genome -reads $reads -dirFasta $breakpointDir/tmp -dirAlign $alignDir -dirVector $vectorDir

 echo "All triplet Alignments finished"


 mkdir -p $pairDir
 mkdir -p $pairDirFinal

echo "Performing segmentation..."

R --slave --vanilla --quiet --no-save  <<MyRScript2

mydir = getwd()
setwd("$srcDIR")
source("get_segmentation.R")
source("get_segmentations_by_DRJpairs.R")
source("merge_segments_on_drj_pairs.R")
setwd(mydir)

segmentASetOfSequences("$seqCoordPairs","$vectorDir","$figureDir","$segmentationResult",zoom=T,margin=20)

getAllValidSegmentationsByDRJPair(segmentationFile="$segmentationResult",drjPairsFile="$drjPairs",pairDir="$pairDir",confirmedDrjPairFile="$confirmedDrjPairs",multipleLog="$multipleLog")

mergeAllDrjs(inputPairDir="$pairDir",outputPairDir="$pairDirFinal")

MyRScript2

echo "All segmentations finished"

echo "Aligning DRJ pairs..."

mkdir -p $alignPairDir
perl $srcDIR/align_drjPairs.pl -coord $drjPairs -bacs $genome -dirFasta $alignPairDir/tmp -dirAlign $alignPairDir/align -dirVector $vectorPairDir

echo "Alignments of DRJ pairs (similarity) finished"

mkdir -p $summaryDir

echo "Drawing summary plots..."

R --slave --vanilla --quiet --no-save  <<MyRScript4

mydir = getwd()
setwd("$srcDIR")
source("draw_summary_scaffolds.R")
setwd(mydir)
drawAllSummary(pairDir="$pairDir",outputDir="$summaryDir",readNames=FALSE,addSimi="$vectorPairDir",merge=TRUE)

MyRScript4

echo "Summary plots finished"
