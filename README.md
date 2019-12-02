# DrjBreakpointFinder

DrjBreakpointFinder is a pipeline adapted from Cassis to discover Direct Repeat Junctions (DRJ) of proviral segments and potential breakpoint in DRJs.

## Reference

Baudet C., Lemaitre C., Dias Z., Gautier C., Tannier E. and Sagot M-F. 2010. Cassis: detection of genomic rearrangement breakpoints. Bioinformatics.

Legeai F., Santos B.F., Robin S., Bretaudeau A., Dikow R.B., Lemaitre C., Jouan V., Ravallec M., Drezen J-M., Tagu D., Gyapay G., Zhou X., Liu Shanlin, Webb B.A., Brady S.G., and Volkoff A-N. 2019. Conserved and specific genomic features of endogenous polydnaviruses revealed by whole genome sequencing of two ichneumonid wasps. Submitted.

## Installation

* git clone git@github.com:stephanierobin/DrjBreakpointFinder.git

## Requirements

* perl
* R
* bash
* BLAST
* clustalw

## Usage

```
sh pipeline.sh -r reads.fa -g genome.fa -i input_directory -o output_directory
```

## User manual


## Input
* a fasta file containing reads from virus sequencing.
* a fasta file containing genome sequence.

## Output

An output directory containing different subdirectories :
* blast
* breakpoint
* drjPairs_alignments
* drjPairs_all_segments
* drjPairs_figures
* drjPairs_merged_segments



## 2 pipelines
