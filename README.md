# DrjBreakpointFinder

DrjBreakpointFinder is a pipeline adapted from Cassis to discover Direct Repeat Junctions (DRJ) of proviral segments and potential breakpoint in DRJs.

## Reference

Baudet C., Lemaitre C., Dias Z., Gautier C., Tannier E. and Sagot M-F. 2010. Cassis: detection of genomic rearrangement breakpoints. Bioinformatics.

Santos B. F., Legeai F., Robin S., Bretaudeau A., Dikow R. B., Lemaitre C., Jouan V., Marc Ravallec M., Drezen J-M., Tagu D., Gyapay G., Zhou X., Liu S., Webb B. A., Brady S. G. and Volkoff A-N. Conserved and specific genomic features of endogenous polydnaviruses revealed by whole genome sequencing of two ichneumonid wasps.

## Installation

*
*


## Requirements

## Usage

```
sh pipeline.sh -r reads.fa -g genome.fa -o output_directory
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
