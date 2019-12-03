# DrjBreakpointFinder

DrjBreakpointFinder is a pipeline adapted from the software [Cassis](http://pbil.univ-lyon1.fr/software/Cassis/) to discover Direct Repeat Junctions (DRJ) of proviral segments and to precisely localize the breakpoint or excision site inside the DRJ sequences.

DrjBreakpointFinder is a [Genscale](http://team.inria.fr/genscale/) and [BIPAA](https://bipaa.genouest.org/is/) tool, developed by:

* Stéphanie Robin
* Claire Lemaitre

## Installation

```
git clone git@github.com:stephanierobin/DrjBreakpointFinder.git
```

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

Example with the small test dataset:

```
sh pipeline.sh -r test/reads.fasta -g test/genome.fasta -i test -o test/drj
```

## User manual


## Input
* `reads.fasta`: a fasta file containing reads from virus sequencing (circular form of the proviral segments).
* `genome.fasta`: a fasta file containing the host genome sequence.

## Output

An output directory containing different subdirectories :
* blast
* breakpoint
* drjPairs_alignments
* drjPairs_all_segments
* drjPairs_figures
* drjPairs_merged_segments



## 2 pipelines



## References

If you use DrjBreakpointFinder, please cite:

Legeai F., Santos B.F., Robin S., Bretaudeau A., Dikow R.B., Lemaitre C., Jouan V., Ravallec M., Drezen J-M., Tagu D., Gyapay G., Zhou X., Liu Shanlin, Webb B.A., Brady S.G., and Volkoff A-N. 2019. Conserved and specific genomic features of endogenous polydnaviruses revealed by whole genome sequencing of two ichneumonid wasps. Preprint on BioRxiv: [https://www.biorxiv.org/content/10.1101/861310v1](https://www.biorxiv.org/content/10.1101/861310v1).

DrjBreakpointFinder is inspired from the software [Cassis](http://pbil.univ-lyon1.fr/software/Cassis/), in particular for the precise identification of the breakpoint (or excision site) location, for more details on the algorithms, see:

Lemaitre C., Tannier E., Gautier C., Sagot M.-F.. [Precise detection of rearrangement breakpoints in mammalian genomes.](http://www.biomedcentral.com/1471-2105/9/286)  *BMC Bioinformatics*, **2008** 9(1):286.

Baudet C., Lemaitre C., Dias Z., Gautier C., Tannier E. and Sagot M-F. 2010. [Cassis: Detection of genomic rearrangement breakpoints.](http://bioinformatics.oxfordjournals.org/cgi/content/short/26/15/1897) *Bioinformatics*, **2010** 26(15):1897-1898.

# Contact

To contact a developer, request help, or for any feedback on DrjBreakpointFinder, please use the issue form of github: https://github.com/stephanierobin/DrjBreakpointFinder/issues