#!/usr/bin/perl
#*****************

#use strict;
use warnings;
use Getopt::Long;


sub Usage
{
    print STDOUT "
-----------------------------------------------------------------------
Create a table with the length of each read, from a fasta file
-----------------------------------------------------------------------
Usage : $0
Mandatory
   -reads    :  input read file (fasta format)
   -out      :  output file 

Details : output format = one line per sequence with its name and length separated by a space
-----------------------------------------------------------------------
\n";
    exit(0);
}


MAIN:
{
    my ($reads,$output_file);

    GetOptions(# lecture de la ligne de commande
	"reads=s"   => \$reads,
	"out=s"   => \$output_file
	);
    &Usage if ((not defined $reads) || (not defined $output_file));

	open(R,"<$reads") or die "cannot open file $reads\n";
	open(outFile,">$output_file") or die "cannot open file $output_file\n";

	my $myLength=0;
	my $name="";

	while(<R>){
		if($_ =~ m/^>(.*)$/){
		    if($name ne ""){
			print outFile $name." ".$myLength."\n";
		    }
		    $name=$1;
		    $myLength=0;
		}
		if($_ =~ m/^(A|T|C|G|-)/){
		    chomp();
		    $myLength+=length($_);
		}
	}
    print outFile $name." ".$myLength."\n";
    
    close(R);
    close (outFile);
}
