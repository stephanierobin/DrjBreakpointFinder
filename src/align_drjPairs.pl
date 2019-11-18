#!/usr/bin/perl
#*****************

use strict;
use warnings;
use Getopt::Long;


sub Usage
{
    print STDOUT "
-----------------------------------------------------------------------
Script to align and create the mismatch vector for each pair of DRJs
-----------------------------------------------------------------------
Usage : $0
Mandatory
   -coord    :  input drjPairs file with the bac name, and the coordinates on the bac
Optionnal
   -bacs     :  path to the all_bacs file (default='../blast/final_bacs_new.fasta')
   -dirFasta :  name of the fasta dir inside the main dir (default='tmp')
   -dirAlign :  name of the output alignment dir inside the main dir (default='align')
   -dirVector :  name of the vector dir inside the main dir (default='vector')
-----------------------------------------------------------------------
\n";
    exit(0);
}


MAIN:
{
    my ( $coord, $bacs, $dirFasta, $dirAlign, $dirVector);

    $bacs="../blast/final_bacs_new.fasta";# chemin des BACs
    $dirFasta="tmp";
    $dirAlign="align";
    $dirVector="vector";

    GetOptions(# lecture de la ligne de commande
	"coord=s"   => \$coord,
	"bacs=s" => \$bacs,
	"dirFasta=s"    => \$dirFasta,
    	"dirAlign=s"    => \$dirAlign,
	"dirVector=s"    => \$dirVector,
	);

    &Usage if not defined $coord;

    my %all_bacs=();
    ## récupère les noms et les ID des séquences à garder
    ## ----------------------------------------
    ## (pour stocker que les séquences d'intéret)
    open(R1,"<$coord") or die "cannot open file $coord\n";
    my $entete=<R1>;
    while(<R1>){
        chomp();# suppressions du '\n'
        my @splity=split();# recupere les info
	my $bac_name=$splity[0];
         $all_bacs{$bac_name}++;
    }
    close(R1);

    my $pathFasta=$dirFasta;
    `mkdir $pathFasta` unless (-d $pathFasta);
    my $pathAlign=$dirAlign;
    `mkdir $pathAlign` unless (-d $pathAlign);
    my $pathVector=$dirVector;
    `mkdir $pathVector` unless (-d $pathVector);


    my %bac_sequence=();

    ## les bacs
    open(BACS,"<$bacs") or die "cannot open file $bacs\n";
    my $current_name="";
    my $ok=0;
    while(<BACS>){
        chomp();
        if (/^>.*$/){
            $ok=0;
            ($current_name)= $_ =~ /^>(.*)$/ ;
            if (defined $all_bacs{$current_name}){
		$ok=1;
		$bac_sequence{$current_name}="";
            }
        }
        if($ok & ! /^>/){
            $bac_sequence{$current_name}.=$_;
        }
    }
    close(BACS);

    ## Pour chaque drjPair, On lit les coordonnées des séquences et on crée les fasta :
    ## --------------------------------------------------------------------------------
    open(R2,"<$coord") or die "cannot open file $coord\n";
    $entete=<R2>;
    while(<R2>){
        chomp();
        my ($bac_name,$inf1,$sup1,$inf2,$sup2,@others)=split();
	my $prefix=join("_",($bac_name,$inf1,$sup1,$inf2,$sup2));

        ## BAC sequence :
        my $full_seq=$bac_sequence{$bac_name};

        my $l1=$sup1-$inf1+1;
        my $seq1=substr($full_seq,$inf1-1,$l1);
        my $name1="bac_".$inf1."-".$sup1;

        my $l2=$sup2-$inf2+1;
        my $seq2=substr($full_seq,$inf2-1,$l2);
        my $name2="bac_".$inf2."-".$sup2;

        my $fastaB1=seq_fasta($seq1,$name1);
        my $fastaB2=seq_fasta($seq2,$name2);

        # On crée le fichier fa
        my $fasta_name=$dirFasta."/".$prefix.".fa";
        open(W4,">$fasta_name")  or die "cannot open file $fasta_name\n";
	print W4 $fastaB1."\n";
        print W4 $fastaB2."\n";
        close(W4);

	# On crée le fichier aln
	my $align=$pathAlign."/".$prefix.".aln";
	my $commande="clustalw -infile=".$fasta_name." -align -type=dna -outfile=". $align." -output=fasta -newtree=tmptree.dnd > clustal.log";
	system($commande);

        # On supprime le fasta déjà utilisé
        unlink($fasta_name) or die "Could not delete the file!\n";
        unlink("tmptree.dnd") or die "Could not delete the file!\n";
	unlink("clustal.log") or die "Could not delete the file!\n";

        # On crée les vecteurs de mismatch
        my $out=$pathVector."/".$prefix.".tab";
	extract_mismatch_vectors_from_drj_ma($align,$out);
    }
    close(R2);
    rmdir($dirFasta);
}


sub seq_fasta {

    my $sequence=shift;
    my $name=shift;
    my $max_char=shift;
    $max_char=60 unless defined $max_char;

    my $fasta=">".$name;
    if($max_char>0){
        my $index=0;
        my $length=length($sequence);
        while ($index<$length){
            $fasta=$fasta."\n".substr($sequence,$index,$max_char);
            $index=$index+$max_char;
        }
    }
    else{
        $fasta=$fasta."\n".$sequence;
    }
    return $fasta;
}


sub extract_mismatch_vectors_from_drj_ma {
     ## to get the differences between drj1 and drj2 along drj1

    my $infile=shift;
    my $outfile=shift;
    
    my %sequences=();
    my $name;
    my $nseq=0;
    
    open(R,"<$infile") or die "cannot open file $infile\n";
    while(<R>){
	chomp();
	 if (/^>.*$/){
	     ($name)= $_ =~ /^>(.*)$/ ;
	     $nseq++;
	     $sequences{$name}="";
	 }
	if(/^(\w|-)/){
	    $sequences{$name}.=$_;
	}
    }
    close(R);

    ## attributing the right name to each sequence
    ## -------------------------------------------
    my %bac_inf=();
    foreach my $seq_name (keys(%sequences)){
	
	if ($seq_name=~/^bac/){
	    my ($inf)= $seq_name=~/^bac_(\d+)-/;
	    $bac_inf{$seq_name}=$inf;
	}
    }
    
    my @bacs = sort { $bac_inf{$a} <=> $bac_inf{$b} } keys(%bac_inf);
    
    ## Finding the mismatches
    ## ----------------------
    my $al_length=length($sequences{$bacs[0]});
    my $b1_letter="";
    my $b2_letter="";
    my @vector1;
    
    # pour compter où on en est sur lesséquences de bacs
    my $compt1=0;
    my @vecPos1;

    for (my $i=0;$i<$al_length;$i++){
	
	$b1_letter=substr($sequences{$bacs[0]},$i,1);
	$b2_letter=substr($sequences{$bacs[1]},$i,1);
	if($b2_letter ne "-"){$compt1++;}
	
	if($b1_letter eq "-"){
	    my $t=@vector1;
	    if($t>0){
		if($b2_letter ne "-") {$vector1[-1]=1;}
	    }
	}
	else{
	    if($b1_letter eq $b2_letter) {push(@vector1,0);} else {push(@vector1,1);}
	    push(@vecPos1,$compt1);
	}

    }
	
    ## writing the results
    open(OUTPUT,">$outfile") or die "cannot open file $outfile\n";
    my $t=@vector1;
    for (my $i=0;$i<$t;$i++){
	print OUTPUT $vector1[$i]." ".$vecPos1[$i]."\n";
	# en plus du vecteur de mismatch le long de b1, on note a quelle position ca correspond dans b2
    }
    close(OUTPUT);
}

