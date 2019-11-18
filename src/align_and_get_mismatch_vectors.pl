#!/usr/bin/perl
#*****************

use strict;
use warnings;
use Getopt::Long;


sub Usage
{
    print STDOUT "
-----------------------------------------------------------------------
Script to align and create the mismatch vector for each read-bac triplet
-----------------------------------------------------------------------
Usage : $0
Mandatory
   -coord    :  input triplet file with the read and bac names, and the coordinates on the bac
Optionnal
   -bacs     :  path to the all_bacs file (default='../blast/final_bacs_new.fasta')
   -reads    :  path to the all_reads file (default='../blast/QZ_A_clean_min100_new.fa')
   -dirFasta :  name of the fasta dir inside the main dir (default='tmp')
   -dirAlign :  name of the output alignment dir inside the main dir (default='align')
   -dirVector :  name of the vector dir inside the main dir (default='vector')
-----------------------------------------------------------------------
\n";
    exit(0);
}


MAIN:
{
    my ( $coord, $bacs, $reads, $dirFasta, $dirAlign, $dirVector);

    $bacs="../blast/final_bacs_new.fasta";# chemin des BACs
    $reads="../blast/QZ_A_clean_min100_new.fa"; # chemin des reads
    $dirFasta="tmp";
    $dirAlign="align";
    $dirVector="vector";

    GetOptions(# lecture de la ligne de commande
	"coord=s"   => \$coord,
	"bacs=s" => \$bacs,
	"reads=s"  => \$reads,
	"dirFasta=s"    => \$dirFasta,
    	"dirAlign=s"    => \$dirAlign,
	"dirVector=s"    => \$dirVector,
	);

    &Usage if not defined $coord;

    my $customDRJ=0;

    my %all_ids=();# création des hash vides
    my %all_bacs=();
    my %all_reads=();
    ## récupère les noms et les ID des séquences à garder
    ## ----------------------------------------
    ## (pour stocker que les séquences d'intéret)
    open(R1,"<$coord") or die "cannot open file $coord\n";
    my $entete=<R1>;
    while(<R1>){
        chomp();# suppressions du '\n'
        my ($id,$read_name,$bac_name,$inf1,$sup1,$inf2,$sup2,$orient)=split();# recupere les info
        $all_ids{$id}++;
        $all_reads{$read_name}++;#affectation des nom et chgt de clé
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
    my %read_sequence=();

    ## On crée les fichiers fasta-------------------------------------------------
    ## les reads
    open(READS,"<$reads") or die "cannot open file $reads\n";
    my $current_name="";
    my $ok=0;
    my $compt=0;
    while(<READS>){
        chomp();
        if (/^>.*$/){
            $ok=0;
            ($current_name)= $_ =~ /^>(.*)$/ ;
            if (defined $all_reads{$current_name}){
            $ok=1;
            $read_sequence{$current_name}="";
            }
        }
        if($ok & ! /^>/){
            $read_sequence{$current_name}.=$_;
        }
    }
    close(READS);

    ## les bacs
    open(BACS,"<$bacs") or die "cannot open file $bacs\n";
    $current_name="";
    $ok=0;
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

    ## Pour chaque triplet, On lit les coordonnées des séquences et on crée les fasta :
    ## --------------------------------------------------------------------------------
    open(R2,"<$coord") or die "cannot open file $coord\n";
    $entete=<R2>;
    while(<R2>){
        chomp();
        my ($id,$read_name,$bac_name,$inf1,$sup1,$inf2,$sup2,$orient)=split();

        ## read sequence
        my $readSeq=$read_sequence{$read_name};
        if ($orient==-1){
            $readSeq=revcomp($readSeq);
        }
        my $fastaR=seq_fasta($readSeq,$read_name);

        ## BAC sequences :
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
        my $fasta_name=$dirFasta."/".$id.".fa";
        open(W4,">$fasta_name")  or die "cannot open file $fasta_name\n";
        print W4 $fastaR."\n";
        print W4 $fastaB1."\n";
        print W4 $fastaB2."\n";
        close(W4);

        # On crée le fichier aln
	my $align=$pathAlign."/".$id.".aln";
	my $commande="clustalw -infile=".$fasta_name." -align -type=dna -outfile=". $align." -output=fasta -newtree=tmptree.dnd > clustal.log";
	system($commande);

        # On supprime le fasta déjà utilisé
        unlink($fasta_name) or die "Could not delete the file!\n";
        unlink("tmptree.dnd") or die "Could not delete the file!\n";
	unlink("clustal.log") or die "Could not delete the file!\n";

        # On crée les vecteurs de mismatch
        #my $infile=$pathAlign."/".$id.".aln";
        my $out=$pathVector."/".$id.".tab";
	extract_mismatch_vectors_from_ma($align,$out,$read_name);
        #$commande="perl ".$dirSource."/extract_mismatch_vectors_from_ma.pl -infile ".$align." -outfile ".$out;
        #system($commande);
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

sub write_fasta {

    my $sequence=shift;
    my $name=shift;
    my $file_name=shift;
    my $max_char=shift;
    $max_char=60 unless defined $max_char;

    open(W,">$file_name")  or die "cannot open file $file_name\n";

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
    print W $fasta."\n";

    close(W);
}

sub revcomp {

    my $seq=shift;
    my $revcomp = reverse($seq);
    $revcomp =~ tr/ACGTacgt/TGCAtgca/;
    return $revcomp;
}

sub extract_mismatch_vectors_from_ma {
     
    my $infile=shift;
    my $outfile=shift;
    my $rname=shift;
    
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
	     if($nseq==1 && not defined $rname){
		 $rname=$name;
	     }
	 }
	if(/^(\w|-)/){
	    $sequences{$name}.=$_;
	}
    }
    close(R);
    
    if(defined $rname){ 
	die "$rname is not a sequence of the input file $infile" if not defined $sequences{$rname};
    }
    
    ## attributing the right name to each sequence
    ## -------------------------------------------
    my %bac_inf=();
    foreach my $seq_name (keys(%sequences)){
	
	if ($seq_name=~/^bac/){
	    my ($inf)= $seq_name=~/^bac_(\d+)-/;
	    $bac_inf{$seq_name}=$inf;
	}
	else{
	    if (not defined $rname){
		$rname=$seq_name;
	    }
	}
    }
    
    my @bacs = sort { $bac_inf{$a} <=> $bac_inf{$b} } keys(%bac_inf);
    
    ## Finding the mismatches
    ## ----------------------
    my $al_length=length($sequences{$rname});
    my $ref_letter="";
    my $b1_letter="";
    my $b2_letter="";
    my @vector1;
    my @vector2;
    
    # pour compter où on en est sur lesséquences de bacs
    my $compt1=0;
    my $compt2=0;
    my @vecPos1;
    my @vecPos2;
    for (my $i=0;$i<$al_length;$i++){
	
	$ref_letter=substr($sequences{$rname},$i,1);
	$b1_letter=substr($sequences{$bacs[0]},$i,1);
	$b2_letter=substr($sequences{$bacs[1]},$i,1);
	if($b1_letter ne "-"){$compt1++;}
	if($b2_letter ne "-"){$compt2++;}
	
	if($ref_letter eq "-"){
	    my $t=@vector1;
	    if($t>0){
		if($b1_letter ne "-") {$vector1[-1]=1;}
		if($b2_letter ne "-") {$vector2[-1]=1;}
	    }
	}
	else{
	    if($b1_letter eq $ref_letter) {push(@vector1,0);} else {push(@vector1,1);}
	    if($b2_letter eq $ref_letter) {push(@vector2,0);} else {push(@vector2,1);}
	    push(@vecPos1,$compt1);
	    push(@vecPos2,$compt2);
	}

    }
	
    ## writing the results
    open(OUTPUT,">$outfile") or die "cannot open file $outfile\n";
    my $t=@vector1;
    for (my $i=0;$i<$t;$i++){
	print OUTPUT $vector1[$i]." ".$vector2[$i]." ".$vecPos1[$i]." ".$vecPos2[$i]."\n";
	# en plus des deux vecteurs de mismatch le long du read, on note a quelle position ca correspond dans les deux séquences de bac
    }
    close(OUTPUT);
}

