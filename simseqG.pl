#!/usr/bin/perl -w
# SimSeqG : SimSeq from Genoscope - A.Gilly
use strict;
    use List::Util qw(min max);
use List::MoreUtils qw(uniq);
use Getopt::Long;
# arguments: bam file, sample=1M, read len=76
##############################################
##############################################
# BEWARE!!!! THIS PROGRAM LOADS THE ENTIRE FASTA IN MEMORY. DO NOT USE ON LARGE GENOMES
##############################################
##############################################
my $VERSION="1.06";
my $READLEN=76;
#my $BAM="/env/cns/proj/projet_AKL/scratch/S/ilmap/test10000h-0.6.1/data.sorted.rmdup.bam";
#my $FASTA="/env/cns/proj/projet_AKL/db/test/Arabidopsis_thaliana.fa";
#my $FASTA="chr1.fa";
my $BAM=""; my $FASTA="";
my $OUTPUT_DIR=".";
my $SSIZE=500000;
my $SAMTOOLS_PATH="/env/cns/src/samtools/samtools-0.1.8/bin/samtools";
my $NUMSIM=60000000;
my $MP_UPLIM=10000;
my $PE_UPLIM=700;
my $OFFSET=0;
# num= samtools view -c, 
$|=1;
my $HELP="";
my $result = GetOptions ( 	'bam=s'		=>\$BAM,
				'out=s'		=>\$OUTPUT_DIR,
				'reference=s'   =>\$FASTA,
				'rlen=s'            =>\$READLEN,
				'sample=i'      =>\$SSIZE,
				'reads=i'       =>\$NUMSIM,
				'mplim=i'       =>\$MP_UPLIM,
				'pelim=i'       =>\$PE_UPLIM,
				'offset=i'      =>\$OFFSET,
				'h|help'             =>\$HELP
				);

usage() if($BAM eq "" || $FASTA eq "" || $HELP ne "");
if (! -d "$OUTPUT_DIR") {
    print "Making directory $OUTPUT_DIR...\n";
    (system("mkdir $OUTPUT_DIR") == 0) or die "Directory creation failed ($!).\n";
}

# print("Calculating sampling step...\n");
# #my $total=`samtools view -c $BAM`;
# #my $step=int(($total-$SSIZE)/1000000);
# my $step=67;
# $step/=2;
our %faindx;
loadFasta();
my %fastadict;
use List::Util qw(sum);

foreach my $key (keys %faindx){
    $fastadict{$key}=length($faindx{$key});
    print("$key ",length($faindx{$key}),"\n");
}
print "INFO:\tFasta sequence loaded in mem (",sum(values(%fastadict)), " bytes).\n";
my %hsample;
my %hpesample;
#print("Sampling 1 in ", $step+1, " reads...\n");
my $name="";
my $sens='';
my $pos=0;
my @qualdist;
my %chrdict;
my $unmapped=0;
open ( IN, "$SAMTOOLS_PATH view -h $BAM|" ) || die "Cannot open $BAM, $!\n";
my $MP_FOR_MAXPE=0;
    while (<IN>){
	my $line=$_;
	if($line =~/^@/){
#	    if(($line =~ /chr1/ || $line =~ /chr4/) && $line =~/^\@SQ\sSN:(.+)\sLN:(\d+)/){
	    if($line =~/^\@SQ\sSN:(.+)\sLN:(\d+)/){
		$chrdict{$1}=$2;
	    }
	    next();
	}
	my @data = split('\t', $line);
	if (($data[1] & 4) == 4){$unmapped++;next;};
	print STDERR "Sampled ", scalar(values %hsample)  ," reads.\r" if( $. % 1000==0);
	if(scalar(values %hpesample)>1000 && $MP_FOR_MAXPE==0){$MP_FOR_MAXPE=scalar(values(%hsample));}
 	last if ( scalar(values %hsample)>$SSIZE && scalar(values %hpesample) > 1000 );
	if($data[6] eq "=" && ( ($data[1] & 4) == 0 ) && ( ($data[1] & 8) == 0 ) && abs($data[8])<10000 && join(' ', @data) =~ /NM:i:[0-2]/){
		if (!(exists($hpesample{$data[0]})) && 
		    ((($data[1] & 16) == 16 && $data[8]>0) ||
		     (($data[1] & 16) == 0 && $data[8]<0))){
		    # This is MP
		    $hsample{$data[0]} = (abs($data[8])+$READLEN);
		    for(my $i=0; $i<$READLEN; $i++){
			push(@{$qualdist[$i]}, ord(substr($data[10], $i, 1)));
		    }
		}
		# This is PE
		$hpesample{$data[0]} = (abs($data[8])+$READLEN) if (!(exists($hsample{$data[0]})) && abs($data[8]) < $PE_UPLIM &&
								 ((($data[1] & 16) == 16 && $data[8]<0) ||
								 (($data[1] & 16) == 0 && $data[8]>0)));
	}
#	my @data = split('\t', $_);
 	# if($data[6] eq "=" && ( ($data[1] & 4) == 0 ) && ( ($data[1] & 8) == 0 ) && abs($data[8])<10000){
	#     if($name eq "" || $name ne $data[0]){
	# 	# first in pair or single
	# 	$name=$data[0];
	# 	$pos=$data[3];
	# 	$sens= (($data[1] & 16) == 16) ? '-': '+';
	#     }else{
	# 	if(($pos<$data[3] && $sens eq '+' && $data[1] & 16 == 16) ||
	# 	   ($pos>$data[3] && $sens eq '-' && $data[1] & 16 == 0)){
	# 	    # This is PE
		    
	# 	    push(@pesample, (abs($data[8])+$READLEN)) if(join(' ', @data) =~ /NM:i:0/);
	# 	}else{
	# 	    print STDERR "Processed ", $., " reads.\r";
	# 	    # This is MP
	# 	    push(@sample, (abs($data[8])+$READLEN));
	# 	}
	#     }
	# }
}
close(IN);
my @sample=sort({$a <=> $b} values %hsample);
my @pesample=sort({$a <=> $b} values %hpesample);

%hsample=();
%hpesample=();

print STDERR "\nDetected a ratio of one PE every ", int($MP_FOR_MAXPE/1000), " MP, ",$unmapped," unmapped total.\n";# " MP. Assuming overestimation of 50%.\n";

#####COMMENT THE FOLLOWING LINE IF NO OVERESTIMATION SHOULD BE ASSUMED#####
#####
##### It is advised only to do this when the mapping is really trustable (i.e. bwasw).
#$MP_FOR_MAXPE*=2;
#####
###########################################################################


# print "x<-c(";
# print(join(',', @sample));
# print(");\n");

# print "y<-c(";
# print(join(',', @pesample));
# print(");\n");

foreach my $i (0..$#qualdist){
#    print "Quality for base $i : there were ", scalar(@{$qualdist[$i]}), "samples, ";
#    print "--------------Base $i\n";
#    print join("\t", sort({$a <=> $b} uniq(@{$qualdist[$i]}))), "\n";
#    print chr(min(uniq(@{$qualdist[$i]}))), "\t",chr(max(uniq(@{$qualdist[$i]}))) , "\n";
    # my @temp=buildCECHisto(sort({$a <=> $b} @{$qualdist[$i]}));
    # $qualdist[$i]=\@temp;
#    print "now histo with ", scalar(@{$qualdist[$i]}), "classes, max, ", max(@{$qualdist[$i]}), ".\n ";
#    print("Position $i\n");
    print STDERR "Computing quality distribution for position $i\r";
    my %histq=histog(@{$qualdist[$i]});
    @{$qualdist[$i]}=[];
    $qualdist[$i]=\%histq;
#    while ( my ($key, $value) = each(%histq) ) { print "$key => $value\n"; }
    
}

print "Computed quality distribution for all positions.       \n";
# build constant-surface histogram for MP/PE
my @mpkeys=buildCECHisto(@sample);
my @pekeys=buildCECHisto(@pesample);
#print("br<-c(0,", join(",", (sort keys %hist)));
#my $maximum=max(sort keys %hist)+500;
my $maximum=max(@sample);
#print ",$maximum);\n";
#print("Sampled ", scalar(@sample), " reads.\nMin was ", min(@sample), ", max was ", max(@sample), ".\n");

my $i=0;

#####################
# MAIN SAMPLING LOOP
#####################

my $sum = sum(values %fastadict);
print $sum;
open (OUT1, ">$OUTPUT_DIR/1.fastq");
open (OUT2, ">$OUTPUT_DIR/2.fastq");
while($i<$NUMSIM){
    if($i % 100 ==0){print STDERR "Simulated $i reads, ",int(($i/$NUMSIM)*100),"% of total.\r";}
    # Choose chromosom & position
    my $position=int(rand($sum+1));
    my $offset=0;
    my $chr="";
    foreach my $key (keys(%fastadict)){
	if($position< ($offset+$fastadict{$key})){
	    $chr=$key;
	    $position-=$offset;
	    last;
	}
	$offset += $fastadict{$key};
    }
    # choose a fragment size
    my $frsize=sampleFromCECH(20, $MP_UPLIM, \@mpkeys);
    $frsize +=sampleFromCECH(20, $PE_UPLIM, \@pekeys);
    if($position+$frsize > $fastadict{$chr}){next;}

    # fetch corresponding sequence
    my $sequence=fetchFromFasta($chr, $position, $frsize);

    # choose random splice width around biotin
    my $splice_l=sampleFromCECH(20, $PE_UPLIM, \@pekeys);
    
    # choose random position on sequence
    ## To have the exact same rate of chimerism as in the original bam, one would have to compute the unmapped/mapped rate and:
    ## - Once in U/M, choose splice_0 at less than 76bp from the end.
    ## - Rest of the time, do it like this:
    my $splice_0=int(rand($splice_l))+($frsize-$splice_l);
    my $left="";
    my $right="";
    # add a PE in the mix once in a while
    if($i % int($MP_FOR_MAXPE/1000)==0){
#	print "PE start, frag size $frsize\n";
	my $plicepos=0;
	my $tries=0;
	while(length($left)<$READLEN && $tries<100){
	    $tries++;
	    $plicepos=sampleFromCECH(20, $PE_UPLIM, \@pekeys);
	    next if(2*$plicepos+2*$READLEN>length($sequence));
	    $left=substr($sequence, $plicepos, $READLEN);
	}
	if ($tries==101){ $i++; next;}
	$tries=0;
	while(length($right)<$READLEN && $tries<100){
	    $tries++;
	    my $rpos=sampleFromCECH(20, $PE_UPLIM, \@pekeys);
	    next if($rpos+$plicepos+$READLEN>length($sequence));
	    $right=Reverse_Complement(substr($sequence, $plicepos+$rpos, $READLEN));	    
	};
	if ($tries==101){ $i++; next;}

#	print "PE end.\n";
	goto writeseq;
	next;
    }

    # apply interval around extremities starting at position offset 
    $left=substr($sequence, $splice_0, $splice_l);
    $right=substr($sequence, 0, $splice_l-(length($sequence)-$splice_0));
    $sequence=$left.$right;
    next if(length($sequence)<2*$READLEN);
    my $lend=length($left);

    # extract mate pairs
    $left=substr($sequence, 0, $READLEN);
    $right=Reverse_Complement(substr($sequence, length($sequence)-$READLEN, $READLEN));
  writeseq: # this label is used only when a PE has been fabricated instead of a MP. It simply jumps over all the MP part.
    if(length($left) != $READLEN || length($right) != $READLEN){$i++;next;};
    my $qual1="";
    my $qual2="";
    for(my $j=0; $j<$READLEN; $j++){
	my $q1=chr(simflathistog($qualdist[$j]));
	my $q2=chr(simflathistog($qualdist[$j]));
	$qual1.=$q1;
	$left=glitchBaseFromQual($left, $q1, $j);
	$qual2.=$q2;
	$right=glitchBaseFromQual($right, $q2, $j);
    }
    my $handle1; my $handle2;
    if(rand()>0.5){$handle1=*OUT1;$handle2=*OUT2;}else{$handle2=*OUT1;$handle1=*OUT2;}
    my $identifiant=$i+$OFFSET;
    print $handle1 "\@SimseqG-$identifiant\n$left\n+SimseqG-$identifiant\n$qual1\n";
    print $handle2 "\@SimseqG-$identifiant\n$right\n+SimseqG-$identifiant\n$qual2\n";
    $i++;
}
close(OUT1);
close(OUT2);
print "\n";

sub fetchFromFasta{
   my ( $chr, $pos, $len )= @_;
   return (substr($faindx{$chr}, $pos, $len));
}

sub loadFasta{
    open ( IN, $FASTA ) || die "Cannot open reference: $!\n";
    my $currchr="";
    my $newchr;
    
    while ( <IN> ){
	chomp;
	my $line=$_;
	if($line =~/>(.*)/){
	    $faindx{$newchr}=$currchr if($currchr ne "");
	    $newchr=$1;
	    $currchr="";
	}
	else{
	    $currchr.=$line;
	}
    }
    $faindx{$newchr}=$currchr;
    undef($currchr);
    close(IN);
    
}

sub Reverse_Complement {
   my ($dna) = @_;
   my $revcomp = reverse($dna);
   $revcomp =~ tr/ACGTacgt/TGCAtgca/;
   return $revcomp;
}

sub buildCECHisto{
    my @tab=@_;
#    @tab=sort(@tab);
    my $eff=scalar(@tab);
    my %hist;
    my $borne=min(@tab);
    my $cureff=0;
    foreach my $val (@tab){
	$cureff+=1;
	if($cureff>($eff/20)){
	    $cureff=0;
	    $borne=$val;
	    $hist{$borne}=1;
	}
	else{
	    $hist{$borne}+=1;
	}
    }
    return sort({$a <=> $b} keys %hist);
    
}

sub sampleFromCECH{
    my($classes, $uplim, $tref)=@_;
    my @skeys=@$tref;
    my $index=int(rand($classes+1));
    if($index==0){
	return int(rand(min(@skeys)+1));
    }
    if($index>($classes-1)){
	return (int(rand($uplim-max(@skeys)+1))+max(@skeys));
    }
    return (int(rand($skeys[$index]-$skeys[$index-1]+1))+$skeys[$index-1]);

}

sub simflathistog{
    my ($href)=@_;
    my %hist=%$href;
    my $qualite=rand();
    foreach my $val (sort({$a <=> $b} keys %hist)){
	if($qualite<=$hist{$val}){return $val;}
    }
    

}

sub histog{
    my @tab=@_;
    @tab=sort({$a <=> $b} @tab);

    my $val=$tab[1];
    my $val0=$val;
    my %freq;
    for(my $i=0; $i<scalar(@tab); $i++){
	if($tab[$i] != $val){$val=$tab[$i];$freq{$val}++;}else{
	    $freq{$val}++;
	}
    }
    my $oldval=0;
    foreach $val (sort({$a <=> $b} keys %freq)){
	if ($val == $val0){$oldval=$freq{$val};}else{$freq{$val}+=$oldval;$oldval=$freq{$val};}
    }
    foreach $val (sort({$a <=> $b} keys %freq)){
	$freq{$val}/=$oldval;
    }
    return %freq;
}

sub glitchBaseFromQual{
    my ($sequence, $qual, $pos)=@_;
    my %alphabet=("A", 1 , "C", 1, "T", 1, "G", 1);
    my $odds=10**((-(ord($qual)-64))/10);
#    print "rand num must be < ", $odds, " for base to be crap. (", (-(ord($qual)-64))/10, ")\n";
    if(rand() < $odds){
	delete $alphabet{substr($sequence, $pos, 1)};
	my $replacement=(keys %alphabet)[int(rand(3))];
#	print STDERR "sequence : ", substr($sequence, 0, $READLEN), "\nposition : ", $pos, "\nreplacement : ", $replacement, "\n";
	substr($sequence, $pos, 1, $replacement);
	return $sequence;
    }else{
	return $sequence;
    }
}

sub usage{
    print("
SimSeqG Mate-Pair sequencing simulator v.$VERSION (2012)

Simulates mate-paired reads that posess similar features to those found in a BAM file: insert size, paired-end proportion and quality distribution for each base are inferred from a read sample. The program automatically generates chimeric reads and sequencing errors. The program will output two FASTQ files.

Usage:

simseqG.pl -bam bamfile -fasta sequence.fasta [OPTIONS...]

Options:

	-bam [filename]	        Mandatory. Specifies the BAM file that the program should use in order to learn read characteristics.
	-reference [filename]	Mandatory. A FASTA sequence that will be used to randomly build reads from. The contigs in this file must match those in the BAM.
	-out	[dirname]	Optional. Output directory. Defaults to current directory (.).
	-rlen [length]		Optional. Integer. Read length. Defaults to 76.
	-sample [size]		Optional. Integer. Minimum number of pairs to be extracted from the BAM for inferring error levels. Defaults to 1 million.
	-reads [number]	        Optional. Number of pairs to simulate. Defaults to 50 million.
	-mplim [maxsize]	Optional. Maximum length for mate-paired insert size. Defaults to 10k. Leave as is if unsure.
	-pelim [maxsize]	Optional. Maximum length for paired-end insert size. Defaults to 700. Leave as is if unsure.
	-offset [number]	Optional. Offset to add to the IDs of reads. Use this if you plan to merge several simseqG outputs.
	-help			This help.


Caveat: this program loads the entire FASTA in memory. Use with caution on large genomes.
Contact: A. Gilly, 3°Et, #319\n");
    exit();

}
