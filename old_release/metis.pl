#!/usr/bin/perl -w
# svfinder
# MA Madoui & JM Aury, Genoscope, 2011

use strict;
use Getopt::Long;
use POSIX qw(floor);
use List::Util qw[min max];
#use diagnostics;

# PARAMETERS
my $VERSION="0.23";
our @RERUNS;
our $TR_ERROR="";
our $REFERENCE="";
our ( $BAM , $INPUT_DIR , $LIB , $EXT , $MAX_COV , $COV_EXTENSION , $READ_SIZE , $SAMTOOLS_PATH , $DOC_MAD , $SPAN_MAD , $SAMPLE_SIZE , $SLC , $MIN_READS, $OUTPUT_DIR, $NOSORT, $NOLUTHER, $MINCOV, $CLSTPRM, $ANNOTATE, $PICARD, $SCORE_DONORS, $BEDTOOLS_PATH, $DONOR_ANNOT, $ACCEPTOR_ANNOT) = ( "", "" , "MP" , 20000 , 1000, 500, 76, "/env/cns/src/samtools/samtools-0.1.8/bin/samtools" , 4 , 4 , 1000000 , "/env/cns/proj/projet_AKL/slc/slclust", 10, "", "", "", 255, "", "", "~agilly/picard-tools-1.65", "", "/env/cns/src/BEDTools/BEDTools-Version-2.16.2/bin", "/env/cns/proj/projet_AKL/db/TAIR10te.bed", "~/TAIR-Complete.bed");
my $result = GetOptions ( 	'discordant=s'	=>\$BAM,
				'in=s'		=>\$INPUT_DIR,
				'lib=s'		=>\$LIB,
				'ext=i'		=>\$EXT,
				'maxcov=i'	=>\$MAX_COV,
				'covext=i'	=>\$COV_EXTENSION,
				'read=i'	=>\$READ_SIZE,
				'samtools=s'	=>\$SAMTOOLS_PATH,
				'dm=f'		=>\$DOC_MAD,
				'sm=f'		=>\$SPAN_MAD,
				'sample=i'	=>\$SAMPLE_SIZE,
				'slc=s'		=>\$SLC,
				'minreads=i'	=>\$MIN_READS,
				'nosort=s'      =>\$NOSORT,
				'annotate'      =>\$ANNOTATE,
				'score_donors'  =>\$SCORE_DONORS,
				'out=s'         =>\$OUTPUT_DIR,
				'bedtools=s'    =>\$BEDTOOLS_PATH,
				'donor_annot=s' =>\$DONOR_ANNOT,
				'acc_annot=s'   =>\$ACCEPTOR_ANNOT,
				'transposition_errors' =>\$TR_ERROR,
				'reference=s'   =>\$REFERENCE);


my $RATIO;
$| = 1;


usage () if ( ($BAM  eq "" && $SCORE_DONORS ne "") || $OUTPUT_DIR eq "" || $LIB !~/^MP$|^PE$/ || $INPUT_DIR eq "" );

my $X;
my $Y;
my $MIN;
my $READS_LENGTH;
my $ISMEDIAN;
my $ISMAD;
my $DOCMEDIAN;
my $DOCMAD;
our $RERUN="NO";




############################################################################
`mkdir -p $OUTPUT_DIR` if ( $OUTPUT_DIR && ! -d $OUTPUT_DIR );
############################################################################

if($TR_ERROR ne "" && $SCORE_DONORS eq ""){
#    die "\nERROR:\t donors should be annotated in order for R-transcriptase error detection to be performed. \n\tEnable with -score_donors.\n\n";
}

print "\nSVfinder Metis module v. $VERSION - Annotation and scoring of translocations\n\n";
#print "\n\"Now Zeus, king of the gods, made Metis his wife first, and she was wisest among gods and mortal men.\"\n";
#printf("%80s", 'Hesiod, Theogony');
print "\n\n";
print "Discordant input:\t$BAM\n";
print "Input dir\t:\t$INPUT_DIR\n";
print "Output dir\t:\t$OUTPUT_DIR\n\n";
print "Donor scoring\tscheduled.\n" if($SCORE_DONORS ne "");
print "Translocation annotation\tscheduled.\n" if ($ANNOTATE ne "");
print "R-transcriptase error detection\tscheduled.\n" if ($TR_ERROR ne "");

print "\n";


die("ERROR:\tDiscordant BAM file does not exist.\n") if (!(-e $BAM));


if($SCORE_DONORS ne "") {
    print("INFO:\tScoring donors...\n");
    score_donors();
}

if($TR_ERROR ne ""){
    print("INFO:\tLooking for R-transcriptase errors...\n");
    rt_detect();
}

if($ANNOTATE ne "") {
$ANNOTATE="/env/cns/proj/projet_AKL/db/TAIR10te.bed";
    print("Annotating...\n");
    die "\nERROR:\tNo annotation file to read from. Goodbye.\n" if ($ANNOTATE eq "");
    die "\nERROR:\tInvalid annotation file. Goodbye.\n" if (! -e $ANNOTATE);
    annotate_TE();
}


sub rt_detect{
    if (!-e "$BAM.bai"){
	`$SAMTOOLS_PATH index $BAM`;
	die "Indexing of $BAM failed, samtools exited abnormally.\n" if($? != 0);
    }
    open (IN, "<:utf8", "$OUTPUT_DIR/sv.formatted") or die "ERROR:\tCould not open $OUTPUT_DIR/sv.formatted for RT-error detection.\n";
    my @SV=<IN>;
    print "INFO:\t$#SV lines read.\n";
    close(IN);
    @SV=grep {$_->[2] eq "TRANS"} map({[split("\t", $_)]} @SV);
    my $lastid=$SV[$#SV]->[1];
    for(my $i=1; $i<$lastid+1; $i++){
	print "INFO:\tGenerating pileup for event $i...\r";
	my @event=grep({$_->[1]==$i} @SV);
	# get donor with max score
	my $mscore;
	if($#event == 0){
	    #there is only one donor, no need to find max.
	    $mscore=$event[0];
	}
	else{
#	    print map{print @$_} @event;
	    my $local_max_score=max(map{$_->[$#{@$_}]} @event);
	    $mscore=(grep({$_->[$#{@$_}]==$local_max_score} @event))[0];
	    die "ERROR:\tNo maximum was found for donor $i. A maximum always exists.\n" if(!defined($mscore));
	}

	# samtools view the reads that map on the acceptor/donor.
	my $donor=$mscore->[7].":".$mscore->[8]."-".$mscore->[9];
	my $dchr=$mscore->[7];
	my $achr=$mscore->[3];
	my $astart=$mscore->[4]-7000;
	my $aend=$mscore->[5]+7000;
	open (IN, "$SAMTOOLS_PATH view $BAM $donor |") or die "ERROR:\tImpossible to view $BAM at region $donor.\n";
	`$SAMTOOLS_PATH view -H $BAM > $OUTPUT_DIR/tmp$i.sam`;
	open (OUT, ">>$OUTPUT_DIR/tmp$i.sam") or die "ERROR:\t opening of $OUTPUT_DIR/tmp$i.sam failed.\n" ;
	while (<IN> ){
	    chomp;
	    my $line=$_;
	    my $matchr= ($dchr eq $achr? "=":$achr);
	    my @data=split('\t', $line);
	    if($data[6] eq $matchr && $data[7] > $astart && $data[7] < $aend){
		print OUT "$line\n";
	    }
	}
	close(IN);
	close(OUT);
	
	# build corresponding pileup, examine.
	open(IN, "$SAMTOOLS_PATH pileup -S -f $REFERENCE $OUTPUT_DIR/tmp$i.sam 2> /dev/null|") or die "ERROR:\tcould not build pileup. Samtools version too old, or missing reference?\n";
	while(<IN>){
	    chomp;
	    my @data=split("\t", $_);
	    my $tmp=$data[4];
	    my $count = ( $tmp =~ tr/[atgcATGC]//);
	    print "\nSuspect position found for event $i at position $data[1]\n" if($count/length($data[4]) > 0.5 && length($data[4])>3);

	}
	close(IN);
    }
	print "INFO:\tDone treating $lastid events...                \n";

}
sub score_donors{
    if (!-e "$BAM.bai"){
	`$SAMTOOLS_PATH index $BAM`;
	die "Sorting of $BAM failed, samtools exited abnormally.\n" if($? != 0);
    }
    open (INSTR, "<:utf8", "$INPUT_DIR/sv.formatted") or die "ERROR:\tCould not open sv.formatted in $INPUT_DIR\n";
    open(OUTSTR, ">:utf8", "$OUTPUT_DIR/sv.tmp") or die "ERROR:\tCould not create sv.tmp in $OUTPUT_DIR\n";
    my $accid=1;
    #$acceptor_ID est dans la col 1.
    #le donneur est 7 8 9
    my %donors;
    my $line;
    my @data;
    while(<INSTR>){
	chomp;
	$line=$_;
	@data=split("\t", $line);
	die("ERROR:\tInput file sv.formatted can't be parsed, please check the format ($#data).\n") if($#data !=13);
#	print "accid=$accid, data1=$data[1]\n";
	if($data[1] == $accid){	$donors{"$data[7]:$data[8]-$data[9]"}=$line;}
	else{
	    if(scalar(keys(%donors))==1) {print OUTSTR $donors{(keys(%donors))[0]}, "\t*\n";%donors=();$donors{"$data[7]:$data[8]-$data[9]"}=$line;$accid++;}else{
		my @d=keys(%donors);
		my @tmp=split('\t', $donors{(keys(%donors))[0]});
		my $acceptor="$tmp[3]:$tmp[4]-$tmp[5]";
		print "ACCEPTOR : $acceptor \n =============";
		my %donor_count=discriminate_donor($BAM,$acceptor ,@d);
		# foreach my $key (sort { (split("\t", $donors{$a}))[3] cmp (split("\t", $donors{$b}))[3] ||
		# 		 (split("\t", $donors{$a}))[4] <=> (split("\t", $donors{$b}))[4]} keys %donors){
		#     print "Donor $key scored $donor_count{$key}.\n";
		#     print OUTSTR $donors{$key}, "\t", $donor_count{$key},"\n";
		    
		# }
		    open (TMP, "<:utf8", "$INPUT_DIR/sv.formatted");
		    while(<TMP>){
			chomp;
			my @td=split("\t", $_);
			foreach my $donor (keys(%donors)){
			    my $score =(defined($donor_count{$donor})?$donor_count{$donor}:0 );
			    if($donor eq "$td[7]:$td[8]-$td[9]" && $td[1] == $accid){
				print "Donor $donor scored $score.\n";
				print OUTSTR $donors{$donor}, "\t", $score,"\n";
			}
		    }
		}
		close(TMP);
		%donors=();
		$donors{"$data[7]:$data[8]-$data[9]"}=$line;
		$accid=$data[1];
	    }
	}
	
    }
	    if(scalar(keys(%donors))==1) {print OUTSTR $donors{(keys(%donors))[0]}, "\t*\n";%donors=();$donors{"$data[7]:$data[8]-$data[9]"}=$line;$accid++;}else{
		my @d=keys(%donors);
		my @tmp=split('\t', $donors{(keys(%donors))[0]});
		my $acceptor="$tmp[3]:$tmp[4]-$tmp[5]";
		print "ACCEPTOR : $acceptor \n =============";
		my %donor_count=discriminate_donor($BAM,$acceptor ,@d);
		# foreach my $key (sort { (split("\t", $donors{$a}))[3] cmp (split("\t", $donors{$b}))[3] ||
		# 		 (split("\t", $donors{$a}))[4] <=> (split("\t", $donors{$b}))[4]} keys %donors){
		#     print "Donor $key scored $donor_count{$key}.\n";
		#     print OUTSTR $donors{$key}, "\t", $donor_count{$key},"\n";
		    
		# }
		    open (TMP, "<:utf8", "$INPUT_DIR/sv.formatted");
		    while(<TMP>){
			chomp;
			my @td=split("\t", $_);
			foreach my $donor (keys(%donors)){
			    my $score =(defined($donor_count{$donor})?$donor_count{$donor}:0 );
			    if($donor eq "$td[7]:$td[8]-$td[9]" && $td[1] == $accid){
				print "Donor $donor scored $score.\n";
				print OUTSTR $donors{$donor}, "\t", $score,"\n";
			}
		    }
		}
		close(TMP);
		%donors=();
		$donors{"$data[7]:$data[8]-$data[9]"}=$line;
		$accid=$data[1];
	    }

    close(IN);
    close(OUT);
    `mv $OUTPUT_DIR/sv.tmp $OUTPUT_DIR/sv.scored.formatted`;
    `cat $OUTPUT_DIR/sv.scored.formatted | cut -f3- > $OUTPUT_DIR/sv.scored.out` if (-e "$OUTPUT_DIR/sv.scored.formatted");
    die "Unable to convert $OUTPUT_DIR/sv.formatted into sv.out.\n" if ($? != 0);
}

    # for each donor
    # see which reads have pairs mapping around acceptor and are uniquely mapped.


# expects an input under the form : chrX:start-end.
# returns a hash containing for each donor the n° of reads matching uniquely to it
sub discriminate_donor{
    my ($file, $acceptor, @donor)=@_;
    my %donor_uniq_count;
    my $achr=(split(':', $acceptor))[0];
    my $astart=(split('-', (split(':', $acceptor))[1]))[0];
    my $aend=(split('-', (split(':', $acceptor))[1]))[0];
    $astart-=7000;
    $aend+=7000;

    foreach my $i (0..$#donor){
	my $count=0;
	print "FOR DONOR :", $donor[$i], "\n";
	my $tmp=$donor[$i];
#	open(OUT, ">$OUTPUT_DIR/$acceptor.$tmp.sam");
	open (IN, "$SAMTOOLS_PATH view $file $donor[$i] |") or die
"Impossible to discriminate donor because there is no file.\n";
	while (<IN> ){
	    chomp;
	    my $line=$_;
	    my $dchr=(split(':', $donor[$i]))[0];
	    my $matchr= ($dchr eq $achr? "=":$achr);
	    my @data=split('\t', $line);
	    my $rlen=length($data[9]);
	    if($data[6] eq $matchr && $data[7] > $astart && $data[7] < $aend){
#		print OUT "$line\n";
#		if($line =~ /MD:Z:$READ_SIZE/ && $line =~ /X0:i:1/){
		if($data[5] =~ /$rlen/ && $line =~ /X0:i:1\s/){
#		print $line, "\n";
		    $donor_uniq_count{$donor[$i]}++;
#	       $line =~ /MD:Z:$READ_SIZE/ ){
		}
	    }
	}
	close(IN);
#	close(OUT);
    }
    return %donor_uniq_count;




}	

# sub discriminate_donor{
#     my ($file, @donor)=@_;
#     my %donor_uniq_count;
#     foreach my $i (0..$#donor){
# 	my $count=0;
# 	open (IN, "$SAMTOOLS_PATH view $file $donor[$i] |") or die
# "Impossible to discriminate donor because there is no file.\n";
# 	while (<IN> ){
# 	    chomp;
# 	    my $line=$_;
#             #BEWARE: Uninitialized value in one of the splits in next line
# 	    my $tmp=(split('XA:Z:', $line))[1];
# 	    if(!defined($tmp) or $tmp eq ""){$count++;next;}
# 	    $tmp=(split ('\t', $tmp))[0];
# 	    my @alternative_hits=split(';', $tmp);

# 	    my $hit=0;
# 	    foreach my $alternative_hit (@alternative_hits){ #pour chaque donneur, on considère chaque read couvrant ce donneur
# 		my @alternative_data=split(',', $alternative_hit);
# 		my $alternative_chr=$alternative_data[0]; # si le read mappe sur un autre donneur, on le jette
# 		my $alternative_start=abs($alternative_data[1]);
# 		for my $j ($i..$#donor){
# 		    if($alternative_chr eq (split(':', $donor[$j]))[0] &&
# $alternative_start < (split('-', (split(':', $donor[$j]))[1]))[1]){
# 			$hit=1;
# 			last;
# 		    }
# 		}
# 		last if($hit==1);
# 	    }
# 	    if($hit==0){
# 		$count++;
# 	    }


# 	}
# 	close(IN);
# 	my $totalcount =`$SAMTOOLS_PATH view -c $file $donor[$i]`;
# 	die "Viewing of $file failed, samtools exited abnormally.\n" if($? != 0);
# 	$donor_uniq_count{$donor[$i]}=($count/$totalcount);
	
#     }
#     return %donor_uniq_count;
# }





#############################################################################
#                           annotate_TE                                                             #
#----------------------------------------------------------------------------                      ###
# Direct copy from SV_FP.pl                                                                       ##### 
# Intersects sv.out with an annotation file to select known TE donor sites.                        ###
# The rest is discarded. Not for use with de-novo discovery of TEs.                                 #
sub annotate_TE {
    my $DIR=$OUTPUT_DIR;
    my $BED_INTER ="$BEDTOOLS_PATH/intersectBed";
    my $F_OUT=((-e "$DIR/sv.scored.out")? "$DIR/sv.scored.out" : "$DIR/sv.out");
    my $F_FOR=((-e "$DIR/sv.scored.formatted")? "$DIR/sv.scored.formatted" : "$DIR/sv.formatted");
    print "Using $F_OUT and $F_FOR";
`cat $F_OUT | grep \"TRANS\" | sort -k6,6 -k7,7n | awk '{OFS=\"\t\";print \$6,\$7,\$8,\$9,\$10, \$11}'> $DIR/TE_donor.tab`;
`cat $F_OUT | grep \"TRANS\" | sort -k6,6 -k7,7n | awk '{OFS=\"\t\";print \$2, \$3, \$4, \$5, \$6, \$7, \$8, \$9, \$10, \$11, \$12}'> $DIR/TE_acceptor.tab`;
`$BED_INTER -a $DIR/TE_acceptor.tab -b $ACCEPTOR_ANNOT -wao > $DIR/TE_acceptor_quesneville.tab`;
`$BED_INTER -a $DIR/TE_donor.tab -b $DONOR_ANNOT -wao > $DIR/TE_donor_quesneville.tab`;

my %donor;

# On traite les TEs qui auraient sauté dans d'autres TEs: Les accepteurs sont contenus dans des régions annotées
open (IN, "$DIR/TE_acceptor_quesneville.tab") || die "$!,\n";
    while ( <IN> ){
	chomp;
	my @data = split /\t/ , $_;
	my $found=0;
	# elsif($donor { "$data[0]:$data[1]..$data[2]" }->[2] != @data[5]){
	#     my $key= "$data[0]:$data[1]..$data[2]" ;
	#     my @chars=('a'..'z');
	#     $key.=$chars[rand @chars];
	#     push
	# }
	
	# generating a list of acceptors/donors
	# we look for an already existing acceptor
        # if there is one, then we check if they have the same donor. If yes, we add the supplementary annotation for that donor
	# if no, it means the acceptor isn't registered yet or it exists but with another donor. Depending on the case, we append annotation or generate a new key.
 	my $key_re="$data[0]:$data[1]..$data[2]";
	foreach my $key (grep(/^$key_re/, keys(%donor))){
	    if($donor{$key}->[1] eq $data[4] && $donor{$key}->[2] == $data[5] && $donor{$key}->[3] == $data[6]){
		$found=1;
		$donor { $key }->[11] .=";".$data[14]; # Si l'accepteur est annoté, on met tous les TEs potentiels dans le champ 11
		last;
	    }
	}
	if($found==0){
	    if ( ! exists $donor { "$data[0]:$data[1]..$data[2]" } ){ 
		push @{$donor { "$data[0]:$data[1]..$data[2]" }} ,  @data[3..14]; # On indexe les donneurs par leurs accepteurs 
	    } 
	    else{
#		$donor { "$data[0]:$data[1]..$data[2]" }->[11] .=";".$data[14]; # Si l'accepteur est annoté, on met tous les TEs potentiels dans le champ 14
		my $key= "$data[0]:$data[1]..$data[2]" ;
		while(exists $donor { $key }){
		    $key="$data[0]:$data[1]..$data[2]" if ($key ne "$data[0]:$data[1]..$data[2]");
		    my @chars=('a'..'z', 'A'..'Z');
		    my $letter=$chars[int rand scalar @chars];
		    $key .=$letter;
		}
		push @{$donor { $key }} ,  @data[3..14];
	    }
	}
	$found=0;
    }
close(IN);

# avant cette étape il existe des keys de %donor pour lesquelles le champ 12 contient 
# On traite les TEs dont le donneur est annoté
    foreach my $site (keys %donor){
	#print $_,"\t", join "\t" , @{$donor{$_}}, "\n";
	open (IN, "$DIR/TE_donor_quesneville.tab") || die "$!,\n";
	while ( my $l = <IN> ){
		chomp $l;
		my @data = split /\t/ , $l;
		if ($donor{$site}->[1] eq $data[0] && $donor{$site}->[2] == $data[1] && $donor{$site}->[3] == $data[2] && $data[9]=~/\w/){ #si le
			$donor{$site}->[12].=";".$data[9];
		}
	}
	close(IN);
}

# Remove incomprehensible duplicates in donor & acceptor annotation
    foreach my $site (keys %donor){
	my %seen = ();
	my @uniq;
	if(defined $donor{$site}->[12]){
	    #$item eq "" || 
	    foreach my $item (split(/;/, $donor{$site}->[12])) {
		push(@uniq, $item) unless ($seen{$item}++);
	    }
	    $donor{$site}->[12]=join(";", grep(!/^$/, @uniq));
	}else{$donor{$site}->[12]='.';}
	@uniq=();
	if(defined $donor{$site}->[11]){
	    foreach my $item (split(/;/, $donor{$site}->[11])) {
		push(@uniq, $item) unless ($seen{$item}++);
	    }
	    $donor{$site}->[11]=join(";", grep(!/^$/, @uniq));
	}else{$donor{$site}->[11]='.';}

    }

# Post traitement pour récupérer les informations de donneurs multiples
    # Ajoute à la fin de chaque élément de %donor un ID unique par accepteur
    `cat $F_FOR | egrep "TRANS" > $DIR/sv.formatted.tmp`;
    open (IN, "<:utf8", "$DIR/sv.formatted.tmp") || die "$!,\n";
    while ( <IN> ){
	chomp;
	my @data = split /\t/ , $_;
	my $key_re=$data[3].":".$data[4]."..".$data[5];
	foreach my $key (grep(/^$key_re/, keys(%donor))){
	    if($donor{$key}->[$#{@{$donor{$key}}}] ne $data[1]){
		push(@{$donor{$key}}, $data[1]);
	    }
	}
#	push(@{$donor{$key}}, $data[0]);
    }
    close(IN);
    `rm $DIR/sv.formatted.tmp`;
    
    
    open (OUT, ">:utf8", "$OUTPUT_DIR/annotated.tmp") || die "Unable to output : $!.\n";
    foreach ( keys %donor ){
	my $key=$_;
#	my $bracket= $donor{$key}->[$#{@{$donor{$key}}}];
	# $bracket=(ord(substr($bracket,0,1))<125)?" ":$bracket;
	# pop(@{$donor{$key}});
	my $id= $donor{$key}->[$#{@{$donor{$key}}}];
	pop(@{$donor{$key}});
	my $keypure=$key;
	chop($keypure) if ($keypure=~m/.+\D$/);
	print OUT $id, "\t", $keypure,"\t", join "\t" , @{$donor{$key}}, "\n";
    }
    
    close(OUT);
    $DIR=$OUTPUT_DIR;
    `sort -n $OUTPUT_DIR/annotated.tmp > $DIR/annotated.1`;
#    `rm $DIR/annotated.tmp`;
#    `cut -f13 $DIR/sv.out > $DIR/tmp`;
    open (IN, "<:utf8", "$OUTPUT_DIR/annotated.1");
    open (OUT,">:utf8", "$OUTPUT_DIR/annotated.out");
    while (<IN>){
	chomp;
	my $line=$_;
	my @data=split('\t', $line);
	my $pattern=join("\\t", @data[2..7]);
	my $corr_line=`grep -P $pattern $F_OUT`;
	open (TMP, "<:utf8","$F_OUT");
	my @lines=grep /$pattern/, <TMP> ;
	$corr_line=$lines[0];
#	print "grep -P \"$pattern\" $DIR/sv.out";
	chomp;
	my @tmp=split('\t', $corr_line);
	my $score=$tmp[$#tmp];
	print OUT $line, "\t", $score;
    }
    close(IN);
    close(OUT);
 #   `paste --delimiters="" $DIR/annotated.1 $DIR/tmp > $DIR/annotated.out`;
   `rm $DIR/annotated.1`;
    die ("ERROR:\t Command \"rm $DIR/annotated.1\" failed.\n") if(($? >> 8) != 0);
    `cut -f1 $F_FOR > $DIR/tmp`;
    die ("ERROR:\t Failed to extract ID column from output file.\n") if(($? >> 8) != 0);
    `paste $DIR/tmp $DIR/annotated.out > $DIR/annotated.f`;
    die ("ERROR:\t Failed to append ID column to output.\n") if(($? >> 8) != 0);
    `cut -f1-11,15- $DIR/annotated.f > $DIR/annotated.formatted`;
    die ("ERROR:\t Failed to remove intersection columns in formatted output file.\n") if(($? >> 8) != 0);
    `mv $DIR/annotated.out $DIR/annotated.o; cut -f1-10,14- $DIR/annotated.o > $DIR/annotated.out`;
    die ("ERROR:\t Failed to remove intersection columns in output file.\n") if(($? >> 8) != 0);
    `rm $DIR/annotated.o $DIR/annotated.f`;
    die ("ERROR:\t Post-run cleaning(1) failed.\n") if(($? >> 8) != 0);
    `rm $DIR/tmp`;
    die ("ERROR:\t Post-run cleaning(2) failed.\n") if(($? >> 8) != 0);

}
