#!/usr/bin/perl -w
# SVfinder Eris Module - A.Gilly
# MA Madoui & JM Aury, Genoscope, 2011

use strict;
use Getopt::Long;
use POSIX qw(floor);
use Switch;
use List::Util qw[min max];
use Devel::Peek;
use FindBin '$Bin';
use Config::Simple;

$|=0;
my %config;
Config::Simple->import_from("$Bin/TE-Tracker.conf", \%config) or die ("\nERROR:\tImpossible to read main configuration file $Bin/TE-Tracker.conf. Please reinstall TE-Tracker. \n\tError message was: \"",Config::Simple->error(),"\"\n");

my $BASEDIR=$config{"GLOBAL.basedir"};
my $VERSION=0.1;
my $PICARD_PATH="/env/cns/src/picard/picard-tools-1.67";
# PARAMETERS
our @RERUNS;
our $NODISCORDANT="";
our ( $BAM , $OUTPUT_DIR , $LIB , $EXT , $MAX_COV , $COV_EXTENSION , $READ_SIZE , $SAMTOOLS_PATH , $DOC_MAD , $SPAN_MAD , $SAMPLE_SIZE , $SLC , $REFERENCE, $NOSORT, $NOLUTHER, $CONCORDANT, $DISCORDANT) = ( "", "" , "MP" , 20000 , 1000, 500, 76, "/env/cns/src/samtools/samtools-0.1.8/bin/samtools" , 4 , 4 , 1000000 , "/env/cns/proj/projet_AKL/slc/slclust", "", "", "", "no", "no");
our @TREAT_BAM=();
our $RUNDATA="";
our $MOCK="";
our $MAXMAP=10000;
our $MAXMIS=1;
our $CHAINLOAD="";
our $FILTER=100000;
my $MAN="";
my @args=@ARGV;
my $result = GetOptions ( 	'bam=s'		=>\$BAM,
				'out=s'		=>\$OUTPUT_DIR,
				'lib=s'		=>\$LIB,
				'ext=i'		=>\$EXT,
				'maxcov=i'	=>\$MAX_COV,
				'covext=i'	=>\$COV_EXTENSION,
				'read=i'	=>\$READ_SIZE,
				'samtools=s'	=>\$SAMTOOLS_PATH,
				'dm=f'		=>\$DOC_MAD,
				'sm=f'		=>\$SPAN_MAD,
				'sample=i'	=>\$SAMPLE_SIZE,
				'nosort=s'        =>\$NOSORT,
				"noluther"      =>\$NOLUTHER,
				'concordant:s'    =>\$CONCORDANT,
				'discordant:s'    =>\$DISCORDANT,
				'treat_bam=s'     =>\@TREAT_BAM,
				'mock'            =>\$MOCK,
				'rundata=s'       =>\$RUNDATA,
				'maxmap=i'        =>\$MAXMAP,
				'maxmis=i'        =>\$MAXMIS,
				'nodiscordant'    =>\$NODISCORDANT,
				'chainload=s'     =>\$CHAINLOAD,
				'f1=s'            =>\$FILTER,
				'picard_path=s'   =>\$PICARD_PATH,
				'man'             =>\$MAN
				);



$| = 1;
man() if ($MAN ne "");
usage() if ( $BAM  eq "" || $OUTPUT_DIR eq "" || $LIB !~/^MP$|^PE$/);
if(scalar @TREAT_BAM > 0 && $TREAT_BAM[0] ne ""){
    foreach my $treat (@TREAT_BAM){
	my @treat_modes=split(":", $treat);
	die("Malformed treat_bam argument: should be [input|discordant|concordant]:i[-j].\n") if (scalar(@treat_modes) !=2);
	if($treat_modes[0] ne "input" and $treat_modes[0] ne "concordant" and $treat_modes[0] ne "discordant"){
	    die("Malformed treat_bam argument, 1° member! Should be [input|concordant|discordant]\n");
	}
	if(!($treat_modes[1] =~ m/^\d(\-\d)?$/)){
	    die("Error in treat_bam parameter, 2° member! Must be of form i or i-j with i, j integers.\n");
	}
    }
}

if ($FILTER =~ /\D/ && $FILTER ne "disabled") {
    die ("ERROR:\tThe f1 argument is a numerical filter. Please use an integer or the string value 'disabled'.\n");
}
die ("ERROR:\t Main BAM-file $BAM not found.\n") if(!-e $BAM);
die ("ERROR:\t Pre-sorted BAM-file $NOSORT not found.\n") if ($NOSORT ne "" && (!-e $NOSORT));

my $X;
my $Y;
my $MIN;
my $READS_LENGTH;
my $ISMEDIAN;
my $ISMAD;
my $DOCMEDIAN;
my $DOCMAD;
our $RERUN="NO";


print "\nSVfinder Eris module v. $VERSION - Detection of discordant reads\n\n";
#print "\n\"So, after all, there was not one Eris (Discord) alone, but all over the earth there are two. As for the one, a man would praise her when he came to understand her.\"\n";
#printf("%80s", 'Hesiod, Works and Days, 11-24');
print("\n\nInput BAM :\t$BAM\n");
print("Sorted input :\t$NOSORT\n") if $NOSORT ne "";
print("Output dir:\t$OUTPUT_DIR\n");
print("Using rundata:\t$RUNDATA\n");
print("Level-1 intrachromosomal filter:\t$FILTER\n");
print("Read size:\t$READ_SIZE bp\n");
print("\n");


if ($NOLUTHER ne ""){
    print("-no-luther mode enabled. Not reforming reads.\n");
}
#print "-$MAXMIS mismatches allowed.\n";
print "-$MAXMAP mappings allowed.\n";
print("\n");

print("Estimation of insert size distribution\tscheduled.\nEstimation of coverage distribution\tscheduled.\n") if($RUNDATA eq "");
print("Detection of over-covered regions\tscheduled.\n") if (!-e "$OUTPUT_DIR/over_cov.txt");
print("Sorting of input file\t\tscheduled.\n") if ($NOSORT eq "");
print("Generation of discordant file\tscheduled.\n") if ($DISCORDANT ne "no");
print("Generation of concordant file\tscheduled.\n") if ($CONCORDANT ne "no");
foreach my $treat (@TREAT_BAM){
    my $TREAT_MODE=(split(":", $treat))[0];
    my $TREAT_ID=(split(":", $treat))[1];
    if(defined($TREAT_MODE) && defined($TREAT_ID)){
	print("Post-treatment of $TREAT_MODE file with $TREAT_ID mismatch filter\tscheduled.\n");
    }
}
print("Chainloading execution to last treated BAM file\tscheduled.\n") if($CHAINLOAD eq "treated");
print("Chainloading execution to concordant BAM file\tscheduled.\n") if($CHAINLOAD eq "concordant");
print("Chainloading execution to discordant BAM file\tscheduled.\n") if($CHAINLOAD eq "discordant");


############################################################################
`mkdir -p $OUTPUT_DIR` if ( $OUTPUT_DIR && ! -d $OUTPUT_DIR );
############################################################################

my @DOC ;
my @SPAN ;
my $root_dir=$OUTPUT_DIR;
my $MAIN_RUN=1;

if ($RUNDATA ne "" && -e "$RUNDATA"){
    open(IN, "$RUNDATA") or die "File $RUNDATA not found.\n";
    my $line=<IN>;
    chomp($line);
    @DOC=split(":", $line);
    print "INFO:\tRead DOC data: $line\n";
    $line=<IN>;
    @SPAN=split(":", $line);
    close(IN);
    print "INFO:\tRead insert size data: $line\n";
}else{
    @DOC=get_DOC_distribution ( $BAM , $DOC_MAD );
    @SPAN= get_insert_size_distribution ( $BAM , $SPAN_MAD );
    
    
    open(OUT, ">$OUTPUT_DIR/rundata.dat");
    print OUT join(":",@DOC), "\n";
    print OUT join(":",@SPAN);
    close(OUT);
}


my(@TRANS, @DEL);
my @CONCORDANTID;
my %OVER_COV = get_high_coverage_positions ( $BAM , $MAX_COV , $COV_EXTENSION , $READ_SIZE );
exit if($MOCK ne "");

our %CURRENT_DISCORDANT_PAIRS;
my $chemin_sorted;


if(scalar(@TREAT_BAM)>1 || !($TREAT_BAM[0]=~/input/)){
    bam2discordant_pairs ( $BAM , $OUTPUT_DIR , $SPAN[2] , $SPAN[3] , $LIB , $EXT );
}


if(lc($DISCORDANT) eq "sorted"){
    print("INFO:\tSorting discordant file...\n");
    `samtools sort $root_dir/discordant.bam $root_dir/discordant.sorted`;
    die("Samtools sort failed on  $root_dir/discordant.bam.\n") if($? != 0);
#    `rm $root_dir/discordant.bam; mv $root_dir/discordant.sorted.bam $root_dir/discordant.bam`;
}

if(lc($CONCORDANT) eq "sorted"){
    print("INFO:\tSorting concordant file...\n");
    `samtools sort $root_dir/concordant.bam $root_dir/discordant.sorted`;
    die("Samtools sort failed on  $root_dir/concordant.bam.\n") if($? != 0);
#    `rm $root_dir/concordant.bam; mv $root_dir/concordant.sorted.bam $root_dir/concordant.bam`;
}

foreach my $treat (@TREAT_BAM){
    my $TREAT_MODE=(split(":", $treat))[0];
    my $TREAT_ID=(split(":", $treat))[1];
    if(defined($TREAT_MODE) && defined($TREAT_ID)){
	switch ($TREAT_MODE) {
	    case "input" {
		if ($NOSORT eq "") {
		    print "INFO:\tSorting by name... (this step can be lengthy)\n";
		    sort_bam_by_name($BAM, "$OUTPUT_DIR/input_sortedby_name");
		    $chemin_sorted="$OUTPUT_DIR/input_sortedby_name.bam";
		}
		else{
		    $chemin_sorted=$NOSORT;
		}
		treat_bam($chemin_sorted, "$OUTPUT_DIR/input.$TREAT_ID.mismatch.ReadOrder.bam", $TREAT_ID);}
	    case "concordant" {treat_bam("$OUTPUT_DIR/concordant.bam", "$OUTPUT_DIR/concordant.$TREAT_ID.mismatch.ReadOrder.bam", $TREAT_ID);}
	    case "discordant" {treat_bam("$OUTPUT_DIR/discordant.bam", "$OUTPUT_DIR/discordant.$TREAT_ID.mismatch.ReadOrder.bam", $TREAT_ID);}
	}
    }
}

if($CHAINLOAD eq "treated"){
    my $TREAT_MODE=(split(":", $TREAT_BAM[$#TREAT_BAM]))[0];
    my $TREAT_ID=(split(":", $TREAT_BAM[$#TREAT_BAM]))[1];
    if(defined($TREAT_MODE) && defined($TREAT_ID)){
	my $mode=$TREAT_MODE;
	print "INFO:\tSorting treated BAM...\n";
	my $PICARD_SORT_TASK=`java -jar $PICARD_PATH/SortSam.jar INPUT=$OUTPUT_DIR/$TREAT_MODE.$TREAT_ID.mismatch.ReadOrder.bam OUTPUT=$OUTPUT_DIR/$TREAT_MODE.$TREAT_ID.mismatch.sorted.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT`;
	($?>>8) and die "ERROR:\tTask exited with errorlevel: ", ($?>>8), ". Exiting.\n";
	$NODISCORDANT="";
	$BAM="$OUTPUT_DIR/$TREAT_MODE.$TREAT_ID.mismatch.sorted.bam";
	$NOSORT="$OUTPUT_DIR/$TREAT_MODE.$TREAT_ID.mismatch.ReadOrder.bam";
	$RUNDATA="$OUTPUT_DIR/rundata.dat";
	if($RUNDATA eq ""){
	    open(IN, "$RUNDATA") or die "File $RUNDATA not found.\n";
	    my $line=<IN>;
	    chomp($line);
	    @DOC=split(":", $line);
	    $line=<IN>;
	    @SPAN=split(":", $line);
	    close(IN);
	}
	bam2discordant_pairs ( $BAM , $OUTPUT_DIR , $SPAN[2] , $SPAN[3] , $LIB , $EXT );

    }
}
#############################################################################
#			get_DOC_distribution
#----------------------------------------------------------------------------
# Take a SAM file, run samtools pileup
# Return an array that insert size distribution 
sub get_DOC_distribution {
    my ( $bam , $n_mad ) = @_;
    # Depth of coverage analysis	
    my @cov;	
    open ( IN, "$SAMTOOLS_PATH pileup $bam |" ) || die "Cannot open $bam, $!\n";
    while ( <IN> ){
	last if ( $. > $SAMPLE_SIZE );
	my @data = split /\t/, $_;
	push ( @cov, $data[3] );
    } 
    close (IN);
    my @doc = get_limits ( \@cov , $n_mad );
    $doc[2] = ( $doc[2] < 1 ) ? 1 : $doc[2] ;	
    print "INFO:\tDepth of coverage analysis\tDOC median: ",$doc[0],"\tDOC MAD: ",$doc[1],"\tDOC inf lim: ",$doc[2],"\tDOC sup lim: ",$doc[3],"\n";
    return @doc;
}
#############################################################################
#			get_insert_size_distribution
#----------------------------------------------------------------------------
# Take a BAM file, run samtools
# Return an array that contains insert size distribution (isd)
sub get_insert_size_distribution {
    my ( $bam , $n_mad ) = @_;
    # Insert size analysis
    my @insert;
    open ( IN, "$SAMTOOLS_PATH view $bam|" ) || die "Cannot open $bam, $!\n";
    while (<IN>){
	last if ( $. > $SAMPLE_SIZE );
	my @data = split /\t/, $_;
	if ( $data[6] eq "=" && $data[8] > 0 ){
	    push ( @insert, $data[8] );
	} 
    } 
    close (IN);
    my @isd = get_limits ( \@insert , $n_mad );
    print "INFO:\tInsert size analysis\tInsert size median: ", $isd[0] ,"\tInsert size MAD: ", $isd[1] ,"\tInsert size inf lim: ", $isd[2] ,"\tInsert size sup lim: " , $isd[3],"\tInsert shift: ", 2 * $n_mad * $isd[1] , "\n"; 
    return @isd;
}
###################################################################################
#			get_limits
#----------------------------------------------------------------------------------
#takes an array of insert sizes
#returns the median, median + n.MAD, median - n.MAD
sub get_limits {
    my ( $data , $n_mad ) = @_;
    my @quartiles = getBoxPlot ( \@$data );
    return	( $quartiles[1] , getMAD ( \@$data , $quartiles[1] ), $quartiles[1]- $n_mad * getMAD ( \@$data , $quartiles[1] ), $quartiles[1] + $n_mad * getMAD ( \@$data , $quartiles[1] ) );
}
###################################################################################
#			getMAD
#----------------------------------------------------------------------------------
# Take an array
# Return the median absolute deviation (MAD)
sub getMAD {
    my ( $array, $med ) = @_;
    my @dev;
    foreach ( @$array ){
	push ( @dev, abs( $_ - $med ) );
    }
    my @MAD = getBoxPlot( \@dev );
    return $MAD[1];
}
##################################################################################
#			getBoxplot
#---------------------------------------------------------------------------------
# Take an array of value 
# Return the quartiles Q1, Q2, Q3, inf and sup limit values 
sub getBoxPlot {
    my ( $array ) = @_;
    my @sorted_array = sort { $a <=> $b } @$array;
    my @results = (quartile ( 1, \@sorted_array ) ,  quartile ( 2, \@sorted_array ), quartile ( 3, \@sorted_array ));
    return ( $results[0] , $results[1] , $results[2] , $results[0] - 1.5 * ( $results[1] - $results[0] ) , $results[2] + 1.5 * ( $results[2] - $results[1] ) );
}

#############################################################################
#			quartile
#----------------------------------------------------------------------------
# Take a value and an integer (1, 2 or 3) 
# Return the quartile value Q1, Q2 or Q3
sub quartile {
    my ( $quart, $array ) = @_ ;
    my $K_quantile = ( ( $quart / 4 ) * ( $#{$array} - 1 ) + 1 );
    my $F_quantile = $K_quantile - POSIX::floor($K_quantile);	#decimal part
    $K_quantile = POSIX::floor($K_quantile);			#entire part
    my $aK_quantile = $array->[ $K_quantile - 1 ];
    return $aK_quantile if ( $F_quantile == 0 );			#if the decimal part is null
    my $aKPlus_quantile = $array->[$K_quantile];
    my $quantile = $aK_quantile + ( $F_quantile * ( $aKPlus_quantile - $aK_quantile ) );#if the decimal part is not null
	return $quantile;
}

###################################################################################
#			get_high_coverage_positions
#----------------------------------------------------------------------------------
# Take a BAM file
# Return a hash of overcovered region
sub get_high_coverage_positions {
    my ( $bam_file , $max_cov, $ext, $read_length ) = @_;
    my %over_cov;
    if (-e "$OUTPUT_DIR/over_cov.txt"){
	print "INFO:\tPrevious run detected. Fetching over-covered regions...";
	open(IN, "$OUTPUT_DIR/over_cov.txt");
	while (<IN>){
	    chomp;
	    my @data=split /\t/;
	    push ( @{$over_cov{$data[0]}}, [$data[1], $data[2]]);
	}
	close(IN);
	my $CMD=`cat $OUTPUT_DIR/over_cov.txt | wc -l`;
	chomp $CMD;
	print "Fetched ", $CMD, " lines.\n";
    }else{
	print "INFO:\tCalculating over-covered ( > $max_cov ) regions...";
	open ( IN , "$SAMTOOLS_PATH pileup $bam_file | awk '{if ( \$4 > $max_cov ){ print \$1\"\\t\"\$2\"\\t\"\$4 } }' |" ) || die "$!\n";
	my $end;
	my ( $ref , $start ) = ( "" , 0 );
	while ( <IN> ){
	    chomp;
	    my @data = split /\t/;
	    if ( $ref eq "" ){
		( $ref , $start , $end ) = ( $data[0] , $data[1] , $data[1] );
	    }
	    elsif ( $ref eq $data[0] ){
		if ( $end + $ext >= $data[1] ) {
		    $end = $data[1];
		}
		else{
		    push ( @{$over_cov{$ref}}, [$start, $end + $read_length] );
#		    print $ref,"\t",$start,"\t",$end,"\n";
		    ( $start , $end ) = ( $data[1] , $data[1] );
		}
	    }
	    else{
		( $ref , $start , $end ) = ( $data[0] , $data[1] , $data[1] );
	    }
	}
	print "Done.\n";
	open ( OV , ">$OUTPUT_DIR/over_cov.txt") || die "Cannot open $OUTPUT_DIR/over_cov.txt, $!\n";
	foreach my $ref ( keys %over_cov ) {
	    foreach my $i ($#{ $over_cov{$ref} }){
		print OV $ref,"\t",$over_cov{$ref}->[$i]->[0],"\t",$over_cov{$ref}->[$i]->[1],"\n";
		
	    }
	}
	close (OV);
    }
    return %over_cov;


}

#############################################################################
#                           sort_BAM_by_name                                        
#----------------------------------------------------------------------------
# Accessory function: sorts the bam file by name using samtools. output in
# global variable NOSORT.

sub sort_bam_by_name{
    my($input, $output)=@_;
    my $cmd=`$SAMTOOLS_PATH sort -n $input $output`;
    die ("samtools exited abnormally: $!\n") if ($? !=0);

}


sub bam2discordant_pairs {
    my ( $bam , $output_dir , $lim_inf , $lim_sup , $lib_type , $extend_overlap ) = @_;
    opendir ( my $open_dir , $output_dir ) || die "Cannot open $output_dir, $!\n";
    my @filelist= grep (/chr.\.[\+\-]\.(del|ins|dup|inv|trans)\.init/, readdir($open_dir));
    if (scalar(@filelist) > 0){
	print "QUERY:\tPrevious run detected. ",scalar(@filelist)," discordant pair files were found. Do you want to reuse them to spare time?(y/N)\n";
	print "\t(Beware, there is no guarantee that the previous run was complete or used the same parameters you did)";
	my $input=<>;
	chomp($input);
	if ($input eq "y" or $input eq "Y"){
	    print "Skipping...\n";
	    $RERUN="YES";
	    return;
	}
    }
    print "INFO:\tDetection of discordant reads\n";
    if ($NOSORT eq "") {
	print "INFO:\tSorting by name... (this step can be lengthy)\n";
	sort_bam_by_name($BAM, "$OUTPUT_DIR/input_sortedby_name");
	$chemin_sorted="$output_dir/input_sortedby_name.bam";
    }
    else{
	$chemin_sorted=$NOSORT;
    }
    my $ref = "";
    my @multiple_hits;
    my $disc_HANDLE;
    my $conc_HANDLE;
    if($CONCORDANT ne "no"){
	open ($conc_HANDLE, "|-", "$SAMTOOLS_PATH view -bS - -o $OUTPUT_DIR/concordant.bam") or die("open: $!");
    }
    if($DISCORDANT ne "no"){
	open ($disc_HANDLE, "|-", "$SAMTOOLS_PATH view -bS - -o $OUTPUT_DIR/discordant.bam") or die("open: $!");
    }
    my $total=`$SAMTOOLS_PATH view -c $chemin_sorted`;
    open ( IN, "$SAMTOOLS_PATH view -h $chemin_sorted |" ) || die "$!\n";
    my $flag="";
    use Time::HiRes qw( gettimeofday tv_interval );
    my $startime=[Time::HiRes::gettimeofday()];
    my $maxsize=0;
    my $count=0;
    my $progress;
    my $remaining;
    my $line="";
    my @data;
    my $memory="";
    while ( <IN>  ){
	chomp;
	$line=$_;
        # Following line to remove duplicates in SAM file
	if($memory eq $line){next;} else {$memory=$line;};
	if($count % 1000 == 0){ 
	    $progress=sprintf('%.2f', ($count/$total)*100);
	    $remaining=$progress !=0 ?((100-$progress)/$progress)*tv_interval($startime):0;
	    $remaining=int($remaining/3600) .':'. int(($remaining-(int($remaining/3600)*3600))/60);
	    print "INFO:\tTreated ",$progress,"% of reads. Estimated time remaining: $remaining.\r";}#, " Size: ", sprintf('%.2f', total_size(\%CURRENT_DISCORDANT_PAIRS)/1000000000) , "Go.\r";}
     
	
	$count++;

	# Following block writes the headers for conc and disc files.
	if($line =~ m/^\@/){
	    print $conc_HANDLE $line, "\n" if $CONCORDANT ne "no";
	    print $disc_HANDLE $line, "\n" if $DISCORDANT ne "no";
	    next;
	}
	@data = split /\t/ , $line;
	next if ( $data[5] eq "*" || $data[6] eq "*" || ( ( $data[3] == $data[7] ) &&  ( $data[6] eq "=" ) ) );#|| $data[4] == 0 ); #on avance jusqu'à la première ligne intéressante suivante
	if ( $ref eq "" ){ #au début on prend la première bonne séquence comme référence et on ajoute tous ses mappings au format perso
	    $ref = $data[0];
	    push ( @{$multiple_hits[0]} ,  get_sam_multiple_records ( \@data ) );
	}
	elsif ( $data[0] eq $ref ){# si une référence existe on ajoute les mappings du read apparié à la suite dans multiple_h
				       push ( @{$multiple_hits[1]} , get_sam_multiple_records ( \@data ) ) ;
				   }
	else{ # une référence existe mais on tombe sur une autre paire valide
	    # on envoie tous les multiple hits déjà trouvés à analyse_one_pair...
	    if(defined($multiple_hits[0]) && defined($multiple_hits[1]) && scalar(@{$multiple_hits[0]}) < $MAXMAP && scalar(@{$multiple_hits[1]})<$MAXMAP) {
		analyse_one_pair_mapping_results ( \@multiple_hits,  $lim_inf , $lim_sup , $lib_type , $output_dir, $extend_overlap, $conc_HANDLE, $disc_HANDLE);
	    }
	    $ref = $data[0];
	    undef(@multiple_hits);
	    push ( @{$multiple_hits[0]} ,  get_sam_multiple_records ( \@data ) ) ;
	}
    }
    analyse_one_pair_mapping_results ( \@multiple_hits,  $lim_inf , $lim_sup , $lib_type , $output_dir, $extend_overlap, $conc_HANDLE, $disc_HANDLE);
    $remaining=tv_interval($startime);
    $remaining=int($remaining/3600) .'h '. int(($remaining-(int($remaining/3600)*3600))/60)."min";
    print "INFO:\tTreated $count reads in $remaining. Writing results to disk...                               \n";
    close($disc_HANDLE) if $DISCORDANT ne "no";
    close($conc_HANDLE) if $CONCORDANT ne "no";
    close (IN);
    print_discordant_pairs (\%CURRENT_DISCORDANT_PAIRS , $output_dir ) if ($NODISCORDANT eq "");
    print "INFO:\tDone.\n";
}


sub DIAGNOSE_printhit{
    my ($pairs)=@_;
    for my $i (0..1){
	print "Mate $i: ", $#{$pairs->[$i]}, " mappings found.";
	print "(Pair ",$pairs->[$i]->[0]->[5], ")\n" if ($#{$pairs->[$i]}>0);
	# for my $j (1..$#{$pairs->[$i]}){
	#     push (@{$clean_matches[$i]}, $pairs->[$i]->[$j]) if($pairs->[$i]->[$j]->[7] <= $MAXMIS);
	#     print "Mapping treated: ", join('.', @{$pairs->[$i]->[$j]}), "\n";
	# }
    }

}

##########################################################################################################
#					get_sam_multiple_records
#---------------------------------------------------------------------------------------------------------
# prend en entrée un tableau de colonnes SAM
# retourne un tableau à un format personnalisé: 
# reference sequence name
# position
# une partie du FLAG (uniqt sens)
# séquence
# read name
# qualité
# ==============
# s'il existe d'autres mappings, pour chaque mapping
# ==============
# le chromosome
# la position
# le sens
# le sens du père
# séquence du père
# read name du père
# qualité du père
# distance
sub get_sam_multiple_records {
    my ( $data ) = @_;
    my @records;
    my $pos = $data->[3];
    my $sens = ( ( $data->[1] & 0x10 )  == 16 ) ? "-" : "+";
    my $cigar = $data->[5];
    push @records , [ $data->[2] , $pos , $sens, $data->[9], $data->[10] , $data->[0] , $data->[4], $data->[1], $data->[5], $data->[6], $data->[7], $data->[8], @{$data}[11..$#{$data}]  ];
    my $line = join "\t" , @$data;

    if ( $line =~/\t(XA\:.+\;)/){
	my @mult = split ( /\;/ , $1);
	my @new_data=();
	for (my $j=0; $j<=$#mult; $j++){
#	foreach my $j ( 0..$#mult ){ 
	    $mult[$j] =~s/.+\://g; #pour chaque autre mapping, on supprime les noms des champs et on parse (chr, pos, cigar, nm)
	    @new_data = split ( /\,/ , $mult[$j] );
	    $sens = $new_data[1]; 
	    $sens =~s/\w*//g;
	    $pos = $new_data[1];
	    $pos =~s/[\+\-]//;
	    $cigar = $new_data[2];
	    push @records ,  [ $new_data[0], $pos, $sens,  $data->[9], $data->[10] , $data->[0] , $data->[4] , $new_data[3], $cigar ];
# Attention le CIGAR n'est jamais utilisé.
	    undef ( @new_data );
	}
	undef ( @mult );
    }
    undef $data;
    return @records;
    undef @records;
}

#############################################################################################
#			analyse_one_pair_mapping_results
#-------------------------------------------------------------------------------------------
# take: an array of array containing one pair mapping features
# load: the selected pair hash of array, keys are pair case (chr_a.sens_a.chr_b.sens_b), values are array of pair features (pair name, pos_a, pos_b)
sub analyse_one_pair_mapping_results {
    my ( $pair_to_analyse , $lim_inf , $lim_sup , $lib , $output_dir , $overlap_limit, $conc_handle, $disc_handle) = @_;
#    print "Records 1 - ", $#{@$pair_to_analyse->[0]}, "\n";
#    print "Records 2 - ", $#{@$pair_to_analyse->[1]}, "\n";
    my $check = is_proper_mapped ( $pair_to_analyse , $lim_inf, $lim_sup , $lib );
#    print "##DIAGNOSE: Size of array: ", sprintf('%.2f', total_size($check)/1000000000) , "Go.\n" if(total_size($check)>100000000);
#    print "Proper mapped  : ", $#{@$check}, "\n";
    
    if ((!defined($check)) || (!defined($check->[0])) || ! @{$check->[0]} ){	# si la paire ne mappe pas (discordante) alors le tableau est vide 
	my %non_overlapping_cases = get_classified_multiple_hits ( $pair_to_analyse, $lim_inf, $lim_sup , $overlap_limit );
#	print "Records 1 - ", join "::", keys(%non_overlapping_cases), "\n";

		if($DISCORDANT ne "no"){
		    print({$disc_handle} join("\t", $pair_to_analyse->[0]->[0]->[5], $pair_to_analyse->[0]->[0]->[7],$pair_to_analyse->[0]->[0]->[0], $pair_to_analyse->[0]->[0]->[1], $pair_to_analyse->[0]->[0]->[6], $pair_to_analyse->[0]->[0]->[8],$pair_to_analyse->[0]->[0]->[9], $pair_to_analyse->[0]->[0]->[10], $pair_to_analyse->[0]->[0]->[11], $pair_to_analyse->[0]->[0]->[3], $pair_to_analyse->[0]->[0]->[4], @{$pair_to_analyse->[0]->[0]}[12..$#{$pair_to_analyse->[0]->[0]}]), "\n");
		    
		    print ({$disc_handle} join("\t", $pair_to_analyse->[1]->[0]->[5], $pair_to_analyse->[1]->[0]->[7], $pair_to_analyse->[1]->[0]->[0], $pair_to_analyse->[1]->[0]->[1], $pair_to_analyse->[1]->[0]->[6], $pair_to_analyse->[1]->[0]->[8], $pair_to_analyse->[1]->[0]->[9], $pair_to_analyse->[1]->[0]->[10], $pair_to_analyse->[1]->[0]->[11], $pair_to_analyse->[1]->[0]->[3], $pair_to_analyse->[1]->[0]->[4], @{$pair_to_analyse->[1]->[0]}[12..$#{$pair_to_analyse->[1]->[0]}]), "\n");# if defined($pair_to_analyse->[1]);
		}

	foreach my $case ( keys %non_overlapping_cases ){
	    foreach my $i ( 0 .. $#{$non_overlapping_cases{$case} } ){
		push @{$CURRENT_DISCORDANT_PAIRS{$case}} , [ $non_overlapping_cases{$case}->[$i]->[0] , $non_overlapping_cases{$case}->[$i]->[1] , $non_overlapping_cases{$case}->[$i]->[2] ];
		
	    }
	}
	undef %non_overlapping_cases;	
    }
    elsif (@{$check->[0]} and scalar(@{$check->[0]})==1 && $CONCORDANT ne "no") # sinon ce ne sont que les MP matches uniques qui nous intéressent
    {
	my @splitted=split('\t', $check->[0]->[0]);
	my $i=$splitted[1];
	my $j=$splitted[2];
	
	if($splitted[1]!=0){ #Si le mapping correct n'est pas le mapping principal, remplacer par les données du XA
	    my $flag= $pair_to_analyse->[0]->[0]->[7];
	    my $length=max($pair_to_analyse->[1]->[$j]->[1], $pair_to_analyse->[0]->[$i]->[1])-min($pair_to_analyse->[1]->[$j]->[1], $pair_to_analyse->[0]->[$i]->[1]); #unsure about this calculation, supposed to be #rightmost_base - #leftmost_base, here max(pos1, 2)-min(pos1, 2)
	    if($pair_to_analyse->[0]->[0]->[1] ne $pair_to_analyse->[0]->[$i]->[1]){
		$flag=($pair_to_analyse->[0]->[$i]->[1] eq "-")? ($flag | 0x10) : ($flag ^ 0x10);
	    }
	    if($NOLUTHER eq ""){
		print({$conc_handle} join("\t", $pair_to_analyse->[0]->[0]->[5], 
					  $flag,
					  $pair_to_analyse->[0]->[$i]->[0], 
					  $pair_to_analyse->[0]->[$i]->[1], 
					  255, 
					  $pair_to_analyse->[0]->[0]->[8],
					  '=', 
					  $pair_to_analyse->[0]->[0]->[10], 
					  $length, 
					  $pair_to_analyse->[0]->[0]->[3], 
					  $pair_to_analyse->[0]->[0]->[4], 
					  "NM:i:$pair_to_analyse->[0]->[$i]->[7]\n"));
	    }	      
	}else {
	    print({$conc_handle} join("\t", $pair_to_analyse->[0]->[0]->[5], $pair_to_analyse->[0]->[0]->[7],$pair_to_analyse->[0]->[0]->[0], $pair_to_analyse->[0]->[0]->[1], $pair_to_analyse->[0]->[0]->[6], $pair_to_analyse->[0]->[0]->[8],$pair_to_analyse->[0]->[0]->[9], $pair_to_analyse->[0]->[0]->[10], $pair_to_analyse->[0]->[0]->[11], $pair_to_analyse->[0]->[0]->[3], $pair_to_analyse->[0]->[0]->[4], @{$pair_to_analyse->[0]->[0]}[12..$#{$pair_to_analyse->[0]->[0]}]), "\n");
	    
	}
	
	if($splitted[2]!=0){#Idem pour le read '2'
				my $flag= $pair_to_analyse->[1]->[0]->[7];
				my $length=max($pair_to_analyse->[1]->[$j]->[1], $pair_to_analyse->[0]->[$i]->[1])-min($pair_to_analyse->[1]->[$j]->[1], $pair_to_analyse->[0]->[$i]->[1]); #unsure about this calculation, supposed to be #rightmost_base - #leftmost_base, here max(pos1, 2)-min(pos1, 2)
				if($pair_to_analyse->[1]->[0]->[1] ne $pair_to_analyse->[1]->[$j]->[1]){
				    $flag=($pair_to_analyse->[1]->[$j]->[1] eq "-")? ($flag | 0x10) : ($flag ^ 0x10);
				}
				if($NOLUTHER eq ""){
				    print({$conc_handle} join("\t", $pair_to_analyse->[0]->[0]->[5], 
							      $flag,
							      $pair_to_analyse->[0]->[$i]->[0], 
							      $pair_to_analyse->[0]->[$i]->[1], 
							      255, 
							      $pair_to_analyse->[0]->[0]->[8],
							      '=', 
							      $pair_to_analyse->[0]->[0]->[10], 
							      $length, 
							      $pair_to_analyse->[0]->[0]->[3], 
							      $pair_to_analyse->[0]->[0]->[4], 
							      "NM:i:$pair_to_analyse->[0]->[$i]->[7]\n"));
				}	      
			    }
	else{
	    
	    print ({$conc_handle} join("\t", $pair_to_analyse->[1]->[0]->[5], $pair_to_analyse->[1]->[0]->[7], $pair_to_analyse->[1]->[0]->[0], $pair_to_analyse->[1]->[0]->[1], $pair_to_analyse->[1]->[0]->[6], $pair_to_analyse->[1]->[0]->[8], $pair_to_analyse->[1]->[0]->[9], $pair_to_analyse->[1]->[0]->[10], $pair_to_analyse->[1]->[0]->[11], $pair_to_analyse->[1]->[0]->[3], $pair_to_analyse->[1]->[0]->[4], @{$pair_to_analyse->[1]->[0]}[12..$#{$pair_to_analyse->[1]->[0]}]), "\n");
	}
    }
    undef $pair_to_analyse;
    undef @$check;
    
    
}



###############################################################################################
#				is_proper_mapped
#-------------------------------------------------------------------------------------------
# take: 1. an array of array containing one pair mapping results (each hit is an array element containing a tab split SAM lines ),
#	2. the lower fragment size limit (bp), 
#	3. the upper fragment size limit (bp), 
#	4. the type of library (MP or PE) 
# return: "" si la paire est discordante, "MP suivi d'un tab-split string:
#
sub is_proper_mapped {
    my ( $pairs, $lim_inf, $lim_sup, $lib ) = @_;
    my @alignments=();
    # nouvel algorithme: retourner uniquement les paires qui possèdent un mapping concordant UNIQUE en MP.
# retourner ici un array au lieu d'une ligne
#plus haut: si l'array MP contient plus d'une ligne, discard
# s'il ne contient qu'une ligne et que cette ligne correspond à l'alignement principal pour 1 : écrire 1 sans XA
# s'il ne contient qu'une ligne et que cette ligne correspond à l'alignement principal pour 2 ! écrire 2 sans XA
# s'il ne contient qu'une ligne et que l'alignement est dans XA pour l'un des deux: récrire l'alignement principal avec les infos de XA et écrire.
    for my $i ( 0 .. $#{$pairs->[0]}){
	for my $j ( 0 .. $#{$pairs->[1]} ) {
	    #test MP - looking for proper mapped pair (MP reverse forward, and size > lim inf , size < lim sup
	    if ( $LIB eq "MP" && $pairs->[0]->[$i]->[0] eq $pairs->[1]->[$j]->[0] && ( ($pairs->[0]->[$i]->[1] > $pairs->[1]->[$j]->[1] && $pairs->[0]->[$i]->[2] eq "+" && $pairs->[1]->[$j]->[2] eq "-" && $pairs->[0]->[$i]->[1] - $pairs->[1]->[$j]->[1] > $lim_inf &&  $pairs->[0]->[$i]->[1] - $pairs->[1]->[$j]->[1] < $lim_sup ) || ($pairs->[0]->[$i]->[1] < $pairs->[1]->[$j]->[1] && $pairs->[0]->[$i]->[2] eq "-" && $pairs->[1]->[$j]->[2] eq "+" && $pairs->[1]->[$j]->[1] - $pairs->[0]->[$i]->[1] >$lim_inf &&  $pairs->[1]->[$j]->[1] - $pairs->[0]->[$i]->[1] < $lim_sup) )  ){
		push(@{$alignments[0]},"MP\t".$i."\t".$j."\t".join ( "\t" , @{$pairs->[0]->[$i]} )."\t".join ("\t", @{$pairs->[1]->[$j]})."\n");	
	    }
	    #test PE - looking for proper mapped pair (PE forward reverse, and size > lim inf , size < lim sup
	    elsif ( $LIB eq "PE" && $pairs->[0]->[$i]->[0] eq $pairs->[1]->[$j]->[0] && ( ($pairs->[0]->[$i]->[1] > $pairs->[1]->[$j]->[1] && $pairs->[0]->[$i]->[2] eq "-" && $pairs->[1]->[$j]->[2] eq "+" && $pairs->[0]->[$i]->[1] - $pairs->[1]->[$j]->[1] > 100 &&  $pairs->[0]->[$i]->[1] - $pairs->[1]->[$j]->[1] < 700 ) || ($pairs->[0]->[$i]->[1] < $pairs->[1]->[$j]->[1] && $pairs->[0]->[$i]->[2] eq "+" && $pairs->[1]->[$j]->[2] eq "-" && $pairs->[1]->[$j]->[1] - $pairs->[0]->[$i]->[1] >100 &&  $pairs->[1]->[$j]->[1] - $pairs->[0]->[$i]->[1] < 700) )  ){
		push(@{$alignments[1]}, "PE\t".join ( "\t" , @{$pairs->[0]->[$i]} )."\t".join ("\t", @{$pairs->[1]->[$j]})."\n");
	    }
	}
    }
    
    return \@alignments;
}
##############################################################################################
#				get_classified_multiple_hits
#------------------------------------------------------------------------------------------
# take: an array of array containing one pair mapping results ( each hit is an array element, each array element is an array containing a tab split SAM lines 
# return: an hash of non overlapped cases of pairs 
sub get_classified_multiple_hits {
    my ( $pairs, $lim_inf, $lim_sup , $overlap_limit ) = @_;
    my %pair_cases=();
    my $read_name = $pairs->[0]->[0]->[5];
    my @clean_matches;
# #    remove all mappings with Edit Distance > $MAXMIS
#     for my $i (0..1){
#     	my $xm=(split("XM:i:", join("\t", @{$pairs->[$i]->[0]})))[1];
#     	$xm=(split("\t", $xm))[0];
#     	push (@{$clean_matches[$i]}, $pairs->[$i]->[0]) if($xm <= $MAXMIS);
#     	for my $j (1..$#{$pairs->[$i]}){
#     	    push (@{$clean_matches[$i]}, $pairs->[$i]->[$j]) if($pairs->[$i]->[$j]->[7] <= $MAXMIS);
#     	}
#     }
#     if(!defined($clean_matches[1]) || !defined($clean_matches[0]) || scalar(@{$clean_matches[0]})==0 || scalar(@{$clean_matches[1]})==0){
#     	return %pair_cases;
#     }
#     $pairs=\@clean_matches;


    # # Reduce the number of mappings by clustering nearby ones
    # my @isolat_array;
    # my $isolats=\@isolat_array;
    # my $REPEAT_THRESHOLD=20000;
    # for my $mate (0..1){
    # 	next if(!defined($pairs->[$mate]));
    # 	my @pair1_sortedchrpos=sort {$a->[0] cmp $b->[0] || $a->[1] <=> $b->[1]} @{$pairs->[$mate]};
    # 	my ($firstplus, $firstminus, $lastplus, $lastminus, $chr);
    # 	$chr="";
    # 	foreach my $alignment (@pair1_sortedchrpos){
    # 	     $chr=$alignment->[0] if($chr eq "");
    # 	    if($chr ne $alignment->[0]){
    # 		push( @{$isolat_array[$mate]}, [ @{$firstplus} ]) if(defined($firstplus));
    # 		push( @{$isolat_array[$mate]}, [ @{$firstminus} ]) if(defined($firstminus));
    # 		undef($firstplus);undef($firstminus);
    # 		if($alignment->[1] eq "+"){$firstplus=$alignment;$lastplus=$alignment;}
    # 		else{$firstminus=$alignment; $lastminus=$alignment;};
    # 		$chr=$alignment->[0];
    # 	    }
    # 	    else{ # même chromosome
    # 		if($alignment->[2] eq "+"){
    # 		    #plus case
    # 		    if(!defined($firstplus)){$firstplus=$alignment; $lastplus=$alignment;next;};
    # 		    if($alignment->[1] > $lastplus->[1]+$REPEAT_THRESHOLD){
    # 			push( @{$isolat_array[$mate]}, [ @{$firstplus} ]);
    # 			$firstplus=$alignment;$lastplus=$alignment;
    # 		    }else{
    # 			$lastplus=$alignment;
    # 		    }
    # 		}
    # 		else{
    # 		    #minus case
    # 		    if(!defined($firstminus)){$firstminus=$alignment; $lastminus=$alignment;next;};
    # 		    if($alignment->[1] > $lastminus->[1]+$REPEAT_THRESHOLD){
    # 			push( @{$isolat_array[$mate]}, [ @{$firstminus} ]);
    # 			$firstminus=$alignment;$lastminus=$alignment;
    # 		    }else{
    # 			$lastminus=$alignment;
    # 		    }
		    
    # 		}
		
    # 	    }
    # 	}
    # 	push( @{$isolat_array[$mate]}, [ @{$firstplus} ]) if(defined($firstplus));
    # 	push( @{$isolat_array[$mate]}, [ @{$firstminus} ]) if(defined($firstminus));

	
    # }
    # $pairs=\@isolat_array;



    for my $i ( 0 .. $#{$pairs->[0]} ){
	for my $j ( 0 .. $#{$pairs->[1]} ) {
	    my ( $chr_a , $chr_b , $pos_a , $pos_b , $sens_a , $sens_b ) = ( $pairs->[0]->[$i]->[0] , $pairs->[1]->[$j]->[0] , $pairs->[0]->[$i]->[1] , $pairs->[1]->[$j]->[1],  $pairs->[0]->[$i]->[2] , $pairs->[1]->[$j]->[2] );
	    if ( $chr_a ne $chr_b ){
		if ( exists $pair_cases{"$chr_b.$sens_b.$chr_a.$sens_a.trans"} ){
		    push @{$pair_cases{"$chr_b.$sens_b.$chr_a.$sens_a.trans"}}, [ $pos_b, $pos_a ] ;
		}
		else{
		    push @{$pair_cases{"$chr_a.$sens_a.$chr_b.$sens_b.trans"}}, [ $pos_a, $pos_b ] ;
		}
	    }
	    else{
		if ( $pos_a > $pos_b ){
		    my ( $chr_c , $sens_c , $pos_c ) = ( $chr_a , $sens_a , $pos_a );
		    ( $chr_a , $sens_a , $pos_a  ) = ( $chr_b , $sens_b , $pos_b );
		    ( $chr_b , $sens_b , $pos_b ) = ( $chr_c , $sens_c , $pos_c );
		    undef($chr_c);undef($sens_c);undef($pos_c);
		}
		if ( $pos_b-$pos_a >= $lim_sup && $sens_a eq "-" && $sens_b eq "+" && $pos_a != $pos_b){
		    push @{$pair_cases{"$chr_a.$sens_a.$chr_b.$sens_b.del"}}, [ $pos_a, $pos_b ] ;
		}
		elsif ( $pos_b-$pos_a <= $lim_sup && $sens_a eq "-" && $sens_b eq "+" && $pos_a != $pos_b ){
		    push @{$pair_cases{"$chr_a.$sens_a.$chr_b.$sens_b.ins"}}, [ $pos_a, $pos_b ] ;
		}
		elsif ( $sens_a eq $sens_b && $pos_a != $pos_b ){
		    push @{$pair_cases{"$chr_a.$sens_a.$chr_b.$sens_b.inv"}}, [ $pos_a, $pos_b ] ;
		}
		elsif ( $sens_a eq "+" && $sens_b eq "-" && $pos_a != $pos_b ){
		    push @{$pair_cases{"$chr_a.$sens_a.$chr_b.$sens_b.dup"}}, [ $pos_a, $pos_b ] ;
		}
	    }
	}
    }
    

    foreach my $key (keys %pair_cases){
	}

    # Get non overlapping pairs
    my %non_overlapped_pair;
    foreach my $case ( keys %pair_cases ){
	my @good_pairs;
	my @sorted_pair_case = sort { $a->[0] <=> $b->[0] } @{$pair_cases{$case}};
	if ( $#{$pair_cases{$case}} > 0 ){
	    @good_pairs = get_pair_cluster ( \@sorted_pair_case, $case , $overlap_limit );
	    foreach my $i ( 0 .. $#good_pairs ){
		push @{$non_overlapped_pair{$case}} , [ $good_pairs[$i]->[0], $good_pairs[$i]->[1] , $read_name ] ;
	    }
	}
	else {
	    push @{$non_overlapped_pair{$case}} , [ $pair_cases{$case}[0]->[0], $pair_cases{$case}[0]->[1] , $read_name ] ; 
	}
	undef ( @good_pairs );
	undef ( @sorted_pair_case);
    }
    undef ( %pair_cases );
    undef ( $pairs );

    
    if($FILTER ne "disabled"){
	# Get likely discordant pairs for low complexity regions
	my $good_case = "";
	my $flag = 0; 
	foreach my $case ( keys %non_overlapped_pair  ) {
	    if ( $case =~/del|inv|dup|ins/ ){
    		foreach my $i ( 0 .. $#{$non_overlapped_pair{$case} } ){
		    if ( ( $non_overlapped_pair{$case}[$i]->[1] - $non_overlapped_pair{$case}[$i]->[0] ) <= $FILTER ){
			$good_case = $case ;
			$flag = 1;
			last;					
		    }
    		}
	    }
	    last if ( $flag == 1 );
	}
	
	# remove artefactual cases
	if  ( $good_case ne "" ){
	    foreach my $case  ( keys %non_overlapped_pair ){
    		delete $non_overlapped_pair{$case} if ( $case ne $good_case && $case =~/del|inv|dup|ins/);
	    }
	}	
    }


    return %non_overlapped_pair;
}

sub sort_for_pairs ($$){
    my ($premier, $deuxieme)=@_;
    if($premier->[0] ne $deuxieme->[0]) {return $premier->[0] cmp $deuxieme->[0];}
    else{return $premier->[1] <=> $deuxieme->[1];}
}

############################################################################################################
#					get_pair_cluster
#-----------------------------------------------------------------------------------------------------------
# take:  array of [pos_a, pos_b], corresponding to a signature
# return:
sub get_pair_cluster{
    my ( $pair , $case , $overlap_limit ) = @_;
    my @clust;
    my $n_clust = 0;
    # initialize cluster
    push @{$clust[$n_clust]} , $pair->[0]; #begin by taking the first (posa, posb)
    for my $k ( 1 .. $#{$pair} ){
	my @overlap = is_overlapping ($pair->[$k]->[0], $pair->[$k]->[0] + $overlap_limit, $pair->[$k-1]->[0], $pair->[$k-1]->[0] + $overlap_limit, 0 );
	if ( $overlap[0] > 0 ){
	    push @{$clust[$n_clust]} , $pair->[$k];
	}
	else{
	    ++$n_clust;
	    push @{$clust[$n_clust]} , $pair->[$k];
	}
    }
    # take a random pair from overlapping region
    my @non_overlapped_pair;
    for my $k ( 0 .. $#clust ){
	my $rand = int ( rand ( $#{$clust[$k]} ) );
	push ( @non_overlapped_pair , [ $clust[$k]->[$rand]->[0], $clust[$k]->[$rand]->[1] ] );
    }
    undef @clust;
    return @non_overlapped_pair;
}

###############################################################################
#			print_discordant_pairs
#------------------------------------------------------------------------------
# take: a hash of discordant pair
# return: print the result in files
sub print_discordant_pairs{
    my ( $pair, $output_dir ) = @_;
    foreach my $case ( keys %$pair ){
	my $A = 0;
	my $B = 1;
	my $different_chr = $case;
	$different_chr =~/[Cc]hr(\w+)\.[\-\+]\.[Cc]hr(\w+)\.[\-\+]\..+/;
	#print $different_chr;
	my $chr1 = $1;
	my $chr2 = $2; 
	my @discordant_pairs = sort { $a->[0] <=> $b->[0] } @{$pair->{$case}};
	#print "chr1 $chr1 chr2 $chr2\n";
	if ( $chr1 ne $chr2 ){
	    if ( $chr1 =~/\d/ && $chr2=~/\d/ && $chr1 > $chr2 ){
		$A = 1;
		$B = 0;				
		$case =~s/[Cc]hr(\w+)\.([\-\+])\.[Cc]hr(\w+)\.([\-\+])\.(.+)/[Cc]hr$3.$4.[Cc]hr$1.$2.$5/;				
	    }
	    if ( $chr1=~/\D/ || $chr2=~/\D/ ){
		my @order = ( $chr1 , $chr2 );
		@order = sort { $a cmp $b } @order;
		if ( $order[0] ne $chr1 ){
		    $A = 1;
		    $B = 0;
		    $case =~s/[Cc]hr(\w+)\.([\-\+])\.[Cc]hr(\w+)\.([\-\+])\.(.+)/[Cc]hr$3.$4.[Cc]hr$1.$2.$5/;					
		}
	    }
	}
	open ( OUT ,">>$output_dir/$case.init" );
	foreach my $i ( 0 .. $#discordant_pairs ){
	    print OUT $discordant_pairs[$i]->[$A] , "\t" , $discordant_pairs[$i]->[$B] , "\t", $discordant_pairs[$i]->[2] , "\n" ;
	}
	close ( OUT );
    }
    # sort and uniq	
    opendir ( my $open_dir , $output_dir ) || die "Cannot open $output_dir, $!\n";
    my @files = readdir ( $open_dir );
    foreach my $file ( @files ){
	if ( $file=~/.\.init$/ ){
	    `sort -k1,1n $output_dir/$file | awk '{if ( vu[\$1"@"\$2]!=1 ){print \$0;vu[\$1"@"\$2]=1}}' > $output_dir/$file.sorted;rm $output_dir/$file; mv $output_dir/$file.sorted $output_dir/$file`; 
	}
    }
    closedir ($open_dir);
}


#############################################################################
#                           treat_BAM                                        
#----------------------------------------------------------------------------
# Accessory function: redefines the bam file so that it contains only perfect
# matches. Lone mates are also removed. Works only on ordered-by-name files.


sub treat_bam{
    my ( $input, $output, $regx ) = @_;
    # preprocesses the bam file to remove all reads that do not map perfectly at least once
    print "INFO:\tTreating file $input with mismatch parameters $regx.\n";
    my $m=$READ_SIZE."M";
    my $newfilename="$output";
    open (OUT, "| $SAMTOOLS_PATH view -Sb - -o $newfilename");
    my %reads;
    my $first="";
    open ( IN, "$SAMTOOLS_PATH view -h $input |" ) || die "Cannot open $input, $!\n";
    my $ct_reads=0;
    while ( <IN> ){
	chomp;
	if ($_ =~ m/^\@/ ) {
	    print (OUT $_, "\n");
	    next;
	}
	my $data=$_;
	my @tab=split('\t', $data);
	if($data =~ m/XT:A:M/){
	    my @multiples=split('XA:Z:', $data);
	    if(scalar(@multiples)>1 && $multiples[1] =~ m/$m,[$regx]/){
		#Sometimes in case of discordant reads the main mapping is quite bad (rescued mate)
		#In this case the secondary mappings contain worthy information
		if($first ne "" && (split('\t', $first))[0] eq $tab[0]){
			$ct_reads++;
		    print (OUT (split('\n', $first))[0], "\n");
		    print (OUT (split('\n', $data))[0], "\n");
		    $first="";
		}else{
		    $first=$data;
		}
		
	    }
	}
	elsif ($tab[5] eq $m && $data =~ m/XM:i:[$regx]/ ) {
	    if($first ne "" && (split('\t', $first))[0] eq $tab[0]){
	    	$ct_reads++;
		print( OUT (split('\n', $first))[0], "\n");
		print(OUT (split('\n', $data))[0], "\n");
		$first="";
	    }else{
		$first=$data;
	    }
	}
    }
    print("\n\n---------- /!\\ -----------\nWARNING:\t No reads remaining after filtering. Are you sure you specified the right read length?\n---------- /!\\ -----------\n\n") if($ct_reads==0);
    close(IN);
    close(OUT);
    print STDERR "Done.\n";
#    `java -jar $PICARD/SamFormatConverter.jar INPUT=$OUTPUT_DIR/perfectmatches.sam OUTPUT=$newfilename 1>&- 2>&`;
#    `rm $OUTPUT_DIR/perfectmatches.sam`;
}

###################################################################################
#			is_overlapping						  #
#----------------------------------------------------------------------------------
# take: an array of monotonously increasing integers (a b c d), see if interval [a:b] and [c:d] are overlapping
# return: an array (intersection_size, union_size, union_min, union_max)
sub is_overlapping{
    use List::Util qw(max min);
    my ( $start_a , $end_a , $start_b , $end_b , $reads_l ) = @_;
    my $overlap_start = ($start_a > $start_b) ? $start_a : $start_b;
    my $overlap_end = ($end_a < $end_b) ? $end_a + $reads_l : $end_b + $reads_l;
    my $overlap_size = ($overlap_start > $overlap_end) ? 0 : $overlap_end - $overlap_start;	
    my $min = min ( $start_a , $end_a + $reads_l , $start_b, $end_b + $reads_l ) ;
    my $max = max ( $start_a , $end_a + $reads_l , $start_b, $end_b + $reads_l ) ; 
    return ( $overlap_size , $max-$min , $min , $max );
}

sub usage{
    print STDERR "\nERROR:\tMandatory parameters were not supplied. Be sure you specify at least -bam and -out.\n";
    print "\t\t (arguments supplied were \'", join "\n\t\t\t\t\t   ", @args, "\'.)\n\n" ;
    print STDERR "Usage :\ttetracker eris -bam [bamfile] -out [out_dir] -lib [PE|MP] -ext [int] -maxcov [int]\n";
    print STDERR "\t\t\t-covext [int] -read [read_size] -samtools [samtools_path] -dm [int]\n";
    print STDERR "\t\t\t-sm [int] -sample [int] -nosort [sorted_bamfile] -concordant [|sorted]\n";
    print STDERR "\t\t\t-discordant [|sorted] -treat_bam [input|concordant|discordant:i-j]\n";
    print STDERR "\t\t\t-mock -rundata [rundata_file] -maxmap [int] -maxmis [int] -picard_path [path]\n";
    print STDERR "\t\t\t-nodiscordant -chainload [treated|concordant|discordant] -f1 [int]\n\n";
    print STDERR "\t Please run \'tetracker eris -man\' for a detailed explanation of command-line options.\n\n";
    exit();
}

sub man{
    system("man $BASEDIR/man/eris 2>/dev/null");
    die("ERROR:\t Impossible to find MAN file in $BASEDIR/man. Please reinstall TE-Tracker.\n") if (($? >>8) !=0);

}
