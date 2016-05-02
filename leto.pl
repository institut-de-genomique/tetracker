#!/usr/bin/perl -w
# svfinder
# MA Madoui & JM Aury, Genoscope, 2011

use strict;
use Getopt::Long;
use POSIX qw(floor);
use List::Util qw[min max];
use IO::Handle;
#use diagnostics;

my $VERSION=0.1;

# PARAMETERS
our @RERUNS;
my $MSIZE=0;
my $HIGHMEM="";
my $HLIMIT=20;
my $RUNDATA="";
my $METHOD="distribution";
our ( $INPUT_DIR , $OUTPUT_DIR , $LIB , $EXT , $MAX_COV , $COV_EXTENSION , $READ_SIZE , $SAMTOOLS_PATH , $DOC_MAD , $SPAN_MAD , $SAMPLE_SIZE , $SLC , $MIN_READS, $NOSORT, $NOLUTHER, $MINCOV, $CLSTPRM, $MODE, $ANNOTATE, $PICARD, $TRIAL) = ( "", "" , "MP" , 20000 , 1000, 500, 76, "/env/cns/src/samtools/samtools-0.1.8/bin/samtools" , 4 , 4 , 1000000 , "/env/cns/proj/projet_AKL/slc/slclust", 10, "");
$CLSTPRM="";
my $MAN="";
my @args=@ARGV;
my $result = GetOptions ( 	'in=s'		=>\$INPUT_DIR,
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
				'slc=s'		=>\$SLC,
				'minreads=i'	=>\$MIN_READS,
				'clustering_parameters=s' =>\$CLSTPRM,
				'a_size=i'      =>\$MSIZE,
				'rundata=s'     =>\$RUNDATA,
				'method=s'      =>\$METHOD,
				'highmem_cmd=s' =>\$HIGHMEM,
				'highmem_lim=i' =>\$HLIMIT,
				'man'         =>\$MAN);


my $RATIO;
$| = 1;
man() if ($MAN ne "");
usage () if ( $INPUT_DIR  eq "" || $OUTPUT_DIR eq "" || $LIB !~/^MP$|^PE$/);


sub usage{
    print STDERR "\nERROR:\tMandatory parameters were not supplied. Be sure you specify at least -in and -out.\n";
    print "\t\t (arguments supplied were \'", join "\n\t\t\t\t\t   ", @args, "\'.)\n\n" ;
    print STDERR "Usage :\ttetracker leto -in [in_dir] -out [out_dir] -lib [PE|MP] -ext [int] -maxcov [int]\n";
    print STDERR "\t\t\t-covext [int] -read [read_size] -samtools [samtools_path] -dm [int]\n";
    print STDERR "\t\t\t-sm [int] -sample [int] -slc [SLC_binary] -minreads [int]\n";
    print STDERR "\t\t\t-clustering_parameters [X:Y] -a_size [int] -rundata [rundata_file]\n";
print STDERR "\t\t\t-method [markov|distribution] -highmem_cmd [cmd] -highmem_lim [int]\n\n";
    print STDERR "\t Please run \'tetracker leto -man\' for a detailed explanation of command-line options.\n\n";
    exit();
}

sub man{
    my $BASEDIR=`pwd`;
    chomp $BASEDIR;
    system("man $BASEDIR/man/letoman 2>/dev/null");
    die("ERROR:\t Impossible to find MAN file in $BASEDIR/man. Please reinstall TE-Tracker.\n") if (($? >>8) !=0);

}



print "\nSVfinder Leto module v. $VERSION - Detection and calling of transposition events\n\n";
#print "\n\"The dark-gowned Leto (the hidden one), always mild, kind to men and to the deathless gods, mild from the beginning, gentlest in all Olympos.\"\n";
#printf("%80s", 'Hesiod, Theogony, 404');
print "\n\n";
print("\n\nInput dir :\t$INPUT_DIR\n");
print("Output dir:\t$OUTPUT_DIR\n");
my @DOC ;
my @SPAN ;

if($RUNDATA eq ""){
open(IN, "$INPUT_DIR/rundata.dat") or die "No rundata.dat file found in input directory. Please run eris.pl again.\n";
}else{
    if(-e $RUNDATA){
	open(IN, $RUNDATA) or die ("ERROR:\t the rundata file $RUNDATA cannot be opened.\n");
    }else{
	die ("ERROR:\t the rundata file $RUNDATA does not exist.\n");
    }
}
my $line=<IN>;
chomp($line);
@DOC=split(":", $line);
print "DOC data:\t$line\n";
$line=<IN>;
@SPAN=split(":", $line);
close(IN);
print "Insert size data:\t$line\n";

if($MSIZE == 0){
    $MSIZE=$SPAN[2];
}
else{
    print "\nINFO:\tWill validate acceptors with size up to $MSIZE.\n";
}
my $X;
my $Y;
my $MIN;
my $READS_LENGTH;
my $ISMEDIAN;
my $ISMAD;
my $DOCMEDIAN;
my $DOCMAD;
our $RERUN="NO";

print "\n";

if ($CLSTPRM ne "") {
    my @params=split (":", $CLSTPRM);
    if (scalar(@params) != 2){
	print "ERROR: If you specify clustering parameters, specify them well, like so : [X:Y].\n";
	exit;
    } else{
	$X=$params[0];
	$Y=$params[1];
	$MIN=$MIN_READS;
	$READS_LENGTH=$READ_SIZE;
	print "-Forcing manual parameters for clustering: X=$X, Y=$Y, MIN=$MIN, READLENGTH=$READS_LENGTH.\n";
    }
}

if($METHOD ne "distr" && $METHOD=~/^markov\,/){
print("INFO:\t Markov-chain parameter evaluation: scheduled.\n");
    # perform Markov-chain evaluation on aligment file
    my $BAM=(split(/\,/, $METHOD))[1];
    die("ERROR:\tfile $BAM could not be found.\n") if( ! -e $BAM);
    $X=get_MC_parameter($BAM);
    $Y=2*$DOC_MAD*$SPAN[1];
print("INFO:\t X set to $X, Y set to $Y following MC calculation.\n");
    $MIN=$MIN_READS;
    $READS_LENGTH=$READ_SIZE;
}




############################################################################
`mkdir -p $OUTPUT_DIR` if ( $OUTPUT_DIR && ! -d $OUTPUT_DIR );
############################################################################


print "\n";



my $root_dir=$OUTPUT_DIR;
#}
my(@TRANS, @DEL);
print(join('.', @TRANS));
my @CONCORDANTID;
`cp $INPUT_DIR/over_cov.txt $OUTPUT_DIR/over_cov.txt` if(!( -e "$OUTPUT_DIR/over_cov.txt"));
my %OVER_COV = get_high_coverage_positions ( "" , $MAX_COV , $COV_EXTENSION , $READ_SIZE );

our %CURRENT_DISCORDANT_PAIRS;


    if (defined($X)){
	print ("Manual parameters set. Generating clusters.\n");
	cluster_disc_pairs ($INPUT_DIR, $X, $Y, $MIN, $READS_LENGTH);
    }
    else{
	print ("Automatic parameters available. Generating clusters. X=",4*$READ_SIZE/$DOC[2],", Y=",2*$DOC_MAD*$SPAN[1] ,".\n");
	cluster_disc_pairs ( $INPUT_DIR , 4*$READ_SIZE/$DOC[2] , 2*$DOC_MAD*$SPAN[1], $MIN_READS , $READ_SIZE);
    }
    sv_calling ($OUTPUT_DIR ,  $DOC[2] , $READ_SIZE, 2 * ( $DOC[1]*$DOC_MAD+$SPAN[0]) , $MSIZE );


#########################################################################################################
#					cluster_disc_pairs
#--------------------------------------------------------------------------------------------------------
#
#
sub cluster_disc_pairs {
    my ( $dir , $x_lim , $y_lim , $min_pair , $read_length ) = @_;
    opendir ( my $open_dir , $INPUT_DIR ) || die "Cannot open $dir, $!\n";
    my @filelist= grep (/chr.\.[\+\-]\.(del|ins|dup|inv|trans)\.init\.clustered/, readdir($open_dir));
    $RERUN="YES";
    if (scalar(@filelist) > 0 && $RERUN eq "NO"){
	print "Previous run detected. ",scalar(@filelist)," discordant pair files were found. Do you want to reuse them to spare time?(y/N)\n";
	print "(Beware, there is no guarantee that the previous run was complete or used the same parameters you did)";
	flush STDIN;
	my $input = <STDIN>;
	chomp($input);
	if ($input eq "y" or $input eq "Y"){
	    print "Skipping...\n";
	    $RERUN="YES";
	}
    }
	


    opendir ( $open_dir , $dir );
    my @files = readdir ( $open_dir );
    print "$#files files were found.\n";
    foreach my $file ( @files ){
	if ( $file =~/chr.*\.[\+\-]\.(del|ins|dup|inv|trans)\.init$/ && !($file =~/(.+)chr(M|C)(.+)/)){
	    if (-e "$OUTPUT_DIR/$file.clustered" && $RERUN eq "YES"){
		print ("Reusing $dir/$file.clustered...\n");
	    }else{
		print "Clustering of $file...\n";
		if ((-s "$dir/$file") > $HLIMIT*1000000){
		    print "Big cluster file found.\n";
print "$HIGHMEM $SLC -i $dir/$file -o $OUTPUT_DIR/$file.clustered -l $x_lim -n $y_lim -m $min_pair -t 0 -r $read_length\n";
			`$HIGHMEM $SLC -i $dir/$file -o $OUTPUT_DIR/$file.clustered -l $x_lim -n $y_lim -m $min_pair -t 0 -r $read_length`;
		}
		else{
		`$SLC -i $dir/$file -o $OUTPUT_DIR/$file.clustered -l $x_lim -n $y_lim -m $min_pair -t 0 -r $read_length`;
		}
	    }
	    if(-e "$OUTPUT_DIR/$file.clustered.merged" && $RERUN eq "YES" ){
		print ("Reusing $dir/$file.clustered.merged...\n");

	    }
	    else{
		print "Merging of $file...\n";
		merge_clusters ( "$OUTPUT_DIR/$file.clustered" , $read_length , $SPAN[0] , $y_lim  , $DOC[2] );
	    }
	    #`rm $dir/$file.clustered` if ( $file =~/chr/ && -s "$dir/$file.clustered" == 0);
	}
    }
    close ( $open_dir );
}
###########################################################################################################
#					merge_clusters
#----------------------------------------------------------------------------------------------------------
#
#
sub merge_clusters { 
	use List::Util qw(max min);	
	my ( $file, $reads_l, $insert, $insert_shift, $min_coverage ) = @_;
	my @clusters = get_clusters ( $file ); # effectuer un split (\t) sur les fichiers .clustered
	my @result;                            # [0] est le n° du cluster
	my @cluster_type = split '\.', $file;  # [1] est #reads, [2] le début de l'accepteur, [3] la fin, [4] et [5] pareil pour le donneur, [6] est une position inconnue, [7] une quantité décimale inconnue, [8] une liste de positions de reads.
	my @is_clust;
	# initialization 
	for my $i ( 0 .. $#clusters ){ $is_clust[$i] = 0 }
	my ( $n, $start_a, $end_a, $start_b, $end_b, $coord ) = ( "", "", "", "", "", "" );	
	my $count=0;
	my $count_trans=0;
#	print "considering file $file\n";
	for my $i ( 0 .. $#clusters ){ # pour chaque cluster dans le fichier
		if( $is_clust[$i] == 0 ){
			( $n , $start_a, $end_a, $start_b, $end_b, $coord ) = ( $clusters[$i][1], $clusters[$i][2], $clusters[$i][3] + $reads_l, $clusters[$i][4], $clusters[$i][5] + $reads_l, $clusters[$i][6] );
			$coord = $clusters[$i][8] if ( $cluster_type[4]=~/del|ins|dup|inv/ ); # coord est le champ [6] si on a une translocation, la liste de positions [8] sinon.
			my $nclust = 1;
#			print "Considering $start_a $end_a $start_b $end_b\n" if $start_b=~/^233/;
			# merging
			for ( my $j = $i+1; $j <= $#clusters; $j++ ){ # pour chaque cluster non encore visité
				last if ( $end_a + ( 2 * $insert_shift ) < $clusters[$j][2] ); #on saute si le début de l'accepteur suivant est trop loin
				                                                               # (on a dépassé les clusters susceptibles de merger)
				next if( $is_clust[$j] == 1 );
				my @overlap_right = is_overlapping ( $start_b , $end_b , $clusters[$j][4] , $clusters[$j][5] , 0 );
                                # return: an array (intersection_size, union_size, union_min, union_max)
				my @overlap_left = is_overlapping ( $start_a , $end_a , $clusters[$j][2] , $clusters[$j][3] , 0 );
				my @set_diffs_right = (abs($clusters[$j][4]-$start_b),abs($clusters[$j][4]-$end_b),abs($clusters[$j][5]-$start_b),abs($clusters[$j][5]-$end_b));
				my @set_diffs_left = (abs($clusters[$j][2]-$start_a),abs($clusters[$j][2]-$end_a),abs($clusters[$j][3]-$start_a),abs($clusters[$j][3]-$end_a));
				my ( $left_is_near, $right_is_near ) = ( 0, 0 );
				foreach my $diff (@set_diffs_right){
					$right_is_near = 1 if ( $diff <= $insert_shift );
				}
				foreach my $diff (@set_diffs_left){
					$left_is_near = 1 if ( $diff <= $insert_shift );
				}
				if ( $left_is_near == 1 && $right_is_near == 1){
					$n += $clusters[$j][1]; 
					 ( $start_a , $end_a , $start_b, $end_b ) = ( $overlap_left[2] , $overlap_left[3] , $overlap_right[2], $overlap_right[3] );
					$coord .= $clusters[$j][8] if ( $cluster_type[4]=~/del|ins|dup|inv/ );
					$is_clust[$j] = 1;
					++$nclust;
					#print ("Cluster ", $clusters[$j][1], "\t", $clusters[$j][2], "\t", $clusters[$j][3], "\t", $clusters[$j][4], "\t", $clusters[$j][5], "\n \thas been merged into $start_a $end_a $start_b $end_b.\n");
				}	
			}
#			print "$nclust clusters\n" if $start_b=~/^233/;
#			 print "$nclust clusters out of $#clusters have been merged into this one.\n";	
			### Get merged or non-merged cluster -------------------------------------
			# filter result for deletion signatures
#			print "Now treating a $cluster_type[4] cluster.";
			if (  $cluster_type[4]=~/del|dup|inv/  && $end_a < $start_b && $end_a-$start_a <= $insert + $insert_shift && $end_b-$start_b <= $insert + $insert_shift){ # si donneurs et accepteurs semblent valides (IS distr)
				my @stat = get_coordonate_stat ( $coord );
				# retourne un tableau (mean, SD) de l'IS
				++$count;
				if ( $start_b - $end_a <=  $stat[0] - $insert + 3 * $stat[1] && $start_b - $end_a >=  $stat[0]-$insert - 3 * $stat[1] ){# ){

					if ( $stat[0] - $insert <= 100000 && $stat[1] < $SPAN[1] * 2 * $SPAN_MAD ){
#					if ($stat[1] < $SPAN[1] * 2 * $SPAN_MAD ){
						if ( $n >= $insert * $min_coverage / ( 2 * $reads_l ) ){
							push @result, [ $count, $n, $start_a, $end_a, $start_b, $end_b, int( $stat[0] - $insert ), $stat[1] ];
						}
					}
					elsif ( $n >= ($end_b-$start_b) * $min_coverage / ( 2* $reads_l ) && $n >= ($end_a-$start_a) * $min_coverage / ( 2 * $reads_l ) ){

						push @result, [ $count, $n, $start_a, $end_a, $start_b, $end_b, int( $stat[0] - $insert ), $stat[1] ];
					}
				    }
			} # result for insertion signatures
			elsif ( $cluster_type[4] eq "ins" && $end_a <= $start_b &&  $end_a-$start_a <= $insert + $insert_shift ){
				my @stat = get_coordonate_stat ( $coord );
				if ( $start_b - $end_a <=  $insert-$stat[0] + 3*$stat[1] && $start_b - $end_a >=  $insert-$stat[0] - 3*$stat[1] && $stat[1] < $SPAN[1] * 2 * $SPAN_MAD ){
					if ( $stat[0] * $min_coverage / $reads_l >= $n ){
						++$count;
						push @result, [ $count, $n, $start_a, $end_a, $start_b, $end_b, int ( $insert - $stat[0] ), $stat[1] ] ;
					}
				}
			} # filter for translocation signatures
			elsif ( $cluster_type[4] eq "trans"){
			    ++$count_trans;
			    if($end_a-$start_a <= $insert + $insert_shift && $end_b-$start_b <= $insert + $insert_shift ){
				++$count;
				push @result, [ $count, $n, $start_a, $end_a, $start_b, $end_b ] ;}
			    else{
#				print "TRANS cluster discarded: acceptor size ", $end_a-$start_a, "donor size: ", $end_b-$start_b, ", expected a maximum of ", $insert + $insert_shift , ".\n"
			    }
			} # filter for inversion signatures
			elsif ( $cluster_type[4] eq "inv"){
			    if( $end_a-$start_a <= $insert + $insert_shift && $end_b-$start_b <= $insert + $insert_shift &&  $end_a < $start_b ){
				++$count;
				push @result, [ $count, $n, $start_a, $end_a, $start_b, $end_b ] ;
			    }#else{print "INV cluster discarded: $start_a $end_a $start_b $end_b, acceptor size ", $end_a-$start_a, " donor size: ", $end_b-$start_b, ", expected a maximum of ", $insert + $insert_shift , ".\n"}
			}									
		}
	}
#	print "In total $count_trans TRANS clusters have been treated.\n";
	# print results
	open ( OUT , ">$file.merged" ) || die "Cannot open $file.merged, $!\n";
	for my $i ( 0 .. $#result){
		print OUT join ("\t", @{$result[$i]}),"\n";
	}
	close ( OUT );
}
#########################################################################################
#				get_clusters
#----------------------------------------------------------------------------------------
# take a file
# return a hash containing the cluster
sub get_clusters {
	my ( $file ) = @_;
	my @clusters;
	open ( IN , $file ) || die "[error]: Cannot open $file, $!\n";
	while ( <IN> ){
		chomp;
		my @data = split '\t', $_;
		push @clusters , [ @data ] ; 
	}
	close ( IN );
#	print "Read ", $#clusters, " clusters from file ", $file, ".\n";
	return @clusters;
}
###################################################################################
#			get_coordonate_stat 					  #
#----------------------------------------------------------------------------------
# take: an array of x y 
# return: mean and standard deviation of the (y-x) values
sub get_coordonate_stat{
	my ( $coord ) = @_;
	$coord =~s/^\((.+)\)$/$1/;
	my @data = split '\)\(', $coord;
	my $mean = fragment_size_mean ( \@data );
	my $sum = 0;
	foreach ( @data ){
		my @xy = split '\,' , $_;
		$sum += ( $xy[1] - $xy[0] - $mean ) * ( $xy[1] - $xy[0] - $mean );
	}
	return ( $mean , sqrt ( $sum / @data ) );
}
###################################################################################
#			fragment_size_mean					  #
#----------------------------------------------------------------------------------
# take: an array of x y
# return: the mean of the (y-x)
sub fragment_size_mean {
	my ( $values ) = @_;
	my $sum = 0;
	foreach  ( @$values ){
		my @xy = split '\,' , $_;
		$sum += $xy[1] - $xy[0];
	}
	return $sum / ( $#{@$values} + 1 );
}


##############################################################################################
#			sv_calling
#---------------------------------------------------------------------------------------------
# 
#
sub sv_calling {
	my ( $in , $cov , $reads , $interval , $lim_inf ) = @_;
	print("Calling structural variants...\n");

	# interchromosomal translocation calling
	my %CLUSTERS_TRANS = get_clusters_data ( $in , "trans" );
	interchromo ( \%CLUSTERS_TRANS , $reads , $interval , $cov, $lim_inf );

	# intrachromosomal translocation calling (sens)
	my %CLUSTERS_DEL = get_clusters_data ( $in , "del" );
	#print(scalar(keys(%CLUSTERS_DEL)), " deletion files read.\n");

	my %CLUSTERS_DUP = get_clusters_data ( $in , "dup" );
	#print(scalar(keys(%CLUSTERS_DUP)), " duplication files read.\n");

	intrachromo ( \%CLUSTERS_DEL , \%CLUSTERS_DUP , $reads , $interval , $cov , "sens" , $lim_inf );

	# intrachromosomal translocation-inversion calling 
	my %CLUSTERS_MINUS_INV = get_clusters_data ( $in , '\-\.inv' );
	my %CLUSTERS_PLUS_INV = get_clusters_data ( $in , '\+\.inv' );
	intrachromo ( \%CLUSTERS_MINUS_INV , \%CLUSTERS_PLUS_INV , $reads , $interval , $cov , "inverse" , $lim_inf );

        #removing previous runs if they exist
#	print "The dir contains ", glob("$OUTPUT_DIR/sv.*"), " ...files.\n";
	`rm $OUTPUT_DIR/sv.*` if (defined(glob("$OUTPUT_DIR/sv.*")) && glob("$OUTPUT_DIR/sv.*") ne "");

	# sortint and selection of SV
	print(scalar(@TRANS), " translocation events called.\n");
	print(scalar(@DEL), " deletion events called.\n");
	sort_select (\@TRANS,"TRANS") if $#TRANS>=0 ;
	sort_select (\@DEL,"DEL") if $#DEL>=0;

        #Remove eyecandy columns from SV file
	if(scalar(@TRANS) !=0 || scalar(@DEL) !=0){
	`cp $OUTPUT_DIR/sv.out $OUTPUT_DIR/sv.formatted`;
	`cut -f3- $OUTPUT_DIR/sv.out > $OUTPUT_DIR/tmp` if ( $? >>8 ==0);
	`mv $OUTPUT_DIR/tmp $OUTPUT_DIR/sv.out 2>&1 /dev/null` if ( $? >>8 ==0);
    }
}
#############################################################################################
#			sort_select
#--------------------------------------------------------------------------------------------
#
#
sub sort_select{
    my ( $called , $type ) = @_;
    my @sorted_trans = sort { $a->[0] cmp $b->[0] || $a->[1] <=> $b->[1] } @$called;
    my @donor = @{$sorted_trans[0]} ;
    my $isnear=0;
    my $count =1;
    my @lines;
    open (OUT , ">>:utf8", "$OUTPUT_DIR/sv.out") || die "Cannot open $OUTPUT_DIR/sv.out, $!,\n";
    foreach my $i ( 0 .. $#sorted_trans ){
	if ( $donor[0] eq $sorted_trans[$i]->[0] ){
	    my @overlap = is_overlapping ( $donor[1], $donor[2] , $sorted_trans[$i]->[1] , $sorted_trans[$i]->[2] , 0 );
	    if ($donor[9] eq "" || $sorted_trans[$i]->[9] eq "") {
		# THIS HAPPENS IN CASE OF DELETIONS
		next;
 	    }
	    # && ( ( $donor[10] / $donor[9]) < ( $sorted_trans[$i]->[10] / $sorted_trans[$i]->[9])) was in the foll. condition
	    if  ( $overlap[0] > 0 ){
		if($isnear==1){
		    # a similar donor already exists
#		    print OUT "\x{251C}";
		    if( test_overcov ( \@donor , \%OVER_COV ) == 0 ){ print OUT "\t", $count, "\t", $type, "\t" , join ( "\t" , @donor ) , "\n";}
		    else { print OUT "\tDonor ignored due to over-coverage.\n";} 
		    
		}
		else{
		    # first such donor
#		print OUT "\x{256D}";
		if( test_overcov ( \@donor , \%OVER_COV ) == 0 ){
		    print OUT "\t", $count, "\t", $type, "\t" , join ( "\t" , @donor ) , "\n";
		} else { print OUT "\tDonor ignored due to over-coverage.\n";} 
		$isnear=1;
		}
		@donor = @{$sorted_trans[$i]};
	    }
	    elsif  ( $overlap[0] == 0 ){
		if($isnear==1){
#		    print OUT "\x{2570}";
		    if( test_overcov ( \@donor , \%OVER_COV ) == 0 ){
			print OUT  "\t", $count, "\t", $type ,"\t" , join ( "\t" , @donor ) , "\n";
		    } else { print OUT "\tDonor ignored due to over-coverage.\n";} 
		    @lines=();
		    $isnear=0;
		}
		else{
		    print OUT " ", "\t", $count, "\t", $type ,"\t" , join ( "\t" , @donor ) , "\n";
		}
		++$count;
		@donor = @{$sorted_trans[$i]};
	    }
	}
	else {
	    if($isnear==1){
#		print OUT "\x{2570}";
		if( test_overcov ( \@donor , \%OVER_COV ) == 0 ){
		    print OUT "\t", $count, "\t", $type ,"\t" , join ( "\t" , @donor ) , "\n";
		} else { print OUT "\tDonor ignored due to over-coverage.\n";} 
		@lines=();
		$isnear=0;
	    }
	    else{
		print OUT "", "\t", $count, "\t", $type ,"\t" , join ( "\t" , @donor ) , "\n" ;
	    }
	    $count++;
	    @donor = @{$sorted_trans[$i]};
	}
    }
 #   print OUT "\x{2570}" if($isnear==1);
    if ( test_overcov ( \@donor , \%OVER_COV ) == 0 ){
	print OUT "\t", $count, "\t", $type ,"\t", join ( "\t" , @donor ) , "\n";
    }else{print OUT "Donor ignored due to over-coverage.\n";} 
    close(OUT);

    if($type eq "TRANS"){
	my %events;
	open (CHECK , "<:utf8", "$OUTPUT_DIR/sv.out") || die "Cannot open $OUTPUT_DIR/sv.out, $!,\n";
	while(<CHECK>){
	    chomp;
	    my $line=$_;
	    next if($line =~/Donor/);
	    my @data=split("\t", $line);
	    push (@{$events{$data[1]}}, $line );
	}
	close(CHECK);
	my @file;
	foreach my $key (sort({$a <=> $b} keys(%events))){
	    my @multiple=@{$events{$key}};
	    next if (scalar(@multiple)==1);
	    my %acceptors;
	    foreach my $line (@multiple){
		my @data = split("\t", $line);
		$acceptors{"$data[8]:$data[9]-$data[10]"}=$line;
		my @tmp=values %acceptors;
		$events{$key}=\@tmp;
		# This ensures no acceptor duplicates in multiple donor clusters
	    }
	}
	# write to file.
	my $id=0;
	my $mem;
	open (OUT , ">:utf8", "$OUTPUT_DIR/sv.out") || die "Cannot open $OUTPUT_DIR/sv.out, $!,\n";
	foreach my $key (sort({$a <=> $b} keys(%events))){
	    print OUT "\x{256D}" if(scalar (@{$events{$key}})>1);
	    print OUT $events{$key}->[0], "\n";
	    if(scalar (@{$events{$key}})>1){
		for (my $i=1; $i<$#{$events{$key}}; $i++){
		    print OUT "\x{251C}", $events{$key}->[$i], "\n";
		}
		print OUT "\x{2570}", $events{$key}->[$#{$events{$key}}], "\n";
	    }
	}
    }
}
#########################################################################################
#				get_clusters_data
#----------------------------------------------------------------------------------------
# take: a directory
# return: a hash containing the clusters data
sub get_clusters_data {
	my ( $dir , $type ) = @_;
	my %clusters;
	opendir ( my $open_dir , $dir ) || die "Cannot open $dir, $!\n";
	my @files = readdir ( $open_dir );
	foreach my $file ( @files ){
		if ( $file=~/.+\.$type\.init\.clustered\.merged$/ ){
			open ( IN , "$dir/$file" )||die"Cannot open $dir/$file,$!\n";
			while ( <IN> ){
				chomp;
				my @data = split '\t', $_;
				push @{$clusters{$file}} , [ @data ] ; 
			}
			close (IN);
		}
	}
	close ( $open_dir );
	return %clusters;
}
###################################################################################
#			interchromo
#----------------------------------------------------------------------------------
# take: a clusters data with interchromosomal signature
# return: interchromosomal translocation
sub interchromo {
	my ( $clusters, $reads_l , $maxint , $min_cov , $lim_inf) = @_;
	my %clusters = %$clusters;
	my @analysed_clusters;
	my $comp = 0;
	my @trans;

	# get all interchromosomal insertion                                     # En entrée, %clusters indice un tableau de clusters par le fichier 
	foreach my $clust (keys %clusters){                                      # dans lequel ils apparaissent. Pour chaque fichier, on boucle sur les autres
		my @cluster_type = split '\.', $clust;                           # (cluster_type et new_cluster_type) et on appelle call_interc_insert, 
		foreach my $new_clust (keys %clusters){	                         # une translocation est possible.
			my @new_cluster_type = split '\.', $new_clust; 
			if ( avoid_duplication ( $clust , $new_clust , \@analysed_clusters ) == 0 ){
				push ( @{$analysed_clusters[$comp]}, $clust, $new_clust );
				++$comp;
				my $sens = is_insertion_possible ( \@cluster_type , \@new_cluster_type );
#				print "Considering $clust versus $new_clust.\n";
				call_interchromosomal_insertion ( $clusters{$clust} , $clusters{$new_clust} , \@cluster_type , \@new_cluster_type ,  $reads_l , $maxint , $sens , $min_cov, $lim_inf ) if ( $sens ne "no" );
			}
		}	
	}	
}

######################################################################################
#			call_interchromosomal_insertion
#-------------------------------------------------------------------------------------
#
#
sub call_interchromosomal_insertion {
	my ( $cluster_data , $new_cluster_data , $cluster_type , $new_cluster_type, $reads_length , $maximal_size_shift , $sens , $min_cov, $lim_inf ) = @_;
	for my $i ( 0 .. $#{$cluster_data} ){
		for my $j ( 0..$#{$new_cluster_data} ){		
			my @insertion = get_insertion ( $cluster_data->[$i] , $new_cluster_data->[$j] , $cluster_type, $new_cluster_type , $sens , $reads_length , $maximal_size_shift );	
			if ( $#insertion > 0){
			    if($insertion[10] >= $insertion[9] * $min_cov / $reads_length ){
				if ( $insertion[9] < $lim_inf ){
				    if ( $insertion[8]/$insertion[9] > 0.8 ){
					push @TRANS, [@insertion] ;
				    }#else{print "Failed test 3 \n";}
                                }
                                elsif($insertion[2]-$insertion[1]<1000){
				    push @TRANS, [@insertion] ;
                                }#else{print "Failed test 2 \n";}
			    }#else{print "Failed test 1 \n";}
			}#else{print "Failed test 0 \n";}
		}
	}
}

##################################################################################
#				get_insertion
#----------------------------------------------------------------------------------
# take 2 arrays of cluster data and their cluster name, the hypothetical sens of insertion, and the read length 
# return an array with insertion features
sub get_insertion {
	my ( $clust_data , $new_clust_data , $clust_type , $new_clust_type , $sens , $reads_l , $maxint ) = @_;
	my @result;
	my @border;
	# Check if the two left clusters are overlapping
	my @ll_overlap = is_overlapping ( $clust_data->[2] , $clust_data->[3] , $new_clust_data->[2] , $new_clust_data->[3] , 0 );	
	if ( $ll_overlap[0] > 100 ){
		@border = getBorders ( $clust_type->[3], $new_clust_type->[3], $clust_data->[4], $clust_data->[5], $new_clust_data->[4], $new_clust_data->[5], $reads_l); # retourne (extrémité1, extrémité2)
		if ( ( $border[1] - $border[0] >= 0 )  && ( ( $border[1] - $border[0]) < $maxint ) ){
			push (@result, $clust_type->[2], $border[0], $border[1], $border[1] - $border[0], $clust_type->[0], $ll_overlap[2], $ll_overlap[3], $sens, $ll_overlap[0], $ll_overlap[1], $clust_data->[1] + $new_clust_data->[1] ); 
		}#else{print "Failed size test 1 with ",$border[1] - $border[0] , " versus $maxint\n";}
	}#else{print "Failed overlap test 1 with ", $ll_overlap[0], "\t",  $clust_data->[2], "\t", $clust_data->[3], "\t",  $new_clust_data->[2] , "\t", $new_clust_data->[3], "\n";}
	# Check if the two right clusters are overlapping
	my @rr_overlap = is_overlapping ( $clust_data->[4], $clust_data->[5], $new_clust_data->[4], $new_clust_data->[5], 0 );
	if ( $rr_overlap[0] > 100 ){
		@border = getBorders ( $clust_type->[1], $new_clust_type->[1], $clust_data->[2], $clust_data->[3], $new_clust_data->[2], $new_clust_data->[3], $reads_l );
		if ( ( $border[1] - $border[0] >= 0 ) && ( ($border[1] - $border[0]) < $maxint ) ){
			push ( @result, $clust_type->[0], $border[0], $border[1] , $border[1] - $border[0], $clust_type->[2], $rr_overlap[2], $rr_overlap[3], $sens, $rr_overlap[0], $rr_overlap[1],$clust_data->[1] + $new_clust_data->[1] );
		}#else{print "Failed size test 2 with ",$border[1] - $border[0] , " versus $maxint\n";}
	}#else{print "Failed overlap test 2 with ", $rr_overlap[0], "\t",  $clust_data->[4], "\t", $clust_data->[5], "\t",  $new_clust_data->[4] , "\t", $new_clust_data->[5], "\n";}
	return @result
}



##################################################################################
#				avoid_duplication
#----------------------------------------------------------------------------------
# take 2 cluster names and an array of cluster name
# return 1 if clusters are the same or if the comparison has already been done
sub avoid_duplication {
	my ($clust, $new_clust, $analysed_clusters) = @_;
	my $done = 0;
	if ($clust eq $new_clust){
		$done = 1;
	}
	elsif ($#{@$analysed_clusters} > 0){
		for my $analysed_cluster (0..$#{@$analysed_clusters}){
			if ( ( ($clust eq $analysed_clusters->[$analysed_cluster]->[0] ) && ($new_clust eq $analysed_clusters->[$analysed_cluster]->[1] ) ) || ( ( $clust eq $analysed_clusters->[$analysed_cluster]->[1] ) && ( $new_clust eq $analysed_clusters->[$analysed_cluster]->[0] ) ) ){
				$done = 1;
			}
		}
	}
	return $done;				
}

###################################################################################
#			intrachromo
#----------------------------------------------------------------------------------
# take: a clusters data with intrachromosomal signature
# return: intrachromosomal translocation
sub intrachromo {
	my ( $clusters_a , $clusters_b , $reads_l , $maxint , $min_cov , $sens , $lim_inf ) = @_;
	my %clusters_a = %$clusters_a;
	my %clusters_b = %$clusters_b;
	my $comp = 0;
	foreach my $clust ( keys %clusters_a ){
		my @cluster_type = split '\.', $clust;
		foreach my $new_clust (keys %clusters_b){	
			my @new_cluster_type = split '\.', $new_clust;
			if ( $cluster_type[0] eq $new_cluster_type[0] ){
				#print "Considering $clust versus $new_clust.\n";
				call_intrachromosomal_insertion ( $clusters_a{$clust} , $clusters_b{$new_clust} , \@cluster_type , \@new_cluster_type ,$reads_l , $maxint , $sens , $min_cov , $lim_inf  );
			}
		}
	}	
}

######################################################################################
#			call_intrachromosomal_insertion
#-------------------------------------------------------------------------------------
# take: data of a deletion and duplication clusters on the same chromosome
# return: intrachromosomal translocation and true deletion
sub call_intrachromosomal_insertion {
	my ( $cluster_data , $new_cluster_data,  $cluster_type , $new_cluster_type, $reads_length , $maximal_size_shift , $sens , $min_cov , $lim_inf ) = @_;
	for my $i ( 0 .. $#{$cluster_data} ){ #for each DEL cluster (ID, numreads, starta, enda, startb, endb, len, ?FLOAT?)
		my $is_linked = 0;
		for my $j ( 0..$#{$new_cluster_data} ){# and for each dup
			my @insertion = get_insertion ( $cluster_data->[$i] , $new_cluster_data->[$j] , $cluster_type, $new_cluster_type , $sens , $reads_length , $maximal_size_shift );	
			if ( $#insertion > 0  && $insertion[10] > ($insertion[9] * $min_cov / $reads_length) && abs ( $insertion[2]-$insertion[5] )> 50000 ){
#			if ( $#insertion > 0  && $insertion[10] > ($insertion[9] * $min_cov / $reads_length)){
				if ( $insertion[9] < $lim_inf ){
					if ( $insertion[8]/$insertion[9] > 0.8 ){
						push @TRANS, [@insertion] ;
						$is_linked = 1;
					    }#else{print "DID NOT PASS COVERAGE TEST\n";}
				}
				elsif($insertion[2]-$insertion[1]<1000){
					push @TRANS, [@insertion] ;
					$is_linked = 1;
				}#else{print "DID NOT PASS TEST 1 \n";}
			}#else{print "DID NOT PASS TEST 0 with $#insertion >0, ",$insertion[10] ,"> ",($insertion[9] * $min_cov / $reads_length) ," and ",( $insertion[2]-$insertion[5] ) , ">50000\n";}
		}
		call_deletion ( $cluster_data->[$i] , $cluster_type , $lim_inf , $min_cov , $reads_length ) if ( $is_linked == 0 );		
	}
}

########################################################################################
#				call_deletion
#--------------------------------------------------------------------------------------
#
#
sub call_deletion {
	my ( $cluster_data , $cluster_type , $min_cluster_size , $min_cov , $reads_length  ) = @_;
	if ( $cluster_type->[4] eq "del" && $cluster_data->[3]-$cluster_data->[2] > $min_cluster_size && $cluster_data->[5]-$cluster_data->[4] > $min_cluster_size && $cluster_data->[1]>$min_cluster_size*$min_cov/$reads_length){
		push @DEL , [$cluster_type->[0] , @$cluster_data[2..3] ,"","","","","","","", $cluster_data->[1] ] ;
	} 
}


###################################################################################
#				is_insertion_possible
#----------------------------------------------------------------------------------
# get two preclusters name
# return true if the overlapping is possible (from chromosome and orientation of the clusters)
sub is_insertion_possible {
	my ( $clust, $new_clust ) = @_;
	if ( ( $clust->[0] eq $new_clust->[0]) && ($clust->[2] eq $new_clust->[2] ) ){
		if ( ( $clust->[1] eq $clust->[3]) && ($new_clust->[1] eq $new_clust->[3] ) ){
			return "inverse";
		}
		elsif ( ( $clust->[1] ne $clust->[3]) && ($new_clust->[1] ne $new_clust->[3] ) ){
			return "sens";
		}
		else{
			return "no";
		}
	}
	else{
		return "no";
	}	
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
###################################################################################
#				getBorders					  #
#----------------------------------------------------------------------------------
# get two intervals position
# return insertion border
sub getBorders {
	my ( $sens_a, $sens_b, $start_a, $end_a, $start_b, $end_b, $reads_l ) = @_; 
	my $left_border;
	my $right_border;
	if ( ( $sens_a eq "+" ) && ( $sens_b eq "-" ) ){
		$right_border = $start_a;
		$left_border = $end_b - $reads_l;
	}
	elsif( ( $sens_a eq "-" ) && ( $sens_b eq "+" ) ){
		$right_border = $start_b;
		$left_border = $end_a - $reads_l;
	}
	my @border = ( $left_border , $right_border );
	return @border;
}
#
#-------------------------------------------------------------------------------------
#
#
sub test_overcov {
	my ( $var , $overcov ) = @_;
	my $donor_is_over = 0;
	my $acceptor_is_over = 0;
	foreach my $over_region (@{$overcov->{$var->[0]}}){
		my @overlap = is_overlapping ($over_region->[0],$over_region->[1],$var->[1],$var->[2],0);
		if ( $overlap[0] > 0 ){
			$acceptor_is_over = 1;
			last;
		}
	}
	foreach my $over_region (@{$overcov->{$var->[4]}}){
		my @overlap = is_overlapping ( $over_region->[0], $over_region->[1], $var->[5],$var->[6],0 );
		if ( $overlap[0] > 0 ){
			#print $over_region->[0],"\t", $over_region->[1],"\t", $var->[5],"\t",$var->[6],"\n";
			$donor_is_over = 1;
			last;
		}
	}
	return 1 if ( $acceptor_is_over == 1 && $donor_is_over == 1 );
	return 0;
}



###################################################################################
#			get_high_coverage_positions
#----------------------------------------------------------------------------------
# Take a BAM file
# Return a hash of overcovered region
sub get_high_coverage_positions {
    my ( $bam_file , $max_cov, $ext, $read_length ) = @_;
    my %over_cov;
    if (-e "$INPUT_DIR/over_cov.txt"){
	print "Previous run detected. Fetching over-covered regions...\n";
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
		    print $ref,"\t",$start,"\t",$end,"\n";
		    ( $start , $end ) = ( $data[1] , $data[1] );
		}
	    }
	    else{
		( $ref , $start , $end ) = ( $data[0] , $data[1] , $data[1] );
	    }
	}
	print "High covered region ( > $max_cov ) detection done\n";
	open ( OV , ">$OUTPUT_DIR/over_cov.txt") || die "Cannot open $OUTPUT_DIR/over_cov.txt, $!\n";
	foreach my $ref ( keys %over_cov ) {
	    foreach my $i ($#{@{$over_cov{$ref}}}){
		print OV $ref,"\t",$over_cov{$ref}->[$i]->[0],"\t",$over_cov{$ref}->[$i]->[1],"\n";
		
	    }
	}
	close (OV);
	close(IN);
    }
    return %over_cov;


}

sub get_MC_parameter {
    my ( $bam) = @_;
    # Depth of coverage analysis	
    my %indexed_cov;	
    open ( IN, "$SAMTOOLS_PATH pileup $bam |" ) || die "Cannot open $bam, $!\n";
    while ( <IN> ){
	last if ( $. > $SAMPLE_SIZE );
	my @data = split /\t/, $_;
	$indexed_cov{$data[1]}= $data[3] ;
    } 
    close (IN);
#    my $cpos=(sort {$a <=> $b} keys %indexed_cov)[0];
    my $cpos=1;
    print "INFO:\tEvaluating optimal X using Markov-Chain procedure v.1.\n";
    print "INFO:\tStarting at pos $cpos.\n";
    my @cov;
    foreach my $pos (sort {$a <=> $b} keys %indexed_cov){
	print "INFO:\tTreating pos $pos \r";
	if($pos>$cpos){
#	    print "\nINFO:\t$cpos\n";
#	    exit;
	    while($pos>$cpos){
		push(@cov, 0);
		$cpos++;
	}
	}else{
	    push @cov, $indexed_cov{$pos};
	}
	$cpos++;

    }
    

    my $nzero=0;
    my $zerosize=0;
    my $totalsize=0;
    my $thissize=0;
    foreach my $cov_elem (@cov){
	if( $cov_elem<1){
	    if( $zerosize==0){
	    # We are starting an empty region
	    $zerosize=1;
	}else{
	    # We are expanding an empty region
	    $zerosize++;
	}
	}elsif($zerosize!=0){
	    # we are leaving a zero region 
	    if($zerosize>10){$nzero++; $totalsize+=$zerosize;}
	    $zerosize=0;
	}
    }
    print "INFO:\t", $nzero, " low-coverage regions were found with total size of $totalsize and average size of ", $totalsize/$nzero, ".\n";
    my $initial_prob=$totalsize/$SAMPLE_SIZE;
    my $p_0_0=($totalsize-$nzero)/$totalsize;
    print("INFO:\t", $nzero, " low-coverage regions were detected. pi_0=$initial_prob, p_0_0=$p_0_0.\n");
    print "INFO:\tComputing probabilities for intervals.\n";
    my $i=1;
    my %massfunction;
    while ($initial_prob*($p_0_0**$i)>0){
#	print "Interval of size $i\tp=", $initial_prob*($p_0_0**$i), "\tlog(p)=", log($initial_prob*($p_0_0**$i)),"\n";
	$massfunction{$i}=$initial_prob*($p_0_0**$i);
	$i++;
    }

    my $numbars=scalar keys %massfunction;
    my $max=$massfunction{1};
    use List::Util qw(sum);
    my $sum=sum(values(%massfunction));
    foreach my $key ( keys %massfunction){
	$massfunction{$key}/=$sum;
    }

    my $cdf=0;
    $i=1;
    while($cdf<0.999999){
	last if(!defined($massfunction{$i}));
	$cdf+=$massfunction{$i};
#	print $cdf, "\n";
#	print "$cdf % were reached at if($cdf*100 % 10==0);
	$i++;
    }
    print "INFO:\tSum was reached at $sum.\n";   
    print "INFO:\tLimit precision percentile reached at interval of size $i.\n";
    return $i;
}
