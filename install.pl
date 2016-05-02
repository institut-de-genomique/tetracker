#!/usr/bin/perl -w
#use lib "$ENV{HOME}/perl_modules/lib/perl5";
use strict; use warnings;
print "\n\nTE-Tracker install script (v 0.3)\n================================\n\n";
print "INFO:\tChecking perl version...\n";
my $PVERSION=$];
print "INFO:\tChecking presence of modules...\n";
my @NMODULES;
eval { require Getopt::Long; };
push @NMODULES, "Getopt::Long" if ($@);
eval { require POSIX; };
push @NMODULES, "POSIX" if ($@);
eval { require Switch; };
push @NMODULES, "Switch" if ($@);
eval { require List::Util; };
push @NMODULES, "List::Util" if ($@);
eval { require Devel::Peek; };
push @NMODULES, "Devel::Peek" if ($@);
eval { require IO::Handle; };
push @NMODULES, "IO::Handle" if ($@);
eval { require Config::Simple; };
push @NMODULES, "Config::Simple" if ($@);
#print "INFO:\tChecking presence and version of samtools...\n";
#my $NULL=`command -v samtools >/dev/null 2>&1`;
#my $SAMTOOLS_RC= ($? >> 8);
print "INFO:\tChecking java...\n";
my $NULL=`command -v java >/dev/null 2>&1`;
my $JAVA_RC= ($? >> 8);


print "\nSUMMARY:\n";
print "\t-\tYour version of Perl is anterior to 5.8.5.\n" if ($PVERSION < 5.008005);
if($#NMODULES>-1){
    print "\t-\t", scalar(@NMODULES), " module(s) were found to be missing (",join(", ", @NMODULES),").\n";
}
#print "\t-\tSamtools is not installed or was not found in your PATH (errorcode $SAMTOOLS_RC).\n" if ($SAMTOOLS_RC != 0);
print "\t-\tJava is not installed or was not found in your PATH (errorcode $JAVA_RC).\n" if ($JAVA_RC != 0);

if($#NMODULES>-1 || $JAVA_RC != 0 || $PVERSION < 5.008005){
    die("ERROR:\tWill not continue with these errors.\n");
}else{
    print("\t-\tEverything seems to be in order. Continuing.\n\n");
}


## PICARD
print "INFO:\tGetting picard...\n";
$NULL=`wget http://sourceforge.net/projects/picard/files/latest/download -O picard.zip >/dev/null 2>&1`;
die("ERROR: \"$!\" while downloading Picard to the current directory. (Are you connected to the Internet?) \n") if (($? >>8) !=0);
print "INFO:\tUncompressing picard...\n";
$NULL=`unzip -o picard.zip`;
die("ERROR: \"$!\" while unzipping picard to the current directory.\n") if (($? >>8) !=0);
$NULL=`rm picard.zip`;
die("ERROR: \"$!\" while removing picard.zip in the current directory.\n") if (($? >>8) !=0);
my $PICARD_PATH=`pwd`;
chomp $PICARD_PATH;
$PICARD_PATH .= "/".(glob("picard*"))[0];
print "INFO:\tPicard path set to $PICARD_PATH\n";


## SAMTOOLS
print "INFO:\tGetting samtools...\n";
$NULL=`wget http://sourceforge.net/projects/samtools/files/samtools/0.1.12/samtools-0.1.12a.tar.bz2/download -O samtools.tar.bz2 >/dev/null 2>&1`;
die("ERROR: \"$!\" while downloading samtools to the current directory.\n") if (($? >>8) !=0);
print "INFO:\tUncompressing samtools...\n";
$NULL=`tar --overwrite -xvjf samtools.tar.bz2`;
die("ERROR: \"$!\" while unzipping samtools to the current directory.\n") if (($? >>8) !=0);
$NULL=`rm samtools.tar.bz2`;
die("ERROR: \"$!\" while removing samtools archive in the current directory.\n") if (($? >>8) !=0);
my $SAM_PATH=`pwd`;
chomp $SAM_PATH;
$SAM_PATH .= "/".(glob("samtools*"))[0];
$SAM_PATH.="/samtools";
print "INFO:\tSamtools path set to $SAM_PATH\n";
print "INFO:\tMaking samtools...\n";
my $PTH=`dirname $SAM_PATH`;
chomp($PTH);
print "\t\t\t$PTH\t\t\n";
$NULL=`cd $PTH ; make`;
die("ERROR: \"$!\" while making samtools. Do you have the required permissions?\n") if (($? >>8) !=0);

## BEDTOOLS
print "INFO:\tGetting BEDtools...\n";
`pwd`;
$NULL=`wget http://bedtools.googlecode.com/files/BEDTools.v2.17.0.tar.gz -O  bedtools.tar.gz >/dev/null 2>&1`;
die("ERROR: \"$!\" while downloading BEDtools to the current directory.\n") if (($? >>8) !=0);
print "INFO:\tUncompressing BEDtools...\n";
$NULL=`tar --overwrite -xvzf bedtools.tar.gz`;
die("ERROR: \"$!\" while unzipping BEDtools to the current directory.\n") if (($? >>8) !=0);
$NULL=`rm bedtools.tar.gz`;
die("ERROR: \"$!\" while removing BEDtools archive in the current directory.\n") if (($? >>8) !=0);
my $BED_PATH=`pwd`;
chomp $BED_PATH;
$BED_PATH .= "/".(glob("bedtools*"))[0];
print "INFO:\tBEDtools path set to $BED_PATH\n";
print "INFO:\tMaking BEDtools...\n";
$NULL=`cd $BED_PATH; make`;
$BED_PATH.="/bin";
die("ERROR: \"$!\" while making BEDtools. Do you have the required permissions?\n") if (($? >>8) !=0);



## WRITING CONFIG

print "INFO:\tPre-writing TE-Tracker.conf...\n";

require Config::Simple;
my $cfg = new Config::Simple(syntax=>'ini');
my $BASEDIR=`pwd`;
chomp $BASEDIR;
$cfg->param("GLOBAL.basedir", $BASEDIR);
$cfg->param("GLOBAL.rms_command", "");
$cfg->param("GLOBAL.rms_args", "");
$cfg->param("GLOBAL.highmem_launch", "\'\'");
$cfg->param("GLOBAL.highmem_thresh", "\'\'");
$cfg->param("GLOBAL.samtools_path", $SAM_PATH);
$cfg->param("GLOBAL.picard_path", $PICARD_PATH);
$cfg->param("GLOBAL.rms_command", "\'\'");
$cfg->param("PREPROCESSING.binary", "eris.pl");
$cfg->param("PREPROCESSING.default_params", "");
$cfg->param("CLUSTERING.binary", "leto.pl");
$cfg->param("CLUSTERING.default_params", "");
$cfg->param("CLUSTERING.slc_path", $BASEDIR."/slclust");
$cfg->param("POSTPROCESSING.binary", "eris.pl");
$cfg->param("POSTPROCESSING.default_params", "");
$cfg->param("POSTPROCESSING.bedtools_path", $BED_PATH);
$cfg->write("$BASEDIR/TE-Tracker.conf");

print "INFO:\tDone.\n\n";
