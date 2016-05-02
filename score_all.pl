#!/usr/bin/perl -w
my $filename=$ARGV[0];
my $bam=$ARGV[1];
open(IN, $filename) or die("Unable to open file $filename.\n");
open(OUT, ">$filename.allscores") or die("Impossible to write to $filename.allscores.\n");
my %donor_uniq_count;
while(<IN>){

    chomp;
    next if ($_=~/^\s\*\s$/);
    my @fields=split(/\t/,$_); 
    next if ($#fields<10);
    if($fields[$#fields] eq '*'){
	print "corrected:";
	# replace by count.
	undef %donor_uniq_count;
	my $file=$bam;
	my $acceptor=$fields[2];
	push(@donor, "$fields[4]:$fields[5]-$fields[6]");
	my $achr=(split(':', $acceptor))[0];
	my $astart=(split(/\.\./, (split(/\:/, $acceptor))[1]))[0];
	my $aend=(split(/\.\./, (split(/\:/, $acceptor))[1]))[1];
	$astart-=7000;
	$aend+=7000;

	my $i=0;
	    my $count=0;
	    my $tmp=$donor[$i];
	    open (MIN, "samtools view $file $donor[$i] |") or die
		"Impossible to discriminate donor because there is no file.\n";
	    $donor_uniq_count{$donor[$i]}=0;
	my $dchr=(split(':', $donor[$i]))[0];
	my $matchr= ($dchr eq $achr? "=":$achr);
#	print "acceptor $acceptor, donneur $donor[$i], dchr $dchr, achr $achr, matchr $matchr, astart $astart, aend $aend\n";
	    while (<MIN> ){
		chomp;
		my $line=$_;
		my @data=split('\t', $line);
		my $rlen=length($data[9]);
		next if(($data[1] & 12) != 0);
		if($data[6] eq $matchr && $data[7] > $astart && $data[7] < $aend){
		    if($data[5] =~ /76/ && $line =~ /X0:i:1\s/){
			$donor_uniq_count{$donor[$i]}++;
		    }
		}
	    }
	    close(MIN);
#	close(OUT);
	$fields[$#fields]=$donor_uniq_count{(keys(%donor_uniq_count))[0]};

    }
    print $fields[$#fields], "\n";
    print OUT join("\t", @fields), "\n";

}
