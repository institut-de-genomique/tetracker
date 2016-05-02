#!/usr/bin/perl -w 
open(IN, "samtools view $ARGV[0] |") or die "Impossible d'ouvrir $ARGV[0].\n";
my %pairdict;
while($.<1000000){
    $_=<IN>;
    chomp;
    my @data=split(/\t/, $_);
    if(($data[1] & 12)==0){
	if(!defined($pairdict{$data[0]}) && $data[8]!=0){
	    # Then we have a mapped pair that was not archived. Store insert size.
	    $pairdict{$data[0]}=$data[8];
	}
    }

}
close(IN);

foreach my $key (keys %pairdict){
    print $pairdict{$key}, "\n";
}
