#!/usr/bin/perl -w

# This script launches an SLC grid on every dir found in the argument file, assuming identical directory structure

open(IN, $ARGV[0]) or die ("Impossible to open $ARGV[0].\n");
my @input_dirs=<IN>;
close(IN);

my $pwd=`pwd`;
chomp $pwd;


foreach my $dir (@input_dirs){
    my $CMD="./tetracker eris $dir/data.sorted.rmdup.bam -out $dir/filter -treat_bam=input:0-1 -nodiscordant -chainload=treated";
    print "$CMD\n";
}
