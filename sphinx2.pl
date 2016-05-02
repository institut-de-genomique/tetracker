#!/usr/bin/perl -w
# Perl Sphinx script - generates transposition events and repeat assembly errors
$| = 1;

# 1. open transposition events file
my @TE;
open ( IN, "$ARGV[0]" ) or die "Impossible to open $ARGV[0] : $!\n";
while ( <IN> ){
    chomp;
    my $line=$_;
    next if ($line =~ /^#/);
    my @data=split('\t', $line);
    push(@TE, \@data);
}
close(IN);
my @sTE=sort {$a->[1] cmp $b->[1] || $a->[2] <=> $b->[2]} @TE;
#my @sTE=@TE;
undef(@TE);

# 2. open fasta reference and construct database
open ( IN, "/env/cns/proj/projet_AKL/db/test/Arabidopsis_thaliana.fa" ) || die "Cannot open reference: $!\n";
my %faindx;
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
#    print %faindx;
# 3. iterate through the transposition event table, modify chromosomes on the fly
foreach my $i (0..$#sTE){
    my $event=$sTE[$i];
    if($event->[0] ne "D"){
	print(join("\t", @{$event}), "\n");
	my $trseq=substr($faindx{$event->[3]}, $event->[4]-1, $event->[5]-$event->[4]);
	$trseq=Reverse_Complement($trseq) if($event->[6] eq "-");
	my $tsd="";
	$tsd=substr($faindx{$event->[1]}, $event->[2]-1,$event->[7]) if ($event->[7]>0);
#	print ("TSD $tsd\n");
	$trseq=$tsd.$trseq;
	substr($faindx{$event->[1]}, $event->[2]-1, 0)=$trseq;
    }

    # boucle d'ajustement pour tous les donneurs/accepteurs suivants
    for(my $j=$i+1; $j<$#sTE+1; $j++){
	my $nextevent=$sTE[$j];
	# dans tous les cas, on pousse donneur et accepteurs suivants s'ils se trouvent en aval de l'accepteur courant
	# (sauf dans le cas de l'accepteur, s'il se trouve à l'intérieur de l'accepteur courant.
	if($nextevent->[1] eq $event->[1] && $nextevent->[2] > $event->[2]+$event->[5]-$event->[4] && $nextevent->[0] ne "D"){	
	    print(join("\t", "BEFORE1: ", @{$nextevent}), "\n");	
	    $nextevent->[2]+= $event->[5]-$event->[4]+1;
	}
	if ($nextevent->[3] eq $event->[1] && $nextevent->[4] > $event->[2]){
	    print(join("\t", "BEFORE2: ", @{$nextevent}), "\n");	
	    $nextevent->[5]+=($event->[5]-$event->[4]+$event->[7]);
	    $nextevent->[4]+=($event->[5]-$event->[4]+$event->[7]); 
	    # écriture de la nouvelle table de transposition
	    print(join("\t", "AFTER: ", @{$nextevent}), "\n");	
	}
	# dans le cas d'une excision, on tire les coordonnées de ce qui se trouve après le donneur
	if($event->[0] eq "X" || $event->[0] eq "D"){
	$nextevent->[2] -= $event->[5]-$event->[4]+1 if($nextevent->[1] eq $event->[3] && $nextevent->[2] > $event->[5] );
	if ($nextevent->[3] eq $event->[3] && $nextevent->[4] > $event->[5]){$nextevent->[4]-=$event->[5]-$event->[4]+1; $nextevent->[5]-=$event->[5]-$event->[4]+1;}
	}
    }
    if($event->[0] eq "X" || $event->[0] eq "D" ){
    	# si on est dans le cas d'une excision, excision du donneur
    	substr($faindx{$event->[3]}, $event->[4], $event->[5]-$event->[4])="";
    }

 
}

# 4. Rebuild the FASTA file
open(OUT, ">./out2.fa") or die "Impossible to write output file: $!\n";
foreach my $chr (sort keys(%faindx)){
    print OUT ">$chr\n";
    my @chrlines = unpack("(A79)*", $faindx{$chr}); ;
    foreach my $line (@chrlines){
	print OUT "$line\n";
    }
}


sub Reverse_Complement {
   my ($dna) = @_;
   my $revcomp = reverse($dna);
   $revcomp =~ tr/ACGTacgt/TGCAtgca/;
   return $revcomp;
}
