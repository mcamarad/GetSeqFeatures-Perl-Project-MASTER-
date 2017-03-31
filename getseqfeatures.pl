#!/usr/bin/perl
use strict;
use warnings;

#Initialization and variables
my ($fa, $i, $j, $k);
my @IDs, my @starts, my @ends, my @strands;
#Opening .gff file
&opening_gff($ARGV[0]);

#Opening .fa file
$fa = &opening_fa($ARGV[1]);


#For-loop to associate each sequence with the features and print the output
print STDERR "#Initializating extraction...\n\n";
for ($i=0; $i < scalar(@IDs); $i++) {
    $j=$starts[$i]-1;
    $k=($ends[$i]-$starts[$i])+1;
    &asso_fa_GFF($strands[$i], $fa, $j, $k, $IDs[$i], $starts[$i], $ends[$i]);
};

print STDERR "#Completed succesfully...\n\n";

#For-loop to GC content of introns
print STDERR "#Initializating GC content of introns...\n\n";
for ($i=0; $i < scalar(@IDs); $i++) {
	&intronizer($ends[$i-1], $starts[$i+1], $IDs[$i], $IDs[$i+2], $IDs[$i+3], $strands[$i], $fa,);
		
};
print STDERR "#Completed succesfully...\n\n";

#####################################################################################################################################		
#####################################################IMPLEMENTING FUNCTIONS CODE#####################################################
#####################################################################################################################################

######################################################Open .fa files
sub opening_fa($$){
	#Initialization of variables
	my ($sequence, $seq);
	#Main block
	open(INPUT, $_[0]) or die "You're introducing wrongly the .fa file. Introduce the .gff file first and then the .fa file/n";
 	while (<INPUT>) {
	   next if /^>/o;		    	#jumps directly to the sequence  
	   chomp $_;
  	   $sequence=substr($_, 0, length($_));
  	   unless (defined $seq){          	#The first time introduces the first sequence 
		 $seq="$sequence";		#and then always else to avoid "use of uninitialized..." issue	
  	   }
  	   else {
		 $seq="$seq"."$sequence";
  	   };

  	};

	close (INPUT);
	print STDERR "#.fa file was read...\n\n";
	return $seq;
}; #opening_fa


######################################################open .gff files 
#(This subroutine works with global arrays and local variables in order to fulfill the task given)
sub opening_gff($$){
	#Initialization of variables
	my ($ID, $group, $start, $end, $strand);
	#Main block
	open(INPUT, $_[0]) or die "You're introducing wrongly the .gff file. Introduce the .gff file first and then the .fa file\n";
	while (<INPUT>) {
	    next if /^\s*$/o;			#"eliminates" blank spaced lines and compiles 
	    chomp $_;
	    $ID =  (split(/\s/, $_))[2];
	    $group = (split(/\s/, $_))[9];
	    $group =~ tr/";//d;			#eliminates trash from the ID of the feature
	    push(@IDs, "$ID."."$group");
	    $start = (split(/\s/, $_))[3];
	    push(@starts, $start);
	    $end = (split(/\s/, $_))[4];
	    push(@ends, $end);
	    $strand = (split(/\s/, $_))[6];
	    push(@strands, $strand);
	};
	
	close (INPUT);
	print STDERR "#.gff file was read...\n\n";
	
}; #opening_gff

######################################################Reverse-complementary
sub com_rev ($$) {
	#Initialization of variables
	my $comrev;
	my ($dna) = @_;
	#Main block
	$dna =~ tr/ACGT/tgca/;		#complementary
	$comrev = reverse($dna);	#reverse
}; #com_rev


######################################################Associate sequences to GFF features
sub asso_fa_GFF ($$) {
	#Initialization of variables
	my ($dnaseq, $comrev, $GC);
	my ($strandsfa, $faseq, $jfa, $kfa, $IDsfa, $startsfa, $endsfa) = @_;
	#Main block
	$dnaseq = substr($faseq, $jfa, $kfa);
	if ($strandsfa eq "+"){
		$GC = &GCcontent($dnaseq);
		print STDOUT "$IDsfa\t$startsfa\t$endsfa\t$GC\t$dnaseq\n";
    	}
    	else{
		$comrev = &com_rev($dnaseq);
		$GC = &GCcontent($comrev);
		print STDOUT "$IDsfa\t$startsfa\t$endsfa\t$GC\t$comrev\n";
    	};
}; #asso_fa_GFF

######################################################Analyzing GC content of introns
#Obtaining intron sequences from coords
sub intronizer ($$){
	#Initialization of variables
	my ($endint, $startint, $IDint, $IDint2, $IDint3, $strandint, $faint) = @_;
	my ($m, $n, $intseq, $gici, $intcomrev);
	#Main block
	if (defined $endint && defined $IDint2){
		$m = $endint-1;				#defines the start pos of intron
		$n = ($startint-$endint)+1;
		$intseq = substr($faint, $m, $n);		#defines the length of the intron	
		if ($IDint eq $IDint2 && $IDint =~ /^exon*/){	#searches for exon... lines	
			if ($strandint eq "+"){		
				$gici= &GCcontent($intseq);
				print STDOUT "intron.$IDint\t$gici\n"; #If you want to see the sequence add \t$intseq
			}
			else{
				$intcomrev = &com_rev($intseq);
				$gici = &GCcontent($intcomrev);
				print STDOUT "intron.$IDint\t$gici\n"; #If you want to see the sequence add \t$intcomrev	
			};
		};
		if (defined $IDint3 && $IDint eq $IDint3 && $IDint2=~ /^start*/){ #Special case for NM_142034
			if ($strandint eq "+"){		
				$gici = &GCcontent($intseq);
				print STDOUT "intron.$IDint\t$gici\n"; #If you want to see the sequence add \t$intseq after $gici
			}
			else{
				$intcomrev = &com_rev($intseq);
				$gici = &GCcontent($intcomrev);
				print STDOUT "intron.$IDint\t$gici\n"; #If you want to see the sequence add \t$intcomrev after $gici	
			};
		};
	};
}; #intronizer

#GC content
sub GCcontent($$){
	# initializing variables
	my ($dna_seq) = @_;
	my ($C,$G);
	$C = $G = 0;

	# Looping through the sequence
	my $length = length($dna_seq);
	for (my $i = 0; $i < $length; $i++) {
	    my $char;
	    $char = uc(substr($dna_seq,$i,1));
	  SWITCH: {
	      $char eq 'G' && ($G++, last SWITCH);
	      $char eq 'C' && ($C++, last SWITCH);
	      $char eq 'g' && ($G++, last SWITCH);
	      $char eq 'c' && ($C++, last SWITCH);
    	  };
   
	};
	my $GC = (($G+$C)/$length);
	return $GC;
}; #GCcontent
