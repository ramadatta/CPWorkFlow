# 15 Mar 2011
# Natascha May
#
# Convert fasta file to one-liner, tab-delimited (go to ####)
# NOTE: Change match pattern to detect start of sequence (go to XXXX)

use warnings;

open (INFILE, $ARGV[0]);
#our ($outfilename) = split(/\./,$ARGV[0]);
open (OUTFILE, ">".$ARGV[0]."_OneLine.fasta");

$first_header = 1;
$first_seq = 1;

while (<INFILE>){
	chomp($_);
	$line = $_;
	$line =~ s/\n|\r//; #**sometimes chomp is not enough**#
	
	if($line =~ /^>/){
		$first_seq = 1;

		if($first_header == 1){
			print OUTFILE "$line";
			$first_header = 0;			
		}		
		else{		
			print OUTFILE "\n$line";
		}
	}
	#if($line !~ /=|GSRW|G7DU|HH3EI6|\||-|_|isotig|Contig|contig|SB|sb|\)/){ #non header line (start of sequence) XXXX
	if($line !~ /^>/){ #non header line (start of sequence) XXXX
		if($first_seq == 1){
			print OUTFILE "\n$line"; #edit this line to change the delimiter for fasta file format ####
			$first_seq = 0;
		}else{
			print OUTFILE "$line";	
		}
	}

}

close(INFILE);
close(OUTFILE);
