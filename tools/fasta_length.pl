# 23 May 2012
# Natascha May
#
# Computes length for multi-fasta file

open (INFILE, $ARGV[0]);
our ($outfilename) = $ARGV[0];
open (OUTFILE, ">".$outfilename."_length");

## Read Fasta File and compute length ###
my $length;
my $totalLength;
my @arr;
$firstentry = 1;

while(<INFILE>){
   chomp;

   # if ">", store accumulated length and clear the length variable
   if(/>(.*)/){ # non-greedy, title stops at first blank space. <------ might need to edit this according to input file **************
	
	$title = $1;

	if($firstentry == 0){ 
		print OUTFILE $length."\n".$title."\t"; # print current length, and next title
		$length = 0;
		next;
	}else{
		print OUTFILE $title."\t"; # print first title
		$firstentry = 0;
		next;	
	}
  }

  # keep adding length of each line while not ">"
  $length += length($_); 
}

print OUTFILE $length."\n"; # print length for final entry

close(INFILE);
close(OUTFILE);

