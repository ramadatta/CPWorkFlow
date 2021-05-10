# This program takes in the invariant sites filename as ARGUMENT 1
# The file which requires the invariant sites to be appended as ARGUMENT 2
# The output file with non-recombinant polymorphic sites + invariant sites

$invarSites_filename="$ARGV[0]";
$fasta2Append="$ARGV[1]";
$outputFasta_prefix="$ARGV[2]";

print "Executing perl $invarSites_filename $fasta2Append $outputFasta_prefix\_InvarSitesAppend.fasta\n";

open FH,"$invarSites_filename";

$str=<FH>;

open HF,"$fasta2Append";

open OUT,">$outputFasta_prefix\_InvarSitesAppend.fasta";
while(<HF>)
{
$header=$_;
$seq=<HF>;
chomp($seq);
print OUT "$header";
print OUT "$seq$str";
}

$str="";


