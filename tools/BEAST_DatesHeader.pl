
$filename="$ARGV[0]";
$output="$ARGV[1]";
open FH,"IsolateDOC.txt";

while(<FH>)
{
chomp($_);

@arr=split(" ",$_);
$dates{$arr[0]}="$arr[1]";

}


open HF,"$filename";
open OUT,">$output\_BEAST_withDates.fasta";
while(<HF>)
{
#print "$_";
	if($_=~/^>(.+)\s*/ && exists($dates{$1}))
	{
	#print "this $1 is replaced by $dates{$1}\n";
	print OUT">$dates{$1}\n";
	}
	else
	{
	print OUT"$_";
	}
}
		


