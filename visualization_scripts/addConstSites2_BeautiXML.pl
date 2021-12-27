# This script takes in XML file from Beauti and the beast constant sites in the order A C G T
# Outputs the XML file with constant sites added into the XML file

($xml, $A, $C, $G, $T) = @ARGV;

open FH,"$xml"; 

print " Adding A C G T constant sites: $A $C $G $T to $xml file\n";

$newxml = $xml =~ s/.xml/_withConstSites.xml/gr; 

open HF,">$newxml"; 

while(<FH>)
{
	#print $_;

	if($_=~/\<data id\=/)
	{
		print HF"\<data id\=\"coreOriginal\" name=\"alignment\"\>";
	}
elsif($_=~/\<\/data\>/)
	{
		print HF"<\/data>\n<data id=\"postGubbins\" spec=\"FilteredAlignment\" filter=\"-\" data=\"\@coreOriginal\" constantSiteWeights=\"$A $C $G $T\"\/>\n"
	}
else
	{
	print HF"$_";
	}

}	
