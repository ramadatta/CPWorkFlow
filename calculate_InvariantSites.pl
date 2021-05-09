# This program takes the faidx stats of reference and core gene alignment 

#!/storage/apps/anaconda3/bin/perl

use Math::Round qw( round );
$Ref_Stats_file = $ARGV[0]; # faidx of reference genome
$CGA_Stats_file = $ARGV[1]; # faidx of core genome alignment


open $info, $Ref_Stats_file or die "Could not open $Ref_Stats_file: $!";

while( $line = <$info>)  {   

	if($line=~/^name/)
	{
	next;
	}
	else
	{
	($name,$start,$end,$A,$T,$C,$G,$N)=split("\t",$line);
	
	print "A:$A,T:$T,C:$C,G:$G,N:$N\n";
	$Aratio=$A/$end;
	$Tratio=$T/$end;
	$Cratio=$C/$end;
	$Gratio=$G/$end;

	print "A:$Aratio,T:$Tratio,C:$Cratio,G:$Gratio\n";
	
	$rounded_Aratio=round($Aratio);
	$rounded_Tratio=round($Tratio);
	$rounded_Cratio=round($Cratio);
	$rounded_Gratio=round($Gratio);

	print "---->Reference Genome nucleotide proportion A:$rounded_Aratio,T:$rounded_Tratio,C:$rounded_Cratio,G:$rounded_Gratio\n";

	}

}

close $info;

open $info2, $CGA_Stats_file or die "Could not open $CGA_Stats_file: $!";

while( $line2 = <$info2>)  {   

	if($line2=~/^name/)
	{
	next;
	}
	else
	{
	($name,$coregenome_start,$coregenome_end,$A,$T,$C,$G,$N)=split("\t",$line2);
	print "core genome end: $coregenome_end\n";
	$diff=$end-$coregenome_end;

	print "Reference genome size (end) is:$end, core genome size (end): $coregenome_end, Difference: $diff\n\n\n";
	$Num_of_As_toAppend=$Aratio*$diff;
	$Num_of_Ts_toAppend=$Tratio*$diff;
	$Num_of_Cs_toAppend=$Cratio*$diff;
	$Num_of_Gs_toAppend=$Gratio*$diff;

	print "A:$Num_of_As_toAppend,T:$Num_of_Ts_toAppend,C:$Num_of_Cs_toAppend,G:$Num_of_Gs_toAppend\n";

	$rounded_Num_of_As_toAppend=round($Num_of_As_toAppend);
	$rounded_Num_of_Ts_toAppend=round($Num_of_Ts_toAppend);
	$rounded_Num_of_Cs_toAppend=round($Num_of_Cs_toAppend);
	$rounded_Num_of_Gs_toAppend=round($Num_of_Gs_toAppend);

	print "A:$rounded_Num_of_As_toAppend,T:$rounded_Num_of_Ts_toAppend,C:$rounded_Num_of_Cs_toAppend,G:$rounded_Num_of_Gs_toAppend\n";

	$roundedSum=$rounded_Num_of_As_toAppend+$rounded_Num_of_Ts_toAppend+$rounded_Num_of_Cs_toAppend+$rounded_Num_of_Gs_toAppend;
	open OUT,">beast2_constsites";
	if($roundedSum==$diff)
	{
	print "$roundedSum==$diff ----> rounding is perfect!\n";
	print OUT"A:$rounded_Num_of_As_toAppend,T:$rounded_Num_of_Ts_toAppend,C:$rounded_Num_of_Cs_toAppend,G:$rounded_Num_of_Gs_toAppend\n";
	$cmd="/storage/apps/SNP_Validation_Scripts/tools/./generateInvariantSites.sh $rounded_Num_of_As_toAppend $rounded_Num_of_Ts_toAppend $rounded_Num_of_Cs_toAppend $rounded_Num_of_Gs_toAppend >invariantSites.txt"; #order of nucleotides as follow: A T C G
system ($cmd); 
	}
	elsif($roundedSum<$diff)
	{
	$toAdd=$diff-$roundedSum;
	print "$roundedSum==$diff ----> rounding is not perfect! Added $toAdd bases to make it equal \n";
	print "A:",($rounded_Num_of_As_toAppend+$toAdd),",T:$rounded_Num_of_Ts_toAppend),C:$rounded_Num_of_Cs_toAppend,G:$rounded_Num_of_Gs_toAppend\n";
	$A_adjusted_added=$rounded_Num_of_As_toAppend+$toAdd;
	print OUT"A:$A_adjusted_added, T:$rounded_Num_of_Ts_toAppend, C:$rounded_Num_of_Cs_toAppend, G:$rounded_Num_of_Gs_toAppend";
	$cmd2="/storage/apps/SNP_Validation_Scripts/tools/./generateInvariantSites.sh $A_adjusted_added $rounded_Num_of_Ts_toAppend $rounded_Num_of_Cs_toAppend $rounded_Num_of_Gs_toAppend >invariantSites.txt"; #order of nucleotides as follow: A T C G
system ($cmd2); 
	}
	elsif($roundedSum>$diff)
	{
	$toSubtract=$roundedSum-$diff;
	print "$roundedSum==$diff ----> rounding is not perfect! Subtracting $toSubtract bases to make it equal \n";
	print "A:",($rounded_Num_of_As_toAppend-$toSubtract),",T:$rounded_Num_of_Ts_toAppend,C:$rounded_Num_of_Cs_toAppend,G:$rounded_Num_of_Gs_toAppend\n";
	$A_adjusted_subtracted=$rounded_Num_of_As_toAppend-$toSubtract;
	print OUT"A:$A_adjusted_subtracted, T:$rounded_Num_of_Ts_toAppend, C:$rounded_Num_of_Cs_toAppend, G:$rounded_Num_of_Gs_toAppend";
	$cmd3="/storage/apps/SNP_Validation_Scripts/tools/./generateInvariantSites.sh $A_adjusted_subtracted $rounded_Num_of_Ts_toAppend $rounded_Num_of_Cs_toAppend $rounded_Num_of_Gs_toAppend >invariantSites.txt"; #order of nucleotides as follow: A T C G

system ($cmd3); 
	}
	last;
	}

}


close $info2;
