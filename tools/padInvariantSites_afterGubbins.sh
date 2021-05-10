#This scripts pads the invariant sites after gubbins step

reference=$1

if [ $# == 0 ]; then
    echo "Usage: $0 <reference fasta> [param2]"
    echo "* reference fasta:  A fasta file considered as reference for calculating invariant sites"
fi

echo "faidx $1" # This is reference
faidx -o "$1"_faidxOut -i nucleotide $1
head "$1"_faidxOut 

# Convert gaps to N - we want to retain this information
sed 's/-/N/g' cleancore.full.recomb.masked.fasta >cleancore.full.recomb.masked_convertedGapstoNs.fasta

# Convert ? to - These are recombination regions want to trim of the column using trimal
sed 's/?/-/g' cleancore.full.recomb.masked_convertedGapstoNs.fasta >cleancore.full.recomb.masked_maskrc_convertedToGaps.fasta

echo "Running Trimal..." 
trimal -in cleancore.full.recomb.masked_maskrc_convertedToGaps.fasta -out cleancore.full.recomb.masked.GapsStripped.trimal.fasta -nogaps


echo "faidx cleancore.full.recomb.masked.fasta"
faidx -o cleancore.full.recomb.masked.GapsStripped.trimal.fasta_faidxOut -i nucleotide cleancore.full.recomb.masked.GapsStripped.trimal.fasta
cat cleancore.full.recomb.masked.GapsStripped.trimal.fasta_faidxOut

## calculate_InvariantSites.pl will find how many invariant sites to append and generates invariantSites.txt using generateInvariantSites.sh
echo "Calculating invariant sites..." 
perl /storage/apps/SNP_Validation_Scripts/tools/calculate_InvariantSites.pl "$1"_faidxOut cleancore.full.recomb.masked.GapsStripped.trimal.fasta_faidxOut # pass the faidx stats of reference and core gene alignment 

echo "Converting Trimal multiline fasta to single line fasta ..." 
## format fasta to oneline
fasta_formatter -i cleancore.full.recomb.masked.GapsStripped.trimal.fasta -o cleancore.full.recomb.masked.GapsStripped.trimal.oneline.fasta

echo "Appending invariant sites" 
## Append the invariant sites to core gene alignment
perl /storage/apps/SNP_Validation_Scripts/tools/AppendInvariantSites.pl invariantSites.txt cleancore.full.recomb.masked.GapsStripped.trimal.oneline.fasta cleancore.full.recomb.masked.GapsStripped.trimal.oneline



echo "faidx after AppendInvariantSites"
faidx -o cleancore.full.recomb.masked.GapsStripped.trimal.oneline_InvarSitesAppend.fasta_faidxOut -i nucleotide cleancore.full.recomb.masked.GapsStripped.trimal.oneline_InvarSitesAppend.fasta
head cleancore.full.recomb.masked.GapsStripped.trimal.oneline_InvarSitesAppend.fasta_faidxOut

