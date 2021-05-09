## Reference genome : https://www.ncbi.nlm.nih.gov/nuccore/NC_011071.1

## 

# Usage:
# ./gubbins_part1.sh list_ecoli_CCgroup24_4samples batch2_21032019_ENT1758_assembly.fasta
: <<'END_COMMENT'
list=`readlink -f $1` # convert provided path to FULL path
basename="${1//list_/}" # modify to set folder name (replace list_ in $1 with nothing
reference=$2 #batch2_21032019_ENT1758_assembly.fasta
END_COMMENT

programname=$0
reference=$1

if [ $# == 0 ]; then
    echo "Usage: $0 <reference fasta> [param2]"
    echo "* reference fasta:  A fasta file considered as reference for calculating invariant sites"
fi

#cd /storage/data/DATA4/analysis/23_Stenotrophomonas_maltophilia_2020/12_Gubbins

#rm core_gene_alignment_woConsensus_editedByGATK_renamed.fasta
perl /storage/apps/SNP_Validation_Scripts/tools/convertFastatoOneLine.pl core_gene_alignment_withConsensus_convertXtoN_convertgaptoLittleN_editedByGATK_renamed.fasta

#Assuming consensus is on top remove consensus sequence from the above fasta file
sed '1,2d' core_gene_alignment_withConsensus_convertXtoN_convertgaptoLittleN_editedByGATK_renamed.fasta_OneLine.fasta >>core_gene_alignment_woConsensus_editedByGATK_renamed.fasta

echo "faidx $1"
faidx -o "$1"_faidxOut -i nucleotide $1
head "$1"_faidxOut 

echo "faidx core_gene_alignment_woConsensus_editedByGATK_renamed.fasta:"
faidx -o core_gene_alignment_woConsensus_editedByGATK_renamed_faidxOut -i nucleotide core_gene_alignment_woConsensus_editedByGATK_renamed.fasta
cat core_gene_alignment_woConsensus_editedByGATK_renamed_faidxOut


## calculate_InvariantSites.pl will find how many invariant sites to append and generates invariantSites.txt using generateInvariantSites.sh

perl /storage/apps/SNP_Validation_Scripts/tools/calculate_InvariantSites.pl "$1"_faidxOut core_gene_alignment_woConsensus_editedByGATK_renamed_faidxOut # pass the faidx stats of reference and core gene alignment 

## Append the invariant sites to core gene alignment
perl /storage/apps/SNP_Validation_Scripts/tools/AppendInvariantSites.pl invariantSites.txt core_gene_alignment_woConsensus_editedByGATK_renamed.fasta core_gene_alignment_woConsensus_editedByGATK_renamed

#faidx -i nucleotide core_gene_alignment_woConsensus_editedByGATK_renamed_InvarSitesAppend.fasta

mkdir gubbins
cd gubbins

time run_gubbins.py --prefix postGubbins --filter_percentage 100 --threads 48 ../core_gene_alignment_woConsensus_editedByGATK_renamed_InvarSitesAppend.fasta --verbose

cd ..

echo "faidx after AppendInvariantSites"
faidx -o core_gene_alignment_woConsensus_editedByGATK_renamed_InvarSitesAppend_faidxOut -i nucleotide core_gene_alignment_woConsensus_editedByGATK_renamed_InvarSitesAppend.fasta
head core_gene_alignment_woConsensus_editedByGATK_renamed_InvarSitesAppend_faidxOut
cd ..
