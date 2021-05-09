#!/usr/bin/bash

list=`readlink -f $1` # convert provided path to FULL path
basename="${1//list_/}" # modify to set folder name (replace list_ in $1 with nothing

echo "$list"
echo $basename

source /storage/apps/anaconda3/etc/profile.d/conda.sh

echo "formatting the alignment fasta (core_gene_alignment.aln) output from roary"
mkdir $basename
fasta_formatter -i core_gene_alignment.aln -o "$basename"/core_gene_alignment_OneLine.fasta -w 0

echo "Entering the SNP validation folder"

echo "Creating, Entering the $basename folder"

cd $basename

echo "Extracting the Core Gene Alignment if there are multiple groups"

for d in `cat $list`; do grep -F -A1 $d core_gene_alignment_OneLine.fasta >>core_gene_alignment_extracted.fasta; done;

conda deactivate # consensus.py and vcf-subset need Python < 3

echo "Creating consensus of the alignment"

(echo ">consensus"; python /storage/apps/SNP_Validation_Scripts/tools/consensus.py core_gene_alignment_extracted.fasta 0 ) >>core_gene_alignment_consensus.fasta

(cat core_gene_alignment_consensus.fasta; cat core_gene_alignment_extracted.fasta) >>core_gene_alignment_withConsensus.fasta

echo "Creating, Entering the snp-sites folder"
mkdir snpsites
cd snpsites

echo "Finding SNPs using snp-sites on the consensus fasta"

snp-sites -mvp -o core_gene_alignment_withConsensus_snpsites ../core_gene_alignment_withConsensus.fasta

echo "Subsetting VCF file without reference sequence"

for sample in `bcftools query -l core_gene_alignment_withConsensus_snpsites.vcf`; do vcf-subset --exclude-ref -c $sample core_gene_alignment_withConsensus_snpsites.vcf > snpsites_${sample}.vcf; done;

conda activate
cd ..

#cd ..

#END

#cd $basename
mkdir snippy
cd snippy


#Edit snippy to use a different vcf-consensus (provide full path)
#make sure the reference has X converted to N
sed 's/X/N/g' ../core_gene_alignment_consensus.fasta >>../core_gene_alignment_consensus_convertXtoN.fasta	


# For ENT samples
for d in `cat $list`; do snippy --cpus 32 --outdir $d --prefix $d --ref ../core_gene_alignment_consensus_convertXtoN.fasta --R1 /storage/data/DATA4/analysis/23_Stenotrophomonas_maltophilia_NovogeneAIT_2020/2_AdapterTrimmed_bbduk_Q30/$d/"$d"*_L2_1_bbmap_adaptertrimmed.fastq --R2 /storage/data/DATA4/analysis/23_Stenotrophomonas_maltophilia_NovogeneAIT_2020/2_AdapterTrimmed_bbduk_Q30/$d/"$d"*_L2_2_bbmap_adaptertrimmed.fastq; done;

cd ../snpsites

#for bcftools isec
for d in `ls -1 snpsites*.vcf`; do sed -i 's/^1/consensus/' $d; done;
for d in `ls -1 snpsites*.vcf`; do bgzip -c $d > $d.gz; done;
for d in `ls -1 snpsites*.gz`; do tabix -p vcf $d; done;

# Get vcf intersect between two snp calling methods (snp-sites based on core genome from prokka-roary and snippy output)

for d in `cat $list`; do bcftools isec -p $d -Ov snpsites_"$d"*.vcf.gz ../snippy/$d/$d.filt.subs.vcf.gz; done;

# This generates 0000.vcf and 0001.vcf...so, 0000.vcf contains unvalidated snps

#for d in `cat /storage/data/DATA4/analysis/23_Stenotrophomonas_maltophilia_2020/10_Roary/list_smaltophilia_group1_21samples`; do bcftools isec -p $d -Ov snpsites_"$d"*.vcf.gz ../snippy/$d/$d.filt.subs.vcf.gz; done;

#END

#cd $basename
#cd snpsites

sed 's/X/N/g' ../core_gene_alignment_withConsensus.fasta >>../core_gene_alignment_withConsensus_convertXtoN.fasta
sed '/^>/! s/-/n/g' ../core_gene_alignment_withConsensus_convertXtoN.fasta >>../core_gene_alignment_withConsensus_convertXtoN_convertgaptoLittleN.fasta

# Creating a coordinate file for GATK

for d in `cat $list`; do grep -v "#" $d/0000.vcf | sed "s/consensus/$d/" >>combined_forGATK.vcf; done;
#for d in `cat /storage/data/DATA4/analysis/23_Stenotrophomonas_maltophilia_2020/10_Roary/list_smaltophilia_group1_21samples`; do grep -v "#" $d/0000.vcf | sed "s/consensus/$d/" >>combined_forGATK.vcf; done;


## Preparing files to run GATk - which masks SNPs those that are not validated.


# (echo "##fileformat=VCFv4.1"; echo -e '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tDUMMY'; awk -v 'OFS=\t' '{$4=$5 ; $5="N" ; print ;}' combined_forGATK.vcf) >>combined_forGATK_reformat.vcf

# Set the ref as a dummy "A" because if set as $4=$5, $5 might have multiple alleles and this causes issues with indexing

(echo "##fileformat=VCFv4.1"; echo -e '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tDUMMY'; awk -v 'OFS=\t' '{$4="A" ; $5="N" ; print ;}' combined_forGATK.vcf) >>combined_forGATK_reformat.vcf  ## for JGH, this was output as type2

sed "s/\*,//" combined_forGATK_reformat.vcf | sed "s/,\*//" >>combined_forGATK_reformat_noAsterix.vcf

gatk IndexFeatureFile -F combined_forGATK_reformat_noAsterix.vcf

cd ..


#(optional) Sometimes the fasta headers dont match between the fasta and vcf. Do a rename after backing up original file:
cp core_gene_alignment_withConsensus_convertXtoN_convertgaptoLittleN.fasta core_gene_alignment_withConsensus_convertXtoN_convertgaptoLittleN.fasta_beforeRename
#sed -i 's/_cutadapt_spades_contigs\.fasta//' core_gene_alignment_withConsensus_convertXtoN_convertgaptoLittleN.fasta
sed -i 's/_spades_contigs\.fasta//' core_gene_alignment_withConsensus_convertXtoN_convertgaptoLittleN.fasta
#(end optional)

samtools faidx core_gene_alignment_withConsensus_convertXtoN_convertgaptoLittleN.fasta

gatk CreateSequenceDictionary -R core_gene_alignment_withConsensus_convertXtoN_convertgaptoLittleN.fasta

## Converts the alternate to Ns

gatk FastaAlternateReferenceMaker -R core_gene_alignment_withConsensus_convertXtoN_convertgaptoLittleN.fasta -O core_gene_alignment_withConsensus_convertXtoN_convertgaptoLittleN_editedByGATK.fasta -V snpsites/combined_forGATK_reformat_noAsterix.vcf

sed 's/^>\S*\s*/>/' core_gene_alignment_withConsensus_convertXtoN_convertgaptoLittleN_editedByGATK.fasta >>core_gene_alignment_withConsensus_convertXtoN_convertgaptoLittleN_editedByGATK_renamed.fasta

# So, in the final finle below, we have an alignment which contains snps which are validated both by core genome approach as we all as mapping based snippy. The unvalidated snps from core genome approach are converted to N to mask the polymorphism.

sed -i 's/:.*$//' core_gene_alignment_withConsensus_convertXtoN_convertgaptoLittleN_editedByGATK_renamed.fasta

cd ..
<<COMMENT1
COMMENT1
