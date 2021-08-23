## This script generates R code for Lollipop plot

cat *snps.filt.vcf | grep -v "#" | grep "TYPE=snp" | awk '{print $2}' | sort | uniq -c | awk '{print $2 "\t" $1}' >snpCount_multiVCF.list

cat *snps.filt.vcf | grep -v "#" | grep "TYPE=snp" | awk '{print $1 "\t" $2}' | sort -u >multiVCF_snps_position.list

awk 'BEGIN{FS=OFS="\t"}{if(NR>1){ print $1,$2-1,$2 }}' multiVCF_snps_position.list | bedtools intersect -wa -wb -a SHG10_pKPC2-km.gff -b stdin | sed 's/ /_/g' | awk '{print $4"\t"$5}'  | sort -nu | awk '{print $1}' | tr '\n' ',' | sed 's/,$//g' >rIranges.int1

awk 'BEGIN{FS=OFS="\t"}{if(NR>1){ print $1,$2-1,$2 }}' multiVCF_snps_position.list | bedtools intersect -wa -wb -a SHG10_pKPC2-km.gff -b stdin | sed 's/ /_/g' | awk '{print $4"\t"$5}'  | sort -nu | awk '{print $2-$1}' | tr '\n' ',' | sed 's/,$//g' >rWidth.int2

awk 'BEGIN{FS=OFS="\t"}{if(NR>1){ print $1,$2-1,$2 }}'  multiVCF_snps_position.list | bedtools intersect -wa -wb -a SHG10_pKPC2-km.gff -b stdin | sed 's/ /_/g' | awk '{print $4"\t"$5"\t"$9}'  | sort -nu | sed -e 's/ID.*;gene=//g' -e 's/ID.*;product=//g' -e 's/;inference=.*//g' -e "s/'//g"  -e "s/[()]/_/g"  -e 's/_$//g' | awk '{print "\""$3"\""}' >rGenes.int3

ncol=`wc -l rGenes.int3 | awk '{print $1}'`;

awk 'BEGIN{FS=OFS="\t"}{if(NR>1){ print $1,$2-1,$2 }}'  multiVCF_snps_position.list | bedtools intersect -wa -wb -a SHG10_pKPC2-km.gff -b stdin | sed 's/ /_/g' | awk '{print $4"\t"$5"\t"$9}'  | sort -nu | sed -e 's/ID.*;gene=//g' -e 's/ID.*;product=//g' -e 's/;inference=.*//g' -e "s/'//g"  -e "s/[()]/_/g"  -e 's/_$//g' | awk '{print "\""$3"\""}' | tr '\n' ',' | sed 's/,$//g' >rGenes.int3


# To declare static Array 
c25=("maroon" "black" "gold1" "brown" "#CAB2D6" "steelblue4" "dodgerblue2" "yellow4" "skyblue2" "green4" "deeppink1" "#FF7F00" "#6A3D9A" "#E31A1C" "blue1" "#FB9A99" "palegreen2" "#FDBF6F" "gray70" "green1" "yellow3" "khaki2" "darkorange4" "orchid1" "darkturquoise")

#echo "number of colors needed are $ncol";

for i in $(seq 1 $ncol)
do
  echo -n "\"${c25[$i]}\", ";
done | sed 's/, $//g' >rColors.int4

#echo "number of gene which need height info is $ncol";

for i in $(seq 1 $ncol)
do
  echo -n "0.08,";
done | sed 's/,$//g' >rHeight.int5

#SNP

awk '{print $2}' multiVCF_snps_position.list | tr '\n' ',' | sed 's/,$//g' >rSNPs.int6

# Create a Rscript
echo "library(trackViewer)"
echo "";

echo "snpCount_df <- read.table(\"snpCount_multiVCF.list\",sep = \"\t\", header = TRUE)";
echo "colnames(snpCount_df) <- c(\"SNP\",\"Count\")";

echo "";

echo -n "features <- GRanges(\"chr1\", IRanges(c("
cat rIranges.int1
echo -n "), width=c("
cat rWidth.int2
echo -n "), names=c("
cat rGenes.int3
echo -n "), fill=c("
cat rColors.int4
echo -n "), height = c("
cat rHeight.int5
echo ")))"
echo "";
echo -n "SNP <- c("
cat rSNPs.int6
echo ")"
echo "";

echo "sample.gr <- GRanges(\"chr1\", IRanges(SNP, width=1, names=paste0(\"snp\", SNP)),
                     color = sample.int(6, length(SNP), replace=TRUE),
                     score = snpCount_df\$Count
                     )
lolliplot(sample.gr, features)"

