
## To draw venn diagram of Transposons, Integrons, AMR genes, Virulence genes from chromsomes and plasmids across multiple species

cd /data02/Analysis/Projects/2_CPE_Transmission/VennDiagram_of_genes_between_species/CPE_Transmission_Data_Analysis


#---------- Identity plasmids and chromosomes using PlasClass

time for j in $(ls *.fasta); do echo $j; time python3.7 ~/sw/PlasClass/classify_fasta.py -f $j -o "$j".plasclass.probs.out -p 60;  done
real	464m57.242s
user	592m31.525s
sys	133m29.767s

mkdir ../1_plasClass_Probs
mv *.plasclass.probs.out ../1_plasClass_Probs
$ awk '{print FILENAME "\t" $0}' *.out >CPE_Trans_PlasClass_Predictions_comb.tab



#---------- Identity Transposons using Bacant (Integron module is not working properly and Resistance genes are used from AMRfinderplus)
time for d in $(ls *.fasta);  do  
echo "$d";
bacant -n "$d" -o "$d"_bacant_out &
done

mkdir ../2_bacant_annot
mv *_bacant_out ../2_bacant_annot
awk '{$13=$14=""; print FILENAME "\t" $1 "\t" $2}' *_bacant_out/transposon.filter.tsv | sed 's/  */\t/g'  | fgrep -v 'QUERY' | sed 's/_bacant_out\/transposon.filter.tsv//g' | awk '{print $1 "\t" $2 "\t" $1"#"$2 "\t" $3}' >CPE_Trans_BacAnt_Transposons.tab


#---------- Identity AMR genes using NCBI AMRfinderplus

# NOTE: Ran this earlier, so just copying the results

mkdir ../3_amrfinderplus_results 

#---------- Identity Virulence factor genes using Abricate

for d in $(ls *.fasta); do echo "abricate -db vfdb $d >"$d".vfdb.tab"; done >abricate_cmds.sh
time parallel --jobs 62 < abricate_cmds.sh &

mkdir ../4_abricate_vf_results 
mv *.vfdb.tab ../4_abricate_vf_results
cd ../4_abricate_vf_results


## Excluding possible contamination samples
-rwxrwxr-x 1 prakki prakki 9.8M Jul  6 14:17 batch12_27072019_ENT1410_assembly_renamed.fasta
-rwxrwxr-x 1 prakki prakki  11M Jul  6 14:17 batch2_21032019_ENT466_assembly_renamed.fasta
-rwxrwxr-x 1 prakki prakki  11M Jul  6 14:17 batch2_21032019_ENT497_assembly_renamed.fasta
-rwxrwxr-x 1 prakki prakki  11M Jul  6 14:18 batch2_21032019_ENT378_assembly_renamed.fasta
-rwxrwxr-x 1 prakki prakki  11M Jul  6 14:18 rbatch2_07082020_ENT792_assembly_renamed.fasta

# These assemblies have abnormal assembly sizes. So excluding them from analysis

## -- Reunning MLST 

for d in $(ls *.fasta); do echo "mlst $d >"$d".mlst.tab"; done >mlst_cmds.sh
time parallel --jobs 62 < mlst_cmds.sh &





























mkdir /home/prakki/anaconda3/pkgs/abricate-1.0.1-ha8f3691_1/db/Bacant_Integron_FinderDB
cp /home/prakki/anaconda3/pkgs/abricate-1.0.1-ha8f3691_1/db/Bacant_Integron_FinderDB/
sudo sed 's/>/>Integron_ABRicateDB~~~/g' Integron.fasta | sed 's/|/~~~/g' | awk -F\~ '{print $0"\~\~\~"$NF}' >sequences
sudo cp sequences /home/prakki/anaconda3/pkgs/abricate-1.0.1-ha8f3691_1/db/Bacant_Integron_FinderDB/
sudo makeblastdb -in sequences -title Bacant_Integron_FinderDB -dbtype nucl -hash_index

% abricate --list
DATABASE  SEQUENCES  DBTYPE  DATE
Bacant_Integron_FinderDB  1094       nucl    2022-Sep-5


% abricate --db tinyamr screen_this.fasta

