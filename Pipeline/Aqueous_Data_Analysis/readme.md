## Molecular epidemiology and genomic analysis of individual patient isolates and environmental sampling reveals transmission of bacterial plasmids harboring carbapenem resistant genes

### Abstract

#### Objectives
The incidence of carbapenem-resistant Enterobacterales (CRE) is found to be increasing since year 2010 in Singapore. To this end, we have used Whole Genome Sequencing (WGS) to characterize and analyze the chromosomal and plasmid molecular epidemiology of CRE from patients and hospital environment in Singapore.

#### Methods
A total of Seventy Enterobacterales isolates (patient n=12; environmental=58;) were sequenced for short-reads from Illumina and long-reads from Oxford Nanopore GridION platforms. The resulting sequence data was assembled and multilocus sequence types (MLST), species assignment, resistome & virulome analysis and plasmid replicon types and phylogenetic relationships between the samples was determined.

#### Results
Majority of the CPE isolates (61.4%; 43/70) were identified as Serratia marcescens followed by Klebsiella Pneumonia (24.2%; 17/70), Enterobacter cloacae (7.14%; 5/70), Escherichia coli (5.7%; 4/70), Klebsiella oxytoca (1.42%; 1/70). blaOXA-48 was found to be the most abundant carbapenemase gene from all the isolates (60/70). Co-carriage of multiple CP genes in samples was not observed. Except for the Enterobacter Species, the blaOXA-48 gene was pre-dominantly found to be present on a 67-kb plasmid in Serratia and Klebsiella species. Overall, X transmission clusters were identified from N=70 isolates. 

#### Conclusions
Combination of short-read and long-read sequencing can effectively help to determine the transmission of carbapenemase gene and in characterizing the specific plasmid carrying the carbapenemase gene. The current study of the isolates further our understanding CPE in the hospital setting. Also, clear associations between patient and environmental isolates were found by studying the transmission of carbapenem resistant plasmids. This reveals that the surrounding environment of patient has certainly plays in CPE transmission dynamics.

[RefPaper1](https://sci-hub.se/https://doi.org/10.1093/jac/dkaa398) , [RefPaper2](https://sci-hub.se/https://doi.org/10.1093/jac/dkz146)

### Introduction:
* Most important problem in healthcare setting
* Which order of bacteria is causing it? What makes it more deadly and powerful to attack the humans?
* CPE - mostly which organisms cause the frequent attacks?
* What causes the Carbapenemase to move? Transmissible mobile genetic elements
* Which sepcies in this study we are concentrating about?
* Are there any prevalent ST clone in our dataset?
* Treatment of antibiotics
* Acquisition of resistance genes through mobile genetic elements
    - specific genomic islands 
    - transposons
    - plasmids
    - integrative and conjugative elements
    - integrons
* Severe infection due to virulence factors
* Toxins of the type III secretion system to circumvent the host immune system and establish infection
* Of the four exotoxins (ExoS, ExoT, ExoU, ExoY), ExoU, a potent phospholipase that disrupts the plasma membrane and leads to rapid cell death, is the most virulent 

### Methods

#### Whole Genome Sequencing and Genome assembly
Sequencing of the 70 isolates was performed using Illumina and Oxford Nanopore GridION platforms. De novo hybrid genome assemblies were generated using the short reads and long reads as an input to Unicycler assembler (Ref). The fastq files are deposited in NCBI SRA database with Bioproject accession (PRJNAXXX).

#### Sequence Typing, Species Identification and Resistome
Multilocus Sequencing Type (MLST) were identified using MLST tool (Ref). Species assignment was done based on Kraken2. The antimicrobial resistant genes were identified using resfinder available CGE website (Ref) and virulent genes by using ABRicate tool with vfdb as database (Ref). Plasmid replicons were assigned using Mob_suite tool (Ref)

#### Bacterial core genome analysis 
For each species, core genome analysis of the bacterial isolates was performed to find the possible transmission clusters. At first the genomes were annotated using prokka (Ref). Then, Roary (Page et al., 2015) was used to generate core genome alignment and SNPs were extracted from the core genome alignment using SNP-sites (https://github.com/sanger-pathogens/snp-sites). For the validation of SNPs generated from SNP-sites tool, another set of SNP sites were generated from mapping reads against the consensus genome sequence using SNIPPY tool. Only SNPs found in both SNP sets were considered as genuine SNPs and the unvalidated SNPs from core genome approach are converted to N to mask the polymorphism. Recombination filtering was performed using Gubbins and the resulting polymorphic snp sites output from gubbins was used to calculate SNP difference between the isolate pairs.

#### Phylogenetic analysis and mutation rate calculation
For the species E.coli and Klebsiella Pneumoniae, mutation rate previously inferred by our group in one of the previous study published (ref) was used. For the Serratia species, to estimate substitution rate, we used SNP alignment generated after eliminating recombined regions by gubbins was imported into BEAUti v2.5 to create an Extensive Markup Language (xml) file for BEAST v2.5 [57]. Using the S. marcescens Db11 (NCBI accession HG326223.1) reference genome the proportion of nucleotide bases of A, T, G, C were calculated and constant sites of A:1040770, T:1029706, C:1522977, G:1520300 were added in xml tag to keep the model representative of S. marcescens Db11 complete genome. We used Strict molecular clock as clock model, bModelTest [58] was used as a substitution model. Three independent runs of BEAST, each with .xml file was given as an input to BEAST to run 10 million MCMC cycles. Mixing and convergence statistics were monitored in Tracer v1.6 (Rambaut et al., 2014). Subsequently, convergence was determined for the given model by checking if the ESS values > 200 for all the parameters. A burn-in of 10 million states was left out from each of the three independent runs of this model. We then combined the results from those runs with the logcombiner program from the BEAST package. We took the age of the root as the height of the root of the Maximum Clade Credibility (MCC) tree reconstructed by combining trees using the tree annotator program from the BEAST package.

### Results

* **Bacterial isolates**
  - From what time to what time these isolates have been collected?
  - Out of the collected, how many have carbapenemases?
  - Where are isolates collected from patients? Which body sites? If environmental, which surrounding sites? Give the count of each site.
  - Out of this, what bacterial species was found and the count of each species.
  - Is there any trend of increasing in the number of isolates by species by year?

  
* **Location of CP genes**
  - Is it chromsomosomal or plasmids?
  - If plasmids what is the size of the plasmids

* **WGS analysis**
   - Genome Characteristics: The average genome size of the de novo assembled genomes for each species - Information about all genome characteristics is summarized in Table. 
   - Resistome analysis - what genes are conferring resistance to aminoglycosides (aacA4 and aadA1), macrolides (mphE and msrE), chloramphenicol (catA1 and catB2), sulphonamides (sul1) and trimethoprim/sulfamethoxazole (dfrB1) in each of the species isolates. In OXA-48-CP isolates, only one gene involved in aminoglycoside antimicrobial resistance was found (aacA4) (Figure S3). Show the facet grid or PCA plot for the resistome
   - Type of plasmid reconstructions showed in each of the 4 species - what is the size of the plasmid with CP gene ? Is it available in NCBI? If not, what is the coverage and identity? Is this the same plasmid in the other species as well carrying the CP gene? What are other plasmids in other species? and their details?
   - Phylogenetic analysis of all the genomes for each species 
   - Typing, plasmid, and genomic transmission analysis - What does MLST analysis reveal? How many groups of STs identified for each species? 
   - Plasmid replicon types were detected, and associations with carbapenemase genes were assessed - Among OXA-48-producing isolates, Inc-- replicon plasmids were common in Serratia marcescens (85.7%) and K. pneumoniae (61.1%), IncL/M plasmids were more common among E. coli and Enterobacter spp. (66.7% and 50.0%, respectively), and IncHI2 plasmids were most common among E. cloacae complex (61.5%). All three types have previously been described to carry the IMP-4 gene in Australia (9, 34). KPC-2-producing K. pneumoniae commonly carried IncFIB (94.6%), ColRNAI (90.2%), and IncX3 (85.7%). In contrast, NDM and OXA isolates carried a wide range of plasmid types. Plasmid replicon sequences were only rarely detected on the same contig as a carbapenemase gene. Fifteen isolates had no plasmid replicon sequences detected at these identity/coverage cutoffs, including all IMI- and SME-producing isolates (chromosomal carbapenemases).
   - Genomic transmission analyses were performed on 33 subgroups of the same species, carbapenemase gene, and ST (where available), including 131 isolates from 119 patients; subgroups included 2 to 14 isolates (median, 3 isolates) (Table S6). Pairwise SNP distances were plotted for each species
   - Plasmid replicon typing results were consistent with transmission analysis results in epidemiologically confirmed clusters, in that all isolates within a cluster carried the same plasmid replicon types, with the exception of the mixed NDM-5/NDM-5 plus OXA-232 cluster (K. pneumoniae ST16), where only the isolates with OXA-232 had the ColKP3-type plasmid detected, suggesting carriage of the OXA-232 gene on this plasmid.
   -  Not much diversity of resistance genes in ST308 isolates - Spread of both of these lineages are more or less heavily dependent on these resistant genes

### Discussion
   - OXA-48-producing Enterobacterales isolates have been reported from hospital and extra-hospital reservoirs. Most of the studies on
this subject have been conducted in developed countries, but the major epicentres of OXA-48-producing Enterobacterales are located in Mediterranean countries such as Algeria, where this resistance mechanism is endemic.
   - Therefore our study encompasses OXA-48 isolates - Plasmid mediated transmission - usual species but serratia - associated with specific plasmids and replicons
   - Along these lines, the presence of additional antimicrobial co-resistance genes could have played a key role in the selection and persistence of this IncL-pVIM-1 variant over the conserved IncL-pOXA-48.

### Conclusions
 
 - In conclusion, how many lineages combined with versions of the same plasmid carrying XX cp genes or OXA-48 have contributed to the dispersion of CPSm in our institution.
 - Moreover, our results also indicate that the IncL-pOXA-48 plasmid has the ability not only to transfer among different Enterobacterales species, including
S. marcescens, but also to exchange genetic material. This fact facilitates the acquisition of other genetic structures carrying multiple resistance genes, including other carbapenemases such as blaVIM-1, thus generating nosocomial pathogens with a higher MDR profile. The transmission of these plasmids among different clinical pathogens including S. marcescens clones with the potential capacity to spread in hospital units with a high risk of infection,
such as ICUs, represents a new challenge and should be considered a public health priority. 
This study reinforces the relevance of species different from K. pneumoniae or E. coli in the CPE landscape and circulating lineages and plasmids in the local
CPE epidemiology.
