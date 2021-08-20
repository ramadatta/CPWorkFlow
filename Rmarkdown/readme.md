bbstats.md - 
# DAT-A-DAT-V - DATa Analysis and DATa Visualization

All the data analysis steps and visualizations are converted into markdown documents and compiled into markdown library

#### **Synopsis**

To keep track of analysis creating multiple markdown documents and preserving the code and steps for the future use.
 
#### Contents 
 
| Document      | Description   | ChartType |
| ------------- |:-------------| :-----|
| [bbstats.md](https://github.com/ramadatta/CPWorkFlow/blob/main/Rmarkdown/bbstats.md) | Data Visualization of Pre and Post Trimmmed Nanopore reads: Statistics from bbtools statswrapper | Dumbbell |
| [readlenDist.md](https://github.com/ramadatta/CPWorkFlow/blob/main/Rmarkdown/readlenDist.md) | Data Visualization of Pre and Post Trimmmed Nanopore reads: Read length distribution | Split Violin|
| [cgeResistome.md](https://github.com/ramadatta/CPWorkFlow/blob/main/Rmarkdown/cgeResistome.md)| Data Visualization of CGE Resistome  |Barplot,Heatmap |
| [Cmp_RawFilt_UnicyclerAssemblies.md](https://github.com/ramadatta/CPWorkFlow/blob/main/Rmarkdown/Cmp_RawFilt_UnicyclerAssemblies.md) | Data Visualization of Raw and Filtered Unicycler Assemblies | Barplot in Facetwrap |
| [plotMob.md](https://github.com/ramadatta/CPWorkFlow/blob/main/Rmarkdown/plotMob.md) | Data Visualization of Replicons from Mobtyper output | Barplot in Facetwrap |
| [proLolli.md](https://github.com/ramadatta/CPWorkFlow/blob/main/Rmarkdown/proLolli.md) | Data Visualization of SNPs across the Genome Sequence | Lollipop Plot |
proLolli.md
| [beasttmrca.md](https://github.com/ramadatta/CPWorkFlow/blob/main/Rmarkdown/cgeResistome.md)| Data Visualization of BEAST Phylotree - in progress | tree |
| [SNPDiffmatrix.md](https://github.com/ramadatta/CPWorkFlow/blob/main/Rmarkdown/cgeResistome.md)| Data Visualization of SNP Difference Matrix - in progress | Heatmap|
| [datesPlot.md](https://github.com/ramadatta/CPWorkFlow/blob/main/Rmarkdown/cgeResistome.md)| Data Visualization of dates between different institutions - in progress | scatterplot |
| [SamplingSite.md](https://github.com/ramadatta/CPWorkFlow/blob/main/Rmarkdown/cgeResistome.md)| Data Visualization of sampling sites between different institutions - in progress| Barplot|
| [WafflePlot.md](https://github.com/ramadatta/CPWorkFlow/blob/main/Rmarkdown/cgeResistome.md)| Data Visualization of how many samples have a gene present and not - in progress  | Waffle|


#### *Notes for authors*

* After knitting the document -  we will successfully generate R Markdown document (Rmd). But this Rmd file cannot be used in github. 

* So we convert Rmd to md by 

``` yaml
    ---
    title: "my project title"
    author: "Rama"
    date: "`r format(Sys.Date())`"
    output:
      html_document:
        keep_md: true
    ---
```

* Then we upload this md document along with the files generated in the folder. Here in this case, for example, I simply dragged and dropped into the browser 1) markdown document 2) bbstats_files folder.

* Without the bbstats_files folder we can view images since images are stored in that folder.


