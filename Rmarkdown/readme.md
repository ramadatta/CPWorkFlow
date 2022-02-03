bbstats.md - 
# DAT-A-DAT-V - DATa Analysis and DATa Visualization

All the data analysis steps and visualizations are converted into markdown documents and compiled into markdown library

#### **Synopsis**

To keep track of analysis creating multiple markdown documents and preserving the code and steps for the future use.
 
#### Contents 
 
|SNo| Document      | Description   | ChartType |
|---| ------------- |:-------------| :-----|
|1| [bbstats.md](https://github.com/ramadatta/CPWorkFlow/blob/main/Rmarkdown/bbstats.md) | Data Visualization of Pre and Post Trimmmed Nanopore reads: Statistics from bbtools statswrapper | Dumbbell |
|2| [readlenDist.md](https://github.com/ramadatta/CPWorkFlow/blob/main/Rmarkdown/readlenDist.md) | Data Visualization of Pre and Post Trimmmed Nanopore reads: Read length distribution | Split Violin|
|3| [cgeResistome.md](https://github.com/ramadatta/CPWorkFlow/blob/main/Rmarkdown/cgeResistome.md)| Data Visualization of CGE Resistome  |Barplot,Heatmap |
|4| [Cmp_RawFilt_UnicyclerAssemblies.md](https://github.com/ramadatta/CPWorkFlow/blob/main/Rmarkdown/Cmp_RawFilt_UnicyclerAssemblies.md) | Data Visualization of Raw and Filtered Unicycler Assemblies | Barplot in Facetwrap |
|5| [plotMob.md](https://github.com/ramadatta/CPWorkFlow/blob/main/Rmarkdown/plotMob.md) | Data Visualization of Replicons from Mobtyper output | Barplot in Facetwrap |
|6| [proLolli.md](https://github.com/ramadatta/CPWorkFlow/blob/main/Rmarkdown/proLolli.md) | Data Visualization of SNPs across the Genome Sequence | Lollipop Plot |
|7| [Comp_shortReads.md](https://github.com/ramadatta/CPWorkFlow/blob/main/Rmarkdown/Comp_shortReads.md) | Comparison of 2 samples at different time points by resistome generated from Raw reads and Assembly | PCA Plot |
|8| [beasttmrca.md](https://github.com/ramadatta/CPWorkFlow/blob/main/Rmarkdown/cgeResistome.md)| Data Visualization of BEAST Phylotree - in progress | tree |
|9| [SNPDiffmatrix.md](https://github.com/ramadatta/CPWorkFlow/blob/main/Rmarkdown/cgeResistome.md)| Data Visualization of SNP Difference Matrix - in progress | Heatmap|
|10| [datesPlot.md](https://github.com/ramadatta/CPWorkFlow/blob/main/Rmarkdown/cgeResistome.md)| Data Visualization of dates between different institutions - in progress | scatterplot |
|11| [SamplingSite.md](https://github.com/ramadatta/CPWorkFlow/blob/main/Rmarkdown/cgeResistome.md)| Data Visualization of sampling sites between different institutions - in progress| Barplot|
|12| [WafflePlot.md](https://github.com/ramadatta/CPWorkFlow/blob/main/Rmarkdown/cgeResistome.md)| Data Visualization of how many samples have a gene present and not - in progress  | Waffle|
|13| [mutRate_BoxPlot.md](https://github.com/ramadatta/tute/blob/main/R/boxplot/mutRate_BoxPlot.md)| Boxplot from pre-calculated values  | Boxplot | 
|14| [SNPAAS.md](https://github.com/ramadatta/tute/blob/main/SNPs/SNPAAS/SNPAAS.md)| Find and Plot SNPs Across All Samples | HeatMap | 
|15| [multiDumbell.md](https://github.com/ramadatta/CPWorkFlow/blob/main/Rmarkdown/multiDumbell.md)| Adding multiple timepoint on dumbbell plot | Dumbell Plot | 



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

* Without the bbstats_files folder we cannot view images since images are stored in that folder.


