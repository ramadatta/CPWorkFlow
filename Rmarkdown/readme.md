bbstats.md - 
# DAT-A-DAT-V - DATa Analysis and DATa Visualization

All the data analysis steps and visualizations are converted into markdown documents and compiled into markdown library

#### **Synopsis**

To keep track of analysis creating multiple markdown documents and preserving the code and steps for the future use.
 
#### Contents 
 
  * [bbstats.md](https://github.com/ramadatta/CPWorkFlow/blob/main/Rmarkdown/bbstats.md) - Data Visualization of Pre and Post Trimmmed Nanopore reads: Statistics from bbtools statswrapper





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


