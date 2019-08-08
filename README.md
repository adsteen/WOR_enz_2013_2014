# Code and data for Steen et al (2019) AEM

This repo contains the code and data required to reproduce [Steen et al (2019) Applied and Environmental Microbiology, DOI: 10.1128/AEM.00102-19](https://aem.asm.org/content/early/2019/07/15/AEM.00102-19). 

To use, clone this project, and create an R project in the directory you have cloned.

Install the following packages:

* tidyverse
* broom
* scales
* cowplot
* plyr
* grid
* gridExtra
* scales

Note that most of the code is written using tidyverse, but some is old code written using plyr. Be careful with this.

In principle, if you run `write_paper.R`, the entire analysis will run on its own. Plots will be printed to the screen but not saved. 

In practice, I recommend running this line-by-line, to see where the errors pop up. 

Note that lines 15-17 of `write_paper.R` set whether plots (and a few tables) will be printed to screen, whether `deep cut' plots (i.e., those that are informative but which didn't end up as figures) will be printed to screen, and whether plots will be saved. Set these appropriately before running.

Do not hesitate to [contact me](asteen1@utk.edu) with questions. 
