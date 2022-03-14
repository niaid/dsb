### Submission March 14 2022 

This is a patch that fixes a markdown html rendering issue.  

The html file for the main vignette in the last CRAN release did not render; the old version of the vignette is currently shown which will be confusing for users. This issue was not apparent in my local builds. As suggested in R dev forums this was caused by an updated knitr header in the vignette--I deleted / renamed the file and changed the DESCRIPTION `VignetteBuilder: knitr, rmarkdown` and this fixed the issue.  

I downloaded the binaries from winbuilder and confirmed that the vignette files now render properly.

I also changed the names of the vignettes in `VignetteIndexEntry` in the R markdown headers so that they match the title of the vignettes, which was my intention but I did not properly specify them before.  


Duration: 1m 12.5s
0 errors ✓ | 0 warnings ✓ | 0 notes ✓  

#### devtools::check_win_release()  
Status: 1 NOTE  
License components with restrictions and base license permitting such:  
  BSD_3_clause + file LICENSE  
**The license is unchanged from the previous accepted version of the package.**  
**The license is been properly specified after guidance from CRAN submission team as in dsb v0.2.0**
