### Submission March 12 2022 

This is a resubmission fixing these two issues below:

*issue 1*  
*Found the following (possibly) invalid URLs: URL: https://mattpm.github.io/dsb/ From: README.md Status: 404 Message: Not Found*  

**I removed the link from the R package so it does not cause any issues in the future. I also fixed the link, https://mattpm.github.io/dsb/ now redirects properly. I noticed some users still use this old link to find the github, but it did not need to be listed in the package.**

*issue 2*  
*Found the following (possibly) invalid DOIs: DOI: doi.org/10.1101/2020.02.24.963603 From: inst/CITATION Message: Invalid DOI indeed, the doi is only 10.1101/2020.02.24.963603*  

**I fixed the citation file as above to only list the DOI number**


#### devtools::check(cran = TRUE)  
Duration: 1m 9.3s  
0 errors ✓ | 0 warnings ✓ | 0 notes ✓  

#### devtools::check_win_release()  
Status: 1 NOTE  
License components with restrictions and base license permitting such:  
  BSD_3_clause + file LICENSE  
**The license is unchanged from the previous accepted version of the package.**  
**The license is been properly specified after guidance from CRAN submission team as in dsb v0.2.0**

