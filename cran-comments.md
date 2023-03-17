### Submission March 17 2023  

R CMD check results ───── dsb 1.0.3 ────  
Duration: 1m 15.3s  
0 errors ✓ | 0 warnings ✓ | 0 notes ✓  

devtools::check_win_release()  
* DONE  
Status: 1 Note  

Package was archived on CRAN  
X-CRAN-Comment: Archived on 2023-03-16 as issues were not corrected in time.  


Found the following (possibly) invalid URLs:  
**These links will be valid as soon as the package is back on CRAN. These links should remain unchanged. These are the correctly formatted links as required by CRAN on previous releases of this package.**   
https://CRAN.R-project.org/package=dsb/vignettes/additional_topics.html  
https://CRAN.R-project.org/package=dsb/vignettes/end_to_end_workflow.html  
https://CRAN.R-project.org/package=dsb/vignettes/no_empty_drops.html  
https://CRAN.R-project.org/package=dsb/vignettes/understanding_dsb.html  



"Possibly misspelled words in DESCRIPTION:"  
**These are not misspelled and are unchanged from the previous CRAN releases of this package.**
  ADT (21:259)  
  Denoise (3:20)  
  Tsang (21:1161)  
  UMI (21:193)  
  denoising (21:79)  
  dsb (21:1041, 21:1093)  
  isotype (21:698)  


Previous submission:  
- Fix issue related to a package startup message only on linux distributions; this was not previously a problem. I was contacted by the CRAN team. The startup message is not at all required for functionality; it is now deleted. 

original message  
*It looks like this package (or a package it requires) has a startup
message which cannot be suppressed: see ?packageStartupMessage.*
