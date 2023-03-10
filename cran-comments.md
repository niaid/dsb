### Submission March 10 2023  

R CMD check results ───── dsb 1.0.3 ────  
Duration: 1m 15.3s  
0 errors ✓ | 0 warnings ✓ | 0 notes ✓  


devtools::check_win_release()  
* DONE  
Status: OK  

- Fix issue related to a package startup message only on linux distributions; this was not previously a problem. I was contacted by the CRAN team. The startup message is not at all required for functionality; it is now deleted. 

original message  
*It looks like this package (or a package it requires) has a startup
message which cannot be suppressed: see ?packageStartupMessage.*
