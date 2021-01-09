## Resubmission Jan 8 2021 

### R CMD check results ── dsb 0.1.0 ───  

` devtools::check(args = c('--as-cran'))`

Duration: 30.6s  

0 errors ✓ | 0 warnings ✓ | 0 notes ✓  

### Last submission from CRAN team: 

#### License components with restrictions and base license permitting such:  
Question about license emailed directly to cran team (Uwe), please let us know if there are additional questions. The license is a BSD3 with additional language making it even more permissive for 3rd party development and there is a 'for research only' (as opposed to being e.g. for patient care) liability statement, these are required by the US National Institutes of Health. 

#### Found the following (possibly) invalid URLs:  
##### URL: https://github.com/niaid/dsb_manuscript

This invalid link is now removed from the packgae vignette.  


##### URL: https://mattpm.github.io/dsb (moved to https://mattpm.github.io/dsb/)  
This link is necessary because it is the only reference to the package in our preprint - it is active and redirects to the current repository. 



## Resubmission

This is a resubmission. See revisions to each point flagged on initial submission below as well as the updated license to comply with my organization. 

## Test environments
* local OS X install, R 3.6.1
* ubuntu 14.04 (on travis-ci), R 3.6.1
* win-builder (devel and release)

## R CMD check results

── R CMD check results ───────────────────────────────────── dsb 0.1.0 ────
Duration: 28.2s

0 errors ✓ | 0 warnings ✓ | 0 notes ✓


On 12/16/20, 9:58 AM, "Uwe Ligges" <ligges@statistik.tu-dortmund.de> wrote:
 
    Thanks, we see:
 
 
       Found the following (possibly) invalid URLs:
         URL: https://github.com/niaid/dsb_manuscript
           From: inst/doc/dsb_normalizing_CITEseq_data.html
                 README.md
           Status: 404
           Message: Not Found


 
         URL: https://mattpm.github.io/dsb (moved to
    https://mattpm.github.io/dsb/)
           From: inst/doc/dsb_normalizing_CITEseq_data.html
                 README.md
           Status: 200
           Message: OK
 
    Please change http --> https, add trailing slashes, or follow moved
    content as appropriate.
  
 
       The Title field should be in title case. Current version is:
       'Normalize and denoise protein data from droplet-based single cell
    profiling (CITE-seq, REAP-seq, Mission Bio Tapestri)'
       In title case that is:
       'Normalize and Denoise Protein Data from Droplet-Based Single Cell
    Profiling (CITE-Seq, REAP-Seq, Mission Bio Tapestri)'
 
       The Description field should not start with the package name,
         'This package' or similar.
 
 
       The Description field contains
         https://www.biorxiv.org/content/10.1101/2020.02.24.963603v1.
       Please enclose URLs in angle brackets (<...>).
 
    Please fix and resubmit.
 
    Best,
    Uwe Ligges

# Revisions for CRAN team: 

Dear Uwe and CRAN team, 

Thanks for your service to the community. Changes to the repository are below: 

1) 404 status of link: https://github.com/niaid/dsb_manuscript
- I removed this unneeded link to another repository that is not yet public.  

2) 200 status of link: https://mattpm.github.io/dsb
-  I would like to keep this link (Status: OK). It redirects to the current repository. This is the only link to the software in our public preprint on BioRxiv. The beta release of this software had a .io website–I removed the website and set up the redirection link after transferring the repository to our organization page github.com/niaid/.  

3) The Title field should be in title case.

Changed title to tile case: 'Normalize and Denoise Protein Data from Droplet-Based Single Cell Profiling (CITE-Seq, REAP-Seq, Mission Bio Tapestri)'
 
4) The Description field should not start with the package name,'This package' or similar.

- Changed to "This lightweight R package provides ..." 
 
5) The Description field contains https://www.biorxiv.org/content/10.1101/2020.02.24.963603v1. Please enclose URLs in angle brackets (<...>).

- enclosed link in brackets. 

*Addition*
Please note I also was notified that my institute uses a templated BSD-3 license. I added file LICENSE and  changed the license field of DESCRIPTION to 'BSD_3_clause + file LICENSE' following the CRAN guideline for BSD-3. https://cran.r-project.org/web/licenses/

*note*
checked win.builder.r.org log files and confirmed potentially misspelled words are correctly spelled. 


```{r}
devtools::check(args = "--as-cran") 
```
#returns ...
── R CMD check results ──────────────── dsb 0.1.0 ────
Duration: 27.2s

0 errors ✓ | 0 warnings ✓ | 0 notes ✓


Have a great holiday. 

-MPM


