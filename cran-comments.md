## Resubmission Jan 23 2021 

This is a resubmission. 

## Test environments
* Winbuilder
  devtools::check_win_release()
 * x86_64-w64-mingw32 (64-bit)
  
* devtools::check_rhub(pkg = ".")
  * Windows Server 2008 R2 SP1, R-devel, 32/64 bit
  * Ubuntu Linux 16.04 LTS, R-release, GCC
  * Fedora Linux, R-devel, clang, gfortran

## R CMD check results
devtools::check(args = c('--as-cran'))  
0 errors ✓ | 0 warnings ✓ | 0 notes ✓  


### Notes on Windows and R builds:   

*possibly mis-spelled words in DESCRIPTION*  
- All are author names or scientific terms. No misspelled words found.  

*The Title field should be in title case. Current version is:*
*'Normalize & Denoise Droplet Single Cell Protein Data (CITE-seq)'*
*In title case that is:*
*'Normalize & Denoise Droplet Single Cell Protein Data (CITE-Seq)'*

- "CITE-seq" is the correct styling of the term; the second half of the name 'seq' is in lower case. ref: https://cite-seq.com/.   

*License components with restrictions and base license permitting such*  
- Please see last submission and emails with Uwe Ligges; this is a standard BSD3 clause license with additional language required by US Government stating it is for research use only and not for patient care.  The license is specified according to guidance from the last submission:  
- BSD_3_clause + file LICENSE | file LICENSE  

## CRAN comments from previous submission

*Please reduce the length of the title to less than 65 characters.*  
- reduced title in description file to < 65 characters.  
 
*Please write references in the description of the DESCRIPTION file in the form authors (year) <doi:...>*   
- changed reference to match styling above.  

*Please always write package names, software names and API ... in single quotes*  
- added single quotes around these elements in the description.  

*Please always explain all acronyms in the description text.*  
- Added acronym explanation for example Unique Molecular Index (UMI). Some acronyms are sequencing methods only referred to by the acronym name in the genomics field-writing what each letter stands may add more confusion than clarity. Checked the CRAN package cinaR for ATAC-seq data analysis and matched the styling it uses by adding single quotes to the method names 'CITE-seq' and 'REAP-seq' matching how the cinaR package does not spell out the acronym for 'ATAC-seq' and adds it in quotes as shown.  

*Please add \value to .Rd files regarding exported methods... write about the structure of the output (class)*  
- added this and small executable example for the package function with note on the output structure and class (R matrix).  

*If a function does not return a value, please document ... Missing Rd-tags: pipe.Rd*    
- Added a value tag for the tidyverse pipe operator '%>%' imported with usethis::use_pipe() which does not return a value.   

*Please add small executable examples in your Rd-files to illustrate the use of the exported function but also enable automatic testing.*  
- Added small executable examples to the DSBNormalizeProtein function. Remaining Rd files are for the example data and pipe operator (see above).  Package also utilizes unit testing for this function in tests/testthat/test_dsb_function.R.  

Thank you! 
Matt  




