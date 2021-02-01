## Resubmission Feb 1 2021 

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


### Cran team comments from 2021/01/24

#### 1) re: "Please only ship the CRAN template for the BSD_3_clause license."  
- Please see last submission and emails with Uwe Ligges. Following instructions sent via email 9 Jan, 2021, I added the CRAN template for BSD-3 and "in addition to the alternate BSD text" as written below. It looks like I have specified license documents in the description exactly as shown below. I have also directly followed the template for CRAN BSD3 clause template  exactly as specified Please do let me know if I misunderstood? https://cran.r-project.org/web/licenses/BSD_3_clause

       ########################### 
       On 1/9/21, 8:53 AM, "Uwe Ligges" <ligges@statistik.tu-dortmund.de> wrote:

           I understood now, but then you need:

           License: BSD_3_clause + file LICENSE | file LICENSE

           and you need in file LICENSE the CRAN template for the BSD_3_clause 
           licene and in addition the alterative BSD text that is currently listed 
           there anyway.
       ###########################

#### 2) The Title field should be in title case. Current version is:  
'Normalize & Denoise Droplet Single Cell Protein Data (CITE-seq)'
In title case that is:
'Normalize & Denoise Droplet Single Cell Protein Data (CITE-Seq)'

Please refer to prev. submitted cran-comments.md. "CITE-seq" with a lower case 's' after the dash in CITE-seq' is the correct styling of the method; please refer to <https://cite-seq.com/>

Thank you! 
-Matt 
