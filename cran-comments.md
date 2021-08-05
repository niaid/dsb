## Submission Aug 1 2021 

Update submission for version 0.2.0

## Test environments
* Winbuilder
  devtools::check_win_release()
 * x86_64-w64-mingw32 (64-bit)
  
* devtools::check_rhub(pkg = ".")  
Bioconductor does not yet build and check packages for R version 4.2 see https://github.com/r-hub/rhub/issues/471  
  * Windows Server 2008 R2 SP1, R-devel, 32/64 bit
  * Ubuntu Linux 16.04 LTS, R-release, GCC
  * Fedora Linux, R-devel, clang, gfortran

## R CMD check results
devtools::check(args = c('--as-cran'))  
0 errors ✓ | 0 warnings ✓ | 0 notes ✓

**Re. custom license used in this package**
Please refer to Cran team comments from 2021/01/24

"Please only ship the CRAN template for the BSD_3_clause license."  
Following instructions sent via email 9 Jan, 2021, the CRAN template for BSD-3 (https://cran.r-project.org/web/licenses/BSD_3_clause) is listed "in addition to the alternate BSD text". License documents are referenced in the description as instructed (see email below).

       ########################### 
       On 1/9/21, 8:53 AM, "Uwe Ligges" <ligges<AT>statistik.tu-dortmund.de> wrote:

           I understood now, but then you need:

           License: BSD_3_clause + file LICENSE | file LICENSE
  
           and you need in file LICENSE the CRAN template for the BSD_3_clause 
           licene and in addition the alterative BSD text that is currently listed 
           there anyway.
       ###########################
  

