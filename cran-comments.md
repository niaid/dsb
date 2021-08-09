## Submission Aug 9 2021 

This is a resubmission -- please read correspondence about license below. I have specified the license as instructed by the CRAN team from my last CRAN release; *the license is unchanged from the last release*. 

Re License components with restrictions and base license permitting such: "Please only ship the CRAN template for the BSD_3_claus license." SEE BELOW.

    ########################### 
    On 1/24/21, 7:55 PM, "Mule, Matthew (NIH/NIAID) [F]" <matthew.mule AT nih.gov> wrote:

    1) re: "Please only ship the CRAN template for the BSD_3_clause license."
    
    Following instructions sent via email 9 Jan, 2021 below to add the CRAN template for BSD-3 which I have added "in addition to the
    alternate BSD text" as written below. Please do let me know if I misunderstood? 
       
    On 1/9/21, 8:53 AM, "Uwe Ligges" <ligges AT statistik.tu-dortmund.de> wrote:

    I understood now, but then you need:

    License: BSD_3_clause + file LICENSE | file LICENSE

    and you need in file LICENSE the CRAN template for the BSD_3_clause licene and in addition the alterative BSD text that is
    currently listed there anyway.

    ###########################

## Submission Aug 7 2021 

Update submission for version 0.2.0

1 NOTE - I (creator / maintainor) changed my email address from: <matthew.mule at nih.gov> to: <mattmule at gmail.com> The latter is a stable long term email address, the former email address is not consistently accessible and will be invalid after 2023. 

## R CMD check results
devtools::check(args = c('--as-cran'))  
0 errors ✓ | 0 warnings ✓ | 0 notes ✓

## Test environments
* Winbuilder
  devtools::check_win_release()
 * x86_64-w64-mingw32 (64-bit)

* devtools::check_rhub(pkg = ".")  
PLEASE NOTE [Bioconductor does not yet build and check packages for R version 4.2](https://github.com/r-hub/rhub/issues/471)
  * Windows Server 2008 R2 SP1, R-devel, 32/64 bit (see note above) 
  * Ubuntu Linux 16.04 LTS, R-release, GCC
  * Fedora Linux, R-devel, clang, gfortran

re: custom license used in this package; this is the same license as in CRAN v0.1.0 Please refer to Cran team comments from 2021/01/24. Following instructions sent via email 9 Jan, 2021, the CRAN template for BSD-3 is listed "in addition to the alternate BSD text". License documents are referenced in the description as instructed (see email below on specification of license file).

   ########################### 
   On 1/9/21, 8:53 AM, "Uwe Ligges" <ligges<AT>statistik.tu-dortmund.de> wrote:

       I understood now, but then you need:

       License: BSD_3_clause + file LICENSE | file LICENSE

       and you need in file LICENSE the CRAN template for the BSD_3_clause 
       licene and in addition the alterative BSD text that is currently listed 
       there anyway.
   ###########################


  
