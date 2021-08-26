dsbmessage <- function()
{
  mesg <-
    cat("dsb package for ADT normalization loaded \n",
        "please cite our paper (under review) at: https://www.biorxiv.org/content/10.1101/2020.02.24.963603v3 \n",
        "questions: https://github.com/niaid/dsb Vignette: https://cran.r-project.org/package=dsb"
    )
  return(mesg)
}

.onAttach <- function(lib, pkg)
{
  # startup message
  msg <- dsbmessage()
  if(!interactive())
  msg[1] <- paste("Package 'dsb'")
  packageStartupMessage(msg)
  invisible()
}
