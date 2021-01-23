dsbmessage <- function()
{
  mesg <-
    cat("dsb package for CITE-seq protein normalization loaded \n",
        "please cite our paper at: https://www.biorxiv.org/content/10.1101/2020.02.24.963603v1 \n",
        "see vignette at https://github.com/niaid/dsb"
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
