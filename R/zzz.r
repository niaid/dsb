dsbmessage <- function()
{
  mesg <-
    cat("loaded dsb package version", as.character(utils::packageVersion("dsb")),"please cite DOI: 10.1101/2020.02.24.963603")
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

