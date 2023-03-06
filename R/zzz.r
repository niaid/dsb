# dsbmessage <- function()
# {
#   mesg <-
#     cat("loaded dsb package version",
#         as.character(utils::packageVersion("dsb")),
#         "please cite https://www.nature.com/articles/s41467-022-29356-8"
#     )
#   return(mesg)
# }

.onAttach <- function(lib, pkg)
{
  # startup message
  #msg <- dsbmessage()

  packageStartupMessage(cat("loaded dsb package version",
                            as.character(utils::packageVersion("dsb")),
                            "please cite https://www.nature.com/articles/s41467-022-29356-8"
  ))
  invisible()
}
