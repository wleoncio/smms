get_functions <- function() {
  functions <- list()
  functions_list <- ls(envir = .GlobalEnv)
  for (func in functions_list) {
    if (is.function(get(func))) {
      functions[[func]] <- get(func)
    }
  }
  return(functions)
}
