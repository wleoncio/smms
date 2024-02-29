get_functions <- function(target_env) {
  functions <- list()
  functions_list <- ls(envir = .GlobalEnv)
  for (func in functions_list) {
    if (is.function(get(func))) {
      functions[[func]] <- get(func)
    }
  }
  list2env(functions, envir = target_env)
}
