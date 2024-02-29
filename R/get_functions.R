#' @title Get functions from the global environment
#' @description Retrieves functions defined in the global environment and
#' returns them to the target environment.
#' @param target_env The environment to which the functions will be returned.
#' @return All functions from the global environment.
#' @author Waldir Leoncio
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
