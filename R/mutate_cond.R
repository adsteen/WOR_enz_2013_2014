## Conditional version of dplyr::mutate
## @details Written by stackoverflow user G. Grothendieck, https://stackoverflow.com/questions/34096162/dplyr-mutate-replace-on-a-subset-of-rows


mutate_cond <- function(.data, condition, ..., envir = parent.frame()) { 
  
  condition <- eval(substitute(condition), .data, envir)
  .data[condition, ] <- .data[condition, ] %>% dplyr::mutate(...)
  .data
}