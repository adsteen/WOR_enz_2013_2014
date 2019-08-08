# Safely get linear models

data_lm <- function(x) {
  m <- tryCatch(
    m <- lm(fl ~ time, data = x),
    error = function(x) NA,
    finally = {}
  )
  m
}
