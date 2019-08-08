# Pull our slope. Could use broom; but I don't feel like it
get_slope <- function(m) {
  slope <- tryCatch(
    slope <- coefficients(m)[2],
    error = function(x) NA,
    finally = {}
  )
}
