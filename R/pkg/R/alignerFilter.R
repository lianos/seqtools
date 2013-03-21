# XxxFilter functions are created as a list of tests to be run against a
# GRanges-like (eg. GappedAlignments) object. Each filter should return a
# logical vector of equal length as the GRanges object which is TRUE
# for the reads to be kept. This `filters` list is then passed to the
# `functionalizeFilters` to collapse these filters into each other by creating
# a function that ANDs them together
functionalizeFilters <- function(filters) {
  if (length(filters) == 0) {
    f <- function(x) x
  } else {
    f <- function(x) {
      keep <- Reduce("&", lapply(filters[-1], function(f) f(x)),
                     filters[[1]](x))
      x[keep]
    }
  }
  f
}

bwaSamseFilter <- function(unique.only=FALSE) {
  filters <- list()

  if (unique.only) {
    filters[[length(filters) + 1L]] <- function(x) {
      stopifnot('X0' %in% names(mcols(x)))
      mcols(x)$X0 == 1
    }
  }

  functionalizeFilters(filters)
}
