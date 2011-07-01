setGeneric("transform", function(`_data`, ...) standardGeneric("transform"),
           useAsDefault=function(`_data`, ...) base::transform)
##' transform function for DataFrame objects that behaves like the one for
##' data.frames
setMethod("transform", "DataFrame", function(`_data`, ...) {
  e <- eval(substitute(list(...)), `_data`, parent.frame())
  tags <- names(e)
  inx <- match(tags, names(`_data`))
  matched <- !is.na(inx)
  if (any(matched)) {
    for (cname in tags[matched]) {
      `_data`[[cname]] <- e[[cname]]
    }
  }
  if (!all(matched)) {
    `_data` <- cbind(`_data`, as(e[!matched], 'DataFrame'))
  }
  `_data`
})

##' @nord
checkVerbose <- function(...) {
  verbose <- list(...)$verbose
  if (is.null(verbose)) verbose <- options()$verbose
  verbose
}

checkTrace <- function(...) {
  do.trace <- list(...)$trace
  if (!is.logical(do.trace)) {
    do.trace <- FALSE
  }
  do.trace
}

##' @nord
dir.exists <- function(path) {
  path <- as.character(path)
  !is.na(file.info(path)$isdir) && file.info(path)$isdir
}

##' Convert NA values in vectors and data.frames to a default value
##'
##' @export
##'
##' @param wut The thing to convert
##' @param to The value convert NA to
##' @return The same type as \code{wut}
convert.na <- function(wut, to=".defaults.") {
  if (is.character(to) && to[1] == ".defaults.") {
    to <- list(logical=FALSE, numeric=0, integer=0L, character="",
               factor="")
  }
  if (is.vector(wut) || is.factor(wut)) {
    wut.type <- is(wut)[1]
    if (is.list(to)) {
      if (!wut.type %in% names(to)) {
        stop("Unknown default conversion value for", wut.type, sep=" ")
      }
      to <- to[[wut.type]]
    }
    if (wut.type == 'factor') {
      levels(wut) <- c(levels(wut), to)
    }
    wut[is.na(wut)] <- to
  } else if (inherits(wut, 'data.frame') || inherits(wut, 'DataFrame')) {
    cols <- 1:ncol(wut)
    if (is(wut, 'data.table')) {
      ## Don't change key columns
      cols <- setdiff(cols, which(colnames(wut) %in% key(wut)))
    }
    for (idx in cols) {
      wut[[idx]] <- convert.na(wut[[idx]], to=to)
    }
  }

  wut
}

