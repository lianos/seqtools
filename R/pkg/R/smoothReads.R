##' Iterates over smother options to clean reads from certain "artifacts"
##'
##' @param reads A GRanges object of short reads
##' @param smooth.by A vector of functions (or names of functions) to use as
##' smoothers.
smoothReads <- function(reads, smooth.by=NULL, ...) {
  verbose <- checkVerbose(...)
  smooth.by <- unique(smooth.by[nchar(smooth.by) > 0L])
  if (!is.null(smooth.by) && length(smooth.by) > 0) {
    if (verbose) {
      cat("smoothing:", paste(smooth.by, collapse=","), "\n")
    }
    for (smoother in smooth.by) {
      if (is.character(smoother)) {
        if (substring(smoother, 1, 9) != 'smoother.') {
          smoother <- paste('smoother', smoother, sep='.')
        }
        smoother <- getFunction(smoother, where=.GlobalEnv)
      }
      if (is.function(smoother)) {
        reads <- smoother(reads, ...)
      } else {
        warning("Unknown smoother `", as.character(smoother), "` ... skipping")
      }
    }
  }
  reads
}

smoother.readsPileup <- function(reads, ...) {
  ## NOTE: For removing read pileups, we assume that all reads are
  ##       the same length
  ## NOTE: For removing pileups -- you might want to get samtools
  ##       to do this for you.
  reads[!duplicated(start(reads)) | !duplicated(end(reads))]
}
