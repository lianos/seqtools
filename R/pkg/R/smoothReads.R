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
  ## TODO: Adapt this functionality to "absorb" the
  ##       filter.singular.positions.by.local.density from SPP
  ## NOTE: For removing read pileups, we assume that all reads are
  ##       the same length
  reads[!duplicated(start(reads)) | !duplicated(end(reads))]
}

## given a tag vector (signed), identify and clean up (either remove or cap) singular positions that exceed local tag density
## vin - tag vector
## cap.fold - maximal fold over enrichment over local density allowed for a single tag position, at which the tag count is capped
## eliminate.fold - max fold enrichment that, when exceeded, results in exclusion of all the tags at that position (e.g. counted as anomaly)
## z.threshold - Z-score used to determine max allowed counts
## filter.singular.positions.by.local.density <- function(tags,window.size=200,cap.fold=4,eliminate.fold=10,z.threshold=3) {
##   ## tabulate tag positions
##   if(length(tags) < 2) {
##     return(tags)
##   }

##   tc <- table(tags)
##   pos <- as.numeric(names(tc))
##   ## storage.mode(pos) <- "double";
##   tc <- as.integer(tc)
##   ## storage.mode(tc) <- "integer";

##   n <- length(pos)
##   whs <- as.integer(floor(window.size / 2))

##   ## storage.mode(n) <- storage.mode(whs) <- "integer";
##   twc <- .Call("cwindow_n_tags_around", pos, tc, pos, whs)
##   twc <- (twc - tc + 1) / window.size # local density

##   pv <- pnorm(z.threshold, lower.tail=FALSE)
##   ## exclude
##   max.counts <- qpois(pv, twc * eliminate.fold, lower.tail=FALSE)
##   tc[tc > max.counts] <- 0

##   ## cap
##   max.counts <- qpois(pv, twc * cap.fold, lower.tail=FALSE)
##   ivi <- which(tc > max.counts)
##   tc[ivi] <- max.counts[ivi] + 1L

##   ## reconstruct tag vector
##   tv <- rep(pos, tc)
##   to <- order(abs(tv))
##   tv <- tv[to]
##   return(tv)
## }
