## There is an S3/S4 mess in here since the foreach and iteratos package
## is S3-ized

## iterate over reads on a chromosome
setOldClass('ireads')
ireads <- function(x, ...) UseMethod('ireads')
ireads.BamFile <- function(x, seqs=seqnames(x), unique.only=FALSE, max.mismatch=0L,
                           smooth.by=NULL, ...) {
  seqs <- as.character(seqs)
  all.seqs <- as.character(seqnames(x))
  stopifnot(all(seqs %in% all.seqs))

  i <- 1L
  n <- length(seqs)

  nextEl <- function() {
    if (i > n) {
      stop("StopIteration")
    }
    these <- getReadsFromSequence(x, seqs[i], unique.only=unique.only,
                                  max.mismatch=max.mismatch, smooth.by=smooth.by,
                                  ...)
    i <<- i + 1L
    these
  }

  obj <- list(nextElem=nextEl)
  class(obj) <- c('ireads', 'abstractiter', 'iter')
  obj
}

setOldClass('ibyseqnames')
ibyseqnames <- function(x, ...) UseMethod('ibyseqnames')
ibyseqnames.GRanges <- function(x, seqs=seqnames(x), ...) {
  seqs <- as.character(seqs)
  all.seqs <- as.character(seqnames(x))
  stopifnot(all(seqs %in% all.seqs))

  i <- 1L
  n <- length(seqs)

  nextEl <- function() {
    if (i > n) {
      stop("StopIteration")
    }
    these <- x[seqnames(x) == seqs[i]]
    i <<- i + 1L
    these
  }

  obj <- list(nextElem=nextEl)
  class(obj) <- c('ibyseqnames', 'abstractiter', 'iter')
  obj
}
ibyseqnames.BamFile <- function(x, seqs=seqnames(x), unique.only=FALSE,
                               max.mismatch=0L, smooth.by=NULL, ...) {
  seqs <- as.character(seqs)
  all.seqs <- as.character(seqnames(x))
  stopifnot(all(seqs %in% all.seqs))

  i <- 1L
  n <- length(seqs)

  nextEl <- function() {
    if (i > n) {
      stop("StopIteration")
    }
    these <- getReadsFromSequence(x, seqs[i], unique.only=unique.only,
                                  max.mismatch=max.mismatch, smooth.by=smooth.by,
                                  ...)
    i <<- i + 1L
    these
  }

  obj <- list(nextElem=nextEl)
  class(obj) <- c('ireads', 'abstractiter', 'iter')
  obj
}

setOldClass('ievent')
##' Considers each "tranascription" unit as an event and creates an iterator
##' over them.
ievent <- function(x, ...) UseMethod('ievent')
ievent.BamFile <- function(x, ...) {

}

ievent.GRanges <- function(x, strand='*', min.count=1L, max.mismatch=0L,
                           type='any', ...) {
  if (strand != '*') {
    x <- x[strand(x) == strand]
  }
  covr <- coverage(ranges(x))
  event.bounds <- slice(covr, min.count, rangesOnly=TRUE)
  n <- length(event.bounds)
  if (n == 0L) {
    stop("StopIteration")
  }

  mm <- data.table(matchMatrix(findOverlaps(ranges(x), event.bounds)))
  key(mm) <- 'subject'
  i <- 1L

  nextEl <- function() {
    if (i > n) {
      stop("StopIteration")
    }
    take <- mm[J(i)]$query
    these <- x[take]
    i <<- i + 1L
    these
  }

  obj <- list(nextElem=nextEl, bounds=event.bounds)
  class(obj) <- c('ievent', 'abstractiter', 'iter')
  obj
}


length.ievent <- function(x) length(ranges(x))
setMethod("ranges", c(x="ievent"),
function(x, ...) {
  x$bounds
})
