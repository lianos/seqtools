## Sequences in BAM file are given in (+) orientation -- even if the sequence
## is aligned to the (-) strand.

all.bam.what <- c("qname", "flag", "rname", "strand", "pos", "qwidth",
                  "mapq",  "cigar","mrnm",  "mpos",   "isize", "seq",
                  "qual")

.default.bam.what <- c('qname', 'strand', 'pos', 'qwidth', 'flag', 'seq',
                       'qual')
.default.bam.what <- all.bam.what

##' "Handy" function to query a bam file -- dumps in default values for
##' \code{scanBamParams}, etc.
##'
##' @exportMethod query
##' @rdname queryBAM-methods
##' @author Steve Lianoglou \email{slianoglou@@gmail.com}
##'
##' @param x A \code{\linkS4class{SeqStore}} object, or a path to a bam file
##' on the filesystem
##' @param what The elements to return for the reads, defaults to
##' \code{'qname', 'strand', 'pos', 'qwidth', 'flag', 'seq', 'qual'}
##' @param which The \code{which} parameter to the
##' \code{\link{ScanBamParam}} function.
##' @param flag The \code{flag} parameter for \code{\link{scasnBamFlag}},
##' defaults to an empty value.
##' @param tag The vars to pass to \code{\link{bamTag}}. Defaults to
##' \code{NM, MD}.
##' @param param A \code{\linkS4class{ScasnBamParam}} object
##' @param max.mismatch Integer indicating the maximum number of mismatches
##' an alignment can have to be returned from the bam query. Defaults to
##' no filtering of reads.
##'
##' @return A list of results, as returned by \code{\link{scanBam}}
setGeneric("query",
function(x, what=.default.bam.what, which=NULL, flag=scanBamFlag(),
         tag=c('NM', 'MD'), param=NULL, max.mismatch=NULL, ...) {
  standardGeneric("query")
})

##' @rdname queryBAM-methods
##' @param unique.only Query the uniquely aligned bam file for the seqstore
setMethod("query", c(x="BamFile"),
function(x, what, which=NULL, flag, tag, param, max.mismatch,
         unique.only=FALSE, ...) {
  if (!isOpen(x)) {
    open(x)
    on.exit(close(x))
  }
  .aligner <- aligner(x)

  if (is.null(param)) {
    param <- ScanBamParam(which=which, flag=flag, what=what)
  }
  if (!is.null(tag) && length(tag) > 0) {
    bamTag(param) <- unique(c(bamTag(param), tag))
  }
  if (is.numeric(max.mismatch)) {
    bamTag(param) <- unique(c('NM', 'CM', bamTag(param)))
  }

  if (.aligner == 'bwa') {
    ## X0 : number of best hits
    ## X1 : number of suboptimal hits
    ## XT : Type: Unique/Repeat/N/Mate-sw
    ## XA : Alternative hits: (chr,pos,CIGAR,NM;)*
    bamTag(param) <- c('Z0', 'X0', 'X1', 'XT', 'XA', bamTag(param))
  }

  bamTag(param) <- unique(c('Z0', bamTag(param)))
  result <- scanBam(x, param=param)
  if (.aligner == 'bwa' && unique.only) {
    result <- filterScanBamByFlag(result, filterBwaUnique, ...)
  }

  if (is.numeric(max.mismatch)) {
    result <- filterScanBamByFlag(result, filterByMismatchTag,
                                  max.mismatch=max.mismatch, ...)
  }

  result
})

filterScanBamByFlag <- function(result, filterf, ...) {
  verbose <- checkVerbose(...)
  if (!is.logical(verbose)) {
    verbose <- options()$verbose
  }
  lapply(result, function(x) {
    if (any(sapply(x, length) == 0L) || !('tag' %in% names(x))) {
      ## Defend against an element of the result being empty
      return(x)
    }
    keep <- filterf(x$tag, ...)
    if (verbose) {
      cat("Removing", 1 - (sum(keep) / length(keep)), "% of reads\n")
    }
    cleaned <- lapply(x[names(x) != 'tag'], '[', keep)
    cleaned$tag <- lapply(x$tag, '[', keep)
    cleaned
  })
}

filterBwaUnique <- function(tag.list, ...) {
  ## Z0 is "my custom" uniqueness flag
  ## flag <- if ('Z0' %in% names(tag.list)) tag.list$Z0 else tag.list$X0
  ## Hack fo rnow -- Z0 tag isn't consistent during testing
  if (!is.null(tag.list$Z0)) {
    tag <- ifelse(is.na(tag.list$Z0), tag.list$X0, tag.list$Z0)
  } else {
    tag <- tag.list$X0
  }
  tag == 1L
}

filterByMismatchTag <- function(tag.list, max.mismatch=1L, ...) {
  ## If BWA was used to align in color space, the mismatch tag seems to be CM
  ## otherwise the mismatch tag is NM
  mm.columns <- c('CM', 'NM')
  use <- sapply(mm.columns, function(x) !is.null(tag.list[[x]]))
  tag.name <- mm.columns[use]
  if (length(tag.name) == 0L) {
    stop("No mismatch information found in BAM file, skpping mismatch filter")
  }
  if (length(tag.name) == 2L) {
    warning("Both CM and NM found in BAM file for mismatch -- using NM")
    tag.name <- 'NM'
  }

  tag.list[[tag.name]] <= max.mismatch
}

.filterByMismatch <- function(result, param, max.mismatch) {
  lapply(result, function(x) {
    keep <- x$tag$NM <= max.mismatch
    cat("... cleaning", length(x$tag$NM) - sum(keep), "mismatched reads\n")
    cleaned <- lapply(x[-which(names(x) == 'tag')], '[', keep)
    cleaned$tag <- lapply(x$tag, '[', keep)
    cleaned
  })
}



