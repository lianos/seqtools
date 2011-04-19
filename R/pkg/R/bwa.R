toBwaMultimapString <- function(granges, nm.tag='tag.NM') {
  ## (chr,pos,CIGAR,NM;)*
  nm <- values(granges)[[nm.tag]]
  if (is.null(nm.tag)) {
    rep(NA_integer_, length(granges))
  }
  nm <- as.integer(nm)

  cigar <- values(granges)$cigar
  if (is.null(cigar)) {
    cigar <- paste(width(granges), 'M', sep="N")
  }
  cigar <- as.character(cigar)

  x <- sprintf('%s,%s%d,%s,%d', as.character(seqnames(granges)),
               as.character(strand(granges)), start(granges), cigar, nm)
  x
}

##' Given a vector of 'multimap strings', returns a GRanges object with
##' the places the read can map to.
##'
##' The \code{values() DataFrame} includes a \code{reference} column which
##' indicates the index of the input \code{mm.string} that the range maps
##' back to.
##'
##' Elements in the input string can be NA, indicating that there is no multimap
##' information. In these cases, the GRanges element in the returned object
##' will simply not have a "back reference" in \code{values()$reference}
##'
##' BWA multimap strings are of the form: \code{(chr,pos,CIGAR,NM;)*}
##'
##' @importFrom ShortRead cigarToQWidth
##'
##' @param mm.string A characteor vector of multimap strings. Each element
##' can have \code{NA}, one, or more multmiapped destinations.
##'
##' @return A \code{GRanges} object with the expanded multimapped regions.
##'
##' Use the \code{values(...)$reference} information to figure out which
##' element in \code{mm.string} each range belongs to.
##'
##' Columns for the edit distance and cigar string for the reads are also
##' included in the \code{elementMetadata} \code{DataFrame}
##' (\code{edit.distance}, and \code{cigar}, respectively).
setGeneric("convertBwaMultimap",
function(x, keep=c('all', 'best'), from=c('sam', 'bam'), to=c('sam', 'bam'), ...) {
  standardGeneric("convertBwaMultimap")
})

setMethod("convertBwaMultimap", c(x="GRanges"),
function(x, keep, from, to,  multimap.column='tag.XA',
         meta.keep=c('qual', 'qname', 'seq'), ...) {
  mm.string <- values(x)[[multimap.column]]
  if (is.null(mm.string)) {
    stop("No remapping information found in ", multimap.column)
  }
  expanded <- convertBwaMultimap(mm.string, keep, from, to, ...)
  if (length(expanded) > 0L && !is.null(meta.keep)) {
    meta <- values(x)
    meta.keep <- meta.keep[meta.keep %in% colnames(meta)]
    xref <- values(expanded)$reference
    for (col in meta.keep) {
      values(expanded)[[col]] <- meta[[col]][xref]
    }
  }
  expanded
})

setMethod("convertBwaMultimap", c(x="character"),
function(x, keep, from, to, ...) {
  if (length(x) == 0L) {
    return(GRanges())
  }
  from <- match.arg(from)
  to <- match.arg(to)
  keep <- match.arg(keep)
  alignments <- strsplit(x, ';', fixed=TRUE)
  n.locs <- sapply(alignments, length)
  n.locs[is.na(alignments)] <- 0L
  pieces <- strsplit(unlist(alignments[n.locs > 0]), ',', fixed=TRUE)
  chrs <- sapply(pieces, '[', 1)
  edit.distance <- as.integer(sapply(pieces, '[', 4))
  strand.start <- sapply(pieces, '[', 2)
  strands <- substring(strand.start, 1, 1)
  starts <- as.integer(substring(strand.start, 2))
  cigar <- sapply(pieces, '[', 3)
  widths <- cigarToQWidth(cigar)
  reference <- rep(seq(n.locs), n.locs)

  gr <- GRanges(chrs, IRanges(starts, width=widths), strands,
                edit.distance=edit.distance, cigar=cigar,
                reference=reference)

  if (keep == 'best') {
    ## Only keep the remappings that have the lowest edit.distance
    is.best <- lapply(split(edit.distance, reference), function(x) x == min(x))
    is.best <- unlist(is.best)
    gr <- gr[is.best]
  }

  values(gr)$flag <- as.integer(ifelse(as.logical(strand(gr) == '-'), 16L, 0L))

  if (from != to) {
    ## SAM and BAM indexing are 1 and 0 based, respectively
    if (from == 'sam' && to == 'bam') {
      ## SAM is 1 based and BAM is 0 based!
      ranges(gr) <- shift(ranges(gr), -1L)
    } else if (from == 'bam' && to == 'sam') {
      ranges(gr) <- shift(ranges(gr), 1L)
    }
  }

  gr
})
