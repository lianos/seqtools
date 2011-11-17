## Fetch reads straight from a BAM file -- no SeqStore or BSgenome.*
## necessary
##
## Some stranded protocols read the library from 3' -> 5' direction.
## If this is the case, set flip.reads=TRUE
setGeneric("getReadsFromSequence",
function(x, seqname, start=NULL, end=NULL, strand=NULL, unique.only=FALSE,
         smooth.by=NULL, with.sequence=FALSE, with.quality=FALSE,
         meta.what=c("flag"), flip.reads=FALSE, ...) {
  standardGeneric("getReadsFromSequence")
})

##' @nord
.noReads <- function() {
  reads <- GRanges(seqnames=character(), ranges=IRanges(), strand=strand())
  values(reads) <- DataFrame(id=integer(), pair.id=integer())
  reads
}

setMethod("getReadsFromSequence", c(x="BamFile"),
function(x, seqname, start, end, strand, unique.only, smooth.by, with.sequence,
         with.quality, meta.what, flip.reads, ...) {
  args <- list(...)
  verbose <- checkVerbose(...)
  trace <- checkTrace(...)
  if (is.null(flip.reads)) {
    flip.reads <- FALSE
  }

  if (!isOpen(x)) {
    open(x)
    on.exit(close(x))
  }

  if (trace) cat("getReadsFromSequence")

  start.time <- as.numeric(Sys.time())

  .seqlengths <- seqlengths(x)

  if (!seqname %in% names(.seqlengths)) {
    stop("Unknown sequence name: ", seqname)
  }

  ## Determine WHICH reads to return
  if (is.null(start) || length(start) == 0L) {
    start <- 1
  }
  if (is.null(end) || length(end) == 0L) {
    end <- .seqlengths[seqname]
  }

  which <- RangesList(IRanges(start=start, end=end))
  names(which) <- seqname

  ## Determine FLAG
  if (!is.null(strand)) {
    if (flip.reads) {
      strand <- swapStrand(strand)
    }
    if (strand %in% c('-', -1)) {
      flag <- scanBamFlag(isMinusStrand=TRUE)
    } else {
      flag <- scanBamFlag(isMinusStrand=FALSE)
    }
  } else {
    flag <- scanBamFlag()
  }

  unique.bam.query <- unique.only && is.null(args$uniqueness.map)
  bam.tag <- unique(c('NM', 'MD', args$tag))

  result <- query(x, what=meta.what, which=which, flag=flag,
                  max.mismatch=args$max.mismatch,
                  unique.only=unique.bam.query, tag=bam.tag)[[1]]

  if (length(result$pos) == 0) {
    return(.noReads())
  }

  ## It seems as if result$pos is.na when the alignments falls off the end of
  ## the chromosome
  end.around <- which(is.na(result$pos))

  if (length(end.around) > 0) {
    is.tag <- which(names(result) == 'tag')
    if (length(is.tag) > 0) {
      tag <- result[[is.tag]]
      result <- result[-is.tag]
    }
    result <- lapply(result, '[', -end.around)
    if (length(is.tag) > 0) {
      tag <- lapply(tag, '[', -end.around)
      result$tag <- tag
    }
  }

  ## id <- result$qname
  if (isPaired(x)) {
    pair.id <- BAMpairId(result$flag, check.paired=FALSE)
  } else {
    pair.id <- rep(1L, length(result$pos))
  }

  strands <- result$strand
  if (flip.reads) {
    strands <- swapStrand(strands)
  }

  reads <- GRanges(seqnames=seqname,
                   ranges=IRanges(start=result$pos, width=result$qwidth),
                   strand=strands, pair.id=pair.id)
                   
  ## Adding more metadata to the reads from BAM file, as requested by caller
  ## (did they want the sequence, too?)
  if (!is.null(meta.what)) {
    if (with.sequence) {
      meta.what <- unique(c(meta.what, 'seq'))
    }
    if (with.quality) {
      meta.what <- unique(c(meta.what, 'qual'))
    }
  } else {
    meta.what <- names(result)
  }

  dont.add <- c('strand', 'tag', 'rname', 'pos', 'qwidth')
  
  meta.what <- meta.what[!meta.what %in% dont.add]
  for (name in meta.what) {
    values(reads)[[name]] <- result[[name]]
  }

  if (flip.reads && ('flag' %in% colnames(values(reads)))) {
    values(reads)$flag <- bitXor(values(reads)$flag, 16L)
  }

  for (tag.name in names(result$tag)) {
    values(reads)[[paste('tag', tag.name, sep=".")]] <- result$tag[[tag.name]]
  }

  ## Apply smoothers to data
  if (!is.null(smooth.by) && length(smooth.by) > 0L) {
    reads <- smoothReads(reads, smooth.by=smooth.by, ...)
  }

  if (verbose) {
    cat("  SeqTools::getReadsFromSequence took",
        as.numeric(Sys.time()) - start.time, "seconds\n")
  }

  reads
})

setMethod("getReadsFromSequence", c(x="character"),
function(x, seqname, start=NULL, end=NULL, strand=NULL, unique.only=TRUE,
         smooth.by=NULL, with.sequence=FALSE, with.quality, meta.what='flag',
         ...) {
  getReadsFromSequence(BamFile(x), seqname, start, end, strand, unique.only,
                       smooth.by, meta.what=meta.what, ...)
})

##' Returns all of the reads from a bamfile
##'
##' @param x BamFile, or path to one
##' @param which.seqnames A character vector indicating which seqnamnes to get
##' reads from, \code{NULL} for all.
setGeneric("getAllReads", function(x, which.seqnames=NULL, ...) {
  standardGeneric("getAllReads")
})

setMethod("getAllReads", c(x="character"),
function(x, which.seqnames, ...) {
  getAllReads(BamFile(x), which.seqnames, ...)
})

setMethod("getAllReads", c(x="BamFile"),
function(x, which.seqnames, ...) {
  if (is.null(which.seqnames)) {
    which.seqnames <- seqlevels(x)
  }
  if (!isOpen(x)) {
    open(x)
    on.exit(close(x))
  }
  all.reads <- lapply(which.seqnames, function(chr) {
    getReadsFromSequence(x, chr, ...)
  })
  
  all.reads <- all.reads[sapply(all.reads, function(x) !is.null(x)&length(x))]
  
  ## make sure all DFs have the same column names
  keep.cols <- colnames(values(all.reads[[1]]))
  if (length(all.reads) > 1) {
    for (i in 2:length(all.reads)) {
      keep.cols <- intersect(keep.cols, colnames(values(all.reads[[i]])))
    }
  }
  
  all.reads <- lapply(all.reads, function(r) {
    meta <- values(r)
    if (ncol(meta) != length(keep.cols)) {
      meta <- meta[, keep.cols]
      values(r) <- meta
    }
    r
  })
  
  suppressWarnings(do.call(c, unname(all.reads)))
})




