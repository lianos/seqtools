## Fetch reads straight from a BAM file -- no SeqStore or BSgenome.*
## necessary
setGeneric("getReadsFromSequence",
function(x, seqname, start=NULL, end=NULL, strand=NULL, unique.only=TRUE,
         smooth.by=NULL, with.sequence=FALSE, with.quality=FALSE,
         meta.what=c("flag"), ...) {
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
         with.quality, meta.what, ...) {
  args <- list(...)
  verbose <- checkVerbose(...)
  trace <- checkTrace(...)

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
    if (strand %in% c('-', -1)) {
      flag <- scanBamFlag(isMinusStrand=TRUE)
    } else {
      flag <- scanBamFlag(isMinusStrand=FALSE)
    }
  } else {
    flag <- scanBamFlag()
  }

  unique.bam.query <- unique.only && is.null(args$uniqueness.map)
  bam.tag <- c('NM', 'MD', args$tag)
  result <- query(x, which=which, flag=flag, max.mismatch=args$max.mismatch,
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

  reads <- GRanges(seqnames=seqname,
                   ranges=IRanges(start=result$pos, width=result$qwidth),
                   strand=result$strand,
                   pair.id=pair.id)

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

  for (tag.name in names(result$tag)) {
    values(reads)[[paste('tag', tag.name, sep=".")]] <- result$tag[[tag.name]]
  }

  ## Apply smoothers to data
  reads <- smoothReads(reads, smooth.by=smooth.by, ...)

  if (verbose) {
    cat("  SeqTools::getReadsFromSequence took",
        as.numeric(Sys.time()) - start.time, "seconds\n")
  }

  reads
})

setMethod("getReadsFromSequence", c(x="character"),
function(x, seqname, start=NULL, end=NULL, strand=NULL, unique.only=TRUE,
         smooth.by=NULL, with.sequence=FALSE, with.quality, meta.what=NULL,
         ...) {
  getReadsFromSequence(BamFile(x), seqname, start, end, strand, unique.only,
                       smooth.by, meta.what=c('flag'), ...)
})

