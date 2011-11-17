################################################################################
## Write GRanges object into SAM records

##' Converts reads sotred in a \code{GRanges} object into a SAM file.
##'
##' This method isn't very smart: all SAM columns need to be in the
##' \code{values()}'s \code{DataFrame} using the same column names that are
##' used in the values returned from \code{\link{scanBam}}. If the columns are
##' missing, then default values will be used.
##'
##' @param x A \code{GRanges} object
##' @param is.paired \code{logical} indicating if the reads are paired.
##' @param tag.prefix columns that who are destined to wind up in the tag fields
##' of the alignments must start with this prefix in order to be processed
##' correctly, eg: \code{tag.NM}, \code{tag.XA}, etc.
toSamTable <- function(x, is.paired=FALSE, tag.prefix="tag.") {
  if (!inherits(x, 'GRanges')) {
    stop("Only work on GRanges objects")
  }
  meta <- values(x)
  if (!'qname' %in% colnames(meta)) {
    stop("Query template name (qname) required")
  }

  seq <- meta$seq
  if (is.null(seq)) {
    seq <- rep('*', nrow(meta))
  } else {
    seq <- ifelse(as.logical(nchar(seq) == 0L), '*', as.character(seq))
  }

  qual <- meta$qual
  if (is.null(qual)) {
    qual <- rep('*', nrow(meta))
  } else {
    qual <- ifelse(as.logical(nchar(qual) == 0L), '*', as.character(qual))
  }

  flag <- meta$flag
  if (is.null(flag)) {
    flag <- ifelse(as.logical(strand(x) == '-'), 16L, 0L)
  } else {
    flag <- as.integer(meta$flag)
  }

  ## If there is no CIGAR info, assume width of string as CIGAR
  cigar <- meta$cigar
  if (is.null(cigar)) {
    cigar <- paste(width(x), 'M', sep="")
  } else {
    cigar <- ifelse(as.logical(is.na(cigar) | nchar(cigar) == 0L),
                    paste(width(x), 'M', sep=""),
                    as.character(cigar))
  }

  ## A value of 255 indicates that mapq score is not available
  mapq <- meta$mapq
  if (is.null(mapq)) {
    mapq <- rep(255L, nrow(meta))
  } else {
    mapq <- as.integer(ifelse(as.logical(is.na(mapq)), 255L, mapq))
  }

  if (!is.paired) {
    mrnm <- rep('*', nrow(meta)) # RNEXT = * if data is not paired end
    mpos <- rep(0L, nrow(meta))  # PNEXT = 0L if data is not paired
    isize <- rep(0L, nrow(meta)) #
  } else {
    mrnm <- meta$mrnm
    mrnm <- as.character(ifelse(as.logical(is.na(mrnm)), '*', mrnm))

    mpos <- meta$mpos
    mpos <- as.character(ifelse(as.logical(is.na(mpos)), 0L, mpos))

    isize <- meta$isize
    isize <- as.integer(ifelse(as.logical(is.na(isize)), 0L, mpos))
  }

  sam <- data.frame(qname=as.character(meta$qname), flag=flag,
                    rname=as.character(seqnames(x)), pos=start(x),
                    mapq=mapq, cigar=cigar, mrnm=mrnm, mpos=mpos, isize=isize,
                    seq=seq, qual=qual)

  sam.tags <- combineIntoSamTagsVector(meta, tag.prefix)
  if (!is.null(sam.tags)) {
    sam$tags <- sam.tags
  }

  sam
}

###############################################################################
## Deal with sam meta tags

##' Parses the type of each column in df to add the appropriate SAM type
##' (eg i,Z, etc.) and returns a character vector of tags that collapse
##' the tags across rows.
##'
##' @param x a \code{list}-like object (a \code{data.frame} will do) with
##' each element representing a tag "column." Each elements of \code{x} must
##' be the same length.
##' have tag info to parse out.
##' @param tag.prefix the prefix used to find the elemnts (columns) of
##' \code{df} that are actually tags.
##' @param sep The separator to collapse multiple tags together per alignment
##' @param .use.c Use C code for combining tags for uber quickness.
combineIntoSamTagsVector <- function(x, tag.prefix="tag.", sep="\t", .use.c=TRUE) {
  ## Ensure that all items in data.frame/list container are the same length
  lengths <- unique(sapply(x, length))
  if (length(lengths) != 1L) {
    stop("The elements in `x` aren't uniform length")
  }
  tag.cols <- grep(tag.prefix, names(x), fixed=TRUE)
  if (length(tag.cols) == 0L) {
    return(NULL)
  }

  tag.names <- gsub(tag.prefix, '', names(x)[tag.cols], fixed=TRUE)
  if (!all(nchar(tag.names) == 2)) {
    stop("tag names can only be two letters long (XA, NM, etc.)")
  }

  ## Tags that are absent from a given alignment will be either NA or NULL
  ## here. These values are converted to 0-length character vectors.
  clean.tags <- lapply(tag.names, function(name) {
    xtags <- x[[paste('tag', name, sep=".")]]
    xtype <- SamTagType(xtags)
    xtags <- as.character(xtags)
    axe <- is.na(xtags) | is.null(xtags)
    xtags <- paste(name, xtype, xtags, sep=":")
    xtags[axe] <- ''
    xtags
  })
  ##names(clean.tags) <- tag.names
  .collapseSamTagStrings(clean.tags, .sep=sep, .use.c=.use.c)
}

##' Returns a character vector as long as the elements included in the
##' sub-lists that will be the tags for each alignment.
##'
##' Tag vectors are passed into \code{...}, which must be the same length.
##' These are merged into one "tag string" per element.
##'
##' Elements in the input vectors that are NULL or NA are not included.
##'
##' @nord
.collapseSamTagStrings <- function(..., .sep="\t", .use.c=TRUE) {
  clean.tags <- list(...)
  if (is.list(clean.tags)) {
    if (length(clean.tags) == 1L) {
      clean.tags <- clean.tags[[1]]
    } else {
      stop("What's going on here?")
    }
  }

  if (.use.c) {
    collapse.tags <- .Call("collapse_sam_tag_list", clean.tags, .sep,
                           PACKAGE="SeqTools")
  } else {
    warning("You are using the R implementation of collapseSamTagStrings",
            immediate.=TRUE)
    collapse.tags <- sapply(seq_along(clean.tags[[1]]), function(i) {
      xtags <- sapply(clean.tags, '[[', i)
      xtags <- xtags[nchar(xtags) > 0]
      do.call(paste, list(xtags,  collapse=.sep))
    })
  }

  collapse.tags
}

##' Determine the SAM "type-label" for this tag.
##'
##' SAM types are:
##'
##' \begin{enumerate}
##'   \item \code{A} printable character
##'   \item \code{Z} printable string (spaces allowed)
##'   \item \code{i} signed 32-bit integer
##'   \item \code{f} single-precision floating number
##'   \item \code{H} hex string (high nybble first)
##' \end
##'
##' If none of the tests for A, Z, i, or f pass, then H is assumed.
##'
##' @param tag.values A character vector of the same type of tags to use
##' inorder to guess the appropriate type.
##' @return \code{character(1)} indicating the inferred type.
SamTagType <- function(tag.values) {
  ## This doesn't correctly detect numeric/integers
  if (is.integer(tag.values)) {
    return('i')
  } else if (is.character(tag.values)) {
    if (all(nchar(tag.values) == 1)) {
      return('A')
    } else {
      return("Z")
    }
  } else if (is.numeric(tag.values)) {
    return('f')
  } else {
    ## HEX? doubtful, but ...
    return("H")
  }
}

################################################################################
## Query SAM files directly

##' Returns the command used for the aligner that created this SAM file
##'
##' @export
##' @seealso \code{link{parseAlignCommandFromBAM}}
##'
##' @param sam.file The path to the SAM file
##' @param n The number of lines to read at the top of the SAM file to
##' find the result
##' @return The command, program, and version used for alignment
setMethod("aligner", c(x="character"),
function(x, ...) {
  lines <- readLines(sam.file, n=1000)
  take <- grep("@PG", lines)
  if (length(take) == 0) {
    cat("No alignment-info in SAM header ")
    return(NULL)
  }
  info <- strsplit(lines[take], "\t")
  if (length(info) > 1) {
    cat("There should only be one line for alignment into in header")
    return(NULL)
  }
  info <- sapply(info[[1]], function(x) gsub("\"", "", x))
  names(info) <- c("PG", "ID", "VN", "CL")
  list(program=gsub("ID:", "", info[['ID']]),
       version=gsub("VN:", "", info[['VN']]),
       command=gsub("CL:", "", info[['CL']]))
})

##' Reads the top \code{n} lines from a SAM file to extract chromosome names.
##'
##' The names of the chromosomes are in the header of the file. These lines
##' start with \code{@SQ}
##'
##' @export
##' @seealso \code{\link{parseChromosomeInfoFromBAM}}
##'
##' @param sam.file The filepath of the SAM file
##' @param n The number of lines to read in from the SAM file. 1000 should be
##' \emph{more} than enough, since these @SQ lines should appear in the top of
##' the file.
##' @return A data.frame, with \code{$name} and \code{$length} columns
setMethod("seqinfo", c(x="character"),
function(x) {
  ext <- substring(x, nchar(x) - 2)
  if ( ext == 'bam') {
    return(seqinfo(BamFile(x)))
  }
  if (ext != 'sam') {
    stop("seqinfo,character is only defined for sam or bam files")
  }
  
  lines <- readLines(x, n=1000)
  take <- grep("@SQ", lines)
  ## @SQ	SN:chr1	LN:247249719
  ## @SQ	SN:chr2	LN:242951149
  ## @SQ	SN:chr3	LN:199501827
  ## @SQ	SN:chr4	LN:191273063
  info <- lapply(strsplit(lines[take], "\t"), function(x) {
    sapply(strsplit(x[-1], ':'), '[', 2)
  })
  df <- as.data.frame(do.call(rbind, info), stringsAsFactors=FALSE)
  colnames(df) <- c('seqnames', 'length')
  df$length <- as.integer(df$length)
  Seqinfo(df$seqnames, df$length)
})

setMethod("seqnames", c(x="character"),
function(x) {
  seqlevels(seqinfo(x))
})

setMethod("seqlevels", c(x="character"),
function(x) {
  seqlevels(seqinfo(x))
})

setMethod("seqlengths", c(x="character"),
function(x) {
  seqlengths(seqinfo(x))
})

################################################################################
## Deprecated (?)


## scanBam2GRanges <- function(result) {
##   defaults <- list(qname='BOGUS', flag=4L, mapq=255L, cigar='1M', mrnm="*",
##                    mpos=0L, isize=0L, qual="I")

##   defaultFill <- function(from, wut, default, n=length(from[[1]])) {
##     x <- from[[wut]]
##     if (is.null(x)[[1]]) {
##       x <- rep(default, n)
##     }
##     x <- ifelse(is.na(x) | is.null(x), default, x)
##     x
##   }

##   x <- GRangesList(lapply(result, function(r) {
##     gr <- GRanges(r$rname, IRanges(r$pos, width=r$qwidth), b$strand)
##     ignore <- c('rname', 'rname', 'pos', 'qwidth', 'strand', 'tag')
##     build <- names(r)[!names(r) %in% ignore]
##     meta <- lapply(build, function(name) {
##       defaultFill(r, name, default[[name]])
##     })
##     names(meta) <- build
##     meta <- DataFrame(meta)
##     flag <- defaultFill
##     meta <- DataFrame(qname=r$qname, flag=r$flag, mapq=r$mapq, cigar=r$cigar,
##                       mrnm=ifelse(is.na(r$mrnm), '*', r$mrnm), mpos=r$mpos,
##                       isize=r$isize,

##   }))
##   names(x) <- names(result)
##   x
## }

