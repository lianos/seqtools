################################################################################
## Write GRanges object into SAM records

## This doesn't try to be fancy. All of the info in the SAM file
## needs to be in the values() of the GRanges object.
##
## Counts aren't expanded.
##

##' Converts reads sotred in a \code{GRanges} into a SAM file.
toSAMTable <- function(x, is.paired=FALSE, tag.prefix="tag.") {
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
    seq <- ifelse(nchar(seq) == 0L, '*', as.character(seq))
  }

  qual <- meta$qual
  if (is.null(qual)) {
    qual <- rep('*', nrow(meta))
  } else {
    qual <- ifelse(nchar(qual) == 0L, '*', as.character(qual))
  }

  flag <- meta$flag
  if (is.null(flag)) {
    flag <- ifelse(strand(x) == '-', 16L, 0L)
  } else {
    flag <- as.integer(meta$flag)
  }

  ## If there is no CIGAR info, assume width of string as CIGAR
  cigar <- meta$cigar
  if (is.null(cigar)) {
    cigar <- paste(width(x), 'M', sep="")
  } else {
    cigar <- ifelse(is.na(cigar) | nchar(cigar) == 0L,
                    paste(width(x), 'M', sep=""),
                    as.character(cigar))
  }

  ## A value of 255 indicates that mapq score is not available
  mapq <- meta$mapq
  if (is.null(mapq)) {
    mapq <- rep(255L, nrow(meta))
  } else {
    mapq <- as.integer(ifelse(is.na(mapq), 255L, mapq))
  }

  if (!is.paired) {
    mrnm <- rep('*', nrow(meta)) # RNEXT = * if data is not paired end
    mpos <- rep(0L, nrow(meta))  # PNEXT = 0L if data is not paired
    isize <- rep(0L, nrow(meta)) #
  } else {
    mrnm <- meta$mrnm
    mrnm <- as.character(ifelse(is.na(mrnm), '*', mrnm))

    mpos <- meta$mpos
    mpos <- as.character(ifelse(is.na(mpos), 0L, mpos))

    isize <- meta$isize
    isize <- as.integer(ifelse(is.na(isize), 0L, mpos))
  }

  ## if (any(colnames(meta) == 'exon.anno')) {
  ##   rename <- which(colnames(meta) == 'exon.anno')
  ##   colnames(meta)[rename] <- 'tag.ZA'
  ## }

  sam <- data.frame(qname=as.character(meta$qname), flag=flag,
                    rname=as.character(seqnames(x)), pos=start(x),
                    mapq=mapq, cigar=cigar, mrnm=mrnm, mpos=mpos, isize=isize,
                    seq=seq, qual=qual)
  
  sam.tags <- parseTagsFromDataFrame(meta, tag.prefix)
  if (!is.null(tags)) {
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
##' @param df a \code{data.frame}-like object. Some of these columns should
##' have tag info to parse out.
##' @param tag.prefix the prefix used to find the elemnts (columns) of
##' \code{df} that are actually tags.
parseTagsFromDataFrame <- function(x, tag.prefix="tag.") {
  tag.cols <- grep(tag.prefix, names(x), fixed=TRUE)
  if (length(tag.cols) == 0L) {
    return(NULL)
  }
  tag.names <- gsub(tag.prefix, '', names(x)[tag.cols], fixed=TRUE)

  ## Tags that are absent from a given alignment will be either NA or NULL
  ## here. These values are converted to 0-length character vectors.
  clean.tags <- lapply(tag.names, function(name) {
    xtags <- x[[paste('tag', name, sep=".")]]
    xtype <- SAMTagType(xtags)
    xtags <- as.character(xtags)
    axe <- is.na(xtags) | is.null(xtags)
    xtags <- paste(name, xtype, xtags, sep=":")
    xtags[axe] <- ''
    xtags
  })
  names(clean.tags) <- tag.names
  collapseSamTagStrings(clean.tags)
}

## Returns a character vector as long as the elements included in the
## sub-lists that will be the tags for each alignment
collapseSamTagStrings <- function(..., .sep="\t", .use.c=TRUE) {
  clean.tags <- list(...)
  if (is.list(clean.tags)) {
    if (length(clean.tags) == 1L) {
      clean.tags <- clean.tags[[1]]
    } else {
      stop("What's going on here?")
    }
  }
  if (is.null(names(clean.tags))) {
    stop("Need names for arguments")
  }
  
  if (.use.c) {
    collapse.tags <- .Call("collapse_sam_tag_list", clean.tags, .sep,
                           PACKAGE="SeqTools")
  } else {
    ## TODO: Speed up tag:concatenation when creating SAM files, this is the
    ## slowest part of this function. Generating the `clean.tags` variable
    ## above isn't so bad.
    collapse.tags <- sapply(seq_along(clean.tags[[1]]), function(i) {
      xtags <- sapply(clean.tags, '[[', i)
      xtags <- xtags[nchar(xtags) > 0]
      do.call(paste, list(xtags,  collapse=.sep))
    })
  }

  collapse.tags
}

SAMTagType <- function(tag.values) {
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

