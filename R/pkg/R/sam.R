################################################################################
## Write GRanges object into SAM records

## This doesn't try to be fancy. All of the info in the SAM file
## needs to be in the values() of the GRanges object.
##
## Counts aren't expanded.
##
toSAMTable <- function(x, is.paired=FALSE) {
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
  tag.cols <- grep('tag.', colnames(meta), fixed=TRUE)
  tag.names <- gsub('tag.', '', colnames(meta)[tag.cols], fixed=TRUE)

  ## Tags that are absent from a given alignment will be either NA or NULL
  ## here. These values are converted to 0-length character vectors.
  clean.tags <- lapply(tag.names, function(name) {
    xtags <- meta[[paste('tag', name, sep=".")]]
    xtype <- SAMTagType(xtags)
    xtags <- as.character(xtags)
    axe <- is.na(xtags) | is.null(xtags)
    xtags <- paste(name, xtype, xtags, sep=":")
    xtags[axe] <- ''
    xtags
  })
  names(clean.tags) <- tag.names

  sam$tags <- collapseSamTagStrings(clean.tags)
  sam
}


## Returns a character vector as long as the elements included in the sub-lists
## that will be the tags for each alignment
collapseSamTagStrings <- function(...) {
  args <- list(...)
  if (is.list(args)) {
    if (length(args) == 1L) {
      args <- args[[1]]
    } else {
      stop("What's going on here?")
    }
  }
  if (is.null(names(args))) {
    stop("Need names for arguments")
  }

  ## TODO: Speed up tag:concatenation when creating SAM files, this is the
  ## slowest part of this function. Generating the `clean.tags` variable
  ## above isn't so bad.
  collapse.tags <- sapply(seq_along(clean.tags[[1]]), function(i) {
    xtags <- sapply(clean.tags, '[[', i)
    xtags <- xtags[nchar(xtags) > 0]
    do.call(paste, list(xtags,  collapse="\t"))
  })

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
  } else if (is.numeric(tag.vdalues)) {
    return('f')
  } else {
    ## HEX? doubtful, but ...
    return("H")
  }
}

## expand.multimap requires the bwa aligner -- it will expand the
## XA flag into SAM records
generateCleanSAM <- function(seqstore, gcache, file=NULL, unique.only=TRUE,
                             max.mismatch=0L, remove=c('antisense'),
                             expand.multimap=TRUE,
                             bam.file=NULL, gene.collapse='longest', flank.up=1000L,
                             flank.down=1000L, stranded=TRUE,
                             path=directory(seqstore),
                             append=FALSE, with.header=!append, chrs=NULL,
                             every=1L) {
  if (is.null(bam.file)) {
    bam.file <- bamFile(seqstore)
  }

  if (is.null(file)) {
    file <- 'clean-alignments.sam'
  }
  file <- file.path(path, file)

  header <- scanBamHeader(bam.file)[[1]]

  if (is.null(chrs)) {
    chrs <- header$targets
  } else {
    bsg <- getBsGenome(gcache)
    chrs <- seqlengths(bsg)[chrs]
  }

  foreach(chr=names(chrs), .packages=c("GenomicFeaturesX", "Rsamtools")) %dopar% {
    cat("===", chr, "===\n")
    ## do.append <- append || with.header || chr != names(chrs)[1]
    do.append <- FALSE

    ac <- tryCatch({
      getAnnotatedChromosome(gcache, chr, gene.collapse, flank.up,
                             flank.down, stranded)
    }, error=function(e) NULL)
    if (is.null(ac)) {
      cat("  ERROR: Getting annotations for", chr, "\n")
      return(NULL)
    }

    chr.reads <- tryCatch({
      getReadsOnChromosome(seqstore, chr, unique.only=unique.only,
                           smooth.by=NULL, max.mismatch=max.mismatch,
                           meta.what=NULL)
    }, error=function(e) NULL)
    if (is.null(chr.reads)) {
      cat("  ERROR: Getting reads for", chr, "\n")
      return(NULL)
    }

    reads <- annotateReads(chr.reads, ac)
    idx <- which(colnames(values(reads)) == 'exon.anno')
    colnames(values(reads))[idx] <- 'tag.ZA'

    if (!is.null(remove)) {
      keep <- !values(reads)$tag.ZA %in% remove
      reads <- reads[keep]
    }

    if (length(reads) == 0) {
      cat("  NOTE: No reads left after removal\n")
      return(NULL)
    }

    sam <- toSAMTable(reads)
    fn <- paste(file, chr, sep=".")
    cat(fn, "\n")
    write.table(sam, file=fn, append=do.append, sep="\t",
                col.names=FALSE, row.names=FALSE, quote=FALSE)
    chr
  }

  cat("Combining files\n")
  if (with.header) {
    df <- data.frame(names(header$text), sapply(header$text, '[[', 1),
                     sapply(header$text, '[[', 2))
    write.table(df, file, append=FALSE, quote=FALSE, sep="\t",
                row.names=FALSE, col.names=FALSE)
  }

  for (chr in names(chrs)) {
    cat("...", chr, "...\n")
    this.file <- paste(file, chr, sep=".")
    if (file.exists(this.file)) {
      cmd <- paste("cat", this.file, ">>", file)
      cat("   ", cmd, "\n")
      code <- system(cmd)
      if (code == 0) {
        unlink(this.file)
      }
    }
  }
}

################################################################################
## Deprecated (?)

##' Converts a GRanges object into a SAM record.
##'
##' This is really meant to turn a CompressedRanges object into a SAM/BAM file
##'
##' The SAM flag only tells us if the read maps to the +/- strand.
##' The quality/other meta info is essentially bogus
asSAM <- function(granges, file, genome='hg19', append=FALSE,
                  with.header=!append, path=".", save=TRUE) {
  file <- file.path(path, file)
  bsg <- getBsGenome(genome)
  if (with.header) {
    header <- data.frame(tag="@SQ",
                         chr=paste("SN", names(seqlengths(bsg)), sep=":"),
                         length=paste("LN", seqlengths(bsg), sep=":"))
    write.table(header, file, append=append, quote=FALSE, sep="\t",
                row.names=FALSE, col.names=FALSE)
  }

  chrs <- unique(as.character(seqnames(granges)))
  for (chr in chrs) {
    cat("===", chr, "===\n")
    do.append <- append || with.header || chr != chrs[1]
    bsg.chr <- unmasked(bsg[[chr]])

    gr <- granges[seqnames(granges) == chr]
    meta <- values(gr)

    n <- if ("count" %in% colnames(meta)) meta$count else rep(1L, nrow(meta))

    flag <- rep(as.integer(ifelse(strand(gr) == '-', 16L, 0L)), n)
    pos <- rep(start(gr), n)
    mapq <- rep(255L, length(pos))
    cigar <- rep(paste(width(gr), "M", sep=""), n)
    mrnm <- rep("*", length(pos))
    mpos <- rep(0L, length(pos))
    isize <- rep(0L, length(pos))
    seq <- rep(as.character(Views(bsg.chr, ranges(gr))), n)
    qual <- sapply(width(gr), function(i) paste(rep("I", i), collapse=""))
    qual <- rep(qual, n)
    sam <- data.frame(qname='BOGUS', flag=flag, rname=chr, pos=pos,
                      mapq=mapq, cigar=cigar, mrnm=mrnm, mpos=mpos, isize=isize,
                      seq=seq, qual=qual)

    write.table(sam, file, append=do.append, sep="\t", col.names=FALSE,
                row.names=FALSE, quote=FALSE)
  }

}


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

