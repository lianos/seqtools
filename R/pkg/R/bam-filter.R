##' Create a new BAM file which only includes defined regions.
##'
##' @param x path to a bam file
##' @param mask A GRanges object that indicates which regions to keep
##' @param destination The name of the bam file you want to output.
##' @param bam.strands \code{'keep', 'flip', 'ignore', 'reassign}.
##' 'keep' means that the strand info in the BAM file is "real".
##' 'flip' will flip the strands before filtering.
##' 'ignore' sets the strands returned from the bamfile to \code{*}.
##' 'reassign' is like ignore, but will reasign the strand of the reads
##' that match the mask to the strand of the mask.
##' @param ... These arguments are passed to the getReadsFromSequence
##' method to tweak which reads you want to retrieve from \code{x}.
##' For instance \code{unique.only=TRUE} would work.
filterBamThrough <- function(x, mask, destination, bam.strands='keep', ...) {
  stop("Export mask to a bed file and use bedtools!")
  
  stopifnot(is.character(x))
  if (!file.exists(x)) {
    stop("Can't read input bam file: ", x)
  }
  stopifnot(inherits(mask, 'GRanges'))
  if (file.exists(destination)) {
    stop(destination, " already exists.")
  }
  bam.strands <- match.arg(bam.strands, c('keep', 'flip', 'ignore', 'reassign'))

  out.path <- dirname(destination)
  out.name <- basename(destination)
  out.name <- gsub("\\.bam", "", out.name)
  out.name <- gsub("\\.sam", "", out.name)
  out.file <- file.path(out.path, paste(out.name, 'sam', sep="."))

  ## write.header
  header <- scanBamHeader(x)[[1]]
  df <- data.frame(names(header$text), sapply(header$text, '[[', 1),
                   sapply(header$text, '[[', 2))
  write.table(df, out.file, append=FALSE, quote=FALSE, sep="\t",
              row.names=FALSE, col.names=FALSE)

  ## bf <- BamFile(x)
  ## open(bf)
  chrs <- as.character(df[[2]])
  chrs <- gsub("SN:", "", chrs[grep("SN:", chrs)])

  out.sams <- mclapply(chrs, mc.preschedule=FALSE, function(chr) {
    cat("===", chr, "===\n")
    chr.mask <- mask[seqnames(mask) == chr]
    reads <- getReadsFromSequence(x, chr, smooth.by=NULL, meta.what=all.bam.what,
                                  with.sequence=TRUE, ...)
    if (length(reads) == 0) {
      return(NULL)
    }

    if (bam.strands %in% c('ignore', 'reasign')) {
      strand(reads) <- '*'
    } else if (bam.strands == 'flip') {
      strand(reads) <- swapStrand(strand(reads))
    }

    freads <- subsetByOverlaps(reads, chr.mask)
    if (length(freads) == 0) {
      return(NULL)
    }

    if (bam.strands == 'reasign') {
      xref <- assignUniqueOverlaps(freads, chr.mask)
      strand(freads) <- strand(chr.mask)[xref]
    }
    out.sam <- file.path(out.path, paste(out.name, chr, 'sam', sep="."))
    cat("Writing chromosome SAM to", out.sam, "\n")
    sam.table <- toSamTable(freads)
    write.table(sam.table, out.sam, sep="\t", col.names=FALSE, row.names=FALSE,
                quote=FALSE, append=FALSE)
    out.sam
  })
  out.sams <- out.sams[!sapply(out.sams, is.null)]

  cat("=== Combining files ===\n")
  for (file in out.sams) {
    cmd <- paste("cat", file, ">>", out.file)
    cat("   ", cmd, "\n")
    result <- system(cmd)
    if (result != 0) {
      stop("error in copying sam")
    }
    unlink(file)
  }

  cat("=== Creating bam file ===\n")
  bam.file <- asBam(out.file, gsub("\\.sam$", "", out.file))
  if (file.exists(paste(bam.file, 'bai', sep="."))) {
    unlink(out.file)
  } else {
    stop("Error creating the bam file")
  }
  invisible(NULL)
}
