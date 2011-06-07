filterBamThrough <- function(x, mask, destination, unique.only=FALSE,
                             max.mismatch=NULL, ...) {
  stopifnot(inherits(mask, 'GRanges'))
  if (missing(destination)) {
    if (inherits(x, 'SeqStore')) {
      out.path <- directory(x)
      out.name <- 'filtered-reads'
    } else if (is.character(x)){
      out.path <- dirname(x)
      out.name <- basename(x)
    }
  } else {
    out.path <- dirname(destination)
    out.name <- basename(destination)
  }
  stopifnot(dir.exists(out.path))

  if (inherits(x, 'SeqStore')) {
    x <- bamFile(x, ...)
  }
  stopifnot(is.character(x))
  if (!file.exists(x)) {
    stop("Can't read input bam file: ", x)
  }

  out.name <- gsub("\\.bam", "", out.name)
  out.name <- gsub("\\.sam", "", out.name)
  out.file <- file.path(out.path, paste(out.name, 'sam', sep="."))

  ## write.header
  header <- scanBamHeader(x)[[1]]
  df <- data.frame(names(header$text), sapply(header$text, '[[', 1),
                   sapply(header$text, '[[', 2))
  write.table(df, out.file, append=FALSE, quote=FALSE, sep="\t",
              row.names=FALSE, col.names=FALSE)

  bf <- BamFile(bam.file)
  chrs <- as.character(unique(seqnames(bf)))
  out.sams <- mclapply(chrs, mc.preschedule=FALSE, function(chr) {
    cat("===", chr, "===\n")
    chr.mask <- mask[seqnames(mask) == chr]
    reads <- getReadsFromSequence(BamFile(bam.file), chr,
                                  unique.only=unique.only,
                                  max.mismatch=max.mismatch,
                                  smooth.by=NULL, meta.what=SeqTools:::all.bam.what,
                                  with.sequence=TRUE, ...)
    freads <- subsetByOverlaps(reads, chr.mask)
    out.sam <- file.path(out.path, paste(out.name, chr, 'sam', sep="."))
    cat("Writing chromosome SAM to", out.sam, "\n")
    sam.table <- toSamTable(freads)
    write.table(sam.table, out.sam, sep="\t", col.names=FALSE, row.names=FALSE,
                quote=FALSE, append=FALSE)
    out.sam
  })

  cat("=== Combining files ===\n")
  for (file in out.sams) {
    cmd <- paste("cat", file, ">>", out.file)
    cat("   ", cmd, "\n")
    system(cmd)

    system(paste("rm", file))
  }

  cat("=== Creating bam file ===\n")
  asBam(out.file, gsub("\\.sam$", "", out.file))

}
