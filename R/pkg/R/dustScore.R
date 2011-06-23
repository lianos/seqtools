
## Core functionality taken from a post by Herve Pages on bioc mailing list
## (TODO: Find URL for post)
setGeneric("dustScore", signature=c('x', 'over'),
function(x, over, k=3, with.sequences=FALSE, ...) {
  standardGeneric("dustScore")
})

setMethod("dustScore", c(x="GenomicRanges", over="BSgenome"),
function(x, over, k, with.sequences, .parallel=TRUE, ...) {
  chrs <- seqlevels(x)

  if (.parallel && length(chrs) > 1L) {
    '%loop%' <- getFunction('%dopar%')
  } else {
    '%loop%' <- getFunction('%do%')
  }

  ds <- foreach(chr=chrs, .packages=c('SeqTools', 'Biostrings'),
                .options.multicore=list(preschedule=FALSE)) %loop% {
    cat("===", chr, "===\n")
    idxs <- which(seqnames(x) == chr)
    if (length(idxs) == 0) {
      return(NULL)
    }
    xx <- x[idxs]
    bsg.chr <- unmasked(over[[chr]])
    df <- dustScore(ranges(xx), bsg.chr, k, with.sequences, ...)
    values(xx) <- as(transform(df, idx=idxs), 'DataFrame')
    xx
  }

  ds <- ds[!sapply(ds, is.null)]
  ds <- do.call(c, unname(ds))
  ds <- ds[values(ds)$idx]
  values(ds) <- values(ds)[, colnames(values(ds)) != 'idx']
  ds
})


setMethod("dustScore", c(x="Ranges", over="XString"),
function(x, over, k, with.sequences, ...) {
  seqs <- Views(over, x)

  s <- .dust.score.core(trinucleotideFrequency(seqs), width(x))

  df <- data.frame(dust.score=s)
  if (with.sequences) {
    df$sequence <- as.character(seqs)
  }

  df
})

setMethod("dustScore", c(x="XStringSet", over="ANY"),
function(x, over, k, with.sequences, ...) {
  s <- .dust.score.core(trinucleotideFrequency(x), width(x))

  df <- data.frame(dust.score=s)
  if (with.sequences) {
    df$sequence <- as.character(seqs)
  }

  df
})

.dust.score.core <- function(tri, widths) {
  stopifnot(is(tri, 'matrix'))
  tri2 <- tri - 1L
  tri2[tri2 < 0] <- 0L

  s <- (tri * tri2 / 2) / (widths - 1L)
  rowSums(s)
}

plotDustScoreDistro <- function(dust.scores, expand=TRUE) {
  ## exon.annos <- unique(dust.scores$exon.anno)
  ## scores <- lapply(exon.annos, function(ea) {
  ##   x <- subset(dust.scores, exon.anno == ea)
  ##   .scores <- x$score
  ##   if (expand) {
  ##     .scores <- rep(.scores, x$count)
  ##   }
  ##   .scores
  ## })
  ## boxplot(scores, names=exon.annos)

  scores <- ddply(dust.scores, .(exon.anno), function(x) {
    .scores <- x$score
    if (expand) {
      .scores <- rep(.scores, x$count)
    }
    data.frame(exon.anno=x$exon.anno[1], score=log2(.scores))
  })

  ggplot(scores, aes(exon.anno, score)) + geom_boxplot()
}
