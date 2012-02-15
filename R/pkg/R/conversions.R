setAs("GRanges", "data.table", function(from) {
  if (length(from) == 0L) {
    return(data.table())
  }
  as.data.table(as(from, 'data.frame'))
})

setAs("GRanges", "data.frame", function(from) {
  x <- lapply(as.data.frame(from), function(xx) {
    if (is.factor(xx)) {
      xx <- as.character(xx)
    }
    xx
  })
  as.data.frame(x, stringsAsFactors=FALSE)
})

setAs("data.table", "GRanges", function(from) {
  as(as.data.frame(from), "GRanges")
})

setAs("data.frame", "GRanges", function(from) {
  if (nrow(from) == 0L || all(is.na(from))) {
    return(GRanges())
  }
  if (!'seqnames' %in% colnames(from)) {
    stop("seqnames required")
  }
  gr.meta.take <- colnames(from) %in% c('seqnames', 'strand')
  gr.meta <- from[, gr.meta.take, drop=FALSE]
  from <- from[, !gr.meta.take, drop=FALSE]

  .ranges <- as(from, 'IRanges')
  DF <- elementMetadata(.ranges)
  elementMetadata(.ranges) <- NULL

  .strand <- if ('strand' %in% colnames(gr.meta)) gr.meta$strand else '*'
  gr <- GRanges(seqnames=gr.meta$seqnames, ranges=.ranges, strand=.strand)
  values(gr) <- DF
  gr
})

setAs("data.table", "IRanges", function(from) {
  as(as.data.frame(from), "IRanges")
})

setAs("data.frame", "IRanges", function(from) {
  if (nrow(from) == 0L || all(is.na(from))) {
    return(IRanges())
  }
  need <- c('start', 'end', 'width')
  have <- colnames(from)[colnames(from) %in% need]
  if (length(have) < 2) {
    stop("Need any two of 'start', 'end', or 'width'")
  }
  ## prefer start/end
  if (all(c('start', 'end') %in% have)) {
    iranges <- IRanges(from$start, from$end)
  } else {
    iranges <- do.call(IRanges, as.list(from[, have]))
  }

  meta.cols <- setdiff(colnames(from), have)

  if (length(meta.cols) > 0) {
    DF <- DataFrame(from[, meta.cols, drop=FALSE])
    colnames(DF) <- meta.cols
    values(iranges) <- DF
  }

  iranges
})
