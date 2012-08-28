
reduceSeqinfo <- function(x, y) {
  stopifnot(is(x, 'Seqinfo'))
  stopifnot(is(y, 'Seqinfo'))

  y[seqlevels(y)[seqlevels(y) %in% seqlevels(x)]]
}


extractSeqinfo <- function(x, y) {
  if (inherits(try(seqinfo(x)), "try-error")) {
    stop("An object with `seqinfo` is required for `x`")
  }
  if (inherits(try(seqinfo(y)), "try-error")) {
    stop("An object with `seqinfo` is required for `x`")
  }

  reduceSeqinfo(seqinfo(x), seqinfo(y))
}

rematchSeqinfo <- function(x, y) {
  if (inherits(try(seqinfo(x)), "try-error")) {
    stop("An object with `seqinfo` is required for `x`")
  }
  if (inherits(try(seqinfo(y)), "try-error")) {
    stop("An object with `seqinfo` is required for `x`")
  }

  si <- extractSeqinfo(x, y)

  si.x <- seqinfo(x)
  si.missed <- names(si.x)[!names(si.x) %in% names(si)]
  si.missed <- si.x[si.missed]

  suppressWarnings({
    si.new <- merge(si, si.missed)
  })

  new2old <- match(seqlevels(si.new), seqlevels(x))
  seqinfo(x, new2old) <- si.new
  x
}

trimSeqinfo <- function(x) {
  si <- tryCatch(seqinfo(x), error=function(e) NULL)
  if (is.null(si)) {
    stop("An object with `seqinfo` is required for `x`")
  }
  si <- si[seqnames(si)[seqnames(si) %in% as.character(seqnames(x))]]
  new2old <- match(seqlevels(si), seqlevels(seqinfo(x)))
  seqinfo(x, new2old) <- si
  x
}
