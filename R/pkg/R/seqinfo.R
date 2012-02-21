extractSeqinfo <- function(x, y) {
  if (inherits(try(seqinfo(x)), "try-error")) {
    stop("An object with `seqinfo` is required for `x`")
  }
  if (inherits(try(seqinfo(y)), "try-error")) {
    stop("An object with `seqinfo` is required for `x`")
  }
  
  si.y <- seqinfo(y)
  si.y[seqlevels(si.y)[seqlevels(si.y) %in% seqlevels(x)]]
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
