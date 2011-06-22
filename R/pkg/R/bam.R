setMethod("header", c(x="BamFile"),
function(x, as.df=FALSE, ...) {
  h <- scanBamHeader(path(x))
  ## if (as.df) {
  ##   dfs <- lapply(h)
  ## }
  h
})

setMethod("isPaired", c(x="BamFile"),
function(x, ...) {
  ## TODO: Implement a smart way to look at a bam file to detect if it's paired
  ## x@is.paired
  FALSE
})

##' @importFrom GenomicRanges seqnames
setMethod("seqnames", c(x="BamFile"),
function(x) {
  names(seqlengths(x))
})


##' @importFrom GenomicRanges seqlevels
setMethod("seqlevels", c(x="BamFile"),
function(x) {
  names(seqlengths(x))
})

setMethod("seqlengths", c(x="BamFile"),
function(x) {
  targets <- lapply(header(x), '[[', 'targets')
  ## TODO: Merge all values across the targets list and check that
  ##       duplicate values have the same specified length
  targets[[1]]
})

setMethod("aligner", c(x="BamFile"),
function(x, ...) {
  ## @PGID:bwaPN:bwaVN:0.5.9rc1-r11
  headr <- header(x)
  info <- lapply(headr, function(.x) .x$text[['@PG']])
  no.info <- sapply(info, is.null)
  if (any(no.info)) {
    warning("No aligner info for ",
            paste(basename(names(no.info)[no.info], collapse=",")))
    info <- info[!no.info]
  }
  if (length(info) == 0) {
    warning("No alignment info found", immediate.=TRUE)
    return(NULL)
  }

  getVal <- function(char.vector, prefix) {
    idx <- grep(tolower(prefix), tolower(char.vector), fixed=TRUE)
    if (length(idx) == 0) {
      ""
    } else {
      char.vector[idx[1]]
    }
  }

  programs <- sapply(info, getVal, 'ID')
  commands <- sapply(info, getVal, 'CL')
  versions <- sapply(info, getVal, 'VN')
  ## list(program=programs, version=versions, command=commands)
  sapply(programs, function(p) gsub("ID:", "", p, ignore.case=TRUE))
})
