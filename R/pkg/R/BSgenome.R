setMethod("getBsGenome", c(x="character"),
function(x, organism=NULL, anno.source='UCSC', ...) {
  lib.name <- 'BSgenome.:organism:.:anno.source:.:genome:'
  if (is.null(organism)) {
    organism <- switch(substring(x, 1, 2),
                       hg='Hsapiens',
                       mm='Mmusculus',
                       sa='Scerevisiae',
                       dm='Dmelanogaster',
                       rn='Rnorvegicus',
                       ce='Celegans',
                       stop("Unknown genome", x, sep=" "))
  }

  lib.name <- gsub(':organism:', organism, lib.name)
  lib.name <- gsub(':anno.source:', anno.source, lib.name)
  lib.name <- gsub(':genome:', x, lib.name)

  suppressPackageStartupMessages({
    found <- require(lib.name, character.only=TRUE)
  })

  if (!found) {
    stop(lib.name, " package required.")
  }

  get(organism, pos=paste('package', lib.name, sep=":"))
})

setMethod("getBsGenome", c(x="GenomicRanges"), function(x, ...) {
  g <- unique(genome(x))
  if (!is.character(g) || length(g) != 1L || nchar(g) < 1) {
    stop("Expected a single, valid genome identifier in seqinfo slot")
  }
  getBsGenome(g, ...)
})

bsgSamHeader <- function(bsgenome) {
  lvls <- seqlengths(bsgenome)
  sprintf("@SQ\tSN:%s\tLN:%d", names(lvls), unname(lvls))
}
