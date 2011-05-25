setMethod("swapStrand", c(x="GRanges"),
function(x, ...) {
  strand(x) <- swapStrand(strand(x))
  x
})

## for Rle or vector
setMethod("swapStrand", c(x="ANY"),
function(x, ...) {
  strands <- strand(as.character(x))
  consider <- strands != '*'
  to.neg <- strands == '+' & consider
  to.pos <- strands == '-' & consider
  if (sum(to.neg) > 0) {
    strands[to.neg] <- '-'
  }
  if (sum(to.pos) > 0) {
    strands[to.pos] <- '+'
  }
  strands
})
