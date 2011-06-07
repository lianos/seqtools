bsgSamHeader <- function(bsgenome) {
  lvls <- seqlengths(bsgenome)
  sprintf("@SQ\tSN:%s\tLN:%d", names(lvls), unname(lvls))
}
