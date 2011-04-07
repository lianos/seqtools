##' BAM files generated from Bowtie use 0-based index
##' @nord
bam.parseTag.MD <- function(md.tag, from.zero.idx=TRUE) {
  ## 18G13 means :
  ##   * 18 basepairs of good alignment
  ##   * then the REFERENCE has a G instead of the character you see here
  ##   * then 13 bp of good alignment
  ##   * In R, substring(the.read, 18, 18) is G instead of whatever is shown
  ##     here in the BAM file
  ##
  ## 0T28G2
  ##   * 0 bp of gooad alignment: first bp is off and is T in reference
  ##   * then 28 bp of good alignment (28+1 position is hosed)
  ##   * 30th position is G in reference
  ##   * 2 more bp of good alignment
  ##
  ## 28A0G2
  ##   * 28 BP of good alignment
  ##   * 29th BP is A in reference
  ##   * 30th BP is G in reference
  ##   * 2 bp of good alignment
  eg <- c("32", "18G13", "32", "30G1", "21C10", "32", "32", "6C25")
  has.mm <- grep('[ACGT]', md.tag)
  if (length(has.mm) > 0) {

  }
}

##' Determine of the flag for a read from a BAM file indicates whether the
##' read is the first in a mate pair.
##'
##' @export
##' @author Steve Lianoglou \email{slianoglou@@gmail.com}
##'
##' @param flag A vector of flags (integers) from the BAM file
##' @param check.paired If this is \code{FALSE} no pair id is checked,
##' everything is considered to be the first pair. This for convenience
##' when a non-paired-run queries this function just as a matter of course.
##'
##' @return A logical vector as long as \code{flag}, where \code{TRUE}
##' indicates that the corresponding read/flag is the first in a mate pair
BAMisFirstPair <- function(flag, check.paired=TRUE) {
  ## Note this will fail if the reads aren't paired -- this will be 0
  ## if all(bitAnd(bam.reads$flag, 0x01) == 0), then these reads are unpaired
  if (check.paired) {
    if (all(bitAnd(flag, 0x01) == 0)) {
      return(rep(TRUE, length(flag)))
    }
  }
  bitAnd(flag, 0x0040) != 0
}

##' Returns the pair id for each read corresponding to \code{flag}.
##'
##' @export
##'
##' @param flag A vector of flags (integers) from the BAM file
##' @param check.paired Convenience function for non-paired-experiments.
##' If \code{FALSE}, every read will have a pair id of \code{1}.
##'
##' @return An integer vector as long as \code{flag} of \code{1} or
##' \code{2}'s.
BAMpairId <- function(flag, check.paired=TRUE) {
  is.first <- BAMisFirstPair(flag, check.paired=check.paired)
  pair.id <- as.integer(is.first)
  pair.id[!is.first] <- 2L
  pair.id
}
