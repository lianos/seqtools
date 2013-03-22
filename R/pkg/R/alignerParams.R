## Construct aligner-specific ScanBamParams for some common alignment options

bwaSamseParams <- function(max.mismatch=NULL, with.alt.alignments=FALSE,
                           flag=scanBamFlag(isUnmappedQuery=FALSE), ...) {
  p <- ScanBamParam(flag=flag)
  tags <- c("X0", "X1", "XM")
  if (with.alt.alignments) {
    tags <- c(tags, "XA")
  }
  bamTag(p) <- tags
  p
}
