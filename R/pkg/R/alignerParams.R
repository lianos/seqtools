## Construct aligner-specific ScanBamParams for some common alignment options

bwaSamseParams <- function(max.mismatch=NULL, with.alt.alignments=FALSE, ...) {
  p <- ScanBamParam()
  tags <- c("X0", "X1", "XM")
  if (with.alt.alignments) {
    tags <- c(tags, "XA")
  }
  bamTag(p) <- tags
  p
}
