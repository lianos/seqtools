context("AREScore")

test_that("AREscores are close", {
  load(file.path(system.file("extdata", "AREscore.verify.rda", package="SeqTools")))
  ss <- unname(AREscore(are.seqs))
  web <- values(are.seqs)$Score
  expect_equal(which(ss == 0), which(web == 0))
  ss <- ss[ss > 0]
  web <- web[web > 0]
  expect_true(all(abs(ss - web) / web < .1))
})

################################################################################
## Utility code to get the testing data
parseResults <- function(fpath) {
  result <- read.table(fpath, header=TRUE, sep="\t")
  pent.pos <- strsplit(result$Pentamer.position, " ")
  pent.pos <- lapply(pent.pos, function(x) {
    if (length(x) == 0) {
      return(IRanges())
    }
    IRanges(as.integer(x), width=5)
  })
  pent.pos <- IRangesList(pent.pos)
  au.blocks <- strsplit(result$AU.block, " ")
  aub <- lapply(au.blocks, function(x) {
    if (length(x) == 0) {
      return(IRanges())
    }
    starts <- as.integer(gsub("\\.+\\d+$", "", x))
    ends <- as.integer(gsub("^\\d+\\.+", "", x))
    IRanges(starts, ends)
  })
  result <- as(result, "DataFrame")
  result$Pentamer.position <- pent.pos
  result$AU.block.position <- IRangesList(aub)
  result
}

if (FALSE) {
  results <- parseResults('/Users/stavros/cBio/projects/seqtools/R/pkg/inst/extdata/AREScore-output.txt')
  are.seqs <- readDNAStringSet('/Users/stavros/cBio/projects/seqtools/R/pkg/inst/extdata/AREScore-input.fa')
  values(are.seqs) <- results[match(names(are.seqs), results$Name),]
  stopifnot(all(rownames(values(are.seqs)) == names(are.seqs)))
  save(are.seqs, file='/Users/stavros/cBio/projects/seqtools/R/pkg/inst/extdata/AREscore.verify.rda')
}
