context("SAM")

test_that("C code to collapse tag produces same result as R", {
  N <- 20
  N.axe <- 5
  tag.list <- lapply(c('a', 'b', 'c'), function(x) {
    x <- sample(letters, N, replace=TRUE)
    x[sample(N, N.axe)] <- ''
    x
  })
  names(tag.list) <- paste('tag', c('XA', 'XB', 'XC'), sep='.')
  
  # .r <- SeqTools:::.collapseSamTagStrings(tag.list, .use.c=FALSE)
  # .c <- SeqTools:::.collapseSamTagStrings(tag.list, .use.c=TRUE)
  suppressWarnings(.r <- combineIntoSamTagsVector(tag.list, .use.c=FALSE))
  .c <- combineIntoSamTagsVector(tag.list, .use.c=TRUE)
  expect_equal(.c, .r)
})

test_that("C code to collapse tags is (much) faster", {
  N <- 3000
  N.axe <- 100
  tag.list <- lapply(c('a', 'b', 'c'), function(x) {
    x <- sample(letters, N, replace=TRUE)
    x[sample(N, N.axe)] <- ''
    x
  })
  names(tag.list) <- paste('tag', c('XA', 'XB', 'XC'), sep=".")
  
  suppressWarnings({
    r.time <- system.time(combineIntoSamTagsVector(tag.list, .use.c=FALSE))
  })
  expect_that(combineIntoSamTagsVector(tag.list, .use.c=TRUE),
              takes_less_than(r.time['elapsed'] / 10))
})

