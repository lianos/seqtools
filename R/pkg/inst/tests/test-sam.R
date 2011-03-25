context("SAM")

test_that("C code to collapse tag produces same result as R", {
  tag.list <- list(a=letters[1:10], b=letters[11:20], c=letters[5:14])
  tag.list$a[c(1, 5, 10)] <- ""
  tag.list$b[c(4, 8, 10)] <- ""
  tag.list$c[c(1, 8, 10)] <- ""
  .r <- collapseSamTagStrings(tag.list, .use.c=FALSE)
  .c <- collapseSamTagStrings(tag.list, .use.c=TRUE)
  expect_equal(.c, .r)
})
