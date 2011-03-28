context("converting NAs in objects")

test_that("NA's in logical vectors are turned to FALSE", {
  expect_that(SeqTools:::convert.na(c(TRUE, TRUE, NA, TRUE, FALSE)),
              is_identical_to(c(TRUE, TRUE, FALSE, TRUE, FALSE)))
})

