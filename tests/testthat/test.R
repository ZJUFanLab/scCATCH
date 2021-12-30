library(scCATCH)
context('scCATCH')

test_that('scCATCH', {
  expect_equal(colnames(demo_marker()), colnames(cellmatch))
})
