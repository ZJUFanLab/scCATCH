library(scCATCH)
context('scCATCH')

test_that('scCATCH', {
  expect_equal(unique(findmarkergenes(mouse_kidney_203_Seurat,'Mouse')[[2]][,1]), c('1','2','3'))
})
