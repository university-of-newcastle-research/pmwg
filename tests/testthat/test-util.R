test_that("Wind+Unwind", {
  test_mat <- matrix(rnorm(16), nrow = 4, ncol = 4)
  # Ensure matrix is positive semi-definite
  test_mat <- test_mat %*% t(test_mat)
  test_vec <- rnorm(10)
  # Now run tests
  # Test wind and unwind are complementary
  expect_equal(test_mat, wind(unwind(test_mat)))
  expect_equal(test_vec, unwind(wind(test_vec)))
  # Test unwind and wind end with correct shape
  expect_identical(length(unwind(test_mat)), 10L)
  expect_identical(dim(wind(test_vec)), c(4L, 4L))
})


test_that("numbers_from_proportion", {
  # Test error checking - now done in run_stage level
  # Test 0 for efficient sampling returning 0 particles to select
  expect_equal(numbers_from_proportion(c(0.5, 0.5, 0))[3], 0)
})


test_that("run_stage", {
  # Test error checking
  # expect_error(
  #   numbers_from_proportion(c(0.3, 0.3, 0.3)),
  #   "The elements of the mix_proportion vector must sum to 1"
  # )
  # expect_error(
  #   numbers_from_proportion(c(0.5, 0.5)),
  #   "mix_proportion vector must have three elements which sum to 1"
  # )
})
