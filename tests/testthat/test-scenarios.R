test_that("scenario inputs are normalized and labeled", {
  scenarios <- small_scenarios()
  expect_s3_class(scenarios, "roed_scenarios")
  expect_equal(dim(scenarios$toxicity), c(2, 2))
  expect_equal(scenarios$dose_names, c("D1", "D2"))
  expect_equal(nrow(as.data.frame(scenarios)), 4)
})

test_that("infeasible Bernoulli correlations are rejected", {
  expect_error(
    roed_scenarios(
      toxicity = matrix(0.90, 1, 1),
      efficacy = matrix(0.10, 1, 1),
      correlation = 1
    ),
    "infeasible"
  )
})

test_that("invalid probability matrices are rejected", {
  expect_error(
    roed_scenarios(c(0.1, 1.2), c(0.2, 0.3)),
    "lie in"
  )
  expect_error(
    roed_scenarios(matrix(0.1, 2, 2), matrix(0.2, 1, 2)),
    "identical dimensions"
  )
})

