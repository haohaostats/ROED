test_that("exact joint probability agrees with independence factorization", {
  scenarios <- roed_scenarios(0.10, 0.50, correlation = 0)
  cache <- ROED:::.build_probability_cache(scenarios, 8)
  observed <- ROED:::.selection_probability(
    cache, n = 8, m_t = 1, m_e = 3,
    p_t = 0.10, p_e = 0.50, q = 0.05
  )
  expected <- pbinom(1, 8, 0.10) *
    pbinom(2, 8, 0.50, lower.tail = FALSE)
  expect_equal(observed, expected, tolerance = 1e-12)
})

test_that("local error is the larger attainable marginal tail", {
  error <- ROED:::.local_error(10, 1, 5, 0.30, 0.25)
  expect_equal(
    unname(error[["a"]]),
    max(pbinom(1, 10, 0.30),
        pbinom(4, 10, 0.25, lower.tail = FALSE))
  )
})

