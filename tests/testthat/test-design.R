test_that("design search returns a feasible deterministic optimum", {
  fit <- small_design()
  expect_s3_class(fit, "roed_design")
  expect_equal(fit$version, "1.0.0")
  expect_equal(sum(fit$design$n), 28)
  expect_lte(fit$operating_characteristics$fwer, 0.20 + 1e-12)
  expect_gte(fit$operating_characteristics$min_g1, 0.50 - 1e-12)
  expect_true(fit$search$single_threaded)

  again <- small_design()
  expect_identical(fit$design, again$design)
})

test_that("operating characteristics and protocol tables are coherent", {
  fit <- small_design()
  oc <- roed_oc(fit)
  protocol <- roed_protocol(fit)
  expect_s3_class(oc, "roed_oc")
  expect_s3_class(protocol, "roed_protocol")
  expect_equal(oc$fwer, fit$operating_characteristics$fwer)
  expect_equal(protocol$target_enrollment, fit$design$n)
  expect_equal(nrow(as.data.frame(oc)), 4)
})

test_that("final selection applies both integer thresholds", {
  fit <- small_design()
  selected <- roed_select(
    fit,
    toxicity = fit$design$m_t,
    efficacy = fit$design$m_e
  )
  expect_true(all(selected$selected))

  failed <- roed_select(
    fit,
    toxicity = fit$design$m_t + 1L,
    efficacy = fit$design$m_e
  )
  expect_false(any(failed$selected))
})

test_that("completed checkpoints return the saved design", {
  path <- tempfile(fileext = ".rds")
  fit <- small_design(checkpoint = path)
  resumed <- roed_resume(path)
  expect_identical(fit$design, resumed$design)
  unlink(path)
})

