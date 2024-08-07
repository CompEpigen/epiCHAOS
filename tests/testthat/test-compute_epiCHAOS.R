
# some sample data
counts <- rbinom(n=2500, size=1, prob=0.2) %>% matrix(50,500)
rownames(counts) <- 1:nrow(counts)
colnames(counts) <- paste0("cell", 1:ncol(counts))
meta <- data.frame(row.names = colnames(counts), cluster=sample(c("group1", "group2"), 500, replace = T))

matrices <- create.group.matrices(counts=counts, meta = meta)

test_that("creating group matrices works", {
  expect_no_error(create.group.matrices(counts=counts, meta = meta))
})

test_that("creating group matrices returns a list", {
  expect_type(create.group.matrices(counts=counts, meta = meta), "list")
})

test_that("computing heterogeneity function works without error", {
  expect_no_error(compute_eITH(x = matrices))
})

test_that("epiCHAOS function works without error", {
  expect_no_error(epiCHAOS(counts=counts, meta = meta, n=10, subsample = 3))
})


