#!/usr/bin/env Rscript
# Test that all required packages load correctly

pass <- TRUE

check_package <- function(pkg, load_fn = library) {
  tryCatch({
    load_fn(pkg, character.only = TRUE)
    cat(sprintf("  [PASS] %s\n", pkg))
    TRUE
  }, error = function(e) {
    cat(sprintf("  [FAIL] %s: %s\n", pkg, conditionMessage(e)))
    FALSE
  })
}

cat("=== Package load tests ===\n")
pkgs <- c("tidyverse", "vroom", "dplyr", "glue", "logger", "argparse", "arrow", "box", "hyprcoloc")
results <- vapply(pkgs, check_package, logical(1))
pass <- all(results)

cat("\n=== Functional tests ===\n")

# Test HyPrColoc with minimal example
tryCatch({
  betas <- matrix(c(0.5, 0.3, 0.4, 0.6, 0.2, 0.3), nrow = 3)
  ses   <- matrix(c(0.1, 0.1, 0.1, 0.1, 0.1, 0.1), nrow = 3)
  ld    <- diag(3)
  res <- hyprcoloc::hyprcoloc(betas, ses, trait.names = c("t1", "t2"), snp.id = c("s1", "s2", "s3"), ld.matrix = ld)
  cat("  [PASS] hyprcoloc() runs without error\n")
}, error = function(e) {
  cat(sprintf("  [FAIL] hyprcoloc() functional test: %s\n", conditionMessage(e)))
  pass <<- FALSE
})

# Test arrow
tryCatch({
  tmp <- tempfile(fileext = ".parquet")
  arrow::write_parquet(data.frame(x = 1:3, y = letters[1:3]), tmp)
  df <- arrow::read_parquet(tmp)
  stopifnot(nrow(df) == 3)
  cat("  [PASS] arrow read/write parquet\n")
}, error = function(e) {
  cat(sprintf("  [FAIL] arrow parquet test: %s\n", conditionMessage(e)))
  pass <<- FALSE
})

# Test box
tryCatch({
  box::use(glue[glue])
  stopifnot(glue("hello {x}", x = "world") == "hello world")
  cat("  [PASS] box::use() works\n")
}, error = function(e) {
  cat(sprintf("  [FAIL] box::use() test: %s\n", conditionMessage(e)))
  pass <<- FALSE
})

cat("\n=== Summary ===\n")
if (pass) {
  cat("All tests PASSED\n")
  quit(status = 0)
} else {
  cat("Some tests FAILED\n")
  quit(status = 1)
}
