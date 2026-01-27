linters <- lintr::linters_with_defaults(
  line_length_linter = lintr::line_length_linter(120),
  object_name_linter = NULL,
  object_length_linter = NULL
)

exclusions <- list(
  "vignettes",
  "man",
  "docs",
  "inst/doc",
  "revdep",
  ".Rproj.user",
  ".beads"
)
