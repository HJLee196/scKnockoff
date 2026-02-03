# dev/dev_history.R
usethis::use_git_ignore(c("*.Rcheck/",
                          "*.tar.gz",
                          "revdep/",
                          ".Rproj.user",
                          ".Rproj.user/",
                          ".Rhistory",
                          ".RData",
                          ".Ruserdata",
                          ".DS_Store",
                          "*.zip",
                          "*.rds",
                          "*.rda"
                          ))
usethis::use_build_ignore("scKnockoff\\.Rproj")


# 1) deps
attachment::att_amend_desc(
  extra.suggests = c("Seurat", "MAST", "glmnet"),
  update.config = TRUE
)

# 2) docs
devtools::document()

# 3) tests
devtools::test()

# 4) examples
devtools::run_examples()

# 5) as-cran check
withr::with_options(
  list(repos = c(CRAN = "https://cloud.r-project.org/")),
  {
    callr::default_repos()
    rcmdcheck::rcmdcheck(args = c("--no-manual", "--as-cran"))
  }
)

# 6) tags + cleanliness
checkhelper::find_missing_tags()
checkhelper::check_clean_userspace()

# 7) spell + urls
spelling::spell_check_package()
urlchecker::url_check()
