if (requireNamespace("tinytest", quietly = TRUE) &&
    utils::packageVersion("tinytest") >= "1.2.2") {

    ## Set a seed to make the test deterministic
    set.seed(808)
    is_at_home <- function() {
        identical(tolower(Sys.getenv("TT_AT_HOME")), "true")
    }
    ## only run the following tests if at home (not by cran)
    tinytest::test_package("intsurv",
                           ncpu = NULL,
                           test_dir = "rcpp-tests",
                           at_home = is_at_home())

    ## run anyway
    ## tinytest::test_package("intsurv",
    ##                        ncup = getOption("Ncpus", 1),
    ##                        side_effects = TRUE)
}
