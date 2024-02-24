
if ( requireNamespace("tinytest", quietly=TRUE) ){
  home <- length(unclass(packageVersion("mig"))[[1]]) == 4
  tinytest::test_package("mig", at_home = home)
}

