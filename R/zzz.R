## .onAttach and .onLoad functions go here:
## http://r-pkgs.had.co.nz/r.html

## Ben's question about how to get drule additions to work in a pkg:
## https://github.com/sgsokol/Deriv/issues/13
## cf. fitode/R/zzz.R
## You need at least one "@importFrom Deriv drule" somewhere
## in order for this to work.
## Note that devtools::document() fails.  You need to
## roxygen2::roxygenise() at least once.

.onLoad <- function(libname, pkgname) {
    ## differentiation rules for Deriv package
    drule[["gerf"]] <- alist(x=dgerf(x,n), n=NULL)
    drule[["dgerf"]] <- alist(x=d2gerf(x,n), n=NULL)
}

.onAttach <- function(libname, pkgname) {
  packageStartupMessage(paste0(
    "Welcome to evolvr\n"
    , "Version: ", packageVersion("evolvr"), "\n"
    , "Authors:  \n  ", get_author_email()[1], "\n  ", get_author_email()[2], "\n  ", get_author_email()[3]
    ))
}

get_author_email <- function(pkg=packageName()) {
  return(eval(parse(text=packageDescription(pkg,fields="Authors@R"))))
}

