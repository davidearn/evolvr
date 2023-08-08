## Functions for saving and restoring the random number generator (RNG) state.
## Taken from Ben Bolker's answer here: 
## http://stackoverflow.com/questions/13997444/print-the-current-random-seed-so-that-i-can-enter-it-with-set-seed-later

#' Save random number generator seed to file
#'
#' @note Taken from \href{http://stackoverflow.com/questions/13997444/print-the-current-random-seed-so-that-i-can-enter-it-with-set-seed-later}{Ben Bolker's answer on stackoverflow}.
#' @export
save_rng <- function(savefile=tempfile()) {
  if (exists(".Random.seed"))  {
    oldseed <- get(".Random.seed", .GlobalEnv)
  } else stop("don't know how to save before set.seed() or r*** call")
  oldRNGkind <- RNGkind()
  save("oldseed","oldRNGkind",file=savefile)
  invisible(savefile)
}

#' Restore random number generator seed from file
#'
#' @note Taken from \href{http://stackoverflow.com/questions/13997444/print-the-current-random-seed-so-that-i-can-enter-it-with-set-seed-later}{Ben Bolker's answer on stackoverflow}.
#' @export
restore_rng <- function(savefile) {
  load(savefile)
  do.call("RNGkind",as.list(oldRNGkind))  ## must be first!
  assign(".Random.seed", oldseed, .GlobalEnv)
}

