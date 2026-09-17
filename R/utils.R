#' Resolve the number of threads for natcpp's parallel functions
#'
#' Internal helper implementing the package-wide default thread policy. A
#' function that accepts a \code{threads} argument should default it to
#' \code{NULL} and pass it through here.
#'
#' The policy follows the advice given by Dirk Eddelbuettel in
#' \url{https://github.com/Rdatatable/data.table/issues/5658}: an explicit
#' user value always wins; otherwise fall back to \code{getOption("Ncpus")}
#' (which CRAN sets for its check farm), then the \code{OMP_THREAD_LIMIT}
#' environment variable, and finally a conservative default of 2. Deliberately
#' avoids grabbing every available core by default, and needs no dependency on
#' \pkg{parallel}. To use all cores, pass \code{threads = 0} (the C++ back ends
#' interpret 0 as \code{hardware_concurrency()}).
#'
#' @param threads User-supplied thread count, or \code{NULL} to apply the
#'   default policy. A supplied value (including 0) is returned unchanged.
#' @return An integer thread count.
#' @noRd
natcpp_threads <- function(threads = NULL) {
  if (!is.null(threads)) return(as.integer(threads)[1])
  n <- getOption("Ncpus")
  if (is.null(n)) {
    tl <- Sys.getenv("OMP_THREAD_LIMIT")
    if (nzchar(tl)) n <- tl
  }
  n <- suppressWarnings(as.integer(n))[1]
  if (length(n) == 0L || is.na(n) || n < 1L) n <- 2L
  n
}
