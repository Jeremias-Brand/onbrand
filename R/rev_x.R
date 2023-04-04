#' Reverse a numeric value within a specified range
#'
#' This function takes a numeric value `x` and an optional range `xrange` as input. It then reverses
#' the position of `x` within the range, returning the new value. If `xrange` is not provided, it defaults
#' to the range specified in `CELL_META$xlim`.
#'
#' @param x A numeric value to be reversed within the specified range.
#' @param xrange A numeric vector of length 2 defining the range within which `x` should be reversed.
#'   Defaults to `CELL_META$xlim`.
#' @return A numeric value representing the reversed position of `x` within the specified range.
#' @export
#' @examples
#' rev_x(5, c(0, 10))  # Should return 5 (5 is reversed to 5 within the range 0 to 10)
#' rev_x(2, c(0, 10))  # Should return 8 (2 is reversed to 8 within the range 0 to 10)
#' rev_x(6)            # Assumes CELL_META$xlim is set and reverses 6 within that range
rev_x <- function(x, xrange = CELL_META$xlim) {
  xrange[2] - x + xrange[1]
}
