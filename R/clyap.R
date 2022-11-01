#' Continuous Lyapunov Equation Solver
#'
#' Computes the solution of \eqn{XA - A^TX - C = 0} using the Bartels & Stewart 1972 approach.
#' This implementation is equivalent to the SciPy implementation solve_continous_lyapunov.
#'
#' @param A \eqn{n x n} numeric or complex matrix.
#' @param C \eqn{n x n} numeric or complex matrix.
#'
#' @return The solution to the above Lyapunov equation.
#' @export
#'
#' @references
#' TODO
#'
#' @examples
#' A <- matrix(c(-3, -2, 0, -1, -1, 0, 0, -5, -1), ncol = 3, byrow = TRUE)
#' C <- matrix(c(7, 9, -1, 4, 6, 5, 1, 2, 3), ncol = 3, byrow = TRUE)
#' X <- clyap(A, C)
clyap <- function(A, C) {
    return(clyap_internal(A, C))
}
