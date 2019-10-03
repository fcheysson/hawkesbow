#' Simple matrix multiplication
#'
#' @param A Matrix
#' @param B Matrix
#'
#' @return Product of matrices
#' @export
#'
#' @examples
#' A <- matrix(1:9, 3, 3)
#' B <- matrix(11:19, 3, 3)
#' matrix_mult(A, B)
matrix_mult <- function(A, B) {
  return (A %*% B)
}
