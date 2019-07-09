#' @useDynLib Jeffs
#' @importFrom Rcpp sourceCpp
NULL

#' @title Characteristic Polynomial
#' @rdname charpoly
#' @param A A square matrix
#' @details Compute the characteristic polynomial
#' @return The coefficients of the characteristic polynomial of A as a vector
#' @export
#' @examples
#' n <- 13
#' A <- matrix(c(-1, 0, 0, 1, 1, 0, 0, -1, -1, -1, 0, 0, -1, 0, -1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0, -1, 1, 0, -1, 1, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0, 0, 1, 0, 0, -1, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 1, 1, -1, -1, -1, -1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, -1, 1, 0, 0, 0, 0,1, 0, 0, 0, 0, 0, 1, -1, -1, 0, 1, 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, -1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1),n ,n ,byrow=TRUE)
#' charpoly(A)

charpoly <- function (A) {
  n <- nrow(A)
  B <- diag(1, n, n)
  p <- rep(1, n + 1)
  for (k in seq_len(n)) {
    B <- A %*% B
    p[k + 1] <- -sum(diag(B))/k
    diag(B) <- diag(B) + p[k + 1]
  }
  p
}


#' @title Adjoint matrix
#' @rdname adjoint
#' @param A A square matrix
#' @details Compute the adjoint matrix
#' @return The adjoint matrix of A
#' @export
#' @examples
#' n <- 13
#' A <- matrix(c(-1, 0, 0, 1, 1, 0, 0, -1, -1, -1, 0, 0, -1, 0, -1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0, -1, 1, 0, -1, 1, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0, 0, 1, 0, 0, -1, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 1, 1, -1, -1, -1, -1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, -1, 1, 0, 0, 0, 0,1, 0, 0, 0, 0, 0, 1, -1, -1, 0, 1, 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, -1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1),n ,n ,byrow=TRUE)
#' adjoint(A)
adjoint <- function(A) round(det(A)*solve(A))


#' @title Permanent of a matrix
#' @rdname permanent
#' @param A A square matrix
#' @details Compute the permanent of matrix
#' @return The permanent of matrix A
#' @export
#' @examples
#' n <- 13
#' A <- matrix(c(-1, 0, 0, 1, 1, 0, 0, -1, -1, -1, 0, 0, -1, 0, -1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0, -1, 1, 0, -1, 1, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0, 0, 1, 0, 0, -1, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 1, 1, -1, -1, -1, -1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, -1, 1, 0, 0, 0, 0,1, 0, 0, 0, 0, 0, 1, -1, -1, 0, 1, 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, -1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1),n ,n ,byrow=TRUE)
#' permanent(A)
permanent <- function(A) Jeffs:::permanent_cpp(A)

#' @title Absoulte feedback (T)
#' @rdname Tmat
#' @param A A square matrix
#' @details Compute the absoulte feedback of a matrix
#' @return The absolute feedback of a matrix A
#' @export
#' @examples
#' n <- 13
#' A <- matrix(c(-1, 0, 0, 1, 1, 0, 0, -1, -1, -1, 0, 0, -1, 0, -1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0, -1, 1, 0, -1, 1, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0, 0, 1, 0, 0, -1, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 1, 1, -1, -1, -1, -1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, -1, 1, 0, 0, 0, 0,1, 0, 0, 0, 0, 0, 1, -1, -1, 0, 1, 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, -1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1),n ,n ,byrow=TRUE)
#' Tmat(A)
Tmat <- function(A){

  if(nrow(A)!=ncol(A))stop('A must be a symetrical matrix')
  n <- ncol(A)
  Tmat_out <- matrix(NA,n,n)
  for(i in seq_len(n)){
    for(j in seq_len(n)){
      Tmat_out[j,i] <- permanent(abs(A[-i,-j]))
    }
  }

  return(Tmat_out)

}









