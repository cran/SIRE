# --------------------------------------------- #
#  Decomposition starting from Gamma and Sigma  #
# --------------------------------------------- #

#' Decomposition starting from Gamma and Sigma
#'
#' @description Function to decompose \eqn{\Gamma'} into recursive and
#' interdependent sub-matrices (internal use)
#' @param Gamma the \eqn{\Gamma'} matrix.
#' @param Sigma the \eqn{\Sigma} matrix.
#'
#' @return A list with components \itemize{
#'\item \code{C}: the matrix highlighting the interdependent mechanisms at deterministic level.
#'\item \code{Psi1}: the matrix highlighting the interdependent mechanisms at stochastic level.
#'\item \code{Psi0}: the matrix highlighting the causal mechanisms.
#'\item \code{powers}: a list containing the matrix powers of \eqn{\Gamma'}.
#' }
#' @export


dec.calc = function(Gamma,
                    Sigma){

  if(sum(diag(Gamma))!= 0 || nrow(Gamma)!=ncol(Gamma)){
    stop("Invalid structure for Gamma, must be hollow and square")
  }

  if(!matrixcalc::is.symmetric.matrix(Sigma)){
    stop("invalid structure for Sigma.")
  }

  gamma.b = (Gamma!=0)+0
  L = nrow(Gamma)
  R = matrix(0, nrow = L, ncol = L)
  powers = list()
  LL = L - 1
  for (q in 1:LL) {
    powers[[q]] = matrixcalc::matrix.power(gamma.b,q)
  }

  R = as.matrix((t(Reduce('+', powers)) != 0) + 0)
  U = matrix(1, nrow = L , ncol = L)
  C = matrixcalc::hadamard.prod(gamma.b, R)
  Psi = matrixcalc::hadamard.prod(gamma.b, ((U - R)))
  C.b = as.matrix((C != 0) + 0)
  I = diag(L)
  sigma.b = (Sigma != 0) + 0
  IR = (I + R)
  P = as.matrix((IR * sigma.b) != 0) + 0
  Psi1 = Psi * t(P)
  Psi0 = Psi * (U - t(P))
  return(list(
    C = C , Psi1 = Psi1 , Psi0 = Psi0 , powers = powers
  ))
}
