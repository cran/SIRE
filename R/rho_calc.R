# ----------------- #
#  Rho Calculation  #
# ----------------- #
#
#' Rho Calculation
#'
#' @description Function to calculate rho (internal use)
#'
#' @param l the equation index for which to calculate rho
#' @param Gamma the \eqn{\Gamma'} matrix
#' @param A the \eqn{A} matrix
#' @param Sigma the \eqn{\Sigma} matrix
#'
#' @return A list with components \itemize{
#'\item \code{S0}: the selection matrix for \eqn{p_j}.
#'\item \code{S1}: the selection matrix for \eqn{\Gamma'}.
#'\item \code{S2}: the selection matrix. for \eqn{A}
#' }
#' @export

rho.calc = function(l , Gamma , A, Sigma) {
  gamma12 = t(Gamma[l,-l])

  S1 = matrix(0, nrow = length(gamma12[gamma12 != 0]), ncol = nrow(Gamma) - 1)
  count1 = 0
  for (j in 1:ncol(S1)) {
    if (gamma12[j] != 0) {
      S1[count1 + 1, j] = 1
      count1 = count1 + 1
    }
  }

  a1 = A[l,]

  S2 = matrix(0, nrow = length(a1[a1 != 0]), ncol = ncol(A))
  count2 = 0
  for (j in 1:ncol(S2)) {
    if (a1[j] != 0) {
      count2 = count2 + 1
      S2[count2, j] = 1
    }
  }

  sigma11 = Sigma[l, l]
  sigma12 = t(Sigma[l,-l])
  sigma21 = t(sigma12)
  Sigma22 = Sigma[-l,-l]

  Gamma22 = Gamma[-l,-l]
  gamma11 = Gamma[l, l]
  gamma21 = t(t(Gamma[-l, l]))

  A2 = A[-l,]

  B = A2 + gamma21 %*% a1
  c = gamma21 + sigma21 / sigma11
  direct.pos = which(Gamma[l,] != 0)

  M = diag(ncol(Gamma22)) - gamma21 %*% gamma12 - Gamma22
  Minv = solve(M)
  fr = Minv %*% c
  Omega = Minv %*% (Sigma22 - (sigma21 %*% sigma12) / sigma11) %*% t(Minv)
  gammas = S1 %*% (t(gamma12) * t(t((dec.calc(Gamma, Sigma)$C + dec.calc(Gamma, Sigma)$Psi1)[l, -l])))

  fs = S1 %*% fr

  S0 = matrix(0, nrow = length(fs[fs != 0]), ncol = nrow(fs))
  count0 = 0
  for (j in 1:ncol(S0)) {
    if (fs[j] != 0) {
      count0 = count0 + 1
      S0[count0, j] = 1
    }
  }

  fss = S0 %*% fs

  as = S2 %*% (a1)

  Pi = S1 %*% Minv %*% B
  Omegas = S1 %*% Omega %*% t(S1)

  rho = matrixcalc::hadamard.prod(gammas, fs)
  rho.tbl = data.frame(cbind(direct.pos, rho))
  colnames(rho.tbl) = c("Feedback eqn.", "rho.est")

  return(list(
    S0 = S0,
    S1 = S1,
    S2 = S2,
    rho.tbl = rho.tbl
  ))
}
