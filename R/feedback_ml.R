 # --------------------------------------------------------------- #
#  Testing for Feedback Effects in a Simultaneous Equation Model  #
# --------------------------------------------------------------- #

#' Testing for Feedback Effects in a Simultaneous Equation Model
#'
#' @param data the data frame containing the data
#' @param out.decompose the decomposition object resulting from \code{causal_decompose()}
#' @param eq.id the equation to be tested for feedback effects
#' @param lb lower bound of the parameter space required for \code{gosolnp}
#' @param ub upper bound of the parameter space required for \code{gosolnp}
#' @param nrestarts number of solver restarts (as in \code{gosolnp})
#' @param nsim number of random parameters to generate for every restart of the solver (as in \code{gosolnp})
#' @param seed.in seed number for gosolnp routine
#' @return A list with components \itemize{
#'\item \code{rho.est}: a data frame with the maximum likelihood estimate of \eqn{rho} and the
#'equations with which each element is involved in feedback-like mechanisms
#'\item \code{loglik}: the value of the log-likelihood of the model
#'\item \code{theta.hessian}: the hessian matrix for the estimated parameters
#'\item \code{rho.jacobian}: the Jacobian matrix of \eqn{\rho} with respect to the entire set of parameters
#'\item \code{wald}: the resulting Wald test statistic
#' }
#'
#' @examples
#' \donttest{
#' data("macroIT")
#' eq.system = list(
#'               eq1 = C ~  CP  + I + CP_1,
#'               eq2 = I ~ K + CP_1,
#'               eq3 = WP ~ I + GDP + GDP_1,
#'               eq4 = GDP ~ C + I + GDP_1,
#'               eq5 = CP ~ WP + T,
#'               eq6 = K ~ I + K_1)
#'
#' instruments = ~ T + CP_1 + GDP_1 + K_1
#'
#' c.dec = causal_decompose(data = macroIT,
#'                          eq.system = eq.system,
#'                          resid.est = "noDfCor",
#'                          instruments = instruments)
#'
#' feedback_ml(data = macroIT,
#'               out.decompose = c.dec,
#'               eq.id = 5,
#'               lb = -200,
#'               ub = 200,
#'               nrestarts = 10,
#'               nsim = 20000,
#'               seed.in = 1)
#'}
#' @export

feedback_ml = function(data, out.decompose, eq.id,
                       lb = -200, ub = 200, nrestarts = 10,
                       nsim = 20000, seed.in = 1) {

    eq.system = out.decompose$eq.system

    if (eq.id < 1 || eq.id > length(eq.system)) {
      stop("Argument 'eq.id' must be between 1 and the length of the system of equation.")
    }

    if (eq.id %% 1 != 0) {
      stop("Argument 'eq.id' must be an integer value.")
    }

    if (eq.id %% 1 != 0) {
      stop("Argument 'eq.id' must be an integer value.")
    }

    # -------------------------------------------- #
    #  extract all the variables (no repetitions)  #
    # -------------------------------------------- #

    all.vars = unique(unlist(lapply(eq.system, function(x) {
      all.vars(x)
    })))

    # ----------------------------- #
    #  extract dependent variables  #
    # ----------------------------- #

    y.vars = unlist(lapply(lapply(eq.system, function(x) {
      all.vars(x)
    }),
    function(x) {
      x[1]
    }))

    # ------------------------------------ #
    #  extract purely exogenous variables  #
    # ------------------------------------ #

    x.vars = setdiff(all.vars, y.vars)

    # ----------------------------------- #
    #  reformulate in alphabetical order  #
    # ----------------------------------- #

    sorted.rhs = lapply(lapply(lapply(eq.system, function(x) {
      all.vars(x)
    }),
    function(x) {
      x[-1]
    }), function(x) {
      sort(x)
    })

    sorted.eqn =
      mapply(function(x, y) {
        stats::reformulate(termlabels = x, response = y)
      },
      sorted.rhs, as.list(y.vars))

    Gamma = out.decompose$Gamma
    A = out.decompose$A
    Sigma = out.decompose$Sigma

    l = eq.id

    # ----------------- #
    #  Starting Values  #
    # ----------------- #

    U = chol(Sigma)

    pos.gamma = matrix(as.numeric(Gamma != 0),
                       nrow = nrow(Gamma),
                       ncol = ncol(Gamma))

    pos.U = matrix(as.numeric(U != 0),
                   nrow = nrow(U),
                   ncol = ncol(U))

    dg = sum(pos.gamma)
    ds = nrow(Sigma)*(nrow(Sigma)+1)/2
    D = dg + ds

    L.gamma = nrow(pos.gamma)
    L.sigma = nrow(pos.U)

    initial.values = numeric(D)
    initial.values = c(Gamma[which((pos.gamma != 0))], (U[which(U != 0)]))
    fb.mat = rho_calc(l = l, Gamma, A, Sigma)

    S1 = fb.mat$S1
    S2 = fb.mat$S2
    S0 = fb.mat$S0

    rho = fb.mat$rho.tbl
    if(all(rho[,2]==0)){
      stop("No feedback effect to test for this equation")
    }

    # ------------------------------ #
    #  Data for likelihood function  #
    # ------------------------------ #

    intercepts = unlist(lapply(out.decompose$systemfit$eq, function(x) {
      x$coefficients[1]
    }))

    y1 = t(as.matrix(data[, which(colnames(data) == y.vars[l])])) - c(intercepts[l])
    y2 = data[, which(colnames(data) %in% y.vars[-l])]
    y2 = as.matrix(y2[, y.vars[-l]]) - t(matrix(intercepts[-l], length(intercepts[-l]), length(y1)))
    z = data[, which(colnames(data) %in% x.vars)]
    z = as.matrix(z[, x.vars])

    ys = S1 %*% t(y2)

    w = rbind(y1 , ys , t(z))

    # -------------------------------------------- #
    # Avoid promise already under evaluation error #
    # -------------------------------------------- #

    l0 = l
    L.gamma0 = L.gamma
    L.sigma0 = L.sigma
    A0 = A
    pos.gamma0 = pos.gamma
    pos.U0 = pos.U
    S00 = S0
    S10 = S1
    S20 = S2

    fn = function(param,
                   X,
                   l = l0,
                   L.gamma = L.gamma0,
                   L.sigma = L.sigma0,
                   A = A0,
                   pos.gamma = pos.gamma0,
                   pos.U = pos.U0,
                   S1 = S10,
                   S2 = S20,
                   S0 = S00) {

      # ------- #
      #  gamma  #
      # ------- #

      gamma =
        array(0, c(L.gamma, L.gamma), dimnames = list(
          paste("X.", 1:L.gamma, sep = ""),
          paste("X.", 1:L.gamma, sep = "")
        ))

      contg = sum(pos.gamma)


      gamma[pos.gamma != 0] = param[1:contg]

      # ------- #
      #  sigma  #
      # ------- #

      U =
        array(0, c(L.sigma, L.sigma), dimnames = list(
          paste("X.", 1:L.sigma, sep = ""),
          paste("X.", 1:L.sigma, sep = "")
        ))

      contu = contg + L.sigma*(L.sigma+1)/2

      U[pos.U != 0] = param[(contg + 1):contu]

      sigma = crossprod(U)

      sigma11 = sigma[l, l]
      sigma12 = t(sigma[l,-l])
      sigma21 = t(sigma12)
      Sigma22 = sigma[-l,-l]

      Gamma22 = gamma[-l,-l]
      gamma11 = gamma[l, l]
      gamma12 = t(gamma[l,-l])
      gamma21 = t(t(gamma[-l, l]))

      a1 = A[l,]
      A2 = A[-l,]

      B = A2 + gamma21 %*% a1
      c = gamma21 + sigma21 / sigma11
      M = diag(ncol(Gamma22)) - gamma21 %*% gamma12 - Gamma22
      Minv = MASS::ginv(M)

      fr = Minv %*% c
      Omega = Minv %*% (Sigma22 - (sigma21 %*% sigma12) / sigma11) %*% t(Minv)

      gammas = S1 %*% t(gamma12)

      fs = S1 %*% fr

      as = S2 %*% (a1)

      Pi0 = S1 %*% Minv %*% B

      Omega0 = (S1 %*% Omega %*% t(S1))

      alpha = c(1, -t(gammas), -a1)
      beta =
        (cbind((-fs),
               diag(length(fs)) + ((fs) %*% t(gammas)),
               ((fs) %*% t(a1)) - Pi0))

      # ---------------- #
      #  -loglikelihood  #
      # ---------------- #
      as.numeric(-(-0.5 * ncol(X) * log(sigma11) - 0.5 * ncol(X) * log(det(Omega0)) -
        (1 / (2 * sigma11)) * (t(alpha) %*% (X %*% t(X)) %*% alpha) -
        0.5 * matrixcalc::matrix.trace((t(beta) %*% MASS::ginv(Omega0) %*% beta) %*% (X %*% t(X)))
      ))
    }

    # -------------------------- #
    #  rho function for Jacobian #
    # -------------------------- #

    eqn = function(param,
                    X,
                    l = l0,
                    L.gamma = L.gamma0,
                    L.sigma = L.sigma0,
                    A = A0,
                    pos.gamma = pos.gamma0,
                    pos.U = pos.U0,
                    S1 = S10,
                    S2 = S20,
                    S0 = S00)

    {
      cont = 0

      # ------- #
      #  gamma  #
      # ------- #

      gamma =
        array(0, c(L.gamma, L.gamma), dimnames = list(
          paste("X.", 1:L.gamma, sep = ""),
          paste("X.", 1:L.gamma, sep = "")
        ))

      contg = sum(pos.gamma)
      gamma[pos.gamma != 0] = param[1:contg]


      # ------- #
      #  sigma  #
      # ------- #

      U =
        array(0, c(L.sigma, L.sigma), dimnames = list(
          paste("X.", 1:L.sigma, sep = ""),
          paste("X.", 1:L.sigma, sep = "")
        ))

      contu = contg + sum(pos.U)
      U[pos.U != 0] = param[(contg + 1):contu]

      sigma = crossprod(U)

      sigma11 = sigma[l, l]
      sigma12 = t(sigma[l,-l])
      sigma21 = t(sigma12)
      Sigma22 = sigma[-l,-l]

      Gamma22 = gamma[-l,-l]
      gamma11 = gamma[l, l]
      gamma12 = t(gamma[l,-l])
      gamma21 = t(t(gamma[-l, l]))

      a1 = A[l,]
      A2 = A[-l,]

      B = A2 + gamma21 %*% a1
      c = gamma21 + sigma21 / sigma11
      M = diag(ncol(Gamma22)) - (gamma21 %*% gamma12) - Gamma22
      Minv = solve(M)
      fr = Minv %*% c
      Omega = Minv %*% (Sigma22 - (sigma21 %*% sigma12) / sigma11) %*% t(Minv)

      gammas = S1 %*% t(gamma12)

      fs = S1 %*% fr

      fss = S0 %*% fs
      gammass = S0 %*% gammas

      rho = fs * gammas

      return(c(rho))
    }

    # ----------------------------#
    #  Optimization - multistart  #
    # --------------------------- #

    fit.ur=Rsolnp::gosolnp(fun=fn,
                  LB = rep(lb,length(initial.values)),
                  UB = rep(ub, length(initial.values)),
                  l = l0,
                  X = w,
                  L.gamma = L.gamma0,
                  L.sigma = L.sigma0,
                  A = A0,
                  pos.gamma = pos.gamma0,
                  pos.U = pos.U0,
                  S1 = S10,
                  S2 = S20,
                  S0 = S00,
                  distr = rep(1,length(initial.values)),
                  n.restarts = nrestarts,
                  n.sim = nsim,
                  rseed = seed.in)

    pars = fit.ur$pars
    logl = fit.ur$values[length(fit.ur$values)]

    # Jacobian matrix

    J.mat = numDeriv::jacobian(
      eqn,
      x = pars,
      X = w,
      l = l0,
      L.gamma = L.gamma0,
      A = A0,
      pos.gamma = pos.gamma0,
      pos.U = pos.U0,
      S1 = S10,
      S2 = S20,
      S0 = S00
    )

    # Hessian matrix

    hessian = numDeriv::hessian(
      func = fn,
      x = pars,
      X = w,
      l = l0,
      L.gamma = L.gamma0,
      L.sigma = L.sigma0,
      A = A0,
      pos.gamma = pos.gamma0,
      pos.U = pos.U0,
      S1 = S10,
      S2 = S20,
      S0 = S00
    )

    Gamma.ml = pos.gamma
    U.ml = pos.U

    gamma.ml = pars[1:dg]
    Gamma.ml[Gamma.ml == 1] = gamma.ml

    u.ml = pars[(dg + 1):D]
    U.ml[U.ml == 1] = u.ml
    Sigma.ml = crossprod(U.ml)

    # Rho estimate

    rho.ml = rho_calc(l,
                      Gamma = Gamma.ml,
                      Sigma = Sigma.ml,
                      A = A)$rho.tbl

    wald = t(rho.ml[,2])%*%MASS::ginv((J.mat)%*%MASS::ginv(hessian)%*%t(J.mat))%*%rho.ml[,2]

    return(list(
      rho.est = rho.ml,
      loglik = logl,
      theta.hessian = hessian,
      rho.jacobian = J.mat,
      wald = wald
    ))
  }
