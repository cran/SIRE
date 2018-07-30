# --------------------------------------------------------------- #
#  Estimation and Decomposition of a Simultaneous Equation Model  #
#       into its Recursive and Interdependent Sub-Systems         #
# --------------------------------------------------------------- #
#' Estimation and decomposition of simultaneous equation model
#'
#' @description Estimate and decompose
#' a Simultaneous Equation Model into its recursive
#' and Interdependent sub-systems
#'
#' @param data the data frame containing the data
#' @param eq.system the system of equations (a list of formula objects, e.g. as in pkg systemfit)
#' @param resid.est the estimation methods for the residual covariance matrix (as in pkg systemfit)
#' @param instruments the intruments used to estimate the model via 3-SLS (as in pkg systemfit)
#' @param p.adj logical indication of whether to adjust the significance level
#'  for multiple comparisons, when testing the significance of the elements of \eqn{\Sigma}
#' @param alpha the significance level \eqn{\alpha} for the tests on the elements of \eqn{\Sigma}
#'
#'@return A list with components \itemize{
#'\item \code{output}: a list containing\itemize{
#'\item \code{eq.system}: the system of equations given as input
#'\item \code{Gamma}: the 3-SLS estimate of \eqn{\Gamma'}
#'\item \code{C}: the matrix highlighting the interdependent mechanisms at deterministic level.
#'\item \code{Psi1}: the matrix highlighting the interdependent mechanisms at stochastic level.
#'\item \code{Psi0}: the matrix highlighting the causal mechanisms.
#'\item \code{A}: the 3-SLS estimate of \eqn{A}
#'\item \code{Sigma}: the 3-SLS estimate of \eqn{Sigma}, with non-significant
#'elements set to zero
#'\item \code{systemfit}: the output from the \code{systemfit} function used to
#'estimate the model
#'\item \code{corrtest}: the output from the \code{corr.test} function used to
#'test the significance on the elements of the residual correlation matrix
#'\item \code{all.graph}: the path diagram of the model, using the package \code{igraph}
#'\item \code{dec.graph}: the path diagram of the decomposed model, with color
#'coding for each vertex}
#'\item \code{path}: a data frame in which every row is a path coefficient,
#'along with indication of the nature of the link (recursive/interdependent)
#'and the index of the equation in which it belongs
#'}
#'
#' @examples
#' data("macroIT")
#' eq.system = list(
#'                eq1 = C ~  CP  + I + CP_1,
#'                eq2 = I ~ K + CP_1,
#'                eq3 = WP ~ I + GDP + GDP_1,
#'                eq4 = GDP ~ C + I + GDP_1,
#'                eq5 = CP ~ WP + T,
#'                eq6 = K ~ I + K_1)
#'
#' instruments = ~ T + CP_1 + GDP_1 + K_1
#'
#' causal.decompose(data = macroIT,
#'                eq.system = eq.system,
#'                resid.est = "noDfCor",
#'                instruments = instruments,
#'                p.adj = TRUE,
#'                alpha = 0.05)
#'@export
causal.decompose = function(data,
                            eq.system,
                            resid.est = "noDfCor",
                            instruments,
                            p.adj = TRUE,
                            alpha = 0.05) {

  if (alpha >= 1 || alpha <= 0) {
    stop("invalid argument alpha. Must be in (0,1).")
  }

  if (!is.logical(p.adj)) {
    stop("invalid argument p.adj. Must be logical.")
  }

  # -------------------------------------------- #
  #  extract all the variables (no repetitions)  #
  # -------------------------------------------- #

  all.vars <- unique(unlist(lapply(eq.system, function(x) {
    all.vars(x)
  })))

  # ----------------------------- #
  #  extract dependent variables  #
  # ----------------------------- #

  y.vars <- unlist(lapply(lapply(eq.system, function(x) {
    all.vars(x)
  }),
  function(x) {
    x[1]
  }))

  # ------------------------------------ #
  #  extract purely exogenous variables  #
  # ------------------------------------ #

  x.vars <- setdiff(all.vars, y.vars)

  # ----------------------------------- #
  #  reformulate in alphabetical order  #
  # ----------------------------------- #

  sorted.rhs <-
    lapply(lapply(lapply(eq.system,
                         function(x) {
                           all.vars(x)
                         }),
                  function(x) {
                    x[-1]
                  }),
           function(x) {
             sort(x)
           })

  sorted.eqn <-
    mapply(function(x, y) {
      stats::reformulate(termlabels = x, response = y)
    },
    sorted.rhs, as.list(y.vars))

  # ---------------------------- #
  #  create binary matrix Gamma  #
  # ---------------------------- #

  L = length(y.vars)

  gamma.b = matrix(0, nrow = L , ncol = L)

  for (r in 1:L) {
    rhs.r = all.vars(sorted.eqn[[r]])[-1]
    flag = match(intersect(y.vars, rhs.r), y.vars)
    gamma.b[r, c(flag)] = 1
  }

  rownames(gamma.b) = c(y.vars)
  colnames(gamma.b) = c(y.vars)

  # ------------------------ #
  #  create binary matrix A  #
  # ------------------------ #

  LX = length(x.vars)

  Ax = matrix(0, nrow = L , ncol = LX)

  for (r in 1:L) {
    rhs.rex = all.vars(sorted.eqn[[r]])[-1]
    flag = match(intersect(x.vars, rhs.rex), x.vars)
    Ax[r, c(flag)] = 1
  }

  # --------------- #
  #  overall graph  #
  # --------------- #

  all.graph <- igraph::graph.adjacency(t(gamma.b))

  # ----------------------- #
  #  decomposition C + Psi  #
  # ----------------------- #

  R = matrix(0, nrow = L, ncol = L)
  powers = list()
  LL = L - 1
  for (q in 1:LL) {
    powers[[q]] = matrixcalc::matrix.power(gamma.b, q)
  }

  R = as.matrix((t(Reduce('+', powers)) != 0) + 0)
  U = matrix(1, nrow = L , ncol = L)
  C = matrixcalc::hadamard.prod(gamma.b, R)
  Psi = matrixcalc::hadamard.prod(gamma.b, ((U - R)))
  C.b = as.matrix((C != 0) + 0)
  I = diag(L)
  C + Psi == gamma.b

  method="3SLS"

  sigma = diag(length(eq.system))
  fiteq <-
    systemfit::systemfit(
      data = data,
      formula = eq.system,
      method = method,
      inst = instruments,
      methodResidCov = resid.est,
      method3sls = "GMM",
      residCovWeighted = F
    )

  vcv = stats::cov2cor(fiteq$residCovEst)
  corrtest = psych::corr.p(n = ncol(data), r = vcv)
  if (p.adj == T)
  {
    tests = (crossprod(matrixcalc::upper.triangle(
      matrix(
        as.numeric(corrtest$p < alpha),
        ncol = length(eq.system),
        nrow = length(eq.system)
      )
    )) != 0) + 0
  } else{
    tests = (crossprod(matrixcalc::lower.triangle(
      matrix(
        as.numeric(corrtest$p < alpha),
        ncol = length(eq.system),
        nrow = length(eq.system)
      )
    )) != 0) + 0
  }

  sigma = (tests * ((fiteq$residCovEst)))

  # --------------------------- #
  #  decomposition Psi0 + Psi1  #
  # --------------------------- #

  if (!Matrix::isDiagonal(sigma)) {
    sigma.b = as.matrix((sigma != 0) + 0)
    IR.b = as.matrix(((I + R) != 0) + 0)
    P = as.matrix(IR.b * sigma.b != 0) + 0
    Psi1 = Psi * t(P)
    Psi0 = Psi * (U - t(P))
    rownames(C) = c(y.vars)
    colnames(C) = c(y.vars)
    rownames(Psi1) = c(y.vars)
    colnames(Psi1) = c(y.vars)
    rownames(Psi0) = c(y.vars)
    colnames(Psi0) = c(y.vars)

    C.graph <- igraph::graph.adjacency(t(C))
    Psi1.graph <- igraph::graph.adjacency(t(Psi1))
    Psi0.graph <- igraph::graph.adjacency(t(Psi0))

    # ------------------------------------------- #
    #  create decomposed graph:                   #
    #  red = recursive, green = true causal       #
    #  dashed red = induced recursion from Sigma  #
    # ------------------------------------------- #

    dec.graph = all.graph
    igraph::E(dec.graph)[match(intersect(igraph::as_ids(igraph::E(dec.graph)),
                                         igraph::as_ids(igraph::E(C.graph))),
                               igraph::as_ids(igraph::E(dec.graph)))]$lty = 1
    igraph::E(dec.graph)[match(intersect(igraph::as_ids(igraph::E(dec.graph)),
                                         igraph::as_ids(igraph::E(C.graph))),
                               igraph::as_ids(igraph::E(dec.graph)))]$color = "red"

    igraph::E(dec.graph)[match(intersect(igraph::as_ids(igraph::E(dec.graph)),
                                         igraph::as_ids(igraph::E(Psi1.graph))),
                               igraph::as_ids(igraph::E(dec.graph)))]$lty = 2
    igraph::E(dec.graph)[match(intersect(igraph::as_ids(igraph::E(dec.graph)),
                                         igraph::as_ids(igraph::E(Psi1.graph))),
                               igraph::as_ids(igraph::E(dec.graph)))]$color = "red"

    igraph::E(dec.graph)[match(intersect(igraph::as_ids(igraph::E(dec.graph)),
                                         igraph::as_ids(igraph::E(Psi0.graph))),
                               igraph::as_ids(igraph::E(dec.graph)))]$lty = 1
    igraph::E(dec.graph)[match(intersect(igraph::as_ids(igraph::E(dec.graph)),
                                         igraph::as_ids(igraph::E(Psi0.graph))),
                               igraph::as_ids(igraph::E(dec.graph)))]$color = "blue"

    # ----------------------------------------- #
    #  complete linked pairs,                   #
    #  '~>' means recursive,                    #
    #  '->' means true causal                   #
    #  '=>' means induced recursion from Sigma  #
    # ----------------------------------------- #

    path = sort(c(
      stringr::str_replace(igraph::as_ids(igraph::E(C.graph)), "\\|", " ~> "),
      stringr::str_replace(igraph::as_ids(igraph::E(Psi1.graph)), "\\|", " => "),
      stringr::str_replace(igraph::as_ids(igraph::E(Psi0.graph)), "\\|", " -> ")
    ))

    syst.beta = lapply(fiteq$eq, function(x) {
      x$coefficients[-1]
    })

    syst.gamma = lapply(syst.beta, function(x) {
      x[names(x) %in% y.vars]
    })

    syst.a = lapply(syst.beta, function(x) {
      x[!(names(x) %in% y.vars)]
    })

    lengths.gamma = unlist(lapply(syst.gamma, length))
    lengths.a = unlist(lapply(syst.a, length))

    vals = unlist(syst.gamma)

    non.zero.g = data.frame(which(gamma.b != 0, arr.ind = T),row.names = NULL)
    non.zero.g = non.zero.g[order(non.zero.g[,1],non.zero.g[,2]),]
    non.zero.g = as.matrix(non.zero.g)

    Gamma.m <-
      as.matrix((Matrix::sparseMatrix(
        i = non.zero.g[, 1], j = non.zero.g[, 2], x = vals
      )))

    C.m = matrixcalc::hadamard.prod(Gamma.m, C)

    Psi0.m = matrixcalc::hadamard.prod(Gamma.m, Psi0)

    Psi1.m = matrixcalc::hadamard.prod(Gamma.m, Psi1)

    vals.a = unlist(syst.a)

    non.zero.a = data.frame(which(Ax != 0, arr.ind = T),row.names = NULL)
    non.zero.a[order(non.zero.a[,1],non.zero.a[,2]),]
    non.zero.a = as.matrix(non.zero.a)

    A.m <-
      as.matrix((Matrix::sparseMatrix(
        i = non.zero.a[, 1], j = non.zero.a[, 2], x = vals.a
      )))

    output = list(
      eq.system = eq.system,
      Gamma = Gamma.m,
      C = C.m,
      Psi1 = Psi1.m,
      Psi0 = Psi0.m,
      A = A.m,
      Sigma = as.matrix(sigma),
      systemfit = fiteq,
      corrtest = corrtest,
      all.graph = all.graph,
      dec.graph = dec.graph
    )

  } else{
    # ------------------- #
    #  create sub-graphs  #
    # ------------------- #

    rownames(C) = c(y.vars)
    colnames(C) = c(y.vars)
    rownames(Psi) = c(y.vars)
    colnames(Psi) = c(y.vars)

    C.graph <- igraph::graph.adjacency(t(C))
    Psi.graph <- igraph::graph.adjacency(t(Psi))

    # -------------------------------------- #
    #  create decomposed graph:              #
    #  red = recursive, green = true causal  #
    # -------------------------------------- #

    dec.graph = all.graph
    igraph::E(dec.graph)[match(intersect(igraph::as_ids(igraph::E(dec.graph)),
                                 igraph::as_ids(igraph::E(C.graph))),
                       igraph::as_ids(igraph::E(dec.graph)))]$color = "red"

    igraph::E(dec.graph)[match(intersect(igraph::as_ids(igraph::E(dec.graph)),
                                 igraph::as_ids(igraph::E(Psi.graph))),
                       igraph::as_ids(igraph::E(dec.graph)))]$color = "blue"

    # ------------------------ #
    #  complete linked pairs,  #
    #  '~>' means recursive,   #
    #  '->' means true causal  #
    # ------------------------ #

    path = sort(c(
      stringr::str_replace(igraph::as_ids(igraph::E(C.graph)), "\\|", " ~> "),
      stringr::str_replace(igraph::as_ids(igraph::E(Psi.graph)), "\\|", " -> ")
    ))

    syst.beta = lapply(fiteq$eq, function(x) {
      x$coefficients[-1]
    })

    syst.gamma = lapply(syst.beta, function(x) {
      x[names(x) %in% y.vars]
    })

    syst.a = lapply(syst.beta, function(x) {
      x[!(names(x) %in% y.vars)]
    })

    lengths.gamma = unlist(lapply(syst.gamma, length))

    lengths.a = unlist(lapply(syst.a, length))

    vals = unlist(syst.gamma)

    non.zero.g = data.frame(which(gamma.b != 0, arr.ind = T),row.names = NULL)
    non.zero.g = non.zero.g[order(non.zero.g[,1],non.zero.g[,2]),]
    non.zero.g = as.matrix(non.zero.g)

    Gamma.m <-
      as.matrix((Matrix::sparseMatrix(
        i = non.zero.g[, 1], j = non.zero.g[, 2], x = vals
      )))

    C.m = matrixcalc::hadamard.prod(Gamma.m, C)

    Psi.m = matrixcalc::hadamard.prod(Gamma.m, Psi)

    vals.a = unlist(syst.a)

    non.zero.a = data.frame(which(Ax != 0, arr.ind = T),row.names = NULL)
    non.zero.a = non.zero.a[order(non.zero.a[,1],non.zero.a[,2]),]
    non.zero.a = as.matrix(non.zero.a)

    A.m <-
      as.matrix((Matrix::sparseMatrix(
        i = non.zero.a[, 1], j = non.zero.a[, 2], x = vals.a
      )))

    output = list(
      eq.system = eq.system,
      Gamma = Gamma.m,
      A = A.m,
      Sigma = sigma,
      systemfit = fiteq,
      C = C,
      Psi = Psi,
      all.graph = all.graph,
      dec.graph = dec.graph
    )
  }

  paths = as.data.frame(path)
  paths$type = "Recursive"
  paths$type[grep("-+", as.character(paths$path))] = "Causal"
  paths$type[grep("~+", as.character(paths$path))] = "Interdependent"
  paths$type[grep("=+", as.character(paths$path))] = "Interdependent from Sigma"
  equation.n = rep(0, length(eq.system))
  for (j in 1:length(eq.system)) {
    equation.n[grep(paste(">", y.vars[j]), paths$path)] = j
  }
  paths$equation.n = equation.n

  decomposition = list(output = output, path = paths)
  return(decomposition)
}
