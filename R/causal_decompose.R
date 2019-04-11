# --------------------------------------------------------------- #
#  Estimation and Decomposition of a Simultaneous Equation Model  #
#       into its Recursive and Interdependent Sub-Systems         #
# --------------------------------------------------------------- #
#' Estimation and decomposition of simultaneous equation model
#'
#' @description Estimate and/or decompose
#' a Simultaneous Equation Model into its recursive
#' and Interdependent sub-systems
#'
#' @param data the data frame containing the data
#' @param eq.system the system of equations (a list of formula objects, e.g. as in pkg \code{systemfit})
#' @param resid.est the estimation methods for the residual covariance matrix (as in \code{systemfit})
#' @param instruments the intruments used to estimate the model via 3-SLS (as in \code{systemfit})
#' @param sigma.in the \eqn{\Sigma} matrix, if the user wants to simulate a particular structure at stochastic level.
#' Overrides 3SLS estimation if specified.
#'
#'@return A list with components \itemize{
#'\item \code{eq.system}: the system of equations given as input
#'\item \code{Gamma}: the 3-SLS estimate of \eqn{\Gamma'}
#'\item \code{C}: the matrix highlighting the interdependent mechanisms at deterministic level.
#'\item \code{Psi1}: the matrix highlighting the interdependent mechanisms at stochastic level.
#'\item \code{Psi0}: the matrix highlighting the causal mechanisms.
#'\item \code{A}: the 3-SLS estimate of \eqn{A}
#'\item \code{Sigma}: the 3-SLS estimate of \eqn{Sigma}
#'\item \code{systemfit}: the output from the \code{systemfit} function used to
#'estimate the model
#'\item \code{all.graph}: the path diagram of the model, using the package \code{igraph}
#'\item \code{dec.graph}: the path diagram of the decomposed model, with color
#'coding for each vertex
#'\item \code{type.out}: the type of analysis performed, either 'simulation' or 'empirical'}
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
#' causal_decompose(data = macroIT,
#'                eq.system = eq.system,
#'                resid.est = "noDfCor",
#'                instruments = instruments,
#'                sigma.in = NULL)
#'@export
causal_decompose <- function(data,
                            eq.system,
                            resid.est = "noDfCor",
                            instruments,
                            sigma.in = NULL) {

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
  C = gamma.b * R
  Psi = gamma.b * (U - R)
  C.b = as.matrix((C != 0) + 0)
  I = diag(L)

  if(!is.null(sigma.in)){
    if (!matrixcalc::is.square.matrix(sigma.in)) {
    stop("invalid specification for sigma.in. Must be square.")
  }
    sigma = sigma.in

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

    type.out = "simulated"
    fiteq = NULL
    A = NULL

  }else{
    fiteq <-
      systemfit::systemfit(
        data = data,
        formula = eq.system,
        method = "3SLS",
        inst = instruments,
        methodResidCov = resid.est,
        method3sls = "GMM",
        residCovWeighted = F
      )

    sigma = fiteq$residCovEst
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

    type.out = "empirical"
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

  Gamma <-
    as.matrix((Matrix::sparseMatrix(
      i = non.zero.g[, 1], j = non.zero.g[, 2], x = vals
    )))

  C = Gamma * C

  Psi0 = Gamma * Psi0

  Psi1 = Gamma * Psi1

  vals.a = unlist(syst.a)

  non.zero.a = data.frame(which(Ax != 0, arr.ind = T),row.names = NULL)
  non.zero.a[order(non.zero.a[,1],non.zero.a[,2]),]
  non.zero.a = as.matrix(non.zero.a)

  A <-
    as.matrix((Matrix::sparseMatrix(
      i = non.zero.a[, 1], j = non.zero.a[, 2], x = vals.a
    )))
  }

    # ------------------------------------------- #
    #  create decomposed graph:                   #
    #  red = recursive, black = true causal       #
    #  dashed red = induced recursion from Sigma  #
    # ------------------------------------------- #

    dec.graph <- igraph::make_empty_graph(n = ncol(gamma.b))
    dec.graph <- igraph::set.vertex.attribute(dec.graph,"name", value=colnames(gamma.b))
    dec.graph <- igraph::set.vertex.attribute(dec.graph,"color", value="white")
    dec.graph <- igraph::set.vertex.attribute(dec.graph,"shape", value="rectangle")
    dec.graph <- igraph::add_edges(dec.graph,c(apply(t(which(C!=0,arr.ind = T)),2,rev)))
    dec.graph <- igraph::set_edge_attr(dec.graph,"color", value = "red")
    dec.graph <- igraph::set_edge_attr(dec.graph,"lty", value = 1)
    dec.graph <- igraph::set_edge_attr(dec.graph,"arrow.mode", value = ">")
    dec.graph <- igraph::set_edge_attr(dec.graph,"arrow.size", value = 0.7)
    dec.graph <- igraph::set_edge_attr(dec.graph,"curved", value = 0.2)
    dec.graph <- igraph::add_edges(dec.graph,c(apply(t(which((Psi0+Psi1)!=0,arr.ind = T)),2,rev)),
                color = "black",
                arrow.mode=">",
                lty=1,
                curved=0)
    dec.graph <- igraph::add_edges(dec.graph,c(apply(t(which(Psi1!=0,arr.ind = T)),2,rev)),
                color = "red",
                arrow.mode="<>",
                lty=2,
                curved=-0.3)
    igraph::V(dec.graph)$shape <- "rectangle"

    output = list(
      eq.system = eq.system,
      Gamma = Gamma,
      C = C,
      Psi1 = Psi1,
      Psi0 = Psi0,
      A = A,
      Sigma = as.matrix(sigma),
      systemfit = fiteq,
      all.graph = all.graph,
      dec.graph = dec.graph,
      type.out = type.out)

  return(output)
}
