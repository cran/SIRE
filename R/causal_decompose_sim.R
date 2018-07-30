# ----------------------------------------------------------------- #
#   Decomposition (No Estimation) of a Simultaneous Equation Model  #
#         into its Recursive and Interdependent Sub-Systems         #
#                    (Gamma = C + Psi1 + Psi0)                      #
# ----------------------------------------------------------------- #
#' Decomposition of simultaneous equation model
#'
#' @description Decompose
#' a Simultaneous Equation Model into its recursive
#' and Interdependent sub-systems
#' @param eq.system the system of equations (a list of formula objects)
#' @param Sigma the \eqn{\Sigma} matrix
#'
#'@return A list with components \itemize{
#'\item \code{output}: a list containing\itemize{
#'\item \code{eq.system}: the system of equations given as input
#'\item \code{Gamma}: the binary matrix \eqn{\Gamma^{b'}}
#'\item \code{C}: the binary matrix highlighting the interdependent mechanisms at deterministic level.
#'\item \code{Psi1}: the binary matrix highlighting the interdependent mechanisms at stochastic level.
#'\item \code{Psi0}: the binary matrix highlighting the causal mechanisms.
#'\item \code{Sigma}: the matrix \eqn{\Sigma} given as input
#'\item \code{all.graph}: the path diagram of the model, using the package \code{igraph}
#'\item \code{dec.graph}: the path diagram of the decomposed model, with color
#'coding for each vertex. Blue means causality, red means interdependence,
#'dashed red means interdependence by effect of error correlation}
#'\item \code{path}: a data frame in which every row is a path coefficient,
#'along with indication of the nature of the link (recursive/interdependent)
#'and the index of the equation in which it belongs
#'}
#'
#'@examples
#'eq.system = list(
#'             eq1 = y1 ~ y5 + y7,
#'             eq2 = y2 ~ z,
#'             eq3 = y3 ~ y11,
#'             eq4 = y4 ~ y3,
#'             eq5 = y5 ~ y10,
#'             eq6 = y6 ~ y5 + y9,
#'             eq7 = y7 ~ y6,
#'             eq8 = y8 ~ y12,
#'             eq9 = y9 ~ y7,
#'             eq10 = y10 ~ y5,
#'             eq11 = y11 ~ y12,
#'             eq12 = y12 ~ y4 + y11,
#'             eq13 = y13 ~ y2 + y6)
#'
#'  # indexes of non-null elements of Sigma
#'  sigma.idx = cbind(c(2,1),
#'                    c(1,5),
#'                    c(13,2),
#'                    c(2,13),
#'                    c(5,1),
#'                    c(1,2))
#'
#'  #fictitious Sigma matrix
#'  Sigma = as.matrix(
#'          Matrix::sparseMatrix(i = sigma.idx[1,] , j = sigma.idx[2,], x = 0.1)) +
#'          diag(length(eq.system))
#'
#'  causal.decompose.sim(eq.system , Sigma)
#' @export

causal.decompose.sim = function(eq.system , Sigma = NULL) {

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

  sorted.rhs =
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

  sorted.eqn =
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

  # --------------- #
  #  overall graph  #
  # --------------- #

  all.graph = igraph::graph.adjacency(t(gamma.b))

  # ----------------------- #
  #  decomposition C + Psi  #
  # ----------------------- #

  R = matrix(0, nrow = L, ncol = L)
  powers = list()
  LL = L - 1
  for (p in 1:LL) {
    powers[[p]] = matrixcalc::matrix.power(gamma.b, p)
  }

  R = as.matrix((t(Reduce('+', powers)) != 0) + 0)
  U = matrix(1, nrow = L , ncol = L)
  C = matrixcalc::hadamard.prod(gamma.b, R)
  Psi = matrixcalc::hadamard.prod(gamma.b, ((U - R)))
  C.b = as.matrix((C != 0) + 0)
  I = diag(L)

  # --------------------------- #
  #  decomposition Psi0 + Psi1  #
  # --------------------------- #

  if(is.null(Sigma)) Sigma = diag(length(eq.system))

  if (!matrixcalc::is.symmetric.matrix(Sigma) || !matrixcalc::is.positive.definite(Sigma)) {
    stop("invalid specification for sigma. Must me symmetric and positive definite.")
  }

  if (!Matrix::isDiagonal(Sigma)) {
    sigma.b = as.matrix((Sigma != 0) + 0)
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

    C.graph = igraph::graph.adjacency(t(C))
    Psi1.graph = igraph::graph.adjacency(t(Psi1))
    Psi0.graph = igraph::graph.adjacency(t(Psi0))
    CPsi.graph = igraph::graph.adjacency(t(C + Psi1))

    # ------------------------------------------- #
    #  create decomposed graph:                   #
    #  red = recursive, green = true causal       #
    #  dashed red = induced recursion from Sigma  #
    # ------------------------------------------- #

    igraph::E(all.graph)$color = "black"
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

    output = list(
      eq.system = eq.system,
      Gamma = gamma.b,
      C = C,
      Psi1 = Psi1,
      Psi0 = Psi0,
      Sigma = Sigma,
      dec.graph = dec.graph,
      all.graph = all.graph
    )
  } else{
    # ------------------- #
    #  create sub-graphs  #
    # ------------------- #

    rownames(C) = c(y.vars)
    colnames(C) = c(y.vars)
    rownames(Psi) = c(y.vars)
    colnames(Psi) = c(y.vars)

    C.graph = igraph::graph.adjacency(t(C))
    Psi.graph = igraph::graph.adjacency(t(Psi))

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

    igraph::E(all.graph)$color = "black"

    # ------------------------ #
    #  complete linked pairs,  #
    #  '~>' means recursive,   #
    #  '->' means true causal  #
    # ------------------------ #

    path = sort(c(
      stringr::str_replace(igraph::as_ids(igraph::E(C.graph)), "\\|", " ~> "),
      stringr::str_replace(igraph::as_ids(igraph::E(Psi.graph)), "\\|", " -> ")
    ))

    output = list(
      eq.system = eq.system,
      Gamma = gamma.b,
      C = C,
      Psi = Psi,
      dec.graph = dec.graph,
      all.graph = all.graph
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

  return(list(output = output, path = data.frame(paths)))
}
