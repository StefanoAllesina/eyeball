#' Builds a food web according to the cascade model
#'
#' The cascade model is described in Cohen et al. "Community food webs: Data and theory", 1990.
#' The species are sorted in order from 1...S and each species has probability C of consuming any of the
#' preceding species.
#' @param S The number of species (integer)
#' @param C The desired connectance (between 0 and 1)
#' @return A connected (one piece) food web of connectance close to the desired value.
#' The food web is a list containing:
#' \describe{
#' \item{$links}{The trophic links from resource to consumer (L x 2 matrix)}
#' \item{$S}{The number of species}
#' \item{$L}{The number of links}
#' \item{$C}{The connectance L / (S choose 2)}
#' \item{$Model}{A code for the model used to generate the food web}
#' \item{$Cycles}{Boolean, whether cycles are present}
#' }
#' All the self-loops and double arrows (a eats b, b eats a) are removed. If the food web contains
#' no cycles, then it is sorted that all the coefficients would be in the upper-triangular part
#' of the corresponding adjacency matrix.
#'
#' If the resulting network has connectance too different from the desired one, or it is not connected,
#' the function will call itself recursively. If too many recursions are performed, it will fail.
#' @examples
#' FW <- eye.foodweb.cascade() # S = 100, C = 0.25
#' FW <- eye.foodweb.cascade(30, 0.2)
#' FW <- eye.foodweb.cascade(S = 250, C = 0.314)
eye.foodweb.cascade <- function(S = 100, C = 0.25){
  ## Build the adjacency matrix
  K <- (matrix(runif(S * S), S, S) < C) * 1
  K[lower.tri(K)] <- 0
  diag(K) <- 0
  ## Now check that it passes all the following criteria
  ## a) the number of connections is about right
  L <- sum(K)
  expected.L <- choose(S, 2) * C
  variance.L <- choose(S, 2) * C * (1. - C)
  if (L < (expected.L - 2 * sqrt(variance.L)) |
      L > (expected.L + 2 * sqrt(variance.L))){
        warning("Failed to build foodweb, wrong number of links. Trying again...")
        return(eye.foodweb.cascade(S, C))
      }
  ## b) The network is connected
  g <- graph.adjacency(K, mode = "directed")
  if (is.connected(g, mode = "weak") == FALSE){
    warning("Failed to build foodweb, not connected. Trying again...")
    return(eye.foodweb.cascade(S, C))
  }
  ## If it passes the tests, then build the foodweb
  fw <- list()
  fw$links <- get.edgelist(g)
  fw$S <- S
  fw$L <- dim(fw$links)[1]
  fw$C <- fw$L / (S * (S - 1) * 0.5)
  fw$Model <- "Cascade"
  fw$Cycles <- FALSE
  return(fw)
}

#' Builds a food web according to the niche model
#'
#' The niche model is described in Williams and Martinez, Nature 2000.
#' The species are sorted in order from 1...S and each species consumes a set of
#' adjacent species.
#' @param S The number of species (integer)
#' @param C The desired connectance (between 0 and 1)
#' @return A connected (one piece) food web of connectance close to the desired value.
#' The food web is a list containing:
#' \describe{
#' \item{$links}{The trophic links from resource to consumer (L x 2 matrix)}
#' \item{$S}{The number of species}
#' \item{$L}{The number of links}
#' \item{$C}{The connectance L / (S choose 2)}
#' \item{$Model}{A code for the model used to generate the food web}
#' \item{$Cycles}{Boolean, whether cycles are present}
#' }
#' All the self-loops and double arrows (a eats b, b eats a) are removed. If the food web contains
#' no cycles, then it is sorted that all the coefficients would be in the upper-triangular part
#' of the corresponding adjacency matrix.
#'
#' If the resulting network has connectance too different from the desired one, or it is not connected,
#' the function will call itself recursively. If too many recursions are performed, it will fail.
#' @examples
#' FW <- eye.foodweb.niche() # S = 100, C = 0.25
#' FW <- eye.foodweb.niche(30, 0.2)
#' FW <- eye.foodweb.niche(S = 250, C = 0.314)
eye.foodweb.niche <- function(S = 100, C = 0.25){
  ## Build the adjacency matrix
  K <- matrix(0, S, S)
  ## assign the niche value for all species
  ## Note: the niche values of the species are ordered
  ni <- sort(runif(S, 0, 1))
  ## determine the radius of the food spectrum for each species
  ## the radius is drawn randomly from a ni * Beta(1, beta) distribution
  beta <- 1.0 / C - 1
  ri <- rbeta(S, 1, beta) * ni
  ## Set species with the lowest niche value to be the basal species
  ri[1] <- 0
  ## The center of the food spectrum for each species are
  ## drawn uniformly from an uniform distribution on [ri/2, ni]
  ci <- numeric(S)
  for(i in 1:S){
    ci[i] <- runif(1, ri[i] / 2, min(ni[i], 1 - ri[i] / 2))
    ## Correction as in Allesina et al Science 2008
    ## Determine the boundary for the food spectrum for each species:
    upper <- ci[i] + ri[i] / 2
    lower <- ci[i] - ri[i] / 2
    ## check which species are falling into that interval: [lower, upper]
    ## and set the corresponding value in the ajacency matrix M to be 1,
    ## indicating a link is established between the two species.
    for(j in 1:S){
      if(ni[j] > lower & ni[j] < upper)
        K[j, i] <- 1
        ## remove double edges
        if (K[i, j] == 1) {
          ## remove one at random
          if (runif(1) < 0.5){
            K[i, j] <- 0
          } else {
            K[j, i] <- 0
          }
        }
    }
  }
  diag(K) <- 0
  ## Now check that it passes all the following criteria
  ## a) the number of connections is about right
  L <- sum(K)
  expected.L <- choose(S, 2) * C
  variance.L <- choose(S, 2) * C * (1. - C)
  if (L < (expected.L - 2 * sqrt(variance.L)) |
        L > (expected.L + 2 * sqrt(variance.L))){
    warning("Failed to build foodweb, wrong number of links. Trying again...")
    return(eye.foodweb.niche(S, C))
  }
  ## b) The network is connected
  g <- graph.adjacency(K, mode = "directed")
  if (is.connected(g, mode = "weak") == FALSE){
    warning("Failed to build foodweb, not connected. Trying again...")
    return(eye.foodweb.niche(S, C))
  }
  ## c) Check if it has no cycles. If so, perform a topological sort
  ts <- topological.sort(g)
  if (length(ts) == S){
    K <- K[ts, ts]
    g <- graph.adjacency(K, mode = "directed")
    Cycles <- FALSE
  } else {
    Cycles <- TRUE
  }
  ## If it passes the tests, then build the foodweb
  fw <- list()
  fw$links <- get.edgelist(g)
  fw$S <- S
  fw$L <- dim(fw$links)[1]
  fw$C <- fw$L / (S * (S - 1) * 0.5)
  fw$Model <- "Niche"
  fw$Cycles <- Cycles
  return(fw)
}

#' Build a food web reading a file
#'
#' Read the binary adjacency matrix contained in the file.
#' Remove self-loops and double edges
#' @param filename A text file containing the adjacency matrix (binary, rows = resources, cols = consumers) of a food web
#' The food web is a list containing:
#' \describe{
#' \item{$links}{The trophic links from resource to consumer (L x 2 matrix)}
#' \item{$S}{The number of species}
#' \item{$L}{The number of links}
#' \item{$C}{The connectance L / (S choose 2)}
#' \item{$Model}{A code for the model used to generate the food web}
#' \item{$Cycles}{Boolean, whether cycles are present}
#' }
#' All the self-loops and double arrows (a eats b, b eats a) are removed. If the food web contains
#' no cycles, then it is sorted that all the coefficients would be in the upper-triangular part
#' of the corresponding adjacency matrix.
#' @examples
#' FW <- eye.foodweb.file("myadjacencymat.txt")
eye.foodweb.file <- function(filename){
  ## Build the adjacency matrix
  K <- (as.matrix(read.table(filename, header = FALSE)) > 0) * 1
  colnames(K) <- NULL
  rownames(K) <- NULL
  S <- dim(K)[1]
  diag(K) <- 0
  ## Check double arrows
  KK <- K + t(K)
  KK[upper.tri(KK)] <- 0
  KK <- (KK == 2) * 1
  K <- K - KK
  g <- graph.adjacency(K, mode = "directed")
  ## Check if it has no cycles. If so, perform a topological sort
  ts <- topological.sort(g)
  if (length(ts) == S){
    K <- K[ts, ts]
    g <- graph.adjacency(K, mode = "directed")
    Cycles <- FALSE
  } else {
    Cycles <- TRUE
  }
  ## Build the foodweb
  fw <- list()
  fw$links <- get.edgelist(g)
  fw$S <- S
  fw$L <- dim(fw$links)[1]
  fw$C <- fw$L / (S * (S - 1) * 0.5)
  fw$Model <- "From File"
  fw$Cycles <- Cycles
  return(fw)
}
