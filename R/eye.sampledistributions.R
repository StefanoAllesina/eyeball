#' Sample pairs of interactions from a bivariate normal distribution
#'
#' @param NumPairs The number of pairs to sample
#' @param mux The desired mean for the first marginal distribution
#' @param muy The desired mean for the second marginal distribution
#' @param sigmax The desired standard deviation for the first marginal distribution
#' @param sigmay The desired standard deviation for the second marginal distribution
#' @param rhoxy The desired correlation
#'
#' @export
#'
#' @return A NumPairs x 2 matrix of interaction strengths sampled from the distribution
#' @examples
#' mypairs <- eye.pairs.from.normal() # default values
#' mypairs <- eye.pairs.from.normal(NumPairs = 120)
eye.pairs.from.normal <- function(NumPairs = 10,
                              mux = -1,
                              muy = 0.5,
                              sigmax = 1/4,
                              sigmay = 1/4,
                              rhoxy = -2/3){
  mus <- c(mux, muy)
  covariance.matrix <- matrix(c(sigmax^2,
                                rhoxy * sigmax * sigmay,
                                rhoxy * sigmax * sigmay,
                                sigmay^2),
                              2, 2)
  Pairs <- mvrnorm(NumPairs, mus, covariance.matrix)
  return(Pairs)
}

#' Sample pairs of interactions from the "Four-Corner" distribution described
#' in Allesina et al. 2014
#'
#' @param NumPairs The number of pairs to sample
#' @param mux The desired mean for the first marginal distribution
#' @param muy The desired mean for the second marginal distribution
#' @param sigmax The desired standard deviation for the first marginal distribution
#' @param sigmay The desired standard deviation for the second marginal distribution
#' @param rhoxy The desired correlation
#'
#' @export
#'
#' @return A NumPairs x 2 matrix of interaction strengths sampled from the distribution
#' @examples
#' mypairs <- eye.pairs.from.fourcorner() # default values
#' mypairs <- eye.pairs.from.fourcorner(NumPairs = 120)
eye.pairs.from.fourcorner <- function(NumPairs = 10,
                                  mux = -1,
                                  muy = 0.5,
                                  sigmax = 1/4,
                                  sigmay = 1/4,
                                  rhoxy = -2/3){
  gamma <- (rhoxy + 1) / 2
  PosNeg <- sign(rnorm(NumPairs))
  SwitchSign <- sign(runif(NumPairs, -(1-gamma), gamma))
  Pairs <- cbind(mux + PosNeg * sigmax, muy + PosNeg * SwitchSign * sigmay)
  return(Pairs)
}

#' Sample pairs of interactions from an empirical distribution
#'
#' @param NumPairs The number of pairs to sample
#' @param Empirical.Distribution A Nx2 matrix with the effects of consumers on resources
#'                                in the first column and those of resources on consumers in the second.
#'
#' @export
#'
#' @return A NumPairs x 2 matrix of interaction strengths sampled from the distribution
#' @examples
#' mypairs <- eye.pairs.from.empirical(100, matrix(runif(20), 10, 2))
eye.pairs.from.empirical <- function(NumPairs = 10,
                                 Empirical.Distribution = NULL){
  if (is.matrix(Empirical.Distribution) == TRUE){
    NR <- dim(Empirical.Distribution)[2]
    if (NR == 2){
      Pairs <- matrix(NumPairs, 2)
      SampledPairs <- sample(1:NR, NumPairs, replace = TRUE)
      Pairs <- Empirical.Distribution[SampledPairs,]
      return(Pairs)
    }
  }
  stop("The Empirical.Distribution must be a N x 2 matrix!")
}
