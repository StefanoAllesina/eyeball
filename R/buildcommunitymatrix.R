eye.parameterize.M <- function(FW,
                               distribution.pairs = "Normal",
                               mux = -1,
                               muy = 0.5,
                               sigmax = 1/4,
                               sigmay = 1/4,
                               rhoxy = -2/3,
                               mu.diagonal = 0,
                               sigma.diagonal = 0){
  S <- FW$S
  M <- matrix(0, S, S)
  ## Get off-diagonal coefficients
  NumPairs <- FW$L
  if (is.matrix(distribution.pairs) == TRUE){
    my.pairs <- eye.pairs.from.empirical(NumPairs = NumPairs,
                                         Empirical.Distribution = distribution.pairs)
  } else {
    if (distribution.pairs == "Normal"){
      my.pairs <- eye.pairs.from.normal(NumPairs = NumPairs,
                                        mux = mux,
                                        muy = muy,
                                        sigmax = sigmax,
                                        sigmay = sigmay,
                                        rhoxy = rhoxy)
    }
    if (distribution.pairs == "Normal"){
      my.pairs <- eye.pairs.from.normal(NumPairs = NumPairs,
                                        mux = mux,
                                        muy = muy,
                                        sigmax = sigmax,
                                        sigmay = sigmay,
                                        rhoxy = rhoxy)
    }
  }
  M[FW$links] <- my.pairs[ , 1]
  M[FW$links[,2:1]] <- my.pairs[ , 2]
  ## Optional: set diagonal
  diag(M) <- rnorm(S, mean = mu.diagonal, sd = sigma.diagonal)
  return(M)
}

eye.buildfoodweb.and.parameterize.M <- function(foodweb.model,
                                                S = 100,
                                                C = 0.1,
                                                distribution.pairs = "Normal",
                                                mux = -1,
                                                muy = 0.5,
                                                sigmax = 1/4,
                                                sigmay = 1/4,
                                                rhoxy = -2/3,
                                                mu.diagonal = 0,
                                                sigma.diagonal = 0){
  FW <- NULL
  # if it's a file
  if (file.exists(foodweb.model) == TRUE){
    FW <- eye.foodweb.file(foodweb.model)
  }
  if (foodweb.model == "Cascade"){
    FW <- eye.foodweb.cascade(S, C)
  }
  if (foodweb.model == "Niche"){
    FW <- eye.foodweb.niche(S, C)
  }
  if (is.null(FW)){
    stop("Invalid food web model. Please enter Cascade, Niche, or a file name for the adjacency matrix.")
  }
  M <- eye.parameterize.M(FW, distribution.pairs, mux, muy, sigmax, sigmay, rhoxy, mu.diagonal, sigma.diagonal)
  return(M)
}
