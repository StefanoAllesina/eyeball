#' @export
eye.approximate.ReL1 <- function(M, calculate.eigenvalues = TRUE){
  if (calculate.eigenvalues == TRUE){
    ev <- eigen(M, only.values = TRUE, symmetric = FALSE)$values
  }
  S <- dim(M)[1]
  NOffDiag <- S * (S - 1)
  ## First, deal with diagonal elements
  d <- mean(diag(M))
  diag(M) <- 0
  ## Second, compute stats
  ## For May and Tang et al.
  mu <- sum(M) / NOffDiag
  sigma2 <- sum(M^2) / NOffDiag - mu^2
  sigma <- sqrt(sigma2)
  rho <- (sum(M * t(M)) / NOffDiag - mu^2) / sigma2
  ## For eyeball
  muU <- mean(M[upper.tri(M)])
  sigmaU2 <- sum(M[upper.tri(M)]^2) * 2 / NOffDiag - muU^2
  muL <- mean(M[lower.tri(M)])
  sigmaL2 <- sum(M[lower.tri(M)]^2) * 2 / NOffDiag - muL^2
  sigmaU <- sqrt(sigmaU2)
  sigmaL <- sqrt(sigmaL2)
  rhoUL <- (sum(M * t(M)) / NOffDiag - muU * muL) / (sigmaU * sigmaL)
  ## May's stability criterion
  ReL1.May <- max((S-1) * mu, sqrt((S-1) * sigma2) - mu) + d
  ## Tang et al. stability criterion
  ReL1.TangEtAl <- max((S-1) * mu, sqrt((S-1) * sigma2) * (1 + rho) - mu) + d
  ## Eyeball approximation
  ## radius for the spectrum of A
  coeff <- (-muL / muU)^(1 / S)
  coeff2 <- (-muL / muU)^(2 / S)
  r.a <- (muU - muL) * coeff / (coeff2 - 1)
  ## center for the spectrum of A
  c.a <- (muL - muU * coeff2) / (coeff2 -1)
  Re.eigenA.middle <- r.a + c.a
  x <- -muL / muU
  thetak <- pi / S
  Re.eigenA.extreme <- (1 + (x^(1/S) * cos(thetak) -1) * (1 + x) / (1 + x^(2/ S) - 2 * x^(1/S) * cos(thetak))) * -muU
  ## Now spectrum of B
  ## Approximation of alpha
  alpha <- S * (sigmaU^2 - sigmaL) / log(sigmaL^2 / sigmaU^2)
  r.b <- (alpha + rhoUL * sigmaU * sigmaL * (S-1)) / sqrt(alpha)
  s.b <- (alpha - rhoUL * sigmaU * sigmaL * (S-1)) / sqrt(alpha)
  ReL1.eyeball <- max(Re.eigenA.extreme, r.b + Re.eigenA.middle) + d
  ReL1.observed <- NULL
  M.eigenvalues <- NULL
  ## if desired, calculate the actual ReL1.observed
  if (calculate.eigenvalues == TRUE){
    ReL1.observed <- max(Re(ev))
    M.eigenvalues <- ev
    ## Debug
    #print(plot(ev, xlim = c(min(Re(ev)), 1 + max(c(ReL1.observed, ReL1.May, ReL1.TangEtAl)))))
    #print(abline(v = ReL1.observed, col = "black"))
    #print(abline(v = ReL1.May, col = "red"))
    #print(abline(v = ReL1.TangEtAl, col = "blue"))
    #print(abline(v = ReL1.eyeball, col = "pink"))
    ## End Debug
  }
  return(list(M = M,
              M.eigenvalues = M.eigenvalues,
              ReL1.observed = ReL1.observed,
              ReL1.May = ReL1.May,
              ReL1.TangEtAl = ReL1.TangEtAl,
              ReL1.eyeball = ReL1.eyeball,
              May.stats = list(S = S, V = sigma2, E = mu, d = d),
              TangEtAl.stats = list(S = S, V = sigma2, E = mu, rho = rho, d = d),
              eyeball.stats = list(S = S, muU = muU, muL = muL,
                                   sigmaU = sigmaU, sigmaL = sigmaL,
                                   rhoUL = rhoUL,
                                   d = d,
                                   radius.A = r.a,
                                   center.A = c.a,
                                   radius.B.h = r.b,
                                   radius.B.v = s.b)
              ))
}
