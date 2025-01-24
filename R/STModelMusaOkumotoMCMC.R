#' Title
#'
#' @param data aaaa
#' @param sites aaaa
#' @param prior aaaa
#' @param iteration aaaa
#' @param burnin aaaaaa
#'
#' @return aaaaa
#' @export
#'
#'
STModelMusaOkumotoMCMC <- function(data, sites, prior = list(
                                     aa1 = 0.001, bb1 = 0.001,
                                     c2 = 1e-05,
                                     d2 = 0.001,
                                     c4 = 2.01,
                                     d4 = 1.005, aa2 = 2.01,
                                     bb2 = 1.005,
                                     c3 = (-2 * log(0.05) / max(dist(sites))),
                                     d3 = 0.1
                                   ), iteration, burnin) {
  b <- 1
  v <- 1
  alpha <- 1
  W <- as.matrix(rep(0, ncol(data)))
  # Hiperparametros

  aa1 <- prior$aa1
  bb1 <- prior$bb1

  c2 <- prior$c2
  d2 <- prior$d2

  c4 <- prior$c4
  d4 <- prior$d4

  aa2 <- prior$aa2
  bb2 <- prior$bb2


  c3 <- prior$c3
  d3 <- prior$d3

  SU2 = 1000
  SU3 = 28.01968
  SU5 = 0.001118574

  X <- cbind(as.matrix(rep(1, ncol(data))), as.matrix(sites))
  Psi <- as.matrix(rep(0, ncol(X)))
  A <- as.matrix(rep(0, ncol(X)))
  B <- diag(100, ncol(X))

  V <- diag(100, ncol(X))
  M <- as.matrix(rep(0, ncol(X)))

  n <- ncol(data)
  m <- nrow(data)
  tempdados <- is.na(data)
  nj <- m - apply(tempdados, 2, sum)
  # p=ncol(X)



  Tt <- array(NA, dim = c(1, n))
  for (y in 1:n) {
    Tt[1, y] <- data[nj[y], y]
  }

  # pp=ncol(X)
  Tt <- t(Tt)

  Malpha <- NULL
  MalphaT <- NULL

  MW <- NULL
  MWT <- NULL

  MPsi <- NULL

  Mv <- NULL
  Mb <- NULL
  MbT <- NULL


  #############################
  ## Programa principal
  ############################ 3
  for (j in 1:iteration) {
    if (j <= burnin) {
      SIGMA <- gSigma(b, v, sites)
      DELTA <- gCorr(b, sites)

      temp <- amostrarWMUSA(W, sites, X, Psi, b, v, nj, SU5, Tt, alpha)
      W <- as.matrix(temp[[1]])
      MWT <- c(MWT, temp[[2]])

      temp <- amostraralphaMUSA(alpha, W, Tt, c2, d2, data, SU2)
      alpha <- temp[[1]]
      MalphaT <- c(MalphaT, temp[[2]])

      temp <- amostrarb(W, v, b, sites, c3, d3, X, Psi, SU3)
      b <- temp[[1]]
      MbT <- c(MbT, temp[[2]])


      aa1 <- 0.5 * ncol(data) + c4
      bb1 <- 0.5 * t(W - X %*% Psi) %*% solve(DELTA) %*% (W - X %*% Psi) + d4
      v <- 1 / rgamma(1, aa1, bb1)


      AA <- solve(solve(V) + t(X) %*% solve(SIGMA) %*% X) %*% (solve(V) %*% M + t(X) %*% solve(SIGMA) %*% W)
      BB <- solve(solve(V) + t(X) %*% solve(SIGMA) %*% X)

      Psi <- as.matrix(MASS::mvrnorm(1, AA, BB))


      SU2 <- sintonizarMUSA(burnin, 0.30, SU2, MalphaT, j)
      SU3 <- sintonizarMUSA(burnin, 0.44, SU3, MbT, j)
      SU5 <- sintonizarNMUSA(burnin, 0.25, SU5, MWT, j)
      print(j)
    } else {
      SIGMA <- gSigma(b, v, sites)
      DELTA <- gCorr(b, sites)

      temp <- amostrarWMUSA(W, sites, X, Psi, b, v, nj, SU5, Tt, alpha)
      W <- as.matrix(temp[[1]])
      MW <- rbind(MW, t(W))
      MWT <- c(MWT, temp[[2]])

      temp <- amostraralphaMUSA(alpha, W, Tt, c2, d2, data, SU2)
      alpha <- temp[[1]]
      Malpha <- c(Malpha, alpha)
      MalphaT <- c(MalphaT, temp[[2]])

      temp <- amostrarb(W, v, b, sites, c3, d3, X, Psi, SU3)
      b <- temp[[1]]
      Mb <- c(Mb, b)
      MbT <- c(MbT, temp[[2]])

      aa1 <- 0.5 * ncol(data) + c4
      bb1 <- 0.5 * t(W - X %*% Psi) %*% solve(DELTA) %*% (W - X %*% Psi) + d4
      v <- 1 / rgamma(1, aa1, bb1)
      Mv <- c(Mv, v)


      AA <- solve(solve(V) + t(X) %*% solve(SIGMA) %*% X) %*% (solve(V) %*% M + t(X) %*% solve(SIGMA) %*% W)
      BB <- solve(solve(V) + t(X) %*% solve(SIGMA) %*% X)

      Psi <- MASS::mvrnorm(1, AA, BB)
      MPsi <- rbind(MPsi, t(Psi))




      print(j)
    }
  }
  resul <- list(SIGMA, DELTA, MW, MWT, Malpha, MalphaT, Mb, MbT, Mv, MPsi)
  names(resul) <- c("SIGMA", "DELTA", "MW", "MWT", "Malpha", "MalphaT", "Mb", "MbT", "Mv", "MPsi")
  return(resul)
}
