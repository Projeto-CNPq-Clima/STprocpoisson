######################################################
gCorr <- function(b, def) {
  n <- nrow(def)
  R <- exp(-b * (as.matrix(dist(def))))
  mat <- R
  mat
}
######################################################
gSigma <- function(b, v, def) {
  n <- nrow(def)
  R <- exp(-b * (as.matrix(dist(def))))
  mat <- v * R
  mat
}
######################################################
logvero <- function(aalpha, WW, TT, nn, yT) {
  res <- sum(nn) * log(aalpha) + sum(WW * nn) - sum(exp(WW) * TT^(aalpha)) + aalpha * sum(log(yT), na.rm = T)
  res
}


######################################################
sintonizar <- function(bar, taxa, tau, mat, i) {
  mat <- as.matrix(mat)



  mater <- (1 / 50) * sum(mat[(i - 49):i, 1])

  if (mater >= taxa) {
    delta <- min(0.01, (i / 50 + 1)^(-0.5))
    temp4 <- log(tau) - delta
    temp5 <- exp(temp4)
    return(temp5)
  } else {
    delta <- min(0.01, (i / 50 + 1)^(-0.5))
    temp4 <- log(tau) + delta
    temp5 <- exp(temp4)
    return(temp5)
  }
}
######################################################
sintonizarN <- function(bar, taxa, tau, mat, i) {
  mat <- as.matrix(mat)



  mater <- (1 / 50) * sum(mat[(i - 49):i, 1])

  if (mater >= taxa) {
    delta <- min(0.01, (i / 50 + 1)^(-0.5))
    temp4 <- log(tau) + delta
    temp5 <- exp(temp4)
    return(temp5)
  } else {
    delta <- min(0.01, (i / 50 + 1)^(-0.5))
    temp4 <- log(tau) - delta
    temp5 <- exp(temp4)
    return(temp5)
  }
}
######################################################
amostrarW <- function(WW, MM, S, XX, PPs, bb, vv, nn, TT, ff) {
  n <- nrow(WW)
  WWprop <- as.matrix(MASS::mvrnorm(1, WW, ff * diag(1, n)))


  SSig <- gSigma(bb, vv, S)

  postWW <- sum(WW * nn) - sum(exp(WW) * t(TT)^exp(MM)) - 0.5 * t(WW - XX %*% PPs) %*% solve(SSig) %*% (WW - XX %*% PPs)
  postWWprop <- sum(WWprop * nn) - sum(exp(WWprop) * t(TT)^exp(MM)) - 0.5 * t(WWprop - XX %*% PPs) %*% solve(SSig) %*% (WWprop - XX %*% PPs)

  prob <- min(exp((postWWprop) - (postWW)), 1)


  u <- runif(1, 0, 1)

  if (u < prob) {
    Wprox <- WWprop

    rejei <- 1
  } else {
    Wprox <- WW
    rejei <- 0
  }





  res <- as.matrix(Wprox)
  res <- list(Wprox, rejei)
  res
}
######################################################
# amostrarW(W,M,sites,X,Psi,b,v,nj,Tt,SU2)
amostrarM <- function(WW, MM, S, ZZ, BBta, bb, vv, yT, nn, TT, ff) {
  n <- nrow(MM)
  MMprop <- as.matrix(MASS::mvrnorm(1, MM, ff * diag(1, n)))
  logyT <- log(yT)
  SSig <- gSigma(bb, vv, S)

  sum1 <- 0
  sum2 <- 0

  for (j in 1:ncol(yT)) {
    res1 <- sum(exp(MM[j, ]) * logyT[, j], na.rm = T)
    res2 <- sum(exp(MMprop[j, ]) * logyT[, j], na.rm = T)
    sum1 <- sum1 + res1
    sum2 <- sum2 + res2
  }

  postMM <- sum(MM * nn) - sum(exp(WW) * t(TT)^exp(MM)) + sum1 - 0.5 * t(MM - ZZ %*% BBta) %*% solve(SSig) %*% (MM - ZZ %*% BBta)
  postMMprop <- sum(MMprop * nn) - sum(exp(WW) * t(TT)^exp(MMprop)) + sum2 - 0.5 * t(MMprop - ZZ %*% BBta) %*% solve(SSig) %*% (MMprop - ZZ %*% BBta)

  prob <- min(exp((postMMprop) - (postMM)), 1)


  u <- runif(1, 0, 1)

  if (u < prob) {
    MMprox <- MMprop

    rejei <- 1
  } else {
    MMprox <- MM
    rejei <- 0
  }





  res <- as.matrix(MMprox)
  res <- list(MMprox, rejei)
  res
}
# amostrarM(W,M,sites,Z,Beta,b,v,yT,nj,Tt,SU2)
######################################################
amostrarb <- function(W, v, b, loca, ab, bb, X, Psi, u1) {
  bprop <- rgamma(1, shape = b * u1, rate = u1)

  SSigprop <- gSigma(bprop, v, loca)

  if ((det(SSigprop) == 0) | (bprop < 0.005)) {
    return(list(b, 0))
  }

  SSig <- gSigma(b, v, loca)
  SSigprop <- gSigma(bprop, v, loca)


  logp <- -0.5 * t(W - X %*% Psi) %*% solve(SSig) %*% (W - X %*% Psi) - 0.5 * log(det(SSig)) + (ab - 1) * log(b) - bb * b

  logpprop <- -0.5 * t(W - X %*% Psi) %*% solve(SSigprop) %*% (W - X %*% Psi) - 0.5 * log(det(SSigprop)) + (ab - 1) * log(bprop) - bb * bprop

  logprob <- logpprop + log(dgamma(b, shape = bprop * u1, rate = u1)) - (logp + log(dgamma(bprop, shape = b * u1, rate = u1)))
  prob <- min(c(1, exp(logprob)))

  u <- runif(1, 0, 1)

  if (u < prob) {
    bprox <- bprop

    rejei <- 1
  } else {
    bprox <- b

    rejei <- 0
  }



  res <- list(bprox, rejei)
  res
}

######################################################
logverosa <- function(MM, WW, ddelta, ttheta, yT, TT, f) {
  suma <- 0
  for (i in 1:ncol(yT)) {
    res <- sum(log(exp(WW[i, ] + MM[i, ]) * as.matrix(yT[, i])^(exp(MM[i, ]) - 1) - ddelta * 2 * pi * f * sin(2 * pi * f * yT[, i] + ttheta)), na.rm = T)
    suma <- suma + res
  }

  res1 <- suma - sum(exp(WW) * t(TT)^(exp(MM)) + ddelta * cos(2 * pi * f * t(TT) + ttheta))
  res1
}
######################################################
amostrarWsa <- function(ddelta, ttheta, WW, MM, loca, XX, PPs, bb, vv, nn, TT, yT, ff, f) {
  n <- nrow(WW)
  WWprop <- as.matrix(MASS::mvrnorm(1, WW, ff * diag(1, n)))

  tema1 <- ifelse((exp(WWprop + MM) * t(TT)^(exp(MM) - 1)) <= (2 * pi * f * ddelta), 1, 0)
  tema2 <- ifelse((exp(WWprop + MM) * (yT[1, ])^(exp(MM) - 1)) <= (2 * pi * f * ddelta), 1, 0)



  if ((sum(tema1) + sum(tema2)) >= 1) {
    return(list(WW, 0))
  } else {

  }


  SSig <- gSigma(bb, vv, loca)

  postWW <- logverosa(MM, WW, ddelta, ttheta, yT, TT, f) - 0.5 * t(WW - XX %*% PPs) %*% solve(SSig) %*% (WW - XX %*% PPs)
  postWWprop <- logverosa(MM, WWprop, ddelta, ttheta, yT, TT, f) - 0.5 * t(WWprop - XX %*% PPs) %*% solve(SSig) %*% (WWprop - XX %*% PPs)

  prob <- min(exp((postWWprop) - (postWW)), 1)


  u <- runif(1, 0, 1)

  if (u < prob) {
    Wprox <- WWprop

    rejei <- 1
  } else {
    Wprox <- WW
    rejei <- 0
  }





  res <- as.matrix(Wprox)
  res <- list(Wprox, rejei)
  res
}
# amostrarWsa(delta,theta,W,M,sites,X,Psi,bw,vw,nj,Tt,yT,SU1,f)
######################################################
amostrarMsa <- function(ddelta, ttheta, WW, MM, loca, XX, PPs, bb, vv, nn, TT, yT, ff, f) {
  n <- nrow(MM)
  MMprop <- as.matrix(MASS::mvrnorm(1, MM, ff * diag(1, n)))

  tema1 <- ifelse((exp(WW + MMprop) * t(TT)^(exp(MMprop) - 1)) <= (2 * pi * f * ddelta), 1, 0)
  tema2 <- ifelse((exp(WW + MMprop) * (yT[1, ])^(exp(MMprop) - 1)) <= (2 * pi * f * ddelta), 1, 0)


  if ((sum(tema1) + sum(tema2)) >= 1) {
    return(list(MM, 0))
  } else {

  }


  SSig <- gSigma(bb, vv, loca)

  postMM <- logverosa(MM, WW, ddelta, ttheta, yT, TT, f) - 0.5 * t(MM - XX %*% PPs) %*% solve(SSig) %*% (MM - XX %*% PPs)
  postMMprop <- logverosa(MMprop, WW, ddelta, ttheta, yT, TT, f) - 0.5 * t(MMprop - XX %*% PPs) %*% solve(SSig) %*% (MMprop - XX %*% PPs)

  prob <- min(exp((postMMprop) - (postMM)), 1)


  u <- runif(1, 0, 1)

  if (u < prob) {
    Mprox <- MMprop

    rejei <- 1
  } else {
    Mprox <- MM
    rejei <- 0
  }





  res <- as.matrix(Mprox)
  res <- list(Mprox, rejei)
  res
}
# amostrarMsa(delta,theta,W,M,sites,X,Beta,bm,vm,nj,Tt,yT,SU2,f)
amostrardelta <- function(ttheta, ddelta, WW, MM, yT, nn, TT, u1, f, d) {
  deltaprop <- runif(1, max(0, ddelta - u1), min(ddelta + u1, d))


  tema1 <- ifelse((exp(WW + MM) * t(TT)^(exp(MM) - 1)) <= (2 * pi * f * deltaprop), 1, 0)
  tema2 <- ifelse((exp(WW + MM) * (yT[1, ])^(exp(MM) - 1)) <= (2 * pi * f * deltaprop), 1, 0)


  if ((sum(tema1) + sum(tema2)) >= 1) {
    return(list(ddelta, 0))
  } else {

  }



  logp <- logverosa(MM, WW, ddelta, ttheta, yT, TT, f) - 0.5 * (log(ddelta) + log(d - ddelta))

  logpprop <- logverosa(MM, WW, deltaprop, ttheta, yT, TT, f) - 0.5 * (log(deltaprop) + log(d - deltaprop))

  logprob <- logpprop + log(dunif(ddelta, max(0, deltaprop - u1), min(d, deltaprop + u1))) - (logp + log(dunif(deltaprop, max(0, ddelta - u1), min(d, ddelta + u1))))

  prob <- min(c(1, exp(logprob)))

  u <- runif(1, 0, 1)

  if (u < prob) {
    bprox <- deltaprop

    rejei <- 1
  } else {
    bprox <- ddelta

    rejei <- 0
  }



  res <- list(bprox, rejei)
  res
}
# amostrardelta(theta,delta,W,M,yT,nj,Tt,0.01,f,100)
########################################################################
amostrartheta <- function(ttheta, ddelta, WW, MM, yT, nn, TT, u1, f) {
  thetaprop <- runif(1, max(0, ttheta - u1), min(ttheta + u1, 2 * pi))

  logp <- logverosa(MM, WW, ddelta, ttheta, yT, TT, f) - 0.5 * (log(ttheta) + log(2 * pi - ttheta))

  logpprop <- logverosa(MM, WW, ddelta, thetaprop, yT, TT, f) - 0.5 * (log(thetaprop) + log(2 * pi - thetaprop))

  logprob <- logpprop + log(dunif(ttheta, max(0, thetaprop - u1), min(2 * pi, thetaprop + u1))) - (logp + log(dunif(thetaprop, max(0, ttheta - u1), min(2 * pi, ttheta + u1))))

  prob <- min(c(1, exp(logprob)))

  u <- runif(1, 0, 1)

  if (u < prob) {
    bprox <- thetaprop

    rejei <- 1
  } else {
    bprox <- ttheta

    rejei <- 0
  }



  res <- list(bprox, rejei)
  res
}
########################################################################
amostrarf <- function(ttheta, ddelta, WW, MM, yT, nn, TT, u1, a, b, f) {
  fprop <- runif(1, max(a, f - u1), min(f + u1, b))

  logp <- logverosa(MM, WW, ddelta, ttheta, yT, TT, f) - 0.5 * (log(f - a) + log(b - f))

  logpprop <- logverosa(MM, WW, ddelta, ttheta, yT, TT, fprop) - 0.5 * (log(fprop - a) + log(b - fprop))

  logprob <- logpprop + log(dunif(f, max(a, fprop - u1), min(b, fprop + u1))) - (logp + log(dunif(fprop, max(a, f - u1), min(b, f + u1))))

  prob <- min(c(1, exp(logprob)))

  u <- runif(1, 0, 1)

  if (u < prob) {
    bprox <- fprop

    rejei <- 1
  } else {
    bprox <- f

    rejei <- 0
  }



  res <- list(bprox, rejei)
  res
}
## Amostrador de gamma-beta MH
amostrargamaGOEL <- function(ggama, eeta, x, ZZ, MM, NNj, Tt, nnj, A, B, ff) {
  n <- ncol(x)

  ggamaprop <- as.matrix(MASS::mvrnorm(1, ggama, ff * solve(t(ZZ) %*% ZZ)))
  tempbeta <- exp(ZZ %*% ggama)
  tempbetaprop <- exp(ZZ %*% ggamaprop)
  tempalpha <- exp(MM %*% eeta)

  xptalpha <- x
  xproptalpha <- x

  TTalpha <- Tt

  for (j in 1:n) {
    xptalpha[, j] <- tempbeta[j, 1] * (x[, j]^tempalpha[j, 1]) # dados elevados a alpha
    xproptalpha[, j] <- tempbetaprop[j, 1] * (x[, j]^tempalpha[j, 1])
    TTalpha[j, 1] <- TTalpha[j, 1]^tempalpha[j, 1]
  }


  pgama <- sum((ZZ %*% ggama) * as.matrix(nnj)) - sum(xptalpha, na.rm = T) - sum(t(NNj) * tempbeta * TTalpha) - 0.5 * t(ggama - A) %*% solve(B) %*% (ggama - A)
  pgamaprop <- sum((ZZ %*% ggamaprop) * as.matrix(nnj)) - sum(xproptalpha, na.rm = T) - sum(t(NNj) * tempbetaprop * TTalpha) - 0.5 * t(ggamaprop - A) %*% solve(B) %*% (ggamaprop - A)


  logprob <- pgamaprop - pgama

  probac <- min(c(1, exp(logprob)))

  u <- runif(1)

  if (u < probac) {
    res <- ggamaprop
    rejei <- 1
  } else {
    res <- ggama
    rejei <- 0
  }

  res <- list(res, rejei)
  res
}
## Amostrador de eta MH
amostraretaGOEL <- function(ggama, eeta, x, ZZ, MM, NNj, Tt, nnj, A, B, ff)
#######################################
{
  n <- ncol(x)

  eetaprop <- as.matrix(MASS::mvrnorm(1, eeta, ff * solve(t(MM) %*% MM)))

  tempbeta <- exp(ZZ %*% ggama)
  tempalphaprop <- exp(MM %*% eetaprop)
  tempalpha <- exp(MM %*% eeta)

  xptalpha <- x
  xproptalpha <- x

  xp <- x
  xprop <- x

  TTalpha <- Tt
  TTalphaprop <- Tt

  for (j in 1:n) {
    xp[, j] <- x[, j]^tempalpha[j, 1]
    xprop[, j] <- x[, j]^tempalphaprop[j, 1]
    xptalpha[, j] <- tempbeta[j, 1] * (x[, j]^tempalpha[j, 1]) # dados elevados a alpha
    xproptalpha[, j] <- tempbeta[j, 1] * (x[, j]^tempalphaprop[j, 1])
    TTalpha[j, 1] <- TTalpha[j, 1]^tempalpha[j, 1]
    TTalphaprop[j, 1] <- TTalphaprop[j, 1]^tempalphaprop[j, 1]
  }


  peta <- sum((MM %*% eeta) * as.matrix(nnj)) - sum(xptalpha, na.rm = T) - sum(t(NNj) * tempbeta * TTalpha) - 0.5 * t(eeta - A) %*% solve(B) %*% (eeta - A) + sum(log(xp), na.rm = T)
  petaprop <- sum((MM %*% eetaprop) * as.matrix(nnj)) - sum(xproptalpha, na.rm = T) - sum(t(NNj) * tempbeta * TTalphaprop) - 0.5 * t(eetaprop - A) %*% solve(B) %*% (eetaprop - A) + sum(log(xprop), na.rm = T)


  logprob <- petaprop - peta

  probac <- min(c(1, exp(logprob)))

  u <- runif(1)

  if (u < probac) {
    res <- eetaprop
    rejei <- 1
  } else {
    res <- eeta
    rejei <- 0
  }

  res <- list(res, rejei)
  res
}

################################################################################################# 3
amostrarWgoel <- function(W, loca, X, Psi, b, v, nj, N, u1) {
  n <- nrow(W)
  Wprop <- MASS::mvrnorm(1, W, u1 * diag(1, n))
  SSig <- gSigma(b, v, loca)

  postW <- sum(as.matrix(nj) * W) - sum(exp(W)) + sum(t(N) * W) - 0.5 * t(W - X %*% Psi) %*% solve(SSig) %*% (W - X %*% Psi)
  postWprop <- sum(as.matrix(nj) * Wprop) - sum(exp(Wprop)) + sum(t(N) * Wprop) - 0.5 * t(Wprop - X %*% Psi) %*% solve(SSig) %*% (Wprop - X %*% Psi)

  prob <- min(exp((postWprop) - (postW)), 1)


  u <- runif(1, 0, 1)

  if (u < prob) {
    Wprox <- Wprop

    rejei <- 1
  } else {
    Wprox <- W
    rejei <- 0
  }





  res <- as.matrix(Wprox)
  res <- list(Wprox, rejei)
  res
}
################################################################### 3
amostrarPsiGOEL <- function(W, X, Psi, V, M, u1, loca, b, v) {
  n <- nrow(Psi)
  Psiprop <- MASS::mvrnorm(1, Psi, u1 * solve(t(X) %*% X))
  SSig <- gSigma(b, v, loca)

  postPsi <- -0.5 * t(Psi - M) %*% solve(V) %*% (Psi - M) - 0.5 * t(W - X %*% Psi) %*% solve(SSig) %*% (W - X %*% Psi)

  postPsiprop <- -0.5 * t(Psiprop - M) %*% solve(V) %*% (Psiprop - M) - 0.5 * t(W - X %*% Psiprop) %*% solve(SSig) %*% (W - X %*% Psiprop)


  prob <- min(exp((postPsiprop) - (postPsi)), 1)


  u <- runif(1, 0, 1)

  if (u < prob) {
    Psiprox <- Psiprop

    rejei <- 1
  } else {
    Psiprox <- Psi
    rejei <- 0
  }





  res <- as.matrix(Psiprox)
  res <- list(res, rejei)
  res
}


############################################### MUSA





######################################################################### 3
# Comparar valor com elementos de vetor
localiz <- function(vetor, valor) {
  vet <- as.matrix(vetor)
  n <- nrow(vet)

  for (j in 1:n) {
    if (vet[j, 1] == valor) {
      return(j)
    }
  }
}
########################################
sintonizarMUSA <- function(bar, taxa, tau, mat, i) {
  mat <- as.matrix(mat)
  temp1 <- seq(50, bar, 50)
  temp2 <- temp1 - 49
  temp3 <- ifelse(temp1 == i, 1, 0)


  if (sum(temp3) == 1) {
    indi <- localiz(temp1, i)
    mater <- (1 / 50) * sum(mat[temp2[indi]:temp1[indi], 1])
    if (mater >= taxa) {
      delta <- min(0.01, (indi + 1)^(-0.5))
      temp4 <- log(tau) - delta
      temp5 <- exp(temp4)
      return(temp5)
    } else {
      delta <- min(0.01, (indi + 1)^(-0.5))
      temp4 <- log(tau) + delta
      temp5 <- exp(temp4)
      return(temp5)
    }
  } else {
    return(tau)
  }
}
####################################################################################
sintonizarNMUSA <- function(bar, taxa, tau, mat, i) {
  mat <- as.matrix(mat)
  temp1 <- seq(50, bar, 50)
  temp2 <- temp1 - 49
  temp3 <- ifelse(temp1 == i, 1, 0)


  if (sum(temp3) == 1) {
    indi <- localiz(temp1, i)
    mater <- (1 / 50) * sum(mat[temp2[indi]:temp1[indi], 1])
    if (mater >= taxa) {
      delta <- min(0.01, (indi + 1)^(-0.5))
      temp4 <- log(tau) + delta
      temp5 <- exp(temp4)
      return(temp5)
    } else {
      delta <- min(0.01, (indi + 1)^(-0.5))
      temp4 <- log(tau) - delta
      temp5 <- exp(temp4)
      return(temp5)
    }
  } else {
    return(tau)
  }
}
####################################################################################
## Amostrador de alpha MH
amostraralphaMUSA <- function(alpha, W, Tt, c1, d1, x, ff) {
  n <- ncol(x)

  alphaprop <- rgamma(1, alpha * ff, rate = ff)



  palpha <- (c1 - 1) * log(alpha) - d1 * alpha - sum(exp(W) * log(1 + Tt / alpha)) - sum(log(x + alpha), na.rm = T)

  palphaprop <- (c1 - 1) * log(alphaprop) - d1 * alphaprop - sum(exp(W) * log(1 + Tt / alphaprop)) - sum(log(x + alphaprop), na.rm = T)

  logprob <- palphaprop + log(dgamma(alpha, alphaprop * ff, rate = ff)) - (palpha + log(dgamma(alphaprop, alpha * ff, rate = ff)))

  probac <- min(c(1, exp(logprob)))

  u <- runif(1)

  if (u < probac) {
    res <- alphaprop
    rejei <- 1
  } else {
    res <- alpha
    rejei <- 0
  }

  res <- list(res, rejei)
  res
}
################################################################################################# 3
amostrarWMUSA <- function(W, loca, X, Psi, b, v, nj, u1, Tt, alpha) {
  n <- nrow(W)
  Wprop <- MASS::mvrnorm(1, W, u1 * diag(1, n))
  SSig <- gSigma(b, v, loca)

  postW <- sum(as.matrix(nj) * W) - sum(exp(W) * log(1 + Tt / alpha)) - 0.5 * t(W - X %*% Psi) %*% solve(SSig) %*% (W - X %*% Psi)
  postWprop <- sum(as.matrix(nj) * Wprop) - sum(exp(Wprop) * log(1 + Tt / alpha)) - 0.5 * t(Wprop - X %*% Psi) %*% solve(SSig) %*% (Wprop - X %*% Psi)
  prob <- min(exp((postWprop) - (postW)), 1)


  u <- runif(1, 0, 1)

  if (u < prob) {
    Wprox <- Wprop

    rejei <- 1
  } else {
    Wprox <- W
    rejei <- 0
  }





  res <- as.matrix(Wprox)
  res <- list(Wprox, rejei)
  res
}


# Helper function to calculate the observed cumulative mean function
mtnp <- function(dados) {
  a <- c(0, dados)
  m <- length(dados)
  k <- 1
  corr <- array(NA, dim = c(m, 2))
  for (i in 1:m) {
    corr[i, 1] <- (a[i] + a[i + 1]) / 2
    tt <- (a[i] + a[i + 1]) / 2
    rmt <- (1 / k) * (i - 1) + (tt - a[i]) / (k * (a[i + 1] - a[i]))
    corr[i, 2] <- rmt
  }
  corr
}
##############################
MeanFunction <- function(eta, gama, t) {
  res <- gama * t^eta
  res
}
#############################


############### compute mean surface


mfWEIBULL <- function(Wl, Ml, tau) {
  res <- exp(Wl) * tau^(exp(Ml))
  res
}
################################### 3
mfWEIBULLSA <- function(delta, eta, gama, f, ttheta, t) {
  res <- gama * t^eta + delta * cos(2 * pi * f * t + ttheta)
  res
}
################################################

mfSURFACEWEIBULL <- function(Wl, Ml, tau, deltaM, ffM, thetaM) {
  res <- exp(Wl) * tau^(exp(Ml)) + deltaM * cos(2 * pi * (ffM) * tau + thetaM)
  res
}


######################## GOELSA



sintonizarNGOELSAT <- function(bar, taxa, tau, mat, i) {
  mat <- as.matrix(mat)



  mater <- (1 / 50) * sum(mat[(i - 49):i, 1])

  if (mater >= taxa) {
    delta <- min(0.01, (i / 50 + 1)^(-0.5))
    temp4 <- log(tau) + delta
    temp5 <- exp(temp4)
    return(temp5)
  } else {
    delta <- min(0.01, (i / 50 + 1)^(-0.5))
    temp4 <- log(tau) - delta
    temp5 <- exp(temp4)
    return(temp5)
  }
}
######################################################
sintonizarGOELSAT <- function(bar, taxa, tau, mat, i) {
  mat <- as.matrix(mat)



  mater <- (1 / 50) * sum(mat[(i - 49):i, 1])

  if (mater >= taxa) {
    delta <- min(0.01, (i / 50 + 1)^(-0.5))
    temp4 <- log(tau) - delta
    temp5 <- exp(temp4)
    return(temp5)
  } else {
    delta <- min(0.01, (i / 50 + 1)^(-0.5))
    temp4 <- log(tau) + delta
    temp5 <- exp(temp4)
    return(temp5)
  }
}

#############################################################
#############################################################
GG1GOELSA <- function(ggamma, eeta, ddelta, ff, ttheta, WW, ZZ, FF, yT, TT) {
  ffalpha <- FF %*% eeta
  zzbeta <- ZZ %*% ggamma
  aalpha <- exp(FF %*% eeta)
  bbeta <- exp(ZZ %*% ggamma)

  suma <- 0

  for (i in 1:ncol(yT)) {
    res <- sum(log(exp(WW[i, ] + ffalpha[i, ] + zzbeta[i, ] - bbeta[i, ] * yT[, i]^(aalpha[i, ]) + (aalpha[i, ] - 1) * log(yT[, i])) - ddelta * 2 * pi * ff * sin(2 * pi * ff * yT[, i] + ttheta)), na.rm = T)

    suma <- suma + res
  }

  res1 <- suma - sum(exp(WW) * (1 - exp(-bbeta * TT^(aalpha))))
  res1
}
#############################################################
GG2GOELSA <- function(ggamma, eeta, ddelta, ff, ttheta, WW, ZZ, FF, yT, TT) {
  ffalpha <- FF %*% eeta
  zzbeta <- ZZ %*% ggamma
  aalpha <- exp(FF %*% eeta)
  bbeta <- exp(ZZ %*% ggamma)

  suma <- 0

  for (i in 1:ncol(yT)) {
    res <- sum(log(exp(WW[i, ] + ffalpha[i, ] + zzbeta[i, ] - bbeta[i, ] * yT[, i]^(aalpha[i, ]) + (aalpha[i, ] - 1) * log(yT[, i])) - ddelta * 2 * pi * ff * sin(2 * pi * ff * yT[, i] + ttheta)), na.rm = T)

    suma <- suma + res
  }

  res1 <- suma - sum(ddelta * cos(2 * pi * ff * TT + ttheta))
  res1
}
###########################################
amostrarthetaGOELSAT <- function(ggamma, eeta, ddelta, ff, ttheta, WW, ZZ, FF, yT, TT, u1) {
  thetaprop <- runif(1, max(0, ttheta - u1), min(ttheta + u1, 2 * pi))

  logp <- GG2GOELSA(ggamma, eeta, ddelta, ff, ttheta, WW, ZZ, FF, yT, TT) - 0.5 * (log(ttheta) + log(2 * pi - ttheta))

  logpprop <- GG2GOELSA(ggamma, eeta, ddelta, ff, thetaprop, WW, ZZ, FF, yT, TT) - 0.5 * (log(thetaprop) + log(2 * pi - thetaprop))

  logprob <- logpprop + log(dunif(ttheta, max(0, thetaprop - u1), min(2 * pi, thetaprop + u1))) - (logp + log(dunif(thetaprop, max(0, ttheta - u1), min(2 * pi, ttheta + u1))))

  prob <- min(c(1, exp(logprob)))

  u <- runif(1, 0, 1)

  if (u < prob) {
    bprox <- thetaprop

    rejei <- 1
  } else {
    bprox <- ttheta

    rejei <- 0
  }



  res <- list(bprox, rejei)
  res
}
#############################################################
amostrardeltaGOELSAT <- function(ggamma, eeta, ddelta, ff, ttheta, WW, ZZ, FF, yT, TT, u1, d) {
  deltaprop <- runif(1, max(0, ddelta - u1), min(ddelta + u1, d))
  aalpha <- exp(FF %*% eeta)
  bbeta <- exp(ZZ %*% ggamma)

  tema1 <- exp(WW + FF %*% eeta + ZZ %*% ggamma - bbeta * (TT^aalpha) + (aalpha - 1) * log(TT))
  tema1 <- ifelse(tema1 <= (deltaprop * 2 * pi * ff), 1, 0)
  tema2 <- exp(WW + FF %*% eeta + ZZ %*% ggamma - bbeta * (as.matrix(yT[1, ])^aalpha) + (aalpha - 1) * log(as.matrix(yT[1, ])))
  tema2 <- ifelse(tema2 <= (deltaprop * 2 * pi * ff), 1, 0)


  if (sum(tema1 + tema2) >= 1) {
    return(list(ddelta, 0))
  } else {

  }


  logp <- GG2GOELSA(ggamma, eeta, ddelta, ff, ttheta, WW, ZZ, FF, yT, TT) - 0.5 * (log(ddelta) + log(d - ddelta))

  logpprop <- GG2GOELSA(ggamma, eeta, deltaprop, ff, ttheta, WW, ZZ, FF, yT, TT) - 0.5 * (log(deltaprop) + log(d - deltaprop))

  logprob <- logpprop + log(dunif(ddelta, max(0, deltaprop - u1), min(d, deltaprop + u1))) - (logp + log(dunif(deltaprop, max(0, ddelta - u1), min(d, ddelta + u1))))

  prob <- min(c(1, exp(logprob)))

  u <- runif(1, 0, 1)

  if (u < prob) {
    bprox <- deltaprop

    rejei <- 1
  } else {
    bprox <- ddelta

    rejei <- 0
  }



  res <- list(bprox, rejei)
  res
}
#############################################################
amostrarfGOELSAT <- function(ggamma, eeta, ddelta, ff, ttheta, WW, ZZ, FF, yT, TT, u1, a, b) {
  fprop <- runif(1, max(a, ff - u1), min(ff + u1, b))

  logp <- GG2GOELSA(ggamma, eeta, ddelta, ff, ttheta, WW, ZZ, FF, yT, TT) - 0.5 * (log(ff - a) + log(b - ff))

  logpprop <- GG2GOELSA(ggamma, eeta, ddelta, fprop, ttheta, WW, ZZ, FF, yT, TT) - 0.5 * (log(fprop - a) + log(b - fprop))

  logprob <- logpprop + log(dunif(ff, max(a, fprop - u1), min(b, fprop + u1))) - (logp + log(dunif(fprop, max(a, ff - u1), min(b, ff + u1))))

  prob <- min(c(1, exp(logprob)))

  u <- runif(1, 0, 1)

  if (u < prob) {
    bprox <- fprop

    rejei <- 1
  } else {
    bprox <- ff

    rejei <- 0
  }



  res <- list(bprox, rejei)
  res
}
######################################################
amostrarWGOELSAT <- function(ggamma, eeta, ddelta, ff, ttheta, WW, ZZ, FF, yT, TT, bb, vv, loca, u1, XX, PPs) {
  n <- nrow(WW)
  WWprop <- as.matrix(MASS::mvrnorm(1, WW, u1 * diag(1, n)))
  aalpha <- exp(FF %*% eeta)
  bbeta <- exp(ZZ %*% ggamma)

  tema1 <- exp(WWprop + FF %*% eeta + ZZ %*% ggamma - bbeta * (TT^aalpha) + (aalpha - 1) * log(TT))
  tema1 <- ifelse(tema1 <= (ddelta * 2 * pi * ff), 1, 0)
  tema2 <- exp(WWprop + FF %*% eeta + ZZ %*% ggamma - bbeta * (as.matrix(yT[1, ])^aalpha) + (aalpha - 1) * log(as.matrix(yT[1, ])))
  tema2 <- ifelse(tema2 <= (ddelta * 2 * pi * ff), 1, 0)

  if (sum(tema1 + tema2) >= 1) {
    return(list(WW, 0))
  } else {

  }


  SSig <- gSigma(bb, vv, loca)

  postWW <- GG1GOELSA(ggamma, eeta, ddelta, ff, ttheta, WW, ZZ, FF, yT, TT) - 0.5 * t(WW - XX %*% PPs) %*% solve(SSig) %*% (WW - XX %*% PPs)
  postWWprop <- GG1GOELSA(ggamma, eeta, ddelta, ff, ttheta, WWprop, ZZ, FF, yT, TT) - 0.5 * t(WWprop - XX %*% PPs) %*% solve(SSig) %*% (WWprop - XX %*% PPs)

  prob <- min(exp((postWWprop) - (postWW)), 1)


  u <- runif(1, 0, 1)

  if (u < prob) {
    Wprox <- WWprop

    rejei <- 1
  } else {
    Wprox <- WW
    rejei <- 0
  }





  res <- as.matrix(Wprox)
  res <- list(Wprox, rejei)
  res
}
##############################################
amostrargamaGOELSAT <- function(ggamma, eeta, ddelta, ff, ttheta, WW, ZZ, FF, yT, TT, A, B, loca, u1) {
  ggamaprop <- as.matrix(MASS::mvrnorm(1, ggamma, u1 * solve(t(ZZ) %*% ZZ)))

  aalpha <- exp(FF %*% eeta)
  bbeta <- exp(ZZ %*% ggamaprop)

  tema1 <- exp(WW + FF %*% eeta + ZZ %*% ggamaprop - bbeta * (TT^aalpha) + (aalpha - 1) * log(TT))
  tema2 <- exp(WW + FF %*% eeta + ZZ %*% ggamaprop - bbeta * (as.matrix(yT[1, ])^aalpha) + (aalpha - 1) * log(as.matrix(yT[1, ])))


  tema1 <- ifelse(tema1 <= (ddelta * 2 * pi * ff), 1, 0)
  tema2 <- ifelse(tema2 <= (ddelta * 2 * pi * ff), 1, 0)

  if (sum(tema1 + tema2) >= 1) {
    return(list(ggamma, 0))
  } else {

  }


  pgama <- GG1GOELSA(ggamma, eeta, ddelta, ff, ttheta, WW, ZZ, FF, yT, TT) - 0.5 * t(ggamma - A) %*% solve(B) %*% (ggamma - A)
  pgamaprop <- GG1GOELSA(ggamaprop, eeta, ddelta, ff, ttheta, WW, ZZ, FF, yT, TT) - 0.5 * t(ggamaprop - A) %*% solve(B) %*% (ggamaprop - A)


  logprob <- pgamaprop - pgama

  probac <- min(c(1, exp(logprob)))

  u <- runif(1)

  if (u < probac) {
    res <- ggamaprop
    rejei <- 1
  } else {
    res <- ggamma
    rejei <- 0
  }

  res <- list(res, rejei)
  res
}
#########################
amostraretaGOELSAT <- function(ggamma, eeta, ddelta, ff, ttheta, WW, ZZ, FF, yT, TT, A, B, loca, u1) {
  etaprop <- as.matrix(MASS::mvrnorm(1, eeta, u1 * solve(t(FF) %*% FF)))

  aalpha <- exp(FF %*% etaprop)
  bbeta <- exp(ZZ %*% ggamma)

  tema1 <- exp(WW + FF %*% etaprop + ZZ %*% ggamma - bbeta * (TT^aalpha) + (aalpha - 1) * log(TT))
  tema1 <- ifelse(tema1 <= (ddelta * 2 * pi * ff), 1, 0)
  tema2 <- exp(WW + FF %*% etaprop + ZZ %*% ggamma - bbeta * (as.matrix(yT[1, ])^aalpha) + (aalpha - 1) * log(as.matrix(yT[1, ])))
  tema2 <- ifelse(tema2 <= (ddelta * 2 * pi * ff), 1, 0)

  if (sum(tema1 + tema2) >= 1) {
    return(list(eeta, 0))
  } else {

  }


  pgama <- GG1GOELSA(ggamma, eeta, ddelta, ff, ttheta, WW, ZZ, FF, yT, TT) - 0.5 * t(eeta - A) %*% solve(B) %*% (eeta - A)
  pgamaprop <- GG1GOELSA(ggamma, etaprop, ddelta, ff, ttheta, WW, ZZ, FF, yT, TT) - 0.5 * t(etaprop - A) %*% solve(B) %*% (etaprop - A)


  logprob <- pgamaprop - pgama

  probac <- min(c(1, exp(logprob)))

  u <- runif(1)

  if (u < probac) {
    res <- etaprop
    rejei <- 1
  } else {
    res <- eeta
    rejei <- 0
  }

  res <- list(res, rejei)
  res
}
############################
amostrarbGOELSAT <- function(WW, vv, bb, loca, ab, ba, XX, PPsi, u1) {
  bprop <- rgamma(1, shape = bb * u1, rate = u1)

  SSigprop <- gSigma(bprop, vv, loca)

  if ((det(SSigprop) == 0) | (bprop < 0.005)) {
    return(list(bb, 0))
  }

  SSig <- gSigma(bb, vv, loca)
  SSigprop <- gSigma(bprop, vv, loca)


  logp <- -0.5 * t(WW - XX %*% PPsi) %*% solve(SSig) %*% (WW - XX %*% PPsi) - 0.5 * log(det(SSig)) + (ab - 1) * log(bb) - ba * bb

  logpprop <- -0.5 * t(WW - XX %*% PPsi) %*% solve(SSigprop) %*% (WW - XX %*% PPsi) - 0.5 * log(det(SSigprop)) + (ab - 1) * log(bprop) - ba * bprop

  logprob <- logpprop + log(dgamma(bb, shape = bprop * u1, rate = u1)) - (logp + log(dgamma(bprop, shape = bb * u1, rate = u1)))
  prob <- min(c(1, exp(logprob)))

  u <- runif(1, 0, 1)

  if (u < prob) {
    bprox <- bprop

    rejei <- 1
  } else {
    bprox <- bb

    rejei <- 0
  }



  res <- list(bprox, rejei)
  res
}



############################ MUSA OKUMOTO SA


################################################################################################# 3
amostrarWMUSAsa <- function(aalpha, WW, ddelta, ttheta, data, TT, ff, loca, X, Psi, b, v, u1) {
  n <- nrow(WW)
  Wprop <- as.matrix(MASS::mvrnorm(1, WW, u1 * diag(1, n)))

  tema <- exp(Wprop) / (aalpha + TT)
  tema <- ifelse(tema <= (ddelta * 2 * pi * ff), 1, 0)


  if (sum(tema) >= 1) {
    return(list(WW, 0))
  } else {

  }



  SSig <- gSigma(b, v, loca)

  postW <- logveroMUSAsa(aalpha, WW, ddelta, ttheta, data, TT, ff) - 0.5 * t(WW - X %*% Psi) %*% solve(SSig) %*% (WW - X %*% Psi)
  postWprop <- logveroMUSAsa(aalpha, Wprop, ddelta, ttheta, data, TT, ff) - 0.5 * t(Wprop - X %*% Psi) %*% solve(SSig) %*% (Wprop - X %*% Psi)
  prob <- min(exp((postWprop) - (postW)), 1)


  u <- runif(1, 0, 1)

  if (u < prob) {
    Wprox <- Wprop

    rejei <- 1
  } else {
    Wprox <- WW
    rejei <- 0
  }





  res <- as.matrix(Wprox)
  res <- list(Wprox, rejei)
  res
}
# amostrarWsa(alpha,W,delta,theta,data,Tt,1/360,sites,X,Psi,b,v,SU5)
################
logveroMUSAsa <- function(aalpha, WW, ddelta, ttheta, data, TT, ff) {
  soma <- 0
  for (j in 1:ncol(data)) {
    temp1 <- sum(log(exp(WW[j, 1]) / (data[, j] + aalpha) - ddelta * 2 * pi * ff * sin(2 * pi * ff * data[, j] + ttheta)), na.rm = T)
    soma <- soma + temp1
  }

  res <- soma + sum(-(exp(WW) * log(1 + TT / aalpha) + ddelta * cos(2 * pi * ff * TT + ttheta)))
  res
}
# logverosa(alpha,W,0,theta,data,Tt,1/360)
##############################################
amostraralphaMUSAsa <- function(aalpha, WW, TT, ttheta, ddelta, c1, d1, data, ff, u1) {
  alphaprop <- rgamma(1, aalpha * u1, rate = u1)

  tema <- exp(WW) / (alphaprop + TT)
  tema <- ifelse(tema <= (ddelta * 2 * pi * ff), 1, 0)


  if (sum(tema) >= 1) {
    return(list(aalpha, 0))
  } else {

  }


  palpha <- (c1 - 1) * log(aalpha) - d1 * aalpha + logveroMUSAsa(aalpha, WW, ddelta, ttheta, data, TT, ff)

  palphaprop <- (c1 - 1) * log(alphaprop) - d1 * alphaprop + logveroMUSAsa(alphaprop, WW, ddelta, ttheta, data, TT, ff)

  logprob <- palphaprop + log(dgamma(aalpha, alphaprop * u1, rate = u1)) - (palpha + log(dgamma(alphaprop, aalpha * u1, rate = u1)))

  probac <- min(c(1, exp(logprob)))

  u <- runif(1)

  if (u < probac) {
    res <- alphaprop
    rejei <- 1
  } else {
    res <- aalpha
    rejei <- 0
  }

  res <- list(res, rejei)
  res
}
#############################################################
amostrardeltaMUSAsa <- function(ttheta, ddelta, WW, aalpha, yT, nn, TT, u1, f, d) {
  deltaprop <- runif(1, max(0, ddelta - u1), min(ddelta + u1, d))

  tema <- exp(WW) / (aalpha + TT)
  tema <- ifelse(tema <= (deltaprop * 2 * pi * f), 1, 0)


  if (sum(tema) >= 1) {
    return(list(ddelta, 0))
  } else {

  }


  logp <- logveroMUSAsa(aalpha, WW, ddelta, ttheta, yT, TT, f) - 0.5 * (log(ddelta) + log(d - ddelta))

  logpprop <- logveroMUSAsa(aalpha, WW, deltaprop, ttheta, yT, TT, f) - 0.5 * (log(deltaprop) + log(d - deltaprop))

  logprob <- logpprop + log(dunif(ddelta, max(0, deltaprop - u1), min(d, deltaprop + u1))) - (logp + log(dunif(deltaprop, max(0, ddelta - u1), min(d, ddelta + u1))))

  prob <- min(c(1, exp(logprob)))

  u <- runif(1, 0, 1)

  if (u < prob) {
    bprox <- deltaprop

    rejei <- 1
  } else {
    bprox <- ddelta

    rejei <- 0
  }



  res <- list(bprox, rejei)
  res
}
########################################################################
amostrarthetaMUSAsa <- function(ttheta, ddelta, WW, aalpha, yT, nn, TT, u1, f) {
  thetaprop <- runif(1, max(0, ttheta - u1), min(ttheta + u1, 2 * pi))

  logp <- logveroMUSAsa(aalpha, WW, ddelta, ttheta, yT, TT, f) - 0.5 * (log(ttheta) + log(2 * pi - ttheta))

  logpprop <- logveroMUSAsa(aalpha, WW, ddelta, thetaprop, yT, TT, f) - 0.5 * (log(thetaprop) + log(2 * pi - thetaprop))

  logprob <- logpprop + log(dunif(ttheta, max(0, thetaprop - u1), min(2 * pi, thetaprop + u1))) - (logp + log(dunif(thetaprop, max(0, ttheta - u1), min(2 * pi, ttheta + u1))))

  prob <- min(c(1, exp(logprob)))

  u <- runif(1, 0, 1)

  if (u < prob) {
    bprox <- thetaprop

    rejei <- 1
  } else {
    bprox <- ttheta

    rejei <- 0
  }



  res <- list(bprox, rejei)
  res
}
####################################
amostrarfMUSAsa <- function(ttheta, ddelta, WW, aalpha, yT, nn, TT, u1, a, b, ff) {
  fprop <- runif(1, max(a, ff - u1), min(ff + u1, b))

  logp <- logveroMUSAsa(aalpha, WW, ddelta, ttheta, yT, TT, ff) - 0.5 * (log(ff - a) + log(b - ff))

  logpprop <- logveroMUSAsa(aalpha, WW, ddelta, ttheta, yT, TT, fprop) - 0.5 * (log(fprop - a) + log(b - fprop))

  logprob <- logpprop + log(dunif(ff, max(a, fprop - u1), min(b, fprop + u1))) - (logp + log(dunif(fprop, max(a, ff - u1), min(b, ff + u1))))

  prob <- min(c(1, exp(logprob)))

  u <- runif(1, 0, 1)

  if (u < prob) {
    bprox <- fprop

    rejei <- 1
  } else {
    bprox <- ff

    rejei <- 0
  }



  res <- list(bprox, rejei)
  res
}

######################### compute mean surface musa ######
mfMUSA <- function(W, alpha, t) {
  x <- exp(W) * log(1 + t / alpha)
  return(x)
}



mfMUSASA <- function(W, alpha, t, delta, f, ttheta) {
  x <- exp(W) * log(1 + t / alpha) + delta * cos(2 * pi * f * t + ttheta)
  return(x)
}



#################### funções goel

mfGOEL <- function(W, Beta, Alpha, t) {
  res <- exp(W) * (1 - exp(-Beta * t^Alpha))
  return(res)
}


mfGOELSA <- function(W, Beta, Alpha, t, delta, f, ttheta) {
  res <- exp(W) * (1 - exp(-Beta * t^Alpha)) + delta * cos(2 * pi * f * t + ttheta)
  return(res)
}



######################### valores iniciais theta e delta

njfunction <- function(data, a, b) {
  nnj <- NULL

  for (i in 1:ncol(data)) {
    res <- sum(ifelse((data[, i] > a) & (data[, i] <= b), 1, 0), na.rm = T)
    nnj <- c(nnj, res)
  }

  Mdata <- array(NA, dim = c(max(nnj), ncol(data)))

  for (i in 1:ncol(data)) {
    Mdata[, i] <- c(data[which((data[, i] > a) & (data[, i] <= b)), i], rep(NA, (max(nnj) - length(which((data[, i] > a) & (data[, i] <= b))))))
  }

  res <- list(Mdata, nnj)
  res
}

##########################
logverouniWEIBULL <- function(xx, Tao, nnt, datos) {
  theta <- xx[1]
  alpha <- xx[2]
  res <- -(-theta * Tao^alpha + nnt * log(theta) + nnt * log(alpha) + alpha * sum(log(datos), na.rm = T))
  res
}
##################################
logverouniMUSA <- function(xx, data, Nn, TT) {
  aalpha <- xx[1]
  WW <- xx[2]
  res <- -(sum(Nn * WW) - sum(log(aalpha + data), na.rm = T) - sum(exp(WW) * log(1 + TT / aalpha)))
  res
}
