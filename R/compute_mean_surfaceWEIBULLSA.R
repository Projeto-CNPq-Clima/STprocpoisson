#' Compute Interpolated Mean Surface for Spatiotemporal Model
#'
#' This function interpolates the mean values on a grid of points (`DNO`) for a spatiotemporal
#' nonhomogeneous Poisson model with Weibull intensity with a seasonal component.. The function utilizes the MCMC outputs
#' from `STModelWeibullMCMCSA` and applies a Gaussian process-based interpolation.
#'
#' @param resultsSA A list containing the output from `STModelWeibullMCMCSA`, including:
#'   - `MW`: Samples for parameter W.
#'   - `MWT`: Acceptance indicators for parameter W.
#'   - `MMj`: Samples for parameter M.
#'   - `MMT`: Acceptance indicators for parameter M.
#'   - `Mvw`: Samples for parameter sigma^2_w.
#'   - `Mvm`: Samples for parameter sigma^2_m.
#'   - `MBeta`: Samples for parameter Beta.
#'   - `Mbw`: Samples for parameter phi_w.
#'   - `MbwT`: Acceptance indicators for phi_w.
#'   - `Mbm`: Samples for parameter phi_m.
#'   - `MbmT`: Acceptance indicators for phi_m.
#'   - `MPsi`: Samples for parameter Psi.
#'   - `Mdelta` : Samples of parameter delta obtained during the MCMC procedure (iteration - burnin).
#'   - `MdeltaT` : A binary vector indicating whether each proposed value of delta was accepted (1) or rejected (0). Used to compute the acceptance rate for delta.
#'   - `Mtheta` : Samples of parameter theta obtained during the MCMC procedure (iteration - burnin).
#'   - `MthetaT` : A binary vector indicating whether each proposed value of delta was accepted (1) or rejected (0). Used to compute the acceptance rate for theta.
#'   - `Mf` : Samples of parameter f obtained during the MCMC procedure (iteration - burnin).
#'   - `MfT` : A binary vector indicating whether each proposed value of delta was accepted (1) or rejected (0). Used to compute the acceptance rate for f.
#'
#' @param sites A matrix with geographic coordinates of the monitoring stations.
#' @param DNO A grid of points where interpolation is to be performed.
#' @param CovXNO Covariates for the scale parameter at the grid points.
#' @param CovZNO Covariates for the shape parameter at the grid points.
#' @param tau A vector of temporal points for which the mean surface is computed.
#'
#' @return A list: (`Surface`) containing the interpolated mean values at the grid points
#' for each temporal point in `tau`. The first column contains the mean values at the initial
#' time step, and subsequent columns contain the differences between consecutive time steps.
#' @export

compute_mean_surfaceWEIBULLSA <- function(resultsSA, sites, DNO, CovXNO, CovZNO,tau) {


  jj=nrow(DNO)
  res=rbind(DNO,as.matrix(sites))
  tt=nrow(res)

  Uns=as.matrix(rep(1,tt))
  Xr=as.matrix(cbind(Uns,res))
  Zr=Xr

  MPsi<-resultsSA$MPsi

  MWNO=array(NA,dim=c(nrow(MPsi),jj))
  MMNO=array(NA,dim=c(nrow(MPsi),jj))
  MfmNO=NULL

  MW<-resultsSA$MW
  Mbw<-resultsSA$Mbw
  Mbv<-resultsSA$Mbv
  MMj<-resultsSA$MMj
  MBeta<-resultsSA$MBeta
  Mvw<-resultsSA$Mvw
  Mbm<-resultsSA$Mbm
  Mvm<-resultsSA$Mvm

  for(h in 1:nrow(MW)){

    WM=as.matrix(MW[h,])
    sigma=gSigma(Mbw[h],Mvw[h],res)
    Psih=as.matrix(MPsi[h,])

    MMr=as.matrix(MMj[h,])
    sigmam=gSigma(Mbm[h],Mvm[h],res)
    Betah=as.matrix(MBeta[h,])

    Mpro=Xr%*%as.matrix(Psih)

    A1=as.matrix(Mpro[(jj+1):tt,])
    A2=as.matrix(Mpro[1:jj,])

    SSA1=sigma[(jj+1):tt,(jj+1):tt]
    SSA2=sigma[1:jj,1:jj]
    SSA12=sigma[1:jj,(jj+1):tt]

    A2est=A2+SSA12%*%solve(SSA1)%*%(as.matrix(WM)-A1)
    SSA2est=SSA2-SSA12%*%solve(SSA1)%*%t(SSA12)
    WNO=as.matrix(MASS::mvrnorm(1,A2est,SSA2est))
    MWNO[h,]=t(WNO)

    Mprom=Zr%*%as.matrix(Betah)

    A1m=as.matrix(Mprom[(jj+1):tt,])
    A2m=as.matrix(Mprom[1:jj,])


    SSA1m=sigmam[(jj+1):tt,(jj+1):tt]
    SSA2m=sigmam[1:jj,1:jj]
    SSA12m=sigmam[1:jj,(jj+1):tt]

    A2estm=A2m+SSA12m%*%solve(SSA1m)%*%(as.matrix(MMr)-A1m)
    SSA2estm=SSA2m-SSA12m%*%solve(SSA1m)%*%t(SSA12m)
    MNO=as.matrix(MASS::mvrnorm(1,A2estm,SSA2estm))
    MMNO[h,]=t(MNO)



  }

  #######################################
  Meanmf=NULL
  Mdelta<-resultsSA$Mdelta
  Mtheta<-resultsSA$Mtheta
  Mf<-resultsSA$Mf
  for(i in 1:length(tau)){

    for(h in 1:nrow(MW)){

      MfmNO=rbind(MfmNO,t(as.matrix(mfSURFACEWEIBULL(MWNO[h,],MMNO[h,],tau[i],Mdelta[h],Mf[h],Mtheta[h]))))

    }

    Meanmf=cbind(Meanmf,as.matrix(apply(MfmNO,2,mean)))
    MfmNO=NULL
  }
  ##################################
  ## Interpolação da função m(t,s)
  ###################################3
  Surface=array(NA,dim=c(nrow(DNO),10))
  Surface[,1]=Meanmf[,1]


  for(kk in 2:10){
    minterest=Meanmf[,kk]-Meanmf[,(kk-1)]
    Surface[,kk]=minterest



  }

  output <- list(Surface,MMNO,MWNO)
  names(output) <- c("Surface","MMNO","MWNO")

  return(output)
}



