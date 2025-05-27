#' STModelGoelMCMC: A Bayesian Space-Time Model Using MCMC
#'
#' This function implements a Bayesian space-time model using MCMC for failure time data across geographical locations.
#' It estimates parameters associated with different components of the model, including covariates, spatial dependencies, and prior distributions.
#'
#' @param data A matrix of failure times, where each column represents a station.
#' @param sites A matrix of geographic coordinates for the stations (e.g., longitude and latitude).
#' @param X A matrix of covariates associated with the parameter W. Default is a column of ones and the coordinates from `sites`.
#' @param Z A matrix of covariates associated with the parameter beta. Default is the same as `X`.
#' @param Fj A matrix of covariates associated with the parameter alpha. Default is the same as `Z`.
#' @param prior A list specifying the hyperparameters for the prior distributions:
#'   - `c3`, `d3`: Shape and rate parameters for the prior of phi_w.
#'   - `BB1`, `AA1`: Mean vector and covariance matrix for the prior of gamma.
#'   - `BB2`, `AA2`: Mean vector and covariance matrix for the prior of eta. Defaults are equal to `BB1` and `AA1`.
#'   - `aa1`, `bb1`: Shape and scale parameters for the inverse-Gamma prior of sigma^2_w.
#'   - `V`, `MM1`: Covariance matrix and mean vector for the prior of Psi.
#' @param iteration The total number of MCMC iterations.
#' @param burnin The number of burn-in iterations to discard.
#'
#' @return A list containing the following components:
#'   - `Mgama`: MCMC chain for the parameter gamma.
#'   - `MgamaT`: A vector of zeros and ones indicating acceptance (1) or rejection (0) of gamma proposals in the Metropolis-Hastings algorithm.
#'   - `Meta`: MCMC chain for the parameter eta.
#'   - `MetaT`: A vector of zeros and ones indicating acceptance (1) or rejection (0) of eta proposals in the Metropolis-Hastings algorithm.
#'   - `Mv`: MCMC chain for the parameter sigma^2.
#'   - `Mb`: MCMC chain for the parameter phi.
#'   - `MbT`: A vector of zeros and ones indicating acceptance (1) or rejection (0) of phi proposals in the Metropolis-Hastings algorithm.
#'   - `MW`: MCMC chain for the parameter W.
#'   - `MWT`: A vector of zeros and ones indicating acceptance (1) or rejection (0) of W proposals in the Metropolis-Hastings algorithm.
#'   - `MPsi`: MCMC chain for the parameter Psi.
#'   - `MPsiT`: A vector of zeros and ones indicating acceptance (1) or rejection (0) of Psi proposals in the Metropolis-Hastings algorithm.
#'
#' @export
STModelGoelMCMCSATEST<- function(data, sites,X=cbind(as.matrix(rep(1,ncol(data))),sites),
                             Z=cbind(as.matrix(rep(1,ncol(data))),sites),
                             Fj=Z,prior=list(BB1=diag(100,ncol(Z)),
                                             AA1=as.matrix(rep(0,ncol(Z))),
                                             BB2=diag(100,ncol(Z)),
                                             AA2=as.matrix(rep(0,ncol(Z))),
                                             A=as.matrix(rep(0,ncol(Z))),
                                             B=diag(100,ncol(Z)),
                                             lgama=as.matrix(c(-3.17022118,0.06602419,1.33314203)),
                                             leta=as.matrix(c(1.384896960,0.032667601,0.004153729)),

                                             aa1=2.01,
                                             bb1=1.005,
                                             c3=(-2*log(0.05)/max(dist(sites)))*0.1,
                                             d3=0.1),iteration,burnin){








  theta=3.525654
  delta=0.01
  f=1/ 366.8474
  W=as.matrix(rep(9,ncol(data)))
  b=1
  v=1
  Psi=as.matrix(rep(0,ncol(X)))


  BB1=prior$BB1
  AA1=prior$AA1
  BB2=prior$BB2
  AA2=prior$AA2
  A=prior$A
  B=prior$B
  lgama=prior$lgama
  leta=prior$leta

  aa1=prior$aa1
  bb1=prior$bb1
  c3=prior$c3
  d3=prior$d3


  n=ncol(data)
  m=nrow(data)
  tempdados=is.na(data)
  nj=m-apply(tempdados,2,sum)

  Tt=array(NA,dim=c(1,n))
  for(y in 1:n){
    Tt[1,y]=data[nj[y],y]
  }

  SU1=0.00001
  SU2=0.0001
  SU3=0.0001
  SU4=100

  MW=NULL
  MWT=NULL
  Meta=NULL
  MetaT=NULL
  Mgamma=NULL
  MgammaT=NULL
  MPsi=NULL
  Mv=NULL
  Mb=NULL
  MbT=NULL
  Mf=NULL
  MfT=NULL
  MdeltaT=NULL
  Mdelta=NULL
  MthetaT=NULL
  Mtheta=NULL


  for(j in 1:iteration){

    if(j<=burnin){
      temp=amostrarWGOELSAT(lgama,leta,delta,f,theta,W,Z,Fj,data,t(Tt),b,v,sites,SU1,X,Psi)
      W=temp[[1]]
      MWT=c(MWT,temp[[2]])

      temp=amostrardeltaGOELSAT(lgama,leta,delta,f,theta,W,Z,Fj,data,t(Tt),0.01,100)
      delta=temp[[1]]
      MdeltaT=c(MdeltaT,temp[[2]])

      temp=amostraretaGOELSAT(lgama,leta,delta,f,theta,W,Z,Fj,data,t(Tt),AA1,BB1,sites,SU2)
      leta=temp[[1]]
      MetaT=c(MetaT,temp[[2]])

      temp=amostrargamaGOELSAT(lgama,leta,delta,f,theta,W,Z,Fj,data,t(Tt),AA1,BB1,sites,SU3)
      lgama=temp[[1]]
      MgammaT=c(MgammaT,temp[[2]])

      varPsi=solve(solve(B)+t(X)%*%solve(gSigma(b,v,sites))%*%X)
      medPsi=(t(A)%*%solve(B)+t(W)%*%solve(gSigma(b,v,sites))%*%X)%*%varPsi
      Psi=as.matrix(MASS::mvrnorm(1,medPsi,varPsi))

      RR=gCorr(b,sites)
      aa=(n/2)+aa1
      bb=0.5*t(W-X%*%Psi)%*%solve(RR)%*%(W-X%*%Psi)+bb1

      v=1/rgamma(1,shape=aa, rate = bb)

      temp=amostrarbGOELSAT(W,v,b,sites,c3,d3,X,Psi,SU4)
      b=temp[[1]]
      MbT=c(MbT,temp[[2]])

      temp=amostrarfGOELSAT(lgama,leta,delta,f,theta,W,Z,Fj,data,t(Tt),0.00001/2,1/(365+10),1/(365-10))
      f=temp[[1]]
      MfT=c(MfT,temp[[2]])

      temp=amostrarthetaGOELSAT(lgama,leta,delta,f,theta,W,Z,Fj,data,t(Tt),0.05)
      theta=temp[[1]]
      MthetaT=c(MthetaT,temp[[2]])

      if((j%%50)==0){
        SU1=sintonizarNGOELSAT(burnin,0.20,SU1,MWT,j)
        SU2=sintonizarNGOELSAT(burnin,0.20,SU2,MetaT,j)
        SU3=sintonizarNGOELSAT(burnin,0.20,SU3,MgammaT,j)
        SU4=sintonizarGOELSAT(burnin,0.44,SU4,MbT,j)

      }else{

      }
      print(j)
    }else{

      temp=amostrarWGOELSAT(lgama,leta,delta,f,theta,W,Z,Fj,data,t(Tt),b,v,sites,SU1,X,Psi)
      W=temp[[1]]
      MWT=c(MWT,temp[[2]])
      MW=rbind(MW,t(W))

      temp=amostrardeltaGOELSAT(lgama,leta,delta,f,theta,W,Z,Fj,data,t(Tt),0.01,100)
      delta=temp[[1]]
      MdeltaT=c(MdeltaT,temp[[2]])
      Mdelta=c(Mdelta,delta)

      temp=amostraretaGOELSAT(lgama,leta,delta,f,theta,W,Z,Fj,data,t(Tt),AA1,BB1,sites,SU2)
      leta=temp[[1]]
      MetaT=c(MetaT,temp[[2]])
      Meta=rbind(Meta,t(leta))

      temp=amostrargamaGOELSAT(lgama,leta,delta,f,theta,W,Z,Fj,data,t(Tt),AA1,BB1,sites,SU2)
      lgama=temp[[1]]
      MgammaT=c(MgammaT,temp[[2]])
      Mgamma=rbind(Mgamma,t(lgama))

      varPsi=solve(solve(B)+t(X)%*%solve(gSigma(b,v,sites))%*%X)
      medPsi=(t(A)%*%solve(B)+t(W)%*%solve(gSigma(b,v,sites))%*%X)%*%varPsi
      Psi=as.matrix(MASS::mvrnorm(1,medPsi,varPsi))
      MPsi=rbind(MPsi,t(Psi))

      RR=gCorr(b,sites)
      aa=(n/2)+aa1
      bb=0.5*t(W-X%*%Psi)%*%solve(RR)%*%(W-X%*%Psi)+bb1

      v=1/rgamma(1,shape=aa, rate = bb)
      Mv=c(Mv,v)

      temp=amostrarbGOELSAT(W,v,b,sites,c3,d3,X,Psi,SU4)
      b=temp[[1]]
      MbT=c(MbT,temp[[2]])
      Mb=c(Mb,b)

      temp=amostrarfGOELSAT(lgama,leta,delta,f,theta,W,Z,Fj,data,t(Tt),0.00001/2,1/(365+10),1/(365-10))
      f=temp[[1]]
      MfT=c(MfT,temp[[2]])
      Mf=c(Mf,f)

      temp=amostrarthetaGOELSAT(lgama,leta,delta,f,theta,W,Z,Fj,data,t(Tt),0.05)
      theta=temp[[1]]
      MthetaT=c(MthetaT,temp[[2]])
      Mtheta=c(Mtheta,theta)


      print(j)
    }

  }
  resul<-list(MW,MWT,Meta,MetaT,Mgamma,MgammaT,MPsi,Mv,Mb,MbT,Mf,MfT,MdeltaT,Mdelta,MthetaT,Mtheta)
  names(resul)<-c("MW","MWT","Meta","MetaT","Mgama","MgamaT","MPsi","Mv","Mb","MbT","Mf","MfT","MdeltaT","Mdelta","MthetaT","Mtheta")
  return(resul)
}

