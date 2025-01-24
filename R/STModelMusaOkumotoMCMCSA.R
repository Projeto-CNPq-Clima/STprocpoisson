#' Title
#'
#' @param data aaaa
#' @param sites aaaa
#' @param X aaaa
#' @param prior aaaa
#' @param iteration aaaa
#' @param burnin aaaa
#'
#' @return aaaa
#' @export

STModelMusaOkumotoMCMCSA <- function(data, sites,X=cbind(as.matrix(rep(1,ncol(data))),(1/100)*sites),
                                     prior = list(Psi=as.matrix(rep(0,ncol(X))),
                                                  A=as.matrix(rep(0,ncol(X))),
                                                  B=diag(100,ncol(X)),
                                                  V=diag(100,ncol(X)),
                                                  M=as.matrix(rep(0,ncol(X))),
                                                  c2=1e-05,
                                                  d2=0.001,
                                                  c4=2.01,
                                                  d4=1.005,
                                                  aa2=2.01,
                                                  bb2=1.005,
                                                  aa1=2.01,
                                                  bb1=1.005,
                                                  c3=(-2*log(0.05)/max(dist(sites)))*0.1,
                                                  d3=0.1
                                     ), iteration, burnin) {



b=1
v=1
alpha=5000
theta=pi/2
delta=0.001
f=1/365

W=as.matrix(rep(10,ncol(data)))

#Hiperparametros
aa1=prior$aa1
bb1=prior$bb1

c2=prior$c2
d2=prior$d2

c4=prior$c4
d4=prior$d4

aa2=prior$aa2
bb2=prior$bb2

c3=prior$c3
d3=prior$d3

Psi=prior$Psi
A=prior$A
B=prior$B

V=prior$V
M=prior$M

SU2=100
SU3=28.01968
SU5= 0.001
n=ncol(data)
m=nrow(data)
tempdados=is.na(data)
nj=m-apply(tempdados,2,sum)
#p=ncol(X)



Tt=array(NA,dim=c(1,n))
for(y in 1:n){
  Tt[1,y]=data[nj[y],y]
}

#pp=ncol(X)
Tt=t(Tt)

Malpha=NULL
MalphaT=NULL

Mtheta=NULL
MthetaT=NULL

Mdelta=NULL
MdeltaT=NULL

Mf=NULL
MfT=NULL
MW=NULL
MWT=NULL

MPsi=NULL

Mv=NULL
Mb=NULL
MbT=NULL


iter=410000
bar=400000


#############################
## Programa principal
############################3
for(j in 1:iter){

  if(j<=bar){

    SIGMA=gSigma(b,v,sites)
    DELTA=gCorr(b,sites)

    temp=amostrarWMUSAsa(alpha,W,delta,theta,data,Tt,0.002726011,sites,X,Psi,b,v,SU5)
    W=as.matrix(temp[[1]])
    MWT=c(MWT,temp[[2]])

    temp=amostraralphaMUSAsa(alpha,W,Tt,theta,delta,c2,d2,data,0.002726011,SU2)
    alpha=temp[[1]]
    MalphaT=c(MalphaT,temp[[2]])

    temp=amostrardeltaMUSAsa(theta,delta,W,alpha,data,nj,Tt,0.01,0.002726011,20)
    delta=temp[[1]]
    MdeltaT=c(MdeltaT,temp[[2]])

    temp=amostrarthetaMUSAsa(theta,delta,W,alpha,data,nj,Tt,0.01,0.002726011)
    theta=temp[[1]]
    MthetaT=c(MthetaT,temp[[2]])

    temp=amostrarfMUSAsa(theta,delta,W,alpha,data,nj,Tt,0.00001/2,1/(365+10),1/(365-10),f)
    f=temp[[1]]
    MfT=c(MfT,temp[[2]])

    temp=amostrarb(W,v,b,sites,c3,d3,X,Psi,SU3)
    b=temp[[1]]
    MbT=c(MbT,temp[[2]])


    aa1=0.5*ncol(data)+c4
    bb1=0.5*t(W-X%*%Psi)%*%solve(DELTA)%*%(W-X%*%Psi)+d4
    v=1/rgamma(1,aa1,bb1)


    AA=solve(solve(V)+t(X)%*%solve(SIGMA)%*%X)%*%( solve(V)%*%M+t(X)%*%solve(SIGMA)%*%W)
    BB=solve(solve(V)+t(X)%*%solve(SIGMA)%*%X)

    Psi=as.matrix(MASS::mvrnorm(1,AA,BB))


    SU2=sintonizar(bar,0.30,SU2,MalphaT,j)
    SU3=sintonizar(bar,0.44,SU3,MbT,j)
    SU5=sintonizarN(bar,0.25,SU5,MWT,j)
    print(j)


  }else{

    SIGMA=gSigma(b,v,sites)
    DELTA=gCorr(b,sites)

    temp=amostrarWMUSAsa(alpha,W,delta,theta,data,Tt,0.002726011,sites,X,Psi,b,v,SU5)
    W=as.matrix(temp[[1]])
    MW=rbind(MW,t(W))
    MWT=c(MWT,temp[[2]])

    temp=amostraralphaMUSAsa(alpha,W,Tt,theta,delta,c2,d2,data,0.002726011,SU2)
    alpha=temp[[1]]
    Malpha=c(Malpha,alpha)
    MalphaT=c(MalphaT,temp[[2]])

    temp=amostrardeltaMUSAsa(theta,delta,W,alpha,data,nj,Tt,0.01,0.002726011,20)
    delta=temp[[1]]
    MdeltaT=c(MalphaT,temp[[2]])
    Mdelta=c(Mdelta,delta)

    temp=amostrarthetaMUSAsa(theta,delta,W,alpha,data,nj,Tt,0.01,0.002726011)
    theta=temp[[1]]
    MthetaT=c(MthetaT,temp[[2]])
    Mtheta=c(Mtheta,theta)

    temp=amostrarfMUSAsa(theta,delta,W,alpha,data,nj,Tt,0.00001/2,1/(365+10),1/(365-10),f)
    f=temp[[1]]
    Mf=c(Mf,f)
    MfT=c(MfT,temp[[2]])


    temp=amostrarb(W,v,b,sites,c3,d3,X,Psi,SU3)
    b=temp[[1]]
    Mb=c(Mb,b)
    MbT=c(MbT,temp[[2]])

    aa1=0.5*ncol(data)+c4
    bb1=0.5*t(W-X%*%Psi)%*%solve(DELTA)%*%(W-X%*%Psi)+d4
    v=1/rgamma(1,aa1,bb1)
    Mv=c(Mv,v)


    AA=solve(solve(V)+t(X)%*%solve(SIGMA)%*%X)%*%( solve(V)%*%M+t(X)%*%solve(SIGMA)%*%W)
    BB=solve(solve(V)+t(X)%*%solve(SIGMA)%*%X)

    Psi=MASS::mvrnorm(1,AA,BB)
    MPsi=rbind(MPsi,t(Psi))




    print(j)
  }

}

}

