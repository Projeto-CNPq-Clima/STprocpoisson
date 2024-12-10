STModelGoelMCMC<-function(data,sites,X=cbind(as.matrix(rep(1,ncol(data))),as.matrix(sites)),Z=X,M=Z,
                   prior=list(c3=(-2*log(0.05)/max(dist(sites)))*0.1,
                              d3=0.1,
                              BB1=diag(100,ncol(Z)),
                              AA1=as.matrix(rep(0,ncol(Z))),
                              aa1=2.01,
                              bb1=1.005,
                              V=diag(100,ncol(X)),
                              MM1=as.matrix(rep(0,ncol(X))), SU1=0.01,
                              SU2=0.01,
                              SU3=28.01968,
                              SU5=0.001118574,
                              SU6=0.0553372),
                   iter=410000,bar=400000){


  pp=ncol(X)
  ##################################
  # Valores iniciais
  ####################################
  b=1.4
  v=2.9
  Psi=as.matrix(rep(0,ncol(X)))
  W=as.matrix(rep(0,ncol(data)))
  lgama=as.matrix(rep(0,ncol(Z)))
  leta=as.matrix(rep(0,ncol(Z)))
  # Hiperparametros
  ####################################3
  c3=prior$c3
  d3=prior$d3
  BB1=prior$BB1
  AA1=prior$AA1
  aa1=prior$aa1
  bb1=prior$bb1
  V=prior$V
  MM1=prior$MM1
  #########################################
  # Parametros computacionais
  #########################################
  SU1=prior$SU1
  SU2=prior$SU2
  SU3=prior$SU3
  SU5=prior$SU5
  SU6=prior$SU6

  n=ncol(data)
  m=nrow(data)
  tempdados=is.na(data)
  nj=m-apply(tempdados,2,sum)
  p=ncol(X)

  Tt=array(NA,dim=c(1,n))
  NN=array(NA,dim=c(1,n))
  iter=410000
  bar=400000
  pul=1
  #########################################
  # Parametros computacionais
  #########################################

  Mgama=NULL
  MgamaT=NULL

  Meta=NULL
  MetaT=NULL
  Mv=NULL
  Mb=NULL
  MbT=NULL
  MW=NULL
  MWT=NULL
  MPsi=NULL
  MPsiT=NULL

  for(y in 1:n){
    Tt[1,y]=data[nj[y],y]
  }

  Tt=t(Tt)

  #############################
  ## Programa principal
  ############################3
  for(j in 1:iter){

    if(j<=bar){

      theta=exp(W)
      beta=exp(Z%*%lgama)
      alpha=exp(M%*%leta)
      for(e in 1:n){
        t1<-data[nj[e],e]
        NN[1,e]<-rpois(1,theta[e,1]*(1-pgamma(beta[e,1]*t1^alpha[e,1],1,1)))
      }

      temp=amostrargamaGOEL(lgama,leta,data,Z,M,NN,Tt,nj,AA1,BB1,SU1)
      lgama=temp[[1]]
      MgamaT=c(MgamaT,temp[[2]])

      temp=amostraretaGOEL(lgama,leta,data,Z,M,NN,Tt,nj,AA1,BB1,SU2)

      leta=temp[[1]]
      MetaT=c(MetaT,temp[[2]])

      RR=gCorr(b,sites)
      aa=(n/2)+aa1
      bb=0.5*t(W-X%*%Psi)%*%solve(RR)%*%(W-X%*%Psi)+bb1

      v=1/rgamma(1,shape=aa, rate = bb)

      temp=amostrarb(W,v,b,sites,c3,d3,X,Psi,SU3)
      b=temp[[1]]
      MbT=c(MbT,temp[[2]])


      temp=amostrarWgoel(W,sites,X,Psi,b,v,nj,NN,SU5)
      W=as.matrix(temp[[1]])
      MWT=c(MWT,temp[[2]])


      temp=amostrarPsiGOEL(W,X,Psi,V,MM1,SU6,sites,b,v)

      Psi=as.matrix(temp[[1]])
      MPsiT=rbind(MPsiT,temp[[2]])

      if((j%%50)==0){
        SU5=sintonizarN(bar,0.30,SU5,MWT,j)
        SU1=sintonizarN(bar,0.25,SU1,MgamaT,j)
        SU2=sintonizarN(bar,0.25,SU2,MetaT,j)
        SU3=sintonizar(bar,0.44,SU3,MbT,j)
        SU6=sintonizarN(bar,0.25,SU6,MPsiT,j)

      }else{

      }

      print(j)


    }else{




      theta=exp(W)
      for(e in 1:n){
        t1<-data[nj[e],e]
        NN[1,e]<-rpois(1,theta[e,1]*(1-pgamma(beta*t1^alpha,1,1)))
      }


      temp=amostrargamaGOEL(lgama,leta,data,Z,M,NN,Tt,nj,AA1,BB1,SU1)

      lgama=temp[[1]]
      Mgama=rbind(Mgama,t(lgama))
      MgamaT=c(MgamaT,temp[[2]])


      temp=amostraretaGOEL(lgama,leta,data,Z,M,NN,Tt,nj,AA1,BB1,SU2)

      leta=temp[[1]]
      Meta=rbind(Meta,t(leta))
      MetaT=c(MetaT,temp[[2]])

      RR=gCorr(b,sites)
      aa=(n/2)+aa1
      bb=0.5*t(W-X%*%Psi)%*%solve(RR)%*%(W-X%*%Psi)+bb1

      v=1/rgamma(1,shape=aa, rate = bb)
      Mv=c(Mv,v)

      temp=amostrarb(W,v,b,sites,c3,d3,X,Psi,SU3)
      b=temp[[1]]
      Mb=c(Mb,b)
      MbT=c(MbT,temp[[2]])

      temp=amostrarWgoel(W,sites,X,Psi,b,v,nj,NN,SU5)
      W=as.matrix(temp[[1]])
      MW=rbind(MW,t(W))
      MWT=c(MWT,temp[[2]])


      temp=amostrarPsiGOEL(W,X,Psi,V,MM1,SU6,sites,b,v)

      Psi=as.matrix(temp[[1]])
      MPsiT=c(MPsiT,temp[[2]])
      MPsi=rbind(MPsi,t(Psi))


      print(j)
    }

  }

  return(list(Mgama,MgamaT,Meta,MetaT,Mv,Mb,MbT,MW,MWT,MPsi,MPsiT))

}
