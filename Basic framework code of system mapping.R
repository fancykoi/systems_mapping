
requiredPackages = c("ggplot2","reshape2","RColorBrewer",'patchwork','orthopolynom','cowplot')
for(packages in requiredPackages){
  if(!require(packages,character.only = TRUE)) install.packages(packages)
  require(packages,character.only = TRUE)
}

library("mclust")
library("deSolve")
library("mvtnorm")

# LSE method
s.mle <- function(s.par,s.y,s.t,x0,y0){
  A <- sum((s.y - com.get_mu(s.par,s.t,x0,y0))^2 )
  A
}

# MLE method
curve.mlefunc<-function( par, y, time.std,x0,y0)
{
  len.cov <- 5;
  par.covar <- par[1:len.cov];
  sig <- SAD3.get_mat(par.covar, times=1:length(time.std),traits = 2);
  
  m  <- length(y[1,]);
  n  <- length(y[,1]);
  
  par1 <- par[-c(1:len.cov)]
  mu0 <- com.get_mu(par1,time.std,x0,y0)
  
  fy0<- dmvnorm(y,mu0,sig)
  A <- -sum(log(fy0));
  
  return (A);
}

# Use SAD(2) model to construct covariance matrix
SAD3.get_mat <- function (par0, times, traits = 1, options = list()) {
  
  par <- par0
  if (class(par0) == "list") 
    par <- unlist(par0)
  t_len <- length(times)
  SAD.3 <- array(0, dim = c(t_len * traits, t_len * traits))
  for (i0 in 1:traits) for (i1 in 1:traits) {
    if (i0 == i1) 
      for (k0 in 1:t_len) for (k1 in k0:t_len) {
        if(k0==k1){
          SAD.3[(i0 - 1) * t_len + k0, (i1 - 1) * t_len + 
                  k1] <- abs(par[i0 * 2])^2 * ((1-par[i0 * 2 - 1]^(2*k0))/(1-par[i0 * 2 - 1]^2))
        }
        if(k0<k1){
          SAD.3[(i0 - 1) * t_len + k0, (i1 - 1) * t_len + 
                  k1] <- abs(par[i0 * 2])^2 *par[i0 * 2 - 1]^(k1 - k0)*sqrt((1-par[i0 * 2 - 1]^(2*k0))/(1-par[i0 * 2 - 1]^(2*k1)))
          SAD.3[(i0 - 1) * t_len + k1, (i1 - 1) * t_len + 
                  k0] <- abs(par[i0 * 2])^2 *par[i0 * 2 - 1]^(k1 - k0)* sqrt((1-par[i0 * 2 - 1]^(2*k0))/(1-par[i0 * 2 - 1]^(2*k1)))
        }
      }
    if (i0<i1)
      for (k0 in 1:t_len) for (k1 in k0:t_len) {
        if(k0==k1){
          SAD.3[(i1-1) * t_len + k0, (i0-1 ) * t_len + 
                  k1] <- par[2*traits+1]*par[i0*2]*par[i1*2]
          SAD.3[(i0-1) * t_len + k0, (i1-1) * t_len + 
                  k1] <- par[2*traits+1]*par[i0*2]*par[i1*2]
        }
        if(k0<k1){
          SAD.3[(i1-1) * t_len + k0, (i0-1 ) * t_len + 
                  k1] <- par[i0*2]*par[i1*2]*((par[i1*2-1]^(k1-k0)-par[i0*2-1]^k0*par[i1*2-1]^k1)/(1-par[i0*2-1]*par[i1*2-1]))/sqrt(((1-par[i0*2-1]^(2*k0))/(1-par[i0*2-1]^2))*
                                                                                                                                      ((1-par[i1*2-1]^(2*k1))/(1-par[i1*2-1]^2)))*par[2*traits+1]
          SAD.3[(i1-1) * t_len + k1, (i0-1 ) * t_len + 
                  k0] <- par[i0*2]*par[i1*2]*((par[i1*2-1]^(k1-k0)-par[i0*2-1]^k0*par[i1*2-1]^k1)/(1-par[i0*2-1]*par[i1*2-1]))/sqrt(((1-par[i0*2-1]^(2*k0))/(1-par[i0*2-1]^2))*
                                                                                                                                      ((1-par[i1*2-1]^(2*k1))/(1-par[i1*2-1]^2)))*par[2*traits+1]
          SAD.3[(i0-1) * t_len + k0, (i1-1 ) * t_len + 
                  k1] <- par[i0*2]*par[i1*2]*((par[i1*2-1]^(k1-k0)-par[i0*2-1]^k0*par[i1*2-1]^k1)/(1-par[i0*2-1]*par[i1*2-1]))/sqrt(((1-par[i0*2-1]^(2*k0))/(1-par[i0*2-1]^2))*
                                                                                                                                      ((1-par[i1*2-1]^(2*k1))/(1-par[i1*2-1]^2)))*par[2*traits+1]
          SAD.3[(i0-1) * t_len + k1, (i1-1 ) * t_len + 
                  k0] <- par[i0*2]*par[i1*2]*((par[i1*2-1]^(k1-k0)-par[i0*2-1]^k0*par[i1*2-1]^k1)/(1-par[i0*2-1]*par[i1*2-1]))/sqrt(((1-par[i0*2-1]^(2*k0))/(1-par[i0*2-1]^2))*
                                                                                                                                      ((1-par[i1*2-1]^(2*k1))/(1-par[i1*2-1]^2)))*par[2*traits+1]
        }
      }
  }
  return(SAD.3)
}

# RK4 solves the LV equation
com.get_mu <- function(par, times, x0,y0)
{
  par0 <- par;
  if (class(par0)!="list")
  {
    par0 <- list(
      r1 = par[1],
      a12 = par[2],
      k1 = par[3],
      r2 = par[4],
      a21 = par[5],
      k2 = par[6]);
  }
  
  state0 <- c(X=x0, Y=y0);
  y <- COMP.f( par0, state0, times );
  
  #if (class(y)=="try-error" )
  #  return(rep(NaN, length(times)*2));
  index.time <- 1:length(times)
  return ( c(y[index.time, 2:3] ) );
}

COMP.f <-function( parameters, state, times ){
  Lorenz<-function(t, state, parameters) 
  {
    with( as.list(c(state, parameters)),
          {
            dX <- r1*X*(1-((X+a12*Y)/k1))
            dY <- r2*Y*(1-((a21*X+Y)/k2))
            
            list(c(dX, dY))
          }
          
    ) # end with(as.list ...
  }
  
  out <- try(rk(y = state, times = times, func = Lorenz, parms = parameters,method="rk4") );
  out;
}

# main function of H0 hypothesis test
curve.get_est_param<-function(dat1,pheno1,pheno2)
{
  phenos <- cbind(pheno1,pheno2)
  ym <- as.numeric(colMeans(phenos))
  value <- c()
  allpar <- list()
  for(ii in 1:10){
    s.par <- c(0.1585815,4.5978148,125.2705139,0.1233422,2.8023220,91.2431834)*runif(6,0.85,1.15)
    res <- optim(s.par,s.mle,s.y=ym,s.t=dat1$sample_times,x0=ym[1],y0=ym[15],
                 method="BFGS",control=list(trace=T,maxit=10000))
    value <- c(value,res$value)
    allpar[[ii]] <- res$par
  }
  p1 <- allpar[[which(value==min(value))]]
  parin1 <- c(runif(1),sd(as.matrix(pheno1)),runif(1),sd(as.matrix(pheno2)),runif(1))
  parin <- c(parin1,p1)
  r0<- optim( parin, curve.mlefunc, y = phenos, time.std=dat1$sample_times,x0=ym[1],y0=ym[15],
              method ="BFGS",control=list(maxit=2000,trace=T))
  
  cat("Estimated parameters:", "log(L)=", r0$value, "PAR=", r0$par, "\n");
  return(c(r0$par,r0$value));
}

# System mapping to obtain H0 parameters
H0 <- curve.get_est_param(dat1,pheno1=dat1$ep.p,pheno2=dat1$sp.p)



