
######################################  Normal Topic Model (NTM) code  ##########################

reverse_transform=function(x) 
{
  out=array(0,length(x)-1);
  for(i in 1:(length(x)-1))
  {
    out[i]=log((x[i]+10^(-5))/(x[length(x)]+10^(-5)));
  }
  return(out)
}

# Data


# Log-Likelihood
loglik_norm <- function(u,y,x) sum((y - x %*% u)^2)

# Transform the parameters: we just have
# to find a bijection between R^3 
# and {(a,b,c,d) \in [0,1]^4 : a+b+c+d=1}.

transform <- function(v) 
{
  # Ensure they in [0,1]
  temp =c(exp(v),1);
  out=temp/(1+sum(exp(v)));
  return(out)
}

library(gtools)
library(SQUAREM)

NTMfit_local <- function(p, y, N, G, K)
{
  data=matrix(y,N,G);
  omega_in=matrix(p[1:(N*K)],N,K);
  alpha_in=matrix(p[-(1:(N*K))],K,G);
  
  ###########   Estimating the effect size alpha  ########################
  
  svd_omega=svd(omega_in);
  temp1=t(svd_omega$v)%*%diag(1/svd_omega$d^2,dim(omega_in)[2])%*%svd_omega$v;
  temp2=t(omega_in)%*%data;
  temp1=solve(t(omega_in)%*%omega_in);
  alpha = temp1%*%temp2;
  
  ###########  Estimating the topic proportions ########################
  
  
  omega=matrix(0,dim(data)[1],K);
  for(n in 1:dim(data)[1])
  {
    omega_vec=omega_in[n,];
    data_vec=data[n,];
    res=optim(reverse_transform(omega_vec), function(v) loglik_norm(transform(v),data_vec,t(alpha)) );
    omega[n,]=transform(res$par);
  }
  
  param_vec_omega=matrix(omega,1,N*K);
  param_vec_alpha=matrix(alpha,1,K*G);
  param_vec_out=c(param_vec_omega,param_vec_alpha);
  return(param_vec_out)
}

NTMfit <- function(data, K, omega_start=NULL, alpha_start=NULL, niter=50, scale=1)
{
  N <- dim(data)[1]; G <- dim(data)[2];
  data_vector=matrix(data,1,N*G);
  if(is.null(omega_start))
    omega_start <- matrix(rdirichlet(N,c(scale/K,scale/K,scale/K,scale/K)), nrow=N);
  if(is.null(alpha_start))
    alpha_start <- matrix(rnorm((K)*G,1,1),nrow=(K));
  
  param_start=c(matrix(omega_start,1,N*K),matrix(alpha_start,1,K*G));
  options(warn=-1)
  system.time(res <- squarem(p=param_start,y=data_vector, N=N, G=G, K=K, fixptfn=NTMfit_local, control=list(maxiter = niter, trace = FALSE)));
  
  omega_out=matrix(res$par[1:(N*K)],N,K);
  alpha_out=matrix(res$par[-(1:(N*K))],K,G);
  
  ll <- list("omega_out"=omega_out, "alpha_out"=alpha_out);
  return(ll)
}