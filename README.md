# R-code used for On Latent Trait Model with Bayesian Marginal Likelihood of Rank-Based Estimation
library("LaplacesDemon")
set.seed(12345)
n=100#sample size
groups=11
epsilon=matrix(NA,n,1)
#Create a sample of "time" observations on the independent variables x
beta=matrix(c(-0.9,0.4, 0.6, 0.2, 0.1),ncol=1)
#x_mat=matrix(NA,n,nrow(beta))
#I need to create 11 design matrices with the covariates
x_m=list()
for(i in 1:groups){
  x_m[[i]]=matrix(NA,nrow=n,ncol=length(beta))
  for(j in 1:length(beta)){
    x_m[[i]][,1]=rnorm(n,15,0.1)
    x_m[[i]][,2]=rpois(n, 2)
    x_m[[i]][,3]=rnorm(n,23, 0.2)
    x_m[[i]][,4]=rpois(n, 2)
    x_m[[i]][,5]=rnorm(n, 16, 0.2)
  }
}

for(i in 1:n){
  epsilon[i,1]=rnorm(1,0,1)
}

sigma2_random_effects=random_effects=c()
latents=latents_cat=matrix(NA,n,groups)
sigma2_random_effects=0.83
lamda2=0.5
omega=5
for(j in 1:groups){
  #sigma2_random_effects[j]=rinvchisq(1, df=5, scale=0.5)
  #lamda^2=0.5, omega=5 therefore mean is omega*lamda^2/(omega-2)=(5*0.5)/3=0.83
  random_effects[j]=rnorm(1,0,sigma2_random_effects)
  latents[,j]=x_m[[j]]%*%beta+random_effects[j]+epsilon
}
#They need to be integer categorical data

latents_cat=ceiling(latents)
max(latents_cat)
min(latents_cat)
for(i in 1:nrow(latents_cat)){
  for(j in 1:ncol(latents_cat)){
    if(latents_cat[i,j]>=(-1)&latents_cat[i,j]<(2)){
      latents_cat[i,j]=((-1)+(2))/2
    }
    else if(latents_cat[i,j]>=(2)&latents_cat[i,j]<(5)){
      latents_cat[i,j]=((0)+(3))/2
    }
    else if(latents_cat[i,j]>=(3)&latents_cat[i,j]<(6)){
      latents_cat[i,j]=((3)+(6))/2
    }
    else if(latents_cat[i,j]>=(6)&latents_cat[i,j]<(9)){
      latents_cat[i,j]=((6)+(9))/2
    }
    
  }
}

#ceiling(latent_cat)
#length(unique(latent_cat))
plot(latents_cat[,1],type="l")
m=11
n_obs=length(latents_cat[,1])
Y=matrix(cbind(latents_cat[,1],latents_cat[,2],latents_cat[,3],latents_cat[,4],
               latents_cat[,5],latents_cat[,6],latents_cat[,7],latents_cat[,8],
               latents_cat[,9],latents_cat[,10],latents_cat[,11]),
         nrow=n_obs)

n_obs=n_obs-1
Y=Y[1:n_obs,]
Y_cat=ceiling(Y)
View(X_11)
X_1=x_m[[1]]
X_2=x_m[[2]]
X_3=x_m[[3]]
X_4=x_m[[4]]
X_5=x_m[[5]]
X_6=x_m[[6]]
X_7=x_m[[7]]
X_8=x_m[[8]]
X_9=x_m[[9]]
X_10=x_m[[10]]
X_11=x_m[[11]]

#I need to rescale the covariates
X_1=scale(X_1,center=TRUE,scale=TRUE)
X_2=scale(X_2,center=TRUE,scale=TRUE)
X_3=scale(X_3,center=TRUE,scale=TRUE)
X_4=scale(X_4,center=TRUE,scale=TRUE)
X_5=scale(X_5,center=TRUE,scale=TRUE)
X_6=scale(X_6,center=TRUE,scale=TRUE)
X_7=scale(X_7,center=TRUE,scale=TRUE)
X_8=scale(X_8,center=TRUE,scale=TRUE)
X_9=scale(X_9,center=TRUE,scale=TRUE)
X_10=scale(X_10,center=TRUE,scale=TRUE)
X_11=scale(X_11,center=TRUE,scale=TRUE)
#Out-of-sample. Without the last observation
m=ncol(Y)
set.seed(12345)


X=array(cbind(X_1[1:n_obs,],X_2[1:n_obs,],
              X_3[1:n_obs,],X_4[1:n_obs,],
              X_5[1:n_obs,],X_6[1:n_obs,],X_7[1:n_obs,],X_8[1:n_obs,],X_9[1:n_obs,],
              X_10[1:n_obs,],X_11[1:n_obs,]),c(n_obs,ncol(X_1),m))

X_init=array(cbind(X_1,X_2,
                   X_3,X_4,
                   X_5,X_6,X_7,X_8,X_9,X_10,X_11),c((n_obs+1),ncol(X_1),m))



#prior hyperparameters
lamda2_prior=0.2
#prior degrees of belief to update sigma square of upsilon (random effects) 
omega_upsilon_prior=17
#Initial values for BETA
BETA=matrix(0,ncol(X_1),1)
#Initial values for z
z=matrix(0.5,nrow=n_obs,ncol=m)
a=b=matrix(NA,nrow=n_obs,ncol=m)
#I need this package for scaled inverse chi-squared
library("geoR")
###MCMC
library("mvtnorm")
library(MCMCpack)
p=ncol(X[,,1])
omega_upsilon_posterior=NULL #posterior degrees of freedom of variance of random effect
lamda2_posterior=NULL #posterior scale parameter of sigma of variance of random effect
mean_upsilon_posterior=var_upsilon_posterior=NULL#post values of mean and
#variance of random effects
n_iter=50000
library("msm")

BETA_POST_5=upsilon_POST_5=sigmasq_upsilon_POST_5=NULL
upsilon_posterior=rep(0.1,m)
sigmasq_upsilon_post=NULL#posterior value of sigmasq of upsilon
z_all_5=NULL
#GIBBS SAMPLER
for(i in 1:n_iter){
  for(j in 1:m){
    omega_upsilon_posterior[j]=(omega_upsilon_prior+1) #updating degrees of freedom of variance of random effect
    
    lamda2_posterior[j]=(lamda2_prior+
                           ((upsilon_posterior[j])^2)/omega_upsilon_posterior[j]) #update scale parameter of variance of random effect
    
    sigmasq_upsilon_post[j]=rinvchisq(1,omega_upsilon_posterior[j],
                                      lamda2_posterior[j]) #variance of random effects
    
    
  } 
  
  
  #given variance I can sample from gamma (random effects)
  mean_upsilon_posterior=array(NA,c(n_obs,1,m))
  
  for(t in 1:n_obs){
    for(j in 1:m){
      mean_upsilon_posterior[t,,j]=(z[t,j]-X[t,,j]%*%BETA)
    }
  }  
  for(j in 1:m){
    var_upsilon_posterior[j]=1/(n_obs+(1/sigmasq_upsilon_post[j]))
  }
  #summing over time:
  mean_upsilon_posterior_sum=apply(mean_upsilon_posterior,3,sum)
  
  for(j in 1:m){
    mean_upsilon_posterior_sum[j]=mean_upsilon_posterior_sum[j]/(n_obs+1/sigmasq_upsilon_post[j])
    upsilon_posterior[j]=rnorm(1,mean_upsilon_posterior_sum[j],
                               sqrt(var_upsilon_posterior[j]))
  } #updating random effects
  
  upsilon_rep=array(NA,c(n_obs,1,m)) 
  #As I have only one random effect gamma for each coutry
  for(j in 1:m){
    upsilon_rep[,,j]=rep(upsilon_posterior[j],n_obs)
  }
  
  for(t in 1:n_obs){
    for(j in 1:m){
      zt=z[t,]
      a[t,j]=max(zt[Y[t,]<Y[t,j]])
      b[t,j]=min(zt[Y[t,j]<Y[t,]]) 
      z[t,j]=rtnorm(1,X[t,,j]%*%BETA+upsilon_rep[t,1,j],1,a[t,j],b[t,j])#latent variable
    }
  }
  muj1=array(NA,c(p,p,m))
  muj1_sum=matrix(NA,p,p)
  muj2=array(NA,c(p,1,m))
  muj2_sum=matrix(NA,p,1)
  muj_final=matrix(NA,p,1)
  varj_final=matrix(NA,p,p)
  for(j in 1:m){
    muj1[,,j]=t(X[,,j])%*%X[,,j]
    #before summing over j
    muj1_sum=apply(muj1,c(1:2),sum)#summing over j
    muj2[,,j]=t(X[,,j])%*%(z[,j]-upsilon_rep[,1,j])
    muj2_sum=apply(muj2,c(1:2),sum)#summing over j
    muj_final=solve(muj1_sum)%*%(muj2_sum)
    varj_final=solve(muj1_sum)
  }
  BETA=mvrnorm(n=1,mu=muj_final,
               varj_final)#sampling from a multivariate Normal
  #with updated mean and variance
  
  BETA_POST_5=rbind(BETA_POST_5,BETA)
  upsilon_POST_5=rbind(upsilon_POST_5,upsilon_posterior)
  sigmasq_upsilon_POST_5=rbind(sigmasq_upsilon_POST_5,sigmasq_upsilon_post)
  
  print(paste("iteration",i))
  z_all_5=rbind(z_all_5,z)
}
#Sensitivity##for BETA###############
install.packages("viridis")
library(bayesplot)
library(viridis)
beta1=cbind(BETA_POST_1[1:50000,1],BETA_POST_2[1:50000,1],BETA_POST_3[1:50000,1],BETA_POST_4[1:50000,1],BETA_POST_5[1:50000,1])
beta2=cbind(BETA_POST_1[1:50000,2],BETA_POST_2[1:50000,2],BETA_POST_3[1:50000,2],BETA_POST_4[1:50000,2],BETA_POST_5[1:50000,2])
beta3=cbind(BETA_POST_1[1:50000,3],BETA_POST_2[1:50000,3],BETA_POST_3[1:50000,3],BETA_POST_4[1:50000,3],BETA_POST_5[1:50000,3])
beta4=cbind(BETA_POST_1[1:50000,4],BETA_POST_2[1:50000,4],BETA_POST_3[1:50000,4],BETA_POST_4[1:50000,4],BETA_POST_5[1:50000,4])
beta5=cbind(BETA_POST_1[1:50000,5],BETA_POST_2[1:50000,5],BETA_POST_3[1:50000,5],BETA_POST_4[1:50000,5],BETA_POST_5[1:50000,5])
chains_beta=array(cbind(beta1[1:50000,],beta2[1:50000,],beta3[1:50000,],beta4[1:50000,],beta5[1:50000,]),c(50000,5,5))
new_colnames <- c("chain1", "chain2","chain3","chain4","chain5")
dimnames(chains_beta)[[2]]=new_colnames
dimnames(chains_beta)[[3]] <- c("\u03B21", "\u03B22","\u03B23","\u03B24","\u03B25")
color_scheme_set("mix-green-red")
color_scheme_set("viridis")
mcmc_trace(chains_beta, pars = c("\u03B21", "\u03B22","\u03B23","\u03B24","\u03B25"), facet_args = list(ncol = 1, strip.position = "left"))
upsilon1=cbind(upsilon_POST_1[1:50000,1],upsilon_POST_2[1:50000,1],upsilon_POST_3[1:50000,1],upsilon_POST_4[1:50000,1],upsilon_POST_5[1:50000,1])
upsilon2=cbind(upsilon_POST_1[1:50000,2],upsilon_POST_2[1:50000,2],upsilon_POST_3[1:50000,2],upsilon_POST_4[1:50000,2],upsilon_POST_5[1:50000,2])
upsilon3=cbind(upsilon_POST_1[1:50000,3],upsilon_POST_2[1:50000,3],upsilon_POST_3[1:50000,3],upsilon_POST_4[1:50000,3],upsilon_POST_5[1:50000,3])
upsilon4=cbind(upsilon_POST_1[1:50000,4],upsilon_POST_2[1:50000,4],upsilon_POST_3[1:50000,4],upsilon_POST_4[1:50000,4],upsilon_POST_5[1:50000,4])
upsilon5=cbind(upsilon_POST_1[1:50000,5],upsilon_POST_2[1:50000,5],upsilon_POST_3[1:50000,5],upsilon_POST_4[1:50000,5],upsilon_POST_5[1:50000,5])
upsilon6=cbind(upsilon_POST_1[1:50000,6],upsilon_POST_2[1:50000,6],upsilon_POST_3[1:50000,6],upsilon_POST_4[1:50000,6],upsilon_POST_5[1:50000,6])
upsilon7=cbind(upsilon_POST_1[1:50000,7],upsilon_POST_2[1:50000,7],upsilon_POST_3[1:50000,7],upsilon_POST_4[1:50000,7],upsilon_POST_5[1:50000,7])
upsilon8=cbind(upsilon_POST_1[1:50000,8],upsilon_POST_2[1:50000,8],upsilon_POST_3[1:50000,8],upsilon_POST_4[1:50000,8],upsilon_POST_5[1:50000,8])
upsilon9=cbind(upsilon_POST_1[1:50000,9],upsilon_POST_2[1:50000,9],upsilon_POST_3[1:50000,9],upsilon_POST_4[1:50000,9],upsilon_POST_5[1:50000,9])
upsilon10=cbind(upsilon_POST_1[1:50000,10],upsilon_POST_2[1:50000,10],upsilon_POST_3[1:50000,10],upsilon_POST_4[1:50000,10],upsilon_POST_5[1:50000,10])
upsilon11=cbind(upsilon_POST_1[1:50000,11],upsilon_POST_2[1:50000,11],upsilon_POST_3[1:50000,11],upsilon_POST_4[1:50000,11],upsilon_POST_5[1:50000,11])
chains_upsilons=array(cbind(upsilon1[1:50000,],upsilon2[1:50000,],upsilon3[1:50000,],upsilon4[1:50000,]
                           ,upsilon5[1:50000,],upsilon6[1:50000,],upsilon7[1:50000,],upsilon8[1:50000,],upsilon9[1:50000,],upsilon10[1:50000,],upsilon11[1:50000,]),c(50000,5,11))
new_colnames <- c("chain1", "chain2","chain3","chain4","chain5")
dimnames(chains_upsilons)[[2]]=new_colnames
dimnames(chains_upsilons)[[3]] <- c("\u03C51", "\u03C52", "\u03C53","\u03C54","\u03C55","\u03C56","\u03C57","\u03C58","\u03C59","\u03C510","\u03C511")
mcmc_trace(chains_upsilons, pars = c("\u03C51", "\u03C52", "\u03C53","\u03C54","\u03C55","\u03C56","\u03C57","\u03C58","\u03C59","\u03C510","\u03C511"), facet_args = list(ncol = 2, strip.position = "left"))
sigma1=cbind(sigmasq_upsilon_POST_1[1:50000,1],sigmasq_upsilon_POST_2[1:50000,1],sigmasq_upsilon_POST_3[1:50000,1],sigmasq_upsilon_POST_4[1:50000,1],sigmasq_upsilon_POST_5[1:50000,1])
sigma2=cbind(sigmasq_upsilon_POST_1[1:50000,2],sigmasq_upsilon_POST_2[1:50000,2],sigmasq_upsilon_POST_3[1:50000,2],sigmasq_upsilon_POST_4[1:50000,2],sigmasq_upsilon_POST_5[1:50000,2])
sigma3=cbind(sigmasq_upsilon_POST_1[1:50000,3],sigmasq_upsilon_POST_2[1:50000,3],sigmasq_upsilon_POST_3[1:50000,3],sigmasq_upsilon_POST_4[1:50000,3],sigmasq_upsilon_POST_5[1:50000,3])
sigma4=cbind(sigmasq_upsilon_POST_1[1:50000,4],sigmasq_upsilon_POST_2[1:50000,4],sigmasq_upsilon_POST_3[1:50000,4],sigmasq_upsilon_POST_4[1:50000,4],sigmasq_upsilon_POST_5[1:50000,4])
sigma5=cbind(sigmasq_upsilon_POST_1[1:50000,5],sigmasq_upsilon_POST_2[1:50000,5],sigmasq_upsilon_POST_3[1:50000,5],sigmasq_upsilon_POST_4[1:50000,5],sigmasq_upsilon_POST_5[1:50000,5])
sigma6=cbind(sigmasq_upsilon_POST_1[1:50000,6],sigmasq_upsilon_POST_2[1:50000,6],sigmasq_upsilon_POST_3[1:50000,6],sigmasq_upsilon_POST_4[1:50000,6],sigmasq_upsilon_POST_5[1:50000,6])
sigma7=cbind(sigmasq_upsilon_POST_1[1:50000,7],sigmasq_upsilon_POST_2[1:50000,7],sigmasq_upsilon_POST_3[1:50000,7],sigmasq_upsilon_POST_4[1:50000,7],sigmasq_upsilon_POST_5[1:50000,7])
sigma8=cbind(sigmasq_upsilon_POST_1[1:50000,8],sigmasq_upsilon_POST_2[1:50000,8],sigmasq_upsilon_POST_3[1:50000,8],sigmasq_upsilon_POST_4[1:50000,8],sigmasq_upsilon_POST_5[1:50000,8])
sigma9=cbind(sigmasq_upsilon_POST_1[1:50000,9],sigmasq_upsilon_POST_2[1:50000,9],sigmasq_upsilon_POST_3[1:50000,9],sigmasq_upsilon_POST_4[1:50000,9],sigmasq_upsilon_POST_5[1:50000,9])
sigma10=cbind(sigmasq_upsilon_POST_1[1:50000,10],sigmasq_upsilon_POST_2[1:50000,10],sigmasq_upsilon_POST_3[1:50000,10],sigmasq_upsilon_POST_4[1:50000,10],sigmasq_upsilon_POST_5[1:50000,10])
sigma11=cbind(sigmasq_upsilon_POST_1[1:50000,11],sigmasq_upsilon_POST_2[1:50000,11],sigmasq_upsilon_POST_3[1:50000,11],sigmasq_upsilon_POST_4[1:50000,11],sigmasq_upsilon_POST_5[1:50000,11])
chains_sigmas=array(cbind(sigma1[1:50000,],sigma2[1:50000,],sigma3[1:50000,],sigma4[1:50000,]
                         ,sigma5[1:50000,],sigma6[1:50000,],sigma7[1:50000,],sigma8[1:50000,],sigma9[1:50000,],sigma10[1:50000,],sigma11[1:50000,]),c(50000,5,11))
subscript_1 <- "\u2081"
subscript_2 <- "\u2082"
subscript_3 <- "\u2083"
subscript_4 <- "\u2084"
subscript_5 <- "\u2085"
subscript_6 <- "\u2086"
subscript_7 <- "\u2087"
subscript_8 <- "\u2088"
subscript_9 <- "\u2089"
subscript_10 <- paste("\u2081", "\u2080", sep = "")
subscript_11 <- paste("\u2081", "\u2081", sep = "")
sigma_1_squared <- paste("\u03C3", subscript_1, "²", sep = "")
sigma_2_squared <- paste("\u03C3", subscript_2, "²", sep = "")
sigma_3_squared <- paste("\u03C3", subscript_3, "²", sep = "")
sigma_4_squared <- paste("\u03C3", subscript_4, "²", sep = "")
sigma_5_squared <- paste("\u03C3", subscript_5, "²", sep = "")
sigma_6_squared <- paste("\u03C3", subscript_6, "²", sep = "")
sigma_7_squared <- paste("\u03C3", subscript_7, "²", sep = "")
sigma_8_squared <- paste("\u03C3", subscript_8, "²", sep = "")
sigma_9_squared <- paste("\u03C3", subscript_9, "²", sep = "")
sigma_10_squared <- paste("\u03C3", subscript_10, "\u00B2", sep = "")
sigma_11_squared <- paste("\u03C3", subscript_11, "\u00B2", sep = "")
new_colnames <- c("chain1", "chain2","chain3","chain4","chain5")
dimnames(chains_sigmas)[[2]]=new_colnames
dimnames(chains_sigmas)[[3]] <- c(sigma_1_squared, sigma_2_squared,sigma_3_squared,sigma_4_squared,sigma_5_squared, sigma_6_squared, sigma_7_squared,sigma_8_squared,sigma_9_squared,sigma_10_squared,sigma_11_squared)
mcmc_trace(chains_sigmas, pars = c(sigma_1_squared, sigma_2_squared,sigma_3_squared,sigma_4_squared,sigma_5_squared, sigma_6_squared, sigma_7_squared,sigma_8_squared,sigma_9_squared,sigma_10_squared,sigma_11_squared), facet_args = list(ncol = 2, strip.position = "left"))
a=mean(BETA_POST_5[,1])
b=mean(BETA_POST_5[,2])
c=mean(BETA_POST_5[,3]) 
d=mean(BETA_POST_5[,4])
e=mean(BETA_POST_5[,5])
cbind(a,b,c,d,e)
a=mean(upsilon_POST_5[,1])
b=mean(upsilon_POST_5[,2])
c=mean(upsilon_POST_5[,3]) 
d=mean(upsilon_POST_5[,4])
e=mean(upsilon_POST_5[,5])
f=mean(upsilon_POST_5[,6])
g=mean(upsilon_POST_5[,7])
h=mean(upsilon_POST_5[,8]) 
i=mean(upsilon_POST_5[,9])
j=mean(upsilon_POST_5[,10])
k=mean(upsilon_POST_5[,11])
cbind(a,b,c,d,e,f,g,h,i,j,k)
###########
a=sd(BETA_POST_5[,1])
b=sd(BETA_POST_5[,2])
c=sd(BETA_POST_5[,3]) 
d=sd(BETA_POST_5[,4])
e=sd(BETA_POST_5[,5])
cbind(a,b,c,d,e)
a=sd(upsilon_POST_5[,1])
b=sd(upsilon_POST_5[,2])
c=sd(upsilon_POST_5[,3]) 
d=sd(upsilon_POST_5[,4])
e=sd(upsilon_POST_5[,5])
f=sd(upsilon_POST_5[,6])
g=sd(upsilon_POST_5[,7])
h=sd(upsilon_POST_5[,8]) 
i=sd(upsilon_POST_5[,9])
j=sd(upsilon_POST_5[,10])
k=sd(upsilon_POST_5[,11])
cbind(a,b,c,d,e,f,g,h,i,j,k)
####To compute RMSE & MAE for beta and upsilon##CHAIN 1######
posterior_means_beta_1=matrix(colMeans(BETA_POST_1),5,1)
mse_beta_1=matrix(NA,5,1)
for(j in 1:5){
  mse_beta_1[j]=(posterior_means_beta_1[j,1]-beta[j,1])^2
}
MSE_beta_1=mean(mse_beta_1)
RMSE_beta_1=sqrt(MSE_beta_1)
mae_beta_1=matrix(NA,5,1 )
for(j in 1:5){
  mae_beta_1[j]=abs(posterior_means_beta_1[j,1]-beta[j,])
}
MAE_beta_1=mean(mae_beta_1)
############################for upsilon chain 1####
posterior_means_upsilon_1=matrix(colMeans(upsilon_POST_1),m,1)
mse_upsilon_1=matrix(NA,m,1)
for(j in 1:m){
  mse_upsilon_1[j]=(posterior_means_upsilon_1[j,1]-upsilon[j,1])^2
}
MSE_upsilon_1=mean(mse_upsilon_1)
RMSE_upsilon_1=sqrt(MSE_upsilon_1)
mae_upsilon_1=matrix(NA,m,1 )
for(j in 1:m){
  mae_upsilon_1[j]=abs(posterior_means_upsilon_1[j,1]-upsilon[j,1])
}
MAE_upsilon_1=mean(mae_upsilon_1)

####To compute RMSE & MAE for beta and upsilon##CHAIN 2######
posterior_means_beta_2=matrix(colMeans(BETA_POST_2),5,1)
mse_beta_2=matrix(NA,5,1)
for(j in 1:5){
  mse_beta_2[j]=(posterior_means_beta_2[j,1]-beta[j,1])^2
}
MSE_beta_2=mean(mse_beta_2)
RMSE_beta_2=sqrt(MSE_beta_2)
mae_beta_2=matrix(NA,5,1 )
for(j in 1:5){
  mae_beta_2[j]=abs(posterior_means_beta_2[j,1]-beta[j,])
}
MAE_beta_2=mean(mae_beta_2)
############################for upsilon chain 2####
posterior_means_upsilon_2=matrix(colMeans(upsilon_POST_2),m,1)
mse_upsilon_2=matrix(NA,m,1)
for(j in 1:m){
  mse_upsilon_2[j]=(posterior_means_upsilon_2[j,1]-upsilon[j,1])^2
}
MSE_upsilon_2=mean(mse_upsilon_2)
RMSE_upsilon_2=sqrt(MSE_upsilon_2)
mae_upsilon_2=matrix(NA,m,1 )
for(j in 1:m){
  mae_upsilon_2[j]=abs(posterior_means_upsilon_2[j,1]-upsilon[j,1])
}
MAE_upsilon_2=mean(mae_upsilon_2)
#######################################
####To compute RMSE & MAE for beta and upsilon##CHAIN 3######
posterior_means_beta_3=matrix(colMeans(BETA_POST_3),5,1)
mse_beta_3=matrix(NA,5,1)
for(j in 1:5){
  mse_beta_3[j]=(posterior_means_beta_3[j,1]-beta[j,1])^2
}
MSE_beta_3=mean(mse_beta_3)
RMSE_beta_3=sqrt(MSE_beta_3)
mae_beta_3=matrix(NA,5,1 )
for(j in 1:5){
  mae_beta_3[j]=abs(posterior_means_beta_3[j,1]-beta[j,])
}
MAE_beta_3=mean(mae_beta_3)
############################for upsilon chain 3####
posterior_means_upsilon_3=matrix(colMeans(upsilon_POST_3),m,1)
mse_upsilon_3=matrix(NA,m,1)
for(j in 1:m){
  mse_upsilon_3[j]=(posterior_means_upsilon_3[j,1]-upsilon[j,1])^2
}
MSE_upsilon_3=mean(mse_upsilon_3)
RMSE_upsilon_3=sqrt(MSE_upsilon_3)
mae_upsilon_3=matrix(NA,m,1 )
for(j in 1:m){
  mae_upsilon_3[j]=abs(posterior_means_upsilon_3[j,1]-upsilon[j,1])
}
MAE_upsilon_3=mean(mae_upsilon_3)
################################
####To compute RMSE & MAE for beta and upsilon##CHAIN 4######
posterior_means_beta_4=matrix(colMeans(BETA_POST_4),5,1)
mse_beta_4=matrix(NA,5,1)
for(j in 1:5){
  mse_beta_4[j]=(posterior_means_beta_4[j,1]-beta[j,1])^2
}
MSE_beta_4=mean(mse_beta_4)
RMSE_beta_4=sqrt(MSE_beta_4)
mae_beta_4=matrix(NA,5,1 )
for(j in 1:5){
  mae_beta_4[j]=abs(posterior_means_beta_4[j,1]-beta[j,])
}
MAE_beta_4=mean(mae_beta_4)
############################for upsilon chain 4####
posterior_means_upsilon_4=matrix(colMeans(upsilon_POST_4),m,1)
mse_upsilon_4=matrix(NA,m,1)
for(j in 1:m){
  mse_upsilon_4[j]=(posterior_means_upsilon_4[j,1]-upsilon[j,1])^2
}
MSE_upsilon_4=mean(mse_upsilon_4)
RMSE_upsilon_4=sqrt(MSE_upsilon_4)
mae_upsilon_4=matrix(NA,m,1 )
for(j in 1:m){
  mae_upsilon_4[j]=abs(posterior_means_upsilon_4[j,1]-upsilon[j,1])
}
MAE_upsilon_4=mean(mae_upsilon_4)
###################################
####To compute RMSE & MAE for beta and upsilon##CHAIN 5######
posterior_means_beta_5=matrix(colMeans(BETA_POST_5),5,1)
mse_beta_5=matrix(NA,5,1)
for(j in 1:5){
  mse_beta_5[j]=(posterior_means_beta_5[j,1]-beta[j,1])^2
}
MSE_beta_5=mean(mse_beta_5)
RMSE_beta_5=sqrt(MSE_beta_5)
mae_beta_5=matrix(NA,5,1 )
for(j in 1:5){
  mae_beta_5[j]=abs(posterior_means_beta_5[j,1]-beta[j,])
}
MAE_beta_5=mean(mae_beta_5)
############################for upsilon chain 5####
posterior_means_upsilon_5=matrix(colMeans(upsilon_POST_5),m,1)
mse_upsilon_5=matrix(NA,m,1)
for(j in 1:m){
  mse_upsilon_5[j]=(posterior_means_upsilon_5[j,1]-upsilon[j,1])^2
}
MSE_upsilon_5=mean(mse_upsilon_5)
RMSE_upsilon_5=sqrt(MSE_upsilon_5)
mae_upsilon_5=matrix(NA,m,1 )
for(j in 1:m){
  mae_upsilon_5[j]=abs(posterior_means_upsilon_5[j,1]-upsilon[j,1])
}
MAE_upsilon_5=mean(mae_upsilon_5)
####################To compute RMSE & MAE######################
posterior_means_z_1=colMeans(z_all_1)
posterior_means_z_2=colMeans(z_all_2)
posterior_means_z_3=colMeans(z_all_3)
posterior_means_z_4=colMeans(z_all_4)
posterior_means_z_5=colMeans(z_all_5)
######CHAIN 1 ######
mse_1=matrix(NA,n,1 )
for(j in 1:n){
  mse_1[j]=sum((posterior_means_z_1-latents[j,])^2)/m
}
MSE_1=sum(mse_1)/n
RMSE_1=sqrt(MSE_1)
mae_1=matrix(NA,n,1 )
for(j in 1:n){
  mae_1[j]=sum(abs(posterior_means_z_1-latents[j,]))/m
}
MAE_1=sum(mae_1)/n
##########CHAIN 2#####
mse_2=matrix(NA,n,1 )
for(j in 1:n){
  mse_2[j]=sum((posterior_means_z_2-latents[j,])^2)/m
}
MSE_2=sum(mse_2)/n
RMSE_2=sqrt(MSE_2)
mae_2=matrix(NA,n,1 )
for(j in 1:n){
  mae_2[j]=sum(abs(posterior_means_z_2-latents[j,]))/m
}
MAE_2=sum(mae_2)/n
#######CHAIN 3 #####
mse_3=matrix(NA,n,1 )
for(j in 1:n){
  mse_3[j]=sum((posterior_means_z_3-latents[j,])^2)/m
}
MSE_3=sum(mse_3)/n
RMSE_3=sqrt(MSE_3)
mae_3=matrix(NA,n,1 )
for(j in 1:n){
  mae_3[j]=sum(abs(posterior_means_z_3-latents[j,]))/m
}
MAE_3=sum(mae_3)/n
#####CHAIN 4#######
mse_4=matrix(NA,n,1 )
for(j in 1:n){
  mse_4[j]=sum((posterior_means_z_4-latents[j,])^2)/m
}
MSE_4=sum(mse_4)/n
RMSE_4=sqrt(MSE_4)
mae_4=matrix(NA,n,1 )
for(j in 1:n){
  mae_4[j]=sum(abs(posterior_means_z_4-latents[j,]))/m
}
MAE_4=sum(mae_4)/n
####CHAIN 5 ########
mse_5=matrix(NA,n,1 )
for(j in 1:n){
  mse_5[j]=sum((posterior_means_z_5-latents[j,])^2)/m
}
MSE_5=sum(mse_5)/n
RMSE_5=sqrt(MSE_5)
mae_5=matrix(NA,n,1 )
for(j in 1:n){
  mae_5[j]=sum(abs(posterior_means_z_5-latents[j,]))/m
}
MAE_5=sum(mae_5)/n
########## psrf_test/Gelman-rubin convergancy test ##########
library(stableGR)
library(coda)
library(bayesplot)
PSRF_1_beta_test <- stable.GR(BETA_POST_1)
PSRF_2_beta_test <- stable.GR(BETA_POST_2)
PSRF_3_beta_test <- stable.GR(BETA_POST_3)
PSRF_4_beta_test <- stable.GR(BETA_POST_4)
PSRF_5_beta_test <- stable.GR(BETA_POST_5)
PSRF_1_upsilon_test <- stable.GR(upsilon_POST_1)
PSRF_2_upsilon_test <- stable.GR(upsilon_POST_2)
PSRF_3_upsilon_test <- stable.GR(upsilon_POST_3)
PSRF_4_upsilon_test <- stable.GR(upsilon_POST_4)
PSRF_5_upsilon_test <- stable.GR(upsilon_POST_5)
PSRF_1_sigma_test <- stable.GR(sigmasq_upsilon_POST_1)
PSRF_2_sigma_test <- stable.GR(sigmasq_upsilon_POST_2)
PSRF_3_sigma_test <- stable.GR(sigmasq_upsilon_POST_3)
PSRF_4_sigma_test <- stable.GR(sigmasq_upsilon_POST_4)
PSRF_5_sigma_test <- stable.GR(sigmasq_upsilon_POST_5)
###95% Credible interval &coverage probability of beta chain 1
Beta_mcmc_object_1 <- as.mcmc(BETA_POST_1)
Beta_credible_interval_1 <- HPDinterval(Beta_mcmc_object_1 , prob = 0.95)
coverage_prob_beta_1 <- matrix(NA, 5,1)  
for (j in 1:5 ){
  if (beta[j,1] >= Beta_credible_interval_1[j, 1] && beta[j,1] <= Beta_credible_interval_1[j, 2]) {
    coverage_prob_beta_1[j,1] <- 1  
  } else {
    coverage_prob_beta_1[j,1] <- 0
  }
}
 mean(coverage_prob_beta_1) 
###95% Credible interval &coverage probability of upsilon chain 1
upsilon_mcmc_object_1 <- as.mcmc(upsilon_POST_1)
upsilon_credible_interval_1 <- HPDinterval(upsilon_mcmc_object_1, prob = 0.95)
coverage_prob_upsilon_1 <- matrix(NA, m,1)  
for (j in 1:m ){
  if (upsilon[j,1] >= upsilon_credible_interval_1[j, 1] && upsilon[j,1] <= upsilon_credible_interval_1[j, 2]) {
    coverage_prob_upsilon_1[j,1] <- 1  
  } else {
    coverage_prob_upsilon_1[j,1] <- 0
  }
}
mean(coverage_prob_upsilon_1)
###95% Credible interval &coverage probability of beta chain 2
Beta_mcmc_object_2 <- as.mcmc(BETA_POST_2)
Beta_credible_interval_2 <- HPDinterval(Beta_mcmc_object_2 , prob = 0.95)
coverage_prob_beta_2 <- matrix(NA, 5,1)  
for (j in 1:5 ){
  if (beta[j,1] >= Beta_credible_interval_2[j, 1] && beta[j,1] <= Beta_credible_interval_2[j, 2]) {
    coverage_prob_beta_2[j,1] <- 1  
  } else {
    coverage_prob_beta_2[j,1] <- 0
  }
}
mean(coverage_prob_beta_2) 
###95% Credible interval &coverage probability of upsilon chain 2
upsilon_mcmc_object_2 <- as.mcmc(upsilon_POST_2)
upsilon_credible_interval_2 <- HPDinterval(upsilon_mcmc_object_2, prob = 0.95)
coverage_prob_upsilon_2 <- matrix(NA, m,1)  
for (j in 1:m ){
  if (upsilon[j,1] >= upsilon_credible_interval_2[j, 1] && upsilon[j,1] <= upsilon_credible_interval_2[j, 2]) {
    coverage_prob_upsilon_2[j,1] <- 1  
  } else {
    coverage_prob_upsilon_2[j,1] <- 0
  }
}
mean(coverage_prob_upsilon_2)
#################################################################
###95% Credible interval &coverage probability of beta chain 3
Beta_mcmc_object_3 <- as.mcmc(BETA_POST_3)
Beta_credible_interval_3 <- HPDinterval(Beta_mcmc_object_3 , prob = 0.95)
coverage_prob_beta_3 <- matrix(NA, 5,1)  
for (j in 1:5 ){
  if (beta[j,1] >= Beta_credible_interval_3[j, 1] && beta[j,1] <= Beta_credible_interval_3[j, 2]) {
    coverage_prob_beta_3[j,1] <- 1  
  } else {
    coverage_prob_beta_3[j,1] <- 0
  }
}
mean(coverage_prob_beta_3) 
###95% Credible interval &coverage probability of upsilon chain 3
upsilon_mcmc_object_3 <- as.mcmc(upsilon_POST_3)
upsilon_credible_interval_3 <- HPDinterval(upsilon_mcmc_object_3, prob = 0.95)
coverage_prob_upsilon_3 <- matrix(NA, m,1)  
for (j in 1:m ){
  if (upsilon[j,1] >= upsilon_credible_interval_3[j, 1] && upsilon[j,1] <= upsilon_credible_interval_3[j, 2]) {
    coverage_prob_upsilon_3[j,1] <- 1  
  } else {
    coverage_prob_upsilon_3[j,1] <- 0
  }
}
mean(coverage_prob_upsilon_3)
####################################################
###95% Credible interval &coverage probability of beta chain 4
Beta_mcmc_object_4 <- as.mcmc(BETA_POST_4)
Beta_credible_interval_4 <- HPDinterval(Beta_mcmc_object_4 , prob = 0.95)
coverage_prob_beta_4 <- matrix(NA, 5,1)  
for (j in 1:5 ){
  if (beta[j,1] >= Beta_credible_interval_4[j, 1] && beta[j,1] <= Beta_credible_interval_4[j, 2]) {
    coverage_prob_beta_4[j,1] <- 1  
  } else {
    coverage_prob_beta_4[j,1] <- 0
  }
}
mean(coverage_prob_beta_4) 
###95% Credible interval &coverage probability of upsilon chain 4
upsilon_mcmc_object_4 <- as.mcmc(upsilon_POST_4)
upsilon_credible_interval_4 <- HPDinterval(upsilon_mcmc_object_4, prob = 0.95)
coverage_prob_upsilon_4 <- matrix(NA, m,1)  
for (j in 1:m ){
  if (upsilon[j,1] >= upsilon_credible_interval_4[j, 1] && upsilon[j,1] <= upsilon_credible_interval_4[j, 2]) {
    coverage_prob_upsilon_4[j,1] <- 1  
  } else {
    coverage_prob_upsilon_4[j,1] <- 0
  }
}
mean(coverage_prob_upsilon_4)
####################################################
###95% Credible interval &coverage probability of beta chain 5
Beta_mcmc_object_5 <- as.mcmc(BETA_POST_5)
Beta_credible_interval_5 <- HPDinterval(Beta_mcmc_object_5 , prob = 0.95)
coverage_prob_beta_5 <- matrix(NA, 5,1)  
for (j in 1:5 ){
  if (beta[j,1] >= Beta_credible_interval_5[j, 1] && beta[j,1] <= Beta_credible_interval_5[j, 2]) {
    coverage_prob_beta_5[j,1] <- 1  
  } else {
    coverage_prob_beta_5[j,1] <- 0
  }
}
mean(coverage_prob_beta_5) 
###95% Credible interval &coverage probability of upsilon chain 5
upsilon_mcmc_object_5 <- as.mcmc(upsilon_POST_5)
upsilon_credible_interval_5 <- HPDinterval(upsilon_mcmc_object_5, prob = 0.95)
coverage_prob_upsilon_5 <- matrix(NA, m,1)  
for (j in 1:m ){
  if (upsilon[j,1] >= upsilon_credible_interval_5[j, 1] && upsilon[j,1] <= upsilon_credible_interval_5[j, 2]) {
    coverage_prob_upsilon_5[j,1] <- 1  
  } else {
    coverage_prob_upsilon_5[j,1] <- 0
  }
}
mean(coverage_prob_upsilon_5)


#Checking the convergence with burn
z_gibbis_burn<-z_all_4[seq(1, nrow(z_all_4[10001:4950000,]),100),]
Beta_gibbis_burn<-BETA_POST_4[seq(1, nrow(BETA_POST_4[10001:50000,]),10),]
upsilon_gibbis_burn<-upsilon_POST_4[seq(1, nrow(upsilon_POST_4[10001:50000,]),10),]
sigma_gibbis_burn<- sigmasq_upsilon_POST_4[seq(1, nrow( sigmasq_upsilon_POST_4[10001:50000,]),10),]
####for Beta#####
par(mfrow=c(5,1),mar=c(3,3,2,1))
traceplot(as.mcmc(Beta_gibbis_burn[,1]), main=expression(paste(beta[1])))
traceplot(as.mcmc(Beta_gibbis_burn[,2]),main=expression(paste(beta[2])))
traceplot(as.mcmc(Beta_gibbis_burn[,3]),main=expression(paste(beta[3])))
traceplot(as.mcmc(Beta_gibbis_burn[,4]),main=expression(paste(beta[4])))
traceplot(as.mcmc(Beta_gibbis_burn[,5]),main=expression(paste(beta[5])))
par(mfrow=c(3,2),mar=c(7,3,3,1))
densplot(as.mcmc(Beta_gibbis_burn[,1]), main=expression(paste(beta[1])))
densplot(as.mcmc(Beta_gibbis_burn[,2]),main=expression(paste(beta[2])))
densplot(as.mcmc(Beta_gibbis_burn[,3]),main=expression(paste(beta[3])))
densplot(as.mcmc(Beta_gibbis_burn[,4]), main=expression(paste(beta[4])))
densplot(as.mcmc(Beta_gibbis_burn[,5]), main=expression(paste(beta[5])))
par(mfrow=c(3,2),mar=c(7,4,3,1))
acf(as.mcmc(Beta_gibbis_burn[,1]),ylab="Autocorrelation",ci=F,lag.max=4000, main=expression(paste(beta[1])))
acf(as.mcmc(Beta_gibbis_burn[,2]),ylab="Autocorrelation",ci=F,lag.max=4000,main=expression(paste(beta[2])))
acf(as.mcmc(Beta_gibbis_burn[,3]),ylab="Autocorrelation",ci=F,lag.max=4000,main=expression(paste(beta[3])))
acf(as.mcmc(Beta_gibbis_burn[,4]),ylab="Autocorrelation",ci=F,lag.max=4000,main=expression(paste(beta[4])))
acf(as.mcmc(Beta_gibbis_burn[,5]),ylab="Autocorrelation",ci=F,lag.max=4000,main=expression(paste(beta[5])))
#Posterior means of Beta
mean(Beta_gibbis_burn[,1])
mean(Beta_gibbis_burn[,2])
mean(Beta_gibbis_burn[,3]) 
mean(Beta_gibbis_burn[,4])
mean(Beta_gibbis_burn[,5])
#Posterior sd of Beta
sd(Beta_gibbis_burn[,1])
sd(Beta_gibbis_burn[,2])
sd(Beta_gibbis_burn[,3]) 
sd(Beta_gibbis_burn[,4])
sd(Beta_gibbis_burn[,5])

coverage_probb <- mean(coverage_prob) 
###to compute coverage probability######
latent_credible_interval <- HPDinterval(as.mcmc(z_all_1), prob = 0.95)
coverage_prob <- matrix(NA, 99,m)  
for (j in 1:99 ){
  for(i in 1:m){
    if (z[j,i] >= latent_credible_interval[i, 1] && z[j,i] <= latent_credible_interval[i, 2]) {
      coverage_prob[j,i] <- 1  
    } else {
      coverage_prob[j,i] <- 0
    }
  }
}
coverage_probb <- mean(coverage_prob) 

######For upsilon####
par(mfrow=c(6,1),mar=c(3,3,2,2))
traceplot(as.mcmc(upsilon_gibbis_burn[,1]),main=expression(paste(upsilon[1])))
traceplot(as.mcmc(upsilon_gibbis_burn[,2]),main=expression(paste(upsilon[2])))
traceplot(as.mcmc(upsilon_gibbis_burn[,3]),main=expression(paste(upsilon[3])))
traceplot(as.mcmc(upsilon_gibbis_burn[,4]),main=expression(paste(upsilon[4])))
traceplot(as.mcmc(upsilon_gibbis_burn[,5]),main=expression(paste(upsilon[5])))
traceplot(as.mcmc(upsilon_gibbis_burn[,6]),main=expression(paste(upsilon[6])))
par(mfrow=c(5,1),mar=c(3,3,2,2))
traceplot(as.mcmc(upsilon_gibbis_burn[,7]),main=expression(paste(upsilon[7])))
traceplot(as.mcmc(upsilon_gibbis_burn[,8]),main=expression(paste(upsilon[8])))
traceplot(as.mcmc(upsilon_gibbis_burn[,9]),main=expression(paste(upsilon[9])))
traceplot(as.mcmc(upsilon_gibbis_burn[,10]),main=expression(paste(upsilon[10])))
traceplot(as.mcmc(upsilon_gibbis_burn[,11]),main=expression(paste(upsilon[11])))
par(mfrow=c(3,2),mar=c(7,3,3,1))
densplot(as.mcmc(upsilon_gibbis_burn[,1]), main=expression(paste(upsilon[1])))
densplot(as.mcmc(upsilon_gibbis_burn[,2]),main=expression(paste(upsilon[2])))
densplot(as.mcmc(upsilon_gibbis_burn[,3]),main=expression(paste(upsilon[3])))
densplot(as.mcmc(upsilon_gibbis_burn[,4]), main=expression(paste(upsilon[4])))
densplot(as.mcmc(upsilon_gibbis_burn[,5]), main=expression(paste(upsilon[5])))
densplot(as.mcmc(upsilon_gibbis_burn[,6]),main=expression(paste(upsilon[6])))
par(mfrow=c(3,2),mar=c(7,3,3,1))
densplot(as.mcmc(upsilon_gibbis_burn[,7]),main=expression(paste(upsilon[7])))
densplot(as.mcmc(upsilon_gibbis_burn[,8]), main=expression(paste(upsilon[8])))
densplot(as.mcmc(upsilon_gibbis_burn[,9]), main=expression(paste(upsilon[9])))
densplot(as.mcmc(upsilon_gibbis_burn[,10]),main=expression(paste(upsilon[10])))
densplot(as.mcmc(upsilon_gibbis_burn[,11]),main=expression(paste(upsilon[11])))
par(mfrow=c(3,2),mar=c(7,4,3,1))
acf(as.mcmc(upsilon_gibbis_burn[,1]),ylab="Autocorrelation",ci=F,lag.max=4000, main=expression(paste(upsilon[1])))
acf(as.mcmc(upsilon_gibbis_burn[,2]),ylab="Autocorrelation",ci=F,lag.max=4000,main=expression(paste(upsilon[2])))
acf(as.mcmc(upsilon_gibbis_burn[,3]),ylab="Autocorrelation",ci=F,lag.max=4000,main=expression(paste(upsilon[3])))
acf(as.mcmc(upsilon_gibbis_burn[,4]),ylab="Autocorrelation",ci=F,lag.max=4000,main=expression(paste(upsilon[4])))
acf(as.mcmc(upsilon_gibbis_burn[,5]),ylab="Autocorrelation",ci=F,lag.max=4000, main=expression(paste(upsilon[5])))
acf(as.mcmc(upsilon_gibbis_burn[,6]),ylab="Autocorrelation",ci=F,lag.max=4000,main=expression(paste(upsilon[6])))
par(mfrow=c(3,2),mar=c(7,4,3,1))
acf(as.mcmc(upsilon_gibbis_burn[,7]),ylab="Autocorrelation",ci=F,lag.max=4000,main=expression(paste(upsilon[7])))
acf(as.mcmc(upsilon_gibbis_burn[,8]),ylab="Autocorrelation",ci=F,lag.max=4000,main=expression(paste(upsilon[8])))
acf(as.mcmc(upsilon_gibbis_burn[,9]),ylab="Autocorrelation",ci=F,lag.max=4000, main=expression(paste(upsilon[9])))
acf(as.mcmc(upsilon_gibbis_burn[,10]),ylab="Autocorrelation",ci=F,lag.max=4000,main=expression(paste(upsilon[10])))
acf(as.mcmc(upsilon_gibbis_burn[,11]),ylab="Autocorrelation",ci=F,lag.max=4000,main=expression(paste(upsilon[11])))
#Posterior means of upsilon
mean(upsilon_gibbis_burn[,1])
mean(upsilon_gibbis_burn[,2])
mean(upsilon_gibbis_burn[,3])
mean(upsilon_gibbis_burn[,4])
mean(upsilon_gibbis_burn[,5])
mean(upsilon_gibbis_burn[,6])
mean(upsilon_gibbis_burn[,7])
mean(upsilon_gibbis_burn[,8])
mean(upsilon_gibbis_burn[,9])
mean(upsilon_gibbis_burn[,10])
mean(upsilon_gibbis_burn[,11])
#Posterior sd of upsilon
sd(upsilon_gibbis_burn[,1])
sd(upsilon_gibbis_burn[,2])
sd(upsilon_gibbis_burn[,3])
sd(upsilon_gibbis_burn[,4])
sd(upsilon_gibbis_burn[,5])
sd(upsilon_gibbis_burn[,6])
sd(upsilon_gibbis_burn[,7])
sd(upsilon_gibbis_burn[,8])
sd(upsilon_gibbis_burn[,9])
sd(upsilon_gibbis_burn[,10])
sd(upsilon_gibbis_burn[,11])
####for sigma ###
par(mfrow=c(4,1),mar=c(3,3,2,2))
traceplot(as.mcmc(sigma_gibbis_burn[,1]),main=expression(paste(sigma[1]^2)))
traceplot(as.mcmc(sigma_gibbis_burn[,2]),main=expression(paste(sigma[2]^2)))
traceplot(as.mcmc(sigma_gibbis_burn[,3]),main=expression(paste(sigma[3]^2)))
traceplot(as.mcmc(sigma_gibbis_burn[,4]),main=expression(paste(sigma[4]^2)))
par(mfrow=c(4,1),mar=c(3,3,2,2))
traceplot(as.mcmc(sigma_gibbis_burn[,5]),main=expression(paste(sigma[5]^2)))
traceplot(as.mcmc(sigma_gibbis_burn[,6]),main=expression(paste(sigma[6]^2)))
traceplot(as.mcmc(sigma_gibbis_burn[,7]),main=expression(paste(sigma[7]^2)))
traceplot(as.mcmc(sigma_gibbis_burn[,8]),main=expression(paste(sigma[8]^2)))
par(mfrow=c(3,1),mar=c(3,3,2,2))
traceplot(as.mcmc(sigma_gibbis_burn[,9]),main=expression(paste(sigma[9]^2)))
traceplot(as.mcmc(sigma_gibbis_burn[,10]),main=expression(paste(sigma[10]^2)))
traceplot(as.mcmc(sigma_gibbis_burn[,11]),main=expression(paste(sigma[11]^2)))
par(mfrow=c(3,2),mar=c(7,3,2,2))
densplot(as.mcmc(sigma_gibbis_burn[,1]), main=expression(paste(sigma[1]^2)))
densplot(as.mcmc(sigma_gibbis_burn[,2]),main=expression(paste(sigma[2]^2)))
densplot(as.mcmc(sigma_gibbis_burn[,3]),main=expression(paste(sigma[3]^2)))
densplot(as.mcmc(sigma_gibbis_burn[,4]), main=expression(paste(sigma[4]^2)))
densplot(as.mcmc(sigma_gibbis_burn[,5]), main=expression(paste(sigma[5]^2)))
densplot(as.mcmc(sigma_gibbis_burn[,6]),main=expression(paste(sigma[6]^2)))
par(mfrow=c(3,2),mar=c(7,3,2,2))
densplot(as.mcmc(sigma_gibbis_burn[,7]),main=expression(paste(sigma[7]^2)))
densplot(as.mcmc(sigma_gibbis_burn[,8]), main=expression(paste(sigma[8]^2)))
densplot(as.mcmc(sigma_gibbis_burn[,9]), main=expression(paste(sigma[9]^2)))
densplot(as.mcmc(sigma_gibbis_burn[,10]),main=expression(paste(sigma[10]^2)))
densplot(as.mcmc(sigma_gibbis_burn[,11]),main=expression(paste(sigma[11]^2)))
par(mfrow=c(3,2),mar=c(7,4,3,1))
acf(as.mcmc(sigma_gibbis_burn[,1]),ylab="Autocorrelation",ci=F,lag.max=4000, main=expression(paste(sigma[1]^2)))
acf(as.mcmc(sigma_gibbis_burn[,2]),ylab="Autocorrelation",ci=F,lag.max=4000,main=expression(paste(sigma[2]^2)))
acf(as.mcmc(sigma_gibbis_burn[,3]),ylab="Autocorrelation",ci=F,lag.max=4000,main=expression(paste(sigma[3]^2)))
acf(as.mcmc(sigma_gibbis_burn[,4]),ylab="Autocorrelation",ci=F,lag.max=4000,main=expression(paste(sigma[4]^2)))
acf(as.mcmc(sigma_gibbis_burn[,5]),ylab="Autocorrelation",ci=F,lag.max=4000, main=expression(paste(sigma[5]^2)))
acf(as.mcmc(sigma_gibbis_burn[,6]),ylab="Autocorrelation",ci=F,lag.max=4000,main=expression(paste(sigma[6]^2)))
par(mfrow=c(3,2),mar=c(7,4,3,1))
acf(as.mcmc(sigma_gibbis_burn[,7]),ylab="Autocorrelation",ci=F,lag.max=4000,main=expression(paste(sigma[7]^2)))
acf(as.mcmc(sigma_gibbis_burn[,8]),ylab="Autocorrelation",ci=F,lag.max=4000,main=expression(paste(sigma[8]^2)))
acf(as.mcmc(sigma_gibbis_burn[,9]),ylab="Autocorrelation",ci=F,lag.max=4000, main=expression(paste(sigma[9]^2)))
acf(as.mcmc(sigma_gibbis_burn[,10]),ylab="Autocorrelation",ci=F,lag.max=4000,main=expression(paste(sigma[10]^2)))
acf(as.mcmc(sigma_gibbis_burn[,11]),ylab="Autocorrelation",ci=F,lag.max=4000,main=expression(paste(sigma[11]^2)))
#Posterior means of sigma
mean(sigma_gibbis_burn[,1])
mean(sigma_gibbis_burn[,2])
mean(sigma_gibbis_burn[,3])
mean(sigma_gibbis_burn[,4])
mean(sigma_gibbis_burn[,5])
mean(sigma_gibbis_burn[,6])
mean(sigma_gibbis_burn[,7])
mean(sigma_gibbis_burn[,8])
mean(sigma_gibbis_burn[,9])
mean(sigma_gibbis_burn[,10])
mean(sigma_gibbis_burn[,11])
#Variation between the groups
library(reshape2)
library(viridis)
library(ggplot2)
library(bayesplot)
new_colnames <- c("1","2","3","4","5","6","7","8","9","10","11")
colnames(g) <- new_colnames
g=data.frame(upsilon_POST_4)
data <- var(g[sapply(g,is.numeric)])
data1 <- melt(data)
colnames(data1)=c("Var1","Var2", "Variations")
x_labels <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11")
y_labels <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11")
ggplot(data1, aes(x = Var1, y = Var2, fill = Variations)) +
  geom_tile() +scale_fill_viridis(discrete = FALSE)+geom_tile()+
  labs(title = "",
       x = "Regions",
       y = "Regions")+scale_x_discrete(labels = x_labels)+
  scale_y_discrete(labels = y_labels)
#Posterior predictive check
#The last observation
last_obs_1=array(NA,c(ncol(X),m,1))
for(j in 1:m){
  last_obs_1[,j,1]=X_init[(n_obs+1),,j]
}

z_pred_outsample_1=array(NA,c(n_iter,1,m))
for(i in 1:1){
  for(j in 1:m){
    z_pred_outsample_1[,i,j]=rnorm(n_iter,BETA_POST_5%*%last_obs_1[,j,i]+upsilon_POST_1[,j],
                                   1)
  }
}

Y_init=latent_cat

posterior_meansz1=matrix(NA,1,m)
posterior_meansz1=t(colMeans(z_pred_outsample_1,dims=1))
#print(posterior_means_z1)
#print(latent_cat[100,])
###################
par(mfrow=c(1,2))
plot(latents[100,],type="b",main="Real data",ylab = "Dependent variable")
#lines(Y_init[16,],col="green")
plot(posterior_meansz1[,1],type="b",main="Predicted data",ylab="Latent variable")
#lines(posterior_means_z2[,2],col="red")
library("ggplot2")
d1 = data.frame(x=c(1,2,3,4,5,6,7,8,9,10,11), y=latents[100,]) # real data
d2 = data.frame(x=c(1,2,3,4,5,6,7,8,9,10,11), y=posterior_meansz1[,1]) # predicted data
par(mar=c(5, 4, 4, 4) + 0.3)              # Additional space for second y-axis
plot(d1$x, d1$y, pch = 16, col = "blue",xlab="",ylab="",cex=2,main="Predictions with GLMM using Rank Likelihood")              # Create first plot
par(new = TRUE)                             # Add new plot
plot(d2$x, d2$y, pch = 17, col = "green",cex=2,              # Create second plot without axes
     axes = FALSE, xlab = "Groups", ylab = "Observed variable Y")
axis(side = 4, at = pretty(range(d2$y)))      # Add second axis
mtext("Latent variable z", side = 4, line = 3)             # Add second axis label
##################################################################
###########Estimates for the classical Model###################
mydata=read.csv('C:\\Users\\user\\Desktop\\data.csv')
View(mydata)
attach(mydata)
library(brms)
library(rstan)
require(brms)
require(rstan)
library(StanHeaders)
##############################
model=brm(formula= Response ~ V1+V2+V3+V4+V5,data=mydata, family = cumulative)
summary(model)
posterior_samples <- posterior_samples(model)
point_estimates <- apply(posterior_samples[,9:13], 2, mean)
head(posterior_samples[,9:13])
head(mse2)
#####MSE and MAE of the classical model
x=matrix(point_estimates,5,1)
y=posterior_samples[,9:13]
mse=matrix(NA,4000,5)
for(i in 1:4000){
  for(j in 1:5){
    mse2[i,j]=(y[i,j]-x[j,1])^2
  }
}
mse1=mean(mse)
RMSE=sqrt(mse1)
######
mae=matrix(NA,4000,5)
for(i in 1:4000){
  for(j in 1:5){
    mae[i,j]=abs(y[i,j]-x[j,1])
  }
}
mae1=mean(mae)

p=posterior_summary(fit2)
######Covarage probability of the classical model 
point_estimates=matrix(point_estimates,5,1)
param_intervals <- posterior_interval(model, prob = 0.95)
param_intervals_2=param_intervals[9:13,]
coverage_prob=matrix(NA, 5,1) 
for (j in 1:5 ){
  if (point_estimates[j,1] >= param_intervals_2[j, 1] && point_estimates[j,1] <= param_intervals_2[j, 2]) {
    coverage_prob[j,1] <- 1  
  } else {
    coverage_prob[j,1] <- 0
  }
}
mean(coverage_prob)
