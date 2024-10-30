# R-code used for "On Latent Trait Model with Bayesian Marginal Likelihood of Rank-Based Estimation"
library("LaplacesDemon")
set.seed(12345)
n=300#sample size
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
lamda2_prior=2
#prior degrees of belief to update sigma square of upsilon (random effects) 
omega_upsilon_prior=15
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
chains_upsilons=array(cbind(upsilon1[1:50000,],upsilon2[1:50000,],upsilon3[1:50000,],upsilon4[1:50000,]                    ,upsilon5[1:50000,],upsilon6[1:50000,],upsilon7[1:50000,],upsilon8[1:50000,],upsilon9[1:50000,],upsilon10[1:50000,],upsilon11[1:50000,]),c(50000,5,11))
new_colnames <- c("chain1", "chain2","chain3","chain4","chain5")
dimnames(chains_upsilons)[[2]]=new_colnames
dimnames(chains_upsilons)[[3]] <- c("\u03C51", "\u03C52", "\u03C53","\u03C54","\u03C55","\u03C56","\u03C57","\u03C58","\u03C59","\u03C510","\u03C511")
mcmc_trace(chains_upsilons, pars = c("\u03C51", "\u03C52", "\u03C53","\u03C54","\u03C55","\u03C56","\u03C57","\u03C58","\u03C59","\u03C510","\u03C511"), facet_args = list(ncol = 2, strip.position = "left"))

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
####for upsilon chain 5####
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
Beta_gibbis_burn<-BETA_POST_4[seq(1, nrow(BETA_POST_4[10001:50000,]),10),]
upsilon_gibbis_burn<-upsilon_POST_4[seq(1, nrow(upsilon_POST_4[10001:50000,]),10),]
##############To compute batch mean and SE####
library(coda)
library(ggplot2)
# Custom function to calculate batch means for each parameter in an MCMC matrix
calculate_batch_means_matrix <- function(samples_matrix, batch_size) {
  n <- nrow(samples_matrix)
  n_batches <- floor(n / batch_size)
  # Trim samples to fit batch size
  samples_matrix <- samples_matrix[1:(n_batches * batch_size), ]
  
  # Reshape matrix into 3D array: (batch_size, n_batches, n_parameters)
  reshaped_samples <- array(samples_matrix, dim = c(batch_size, n_batches, ncol(samples_matrix)))
  
  # Calculate batch means across batches for each parameter
  batch_means <- apply(reshaped_samples, 2, colMeans)
  
  return(batch_means)
}

set.seed(123)

batch_sizes <- c(100, 500, 1000, 1500, 2000)

# Data frame to store means and standard errors for each batch size and each parameter
results <- data.frame(batch_size = integer(), parameter = character(), 
                      batch_mean = numeric(), standard_error = numeric())

# Loop over batch sizes and calculate means and standard errors
for (batch_size in batch_sizes) {
  # Calculate batch means for each parameter at this batch size
  batch_means <- calculate_batch_means_matrix(Beta_gibbis_burn, batch_size)
  
  # Calculate standard errors for each parameter based on batch means
  standard_errors <- apply(batch_means, 1, function(x) sd(x) / sqrt(length(x)))
  
  # Store results in data frame
  for (i in seq_along(standard_errors)) {
    results <- rbind(results, data.frame(
      batch_size = batch_size,
      parameter = paste0("\u03B2", i),
      batch_mean = mean(batch_means[i, ]),
      standard_error = standard_errors[i]
    ))
  }
}

# Display results
mydata2=(results)
##################################################
library(ggplot2)
library(cowplot)

# Plot 1: Standard error vs. batch size
plot1 <- ggplot(mydata2, aes(x = batch_size, y = standard_error, color = Parameters)) +
  geom_line() +                 
  geom_point() +                
  labs(
    x = "Sample size",
    y = "Batch standard error"
  ) +
  theme_minimal() +  
  theme(panel.background = element_rect(fill = "grey"), legend.position = "none") + # Remove legend here
  scale_x_continuous(breaks = batch_sizes)

# Plot 2: Batch mean vs. batch size
plot2 <- ggplot(results, aes(x = batch_size, y = batch_mean, color = Parameters)) +
  geom_line() +                 
  geom_point() +                
  labs(
    x = "Sample size",
    y = "Batch mean"
  ) +
  theme_minimal() +              
  theme(panel.background = element_rect(fill = "grey"), legend.position = "none") + # Remove legend here
  scale_x_continuous(breaks = batch_sizes)

# Extract the legend from one of the plots
legend <- get_legend(
  ggplot(results, aes(x = batch_size, y = batch_mean, color = Parameters)) +
    geom_line() +
    theme_minimal() +
    theme(legend.position = "right")  # Adjust legend position if needed
)

# Combine plots and legend
combined_plot <- plot_grid(plot1, plot2, legend, ncol = 3, rel_widths = c(1, 1, 0.2))

# Display the combined plot with a common legend
print(combined_plot)
####################################################################for upsilon########
calculate_batch_means_matrix <- function(samples_matrix, batch_size) {
  n <- nrow(samples_matrix)
  n_batches <- floor(n / batch_size)
  # Trim samples to fit batch size
  samples_matrix <- samples_matrix[1:(n_batches * batch_size), ]
  
  # Reshape matrix into 3D array: (batch_size, n_batches, n_parameters)
  reshaped_samples <- array(samples_matrix, dim = c(batch_size, n_batches, ncol(samples_matrix)))
  
  # Calculate batch means across batches for each parameter
  batch_means <- apply(reshaped_samples, 2, colMeans)
  
  return(batch_means)
}

# Example MCMC samples in matrix form (10000 samples, 3 parameters)
set.seed(123)

# Different batch sizes to test
batch_sizes <- c(100, 500, 1000, 1500, 2000)

# Data frame to store means and standard errors for each batch size and each parameter
results <- data.frame(batch_size = integer(), parameter = character(), 
                      batch_mean = numeric(), standard_error = numeric())

# Loop over batch sizes and calculate means and standard errors
for (batch_size in batch_sizes) {
  # Calculate batch means for each parameter at this batch size
  batch_means <- calculate_batch_means_matrix(upsilon_gibbis_burn, batch_size)
  
  # Calculate standard errors for each parameter based on batch means
  standard_errors <- apply(batch_means, 1, function(x) sd(x) / sqrt(length(x)))
  
  # Store results in data frame
  for (i in seq_along(standard_errors)) {
    results <- rbind(results, data.frame(
      batch_size = batch_size,
      parameter = paste0("\u03C5", i),
      batch_mean = mean(batch_means[i, ]),
      standard_error = standard_errors[i]
    ))
  }
}
upsilondata2=(results)
######################
upsilondata2$Parameters <- factor(upsilondata2$Parameters, levels = paste0("υ", 1:11))

# Define batch sizes for axis breaks
batch_sizes <- unique(upsilondata2$batch_size)
# Plot 1: Standard error vs. batch size
plot1 <- ggplot(upsilondata2, aes(x = batch_size, y = standard_error, color = Parameters)) +
  geom_line() +                 
  geom_point() +                
  labs(
    x = "Sample size",
    y = "Batch standard error"
  ) +
  theme_minimal() +  
  theme(panel.background = element_rect(fill = "grey"), legend.position = "none") + # Remove legend here
  scale_x_continuous(breaks = batch_sizes)

# Plot 2: Batch mean vs. batch size
plot2 <- ggplot(upsilondata2, aes(x = batch_size, y = batch_mean, color = Parameters)) +
  geom_line() +                 
  geom_point() +                
  labs(
    x = "Sample size",
    y = "Batch mean"
  ) +
  theme_minimal() +              
  theme(panel.background = element_rect(fill = "grey"), legend.position = "none") + # Remove legend here
  scale_x_continuous(breaks = batch_sizes)

# Extract the legend from one of the plots
legend <- get_legend(
  ggplot(upsilondata2, aes(x = batch_size, y = batch_mean, color = Parameters)) +
    geom_line() +
    theme_minimal() +
    theme(legend.position = "right")  # Adjust legend position if needed
)

# Combine plots and legend
combined_plot <- plot_grid(plot1, plot2, legend, ncol = 3, rel_widths = c(1, 1, 0.2))

# Display the combined plot with a common legend
print(combined_plot)
#####################
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
###To chech posterior predictive
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
####R Code for real-data##########################################################################################
####################################################################################################
library("LaplacesDemon")
set.seed(12345)
data=read.csv('D:\\Myfiles\\Mythesisresearch\\birthweight.csv')
attach(data)
X_tig=as.matrix(subset(data[,2:6], Regions==1))
X_afa=as.matrix(subset(data[,2:6], Regions==2))
X_amh=as.matrix(subset(data[,2:6], Regions==3))
X_oro=as.matrix(subset(data[,2:6], Regions==4))
X_som=as.matrix(subset(data[,2:6], Regions==5))
X_ben=as.matrix(subset(data[,2:6], Regions==6))
X_snn=as.matrix(subset(data[,2:6], Regions==7))
X_gam=as.matrix(subset(data[,2:6], Regions==8))
X_har=as.matrix(subset(data[,2:6], Regions==9))
X_add=as.matrix(subset(data[,2:6], Regions==10))
X_dir=as.matrix(subset(data[,2:6], Regions==11))
#I need to rescale the covariates
X_tig=scale(X_tig,center=TRUE,scale=TRUE)
X_afa=scale(X_afa,center=TRUE,scale=TRUE)
X_amh=scale(X_amh,center=TRUE,scale=TRUE)
X_oro=scale(X_oro,center=TRUE,scale=TRUE)
X_som=scale(X_som,center=TRUE,scale=TRUE)
X_ben=scale(X_ben,center=TRUE,scale=TRUE)
X_snn=scale(X_snn,center=TRUE,scale=TRUE)
X_gam=scale(X_gam,center=TRUE,scale=TRUE)
X_har=scale(X_har,center=TRUE,scale=TRUE)
X_add=scale(X_add,center=TRUE,scale=TRUE)
X_dir=scale(X_dir,center=TRUE,scale=TRUE)
rownames(X_tig) <- 1:nrow(X_tig)
rownames(X_afa) <- 1:nrow(X_afa)
rownames(X_amh) <- 1:nrow(X_amh)
rownames(X_oro) <- 1:nrow(X_oro)
rownames(X_som) <- 1:nrow(X_som)
rownames(X_ben) <- 1:nrow(X_ben)
rownames(X_snn) <- 1:nrow(X_snn)
rownames(X_gam) <- 1:nrow(X_gam)
rownames(X_har) <- 1:nrow(X_har)
rownames(X_add) <- 1:nrow(X_add)
rownames(X_dir) <- 1:nrow(X_dir)
set.seed(12345)
X=list(X_tig,X_afa,X_amh,X_oro,X_som,X_ben,X_snn,X_gam,X_har,X_add,X_dir)
y_tig=matrix(subset(data[,1], Regions==1),nrow(X_tig),1)
y_afa=matrix(subset(data[,1], Regions==2),nrow(X_afa),1)
y_amh=matrix(subset(data[,1], Regions==3),nrow(X_amh),1)
y_oro=matrix(subset(data[,1], Regions==4),nrow(X_oro),1)
y_som=matrix(subset(data[,1], Regions==5),nrow(X_som),1)
y_ben=matrix(subset(data[,1], Regions==6),nrow(X_ben),1)
y_snn=matrix(subset(data[,1], Regions==7),nrow(X_snn),1)
y_gam=matrix(subset(data[,1], Regions==8),nrow(X_gam),1)
y_har=matrix(subset(data[,1], Regions==9),nrow(X_har),1)
y_add=matrix(subset(data[,1], Regions==10),nrow(X_add),1)
y_dir=matrix(subset(data[,1], Regions==11),nrow(X_dir),1)
Y=list(y_tig,y_afa,y_amh,y_oro,y_som,y_ben,y_snn,y_gam,y_har,y_add,y_dir)
#Iinitial values for z
z_tig=matrix(rep(5,nrow(X_tig)),nrow(X_tig),1)
z_afa=matrix(rep(5,nrow(X_afa)),nrow(X_afa),1)
z_amh=matrix(rep(5,nrow(X_amh)),nrow(X_amh),1)
z_oro=matrix(rep(5,nrow(X_oro)),nrow(X_oro),1)
z_som=matrix(rep(5,nrow(X_som)),nrow(X_som),1)
z_ben=matrix(rep(5,nrow(X_ben)),nrow(X_ben),1)
z_snn=matrix(rep(5,nrow(X_snn)),nrow(X_snn),1)
z_gam=matrix(rep(5,nrow(X_gam)),nrow(X_gam),1)
z_har=matrix(rep(5,nrow(X_har)),nrow(X_har),1)
z_add=matrix(rep(5,nrow(X_add)),nrow(X_add),1)
z_dir=matrix(rep(5,nrow(X_dir)),nrow(X_dir),1)
z=list(z_tig,z_afa,z_amh,z_oro,z_som,z_ben,z_snn,z_gam,z_har,z_add,z_dir)
m=11
n_obs=matrix(c(length(z_tig),length(z_afa),length(z_amh),length(z_oro),length(z_som),length(z_ben),length(z_snn),length(z_gam),length(z_har),length(z_add),length(z_dir)),m,1)
#prior hyperparameters
lamda2_prior=0.3
#prior degrees of belief to update sigma square of upsilon (latent traits) 
omega_upsilon_prior=19
#Initial values for BETA
beta=matrix(c(-1.2,-0.9, -0.8,-0.9, 1))
BETA=matrix(NA,ncol(X_tig),1)
#I need this package for scaled inverse chi-squared
library("geoR")
library("mvtnorm")
library(MCMCpack)
p=ncol(X[[1]])
omega_upsilon_posterior=NULL #posterior degrees of freedom of variance of latent trait
lamda2_posterior=NULL #posterior scale parameter of sigma of variance of random latent trait
mean_upsilon_posterior=var_upsilon_posterior=NULL#post values of mean and
#variance of latent trait
n_iter=30000
library("msm")

BETA_POST=upsilon_POST=sigmasq_upsilon_POST=NULL
upsilon_posterior=rep(0.1,m)
sigmasq_upsilon_post=NULL#posterior value of sigmasq of upsilon
z_all=NULL
#GIBBS SAMPLER
for(i in 1:n_iter){
  for(j in 1:m){
    omega_upsilon_posterior[j]=(omega_upsilon_prior+1) #updating degrees of freedom of variance of latent trait
    
    lamda2_posterior[j]=(lamda2_prior+
                           ((upsilon_posterior[j])^2)/omega_upsilon_posterior[j]) #update scale parameter of variance of latent trait
    
    sigmasq_upsilon_post[j]=rinvchisq(1,omega_upsilon_posterior[j],
                                      lamda2_posterior[j]) #variance of latent trait
    
    
  } 
  
  #given variance I can sample from upsilon (latent traits)
  mean_upsilon_posterior=list()
  mean_upsilon_posterior[[1]]=matrix(NA,length(z_tig),1)
  mean_upsilon_posterior[[2]]=matrix(NA,length(z_afa),1)
  mean_upsilon_posterior[[3]]=matrix(NA,length(z_amh),1)
  mean_upsilon_posterior[[4]]=matrix(NA,length(z_oro),1)
  mean_upsilon_posterior[[5]]=matrix(NA,length(z_som),1)
  mean_upsilon_posterior[[6]]=matrix(NA,length(z_ben),1)
  mean_upsilon_posterior[[7]]=matrix(NA,length(z_snn),1)
  mean_upsilon_posterior[[8]]=matrix(NA,length(z_gam),1)
  mean_upsilon_posterior[[9]]=matrix(NA,length(z_har),1)
  mean_upsilon_posterior[[10]]=matrix(NA,length(z_add),1)
  mean_upsilon_posterior[[11]]=matrix(NA,length(z_dir),1)
  for(t in 1:length(z_tig)){
      mean_upsilon_posterior[[1]][t,]=(z[[1]][t,]-X[[1]][t,]%*%beta)
  }
  for(t in 1:length(z_afa)){
    mean_upsilon_posterior[[2]][t,]=(z[[2]][t,]-X[[2]][t,]%*%beta)
  }
  for(t in 1:length(z_amh)){
    mean_upsilon_posterior[[3]][t,]=(z[[3]][t,]-X[[3]][t,]%*%beta)
  }
  for(t in 1:length(z_oro)){
    mean_upsilon_posterior[[4]][t,]=(z[[4]][t,]-X[[4]][t,]%*%beta)
  }
  for(t in 1:length(z_som)){
    mean_upsilon_posterior[[5]][t,]=(z[[5]][t,]-X[[5]][t,]%*%beta)
  }
  for(t in 1:length(z_ben)){
    mean_upsilon_posterior[[6]][t,]=(z[[6]][t,]-X[[6]][t,]%*%beta)
  }
  for(t in 1:length(z_snn)){
    mean_upsilon_posterior[[7]][t,]=(z[[7]][t,]-X[[7]][t,]%*%beta)
  }
  for(t in 1:length(z_gam)){
    mean_upsilon_posterior[[8]][t,]=(z[[8]][t,]-X[[8]][t,]%*%beta)
  }
  for(t in 1:length(z_har)){
    mean_upsilon_posterior[[9]][t,]=(z[[9]][t,]-X[[9]][t,]%*%beta)
  }
  for(t in 1:length(z_add)){
    mean_upsilon_posterior[[10]][t,]=(z[[10]][t,]-X[[10]][t,]%*%beta)
  }
  for(t in 1:length(z_dir)){
    mean_upsilon_posterior[[11]][t,]=(z[[11]][t,]-X[[11]][t,]%*%beta)
  }
  ###########################
  
  
  for(j in 1:m){
    var_upsilon_posterior[j]=1/(n_obs[j,]+(1/sigmasq_upsilon_post[j]))
  }
  #summing over time:
  mean_upsilon_posterior_sum=c()
  for (matrix in mean_upsilon_posterior) {
    matrix_sum =sum(matrix)
    mean_upsilon_posterior_sum =c(mean_upsilon_posterior_sum, matrix_sum)
  }
  #updating latent traits
  for(j in 1:m){
    mean_upsilon_posterior_sum[j]=mean_upsilon_posterior_sum[j]/(n_obs[j,]+1/sigmasq_upsilon_post[j])
    upsilon_posterior[j]=rnorm(1,mean_upsilon_posterior_sum[j],
                             sqrt(var_upsilon_posterior[j]))
  }
  #As I have only one latent trait upsilon for each region
  upsilon_rep=list()
  for(j in 1:m){
    upsilon_rep[[j]]=rep(upsilon_posterior[j],n_obs[j,])
  }
  ###to convert the list into matrix 
  upsilon_rep=lapply(upsilon_rep, function(x) matrix(x, nrow = length(x), ncol = 1))
  a=b=list()
  a[[1]]=matrix(NA,length(z_tig),1)
  a[[2]]=matrix(NA,length(z_afa),1)
  a[[3]]=matrix(NA,length(z_amh),1)
  a[[4]]=matrix(NA,length(z_oro),1)
  a[[5]]=matrix(NA,length(z_som),1)
  a[[6]]=matrix(NA,length(z_ben),1)
  a[[7]]=matrix(NA,length(z_snn),1)
  a[[8]]=matrix(NA,length(z_gam),1)
  a[[9]]=matrix(NA,length(z_har),1)
  a[[10]]=matrix(NA,length(z_add),1)
  a[[11]]=matrix(NA,length(z_dir),1)
  b[[1]]=matrix(NA,length(z_tig),1)
  b[[2]]=matrix(NA,length(z_afa),1)
  b[[3]]=matrix(NA,length(z_amh),1)
  b[[4]]=matrix(NA,length(z_oro),1)
  b[[5]]=matrix(NA,length(z_som),1)
  b[[6]]=matrix(NA,length(z_ben),1)
  b[[7]]=matrix(NA,length(z_snn),1)
  b[[8]]=matrix(NA,length(z_gam),1)
  b[[9]]=matrix(NA,length(z_har),1)
  b[[10]]=matrix(NA,length(z_add),1)
  b[[11]]=matrix(NA,length(z_dir),1)
  for(t in 1:length(z_tig)){
    zt=z[[1]]
    a[[1]][t,]=max(zt[Y[[1]][,1]<Y[[1]][t,]])
    b[[1]][t,]=min(zt[Y[[1]][t,]<Y[[1]][,1]])
    z[[1]][t,]=rtnorm(1,X[[1]][t,]%*%beta+upsilon_rep[[1]][t,],1,a[[1]][t,],b[[1]][t,])
  }
  for(t in 1:length(z_afa)){
    zt=z[[2]]
    a[[2]][t,]=max(zt[Y[[2]][,1]<Y[[2]][t,]])
    b[[2]][t,]=min(zt[Y[[2]][t,]<Y[[2]][,1]])
    z[[2]][t,]=rtnorm(1,X[[2]][t,]%*%beta+upsilon_rep[[2]][t,],1,a[[2]][t,],b[[2]][t,])
  }
  for(t in 1:length(z_amh)){
    zt=z[[3]]
    a[[3]][t,]=max(zt[Y[[3]][,1]<Y[[3]][t,]])
    b[[3]][t,]=min(zt[Y[[3]][t,]<Y[[3]][,1]])
    z[[3]][t,]=rtnorm(1,X[[3]][t,]%*%beta+upsilon_rep[[3]][t,],1,a[[3]][t,],b[[3]][t,])
  }
  for(t in 1:length(z_oro)){
    zt=z[[4]]
    a[[4]][t,]=max(zt[Y[[4]][,1]<Y[[4]][t,]])
    b[[4]][t,]=min(zt[Y[[4]][t,]<Y[[4]][,1]])
    z[[4]][t,]=rtnorm(1,X[[4]][t,]%*%beta+upsilon_rep[[4]][t,],1,a[[4]][t,],b[[4]][t,])
  } 
  for(t in 1:length(z_som)){
    zt=z[[5]]
    a[[5]][t,]=max(zt[Y[[5]][,1]<Y[[5]][t,]])
    b[[5]][t,]=min(zt[Y[[5]][t,]<Y[[5]][,1]])
    z[[5]][t,]=rtnorm(1,X[[5]][t,]%*%beta+upsilon_rep[[5]][t,],1,a[[5]][t,],b[[5]][t,])
  }
  for(t in 1:length(z_ben)){
    zt=z[[6]]
    a[[6]][t,]=max(zt[Y[[6]][,1]<Y[[6]][t,]])
    b[[6]][t,]=min(zt[Y[[6]][t,]<Y[[6]][,1]])
    z[[6]][t,]=rtnorm(1,X[[6]][t,]%*%beta+upsilon_rep[[6]][t,],1,a[[6]][t,],b[[6]][t,])
  }
  for(t in 1:length(z_snn)){
    zt=z[[7]]
    a[[7]][t,]=max(zt[Y[[7]][,1]<Y[[7]][t,]])
    b[[7]][t,]=min(zt[Y[[7]][t,]<Y[[7]][,1]])
    z[[7]][t,]=rtnorm(1,X[[7]][t,]%*%beta+upsilon_rep[[7]][t,],1,a[[7]][t,],b[[7]][t,])
  }
  for(t in 1:length(z_gam)){
    zt=z[[8]]
    a[[8]][t,]=max(zt[Y[[8]][,1]<Y[[8]][t,]])
    b[[8]][t,]=min(zt[Y[[8]][t,]<Y[[8]][,1]])
    z[[8]][t,]=rtnorm(1,X[[8]][t,]%*%beta+upsilon_rep[[8]][t,],1,a[[8]][t,],b[[8]][t,])
  }
  for(t in 1:length(z_har)){
    zt=z[[9]]
    a[[9]][t,]=max(zt[Y[[9]][,1]<Y[[9]][t,]])
    b[[9]][t,]=min(zt[Y[[9]][t,]<Y[[9]][,1]])
    z[[9]][t,]=rtnorm(1,X[[9]][t,]%*%beta+upsilon_rep[[9]][t,],1,a[[9]][t,],b[[9]][t,])
    
  }
  for(t in 1:length(z_add)){
    zt=z[[10]]
    a[[10]][t,]=max(zt[Y[[10]][,1]<Y[[10]][t,]])
    b[[10]][t,]=min(zt[Y[[10]][t,]<Y[[10]][,1]])
    z[[10]][t,]=rtnorm(1,X[[10]][t,]%*%beta+upsilon_rep[[10]][t,],1,a[[10]][t,],b[[10]][t,])
  }
  for(t in 1:length(z_dir)){
    zt=z[[11]]
    a[[11]][t,]=max(zt[Y[[11]][,1]<Y[[11]][t,]])
    b[[11]][t,]=min(zt[Y[[11]][t,]<Y[[11]][,1]])
    z[[11]][t,]=rtnorm(1,X[[11]][t,]%*%beta+upsilon_rep[[11]][t,],1,a[[11]][t,],b[[11]][t,])
  }
  muj1=array(NA,c(p,p,m))
  muj1_sum=matrix(NA,p,p)
  muj2=array(NA,c(p,1,m))
  muj2_sum=matrix(NA,p,1)
  muj_final=matrix(NA,p,1)
  varj_final=matrix(NA,p,p)
  for(j in 1:m){
    muj1[,,j]=t(X[[j]])%*%X[[j]]
    #before summing over j
    muj1_sum=apply(muj1,c(1:2),sum)#summing over j
    muj2[,,j]=t(X[[j]])%*%(z[[j]]-upsilon_rep[[j]])
    muj2_sum=apply(muj2,c(1:2),sum)#summing over j
    muj_final=solve(muj1_sum)%*%(muj2_sum)
    varj_final=solve(muj1_sum)
  }
  BETA=mvrnorm(n=1,mu=muj_final,
               varj_final)#sampling from a multivariate Normal
  #with updated mean and variance
  
  BETA_POST=rbind(BETA_POST,BETA)
  upsilon_POST=rbind(upsilon_POST,upsilon_posterior)
  sigmasq_upsilon_POST=rbind(sigmasq_upsilon_POST,sigmasq_upsilon_post)
  
  print(paste("iteration",i))
  z_all=rbind(z_all,z)
}
library(coda)
library(bayesplot)
library(stableGR)
Beta_gibbis_burn<-BETA_POST[seq(1, nrow(BETA_POST[5001:30000,]),5),]
upsilon_gibbis_burn<-upsilon_POST[seq(1, nrow(upsilon_POST[5001:30000,]),5),]
mean(Beta_gibbis_burn[,1])
mean(Beta_gibbis_burn[,2])
mean(Beta_gibbis_burn[,3])
mean(Beta_gibbis_burn[,4])
mean(Beta_gibbis_burn[,5])
sd(Beta_gibbis_burn[,1])
sd(Beta_gibbis_burn[,2])
sd(Beta_gibbis_burn[,3])
sd(Beta_gibbis_burn[,4])
sd(Beta_gibbis_burn[,5])
Beta_mcmc_object=as.mcmc(Beta_gibbis_burn)
Beta_credible_interval=HPDinterval(Beta_mcmc_object , prob = 0.95)
PSRF_beta_test=stable.GR(Beta_gibbis_burn)
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
upsilon_mcmc_object=as.mcmc(upsilon_gibbis_burnn)
upsilon_credible_interval=HPDinterval(upsilon_mcmc_object , prob = 0.95)
PSRF_upsilon_test=stable.GR(upsilon_gibbis_burn)
#####Regional variation 
library(reshape2)
library(viridis)
library(ggplot2)
new_colnames <- c("1","2","3","4","5","6","7","8","9","10","11")
colnames(g) <- new_colnames
g=data.frame(sigma_gibbis_burn)
data <- var(g[sapply(g,is.numeric)])
data1 <- melt(data)
data2=read.csv(('D:\\Myfiles\\Mythesisresearch\\Variation.csv'))
data2=data.frame(data2)
colnames(data2)=c("Var1","Var2", "Variations")
x_labels <- c("Tigray", "Afar", "Amhara", "Oromia", "Somale", "Benishangul", "SNN", "Gambela", "Harari", "Addis Abeba", "Dire Dawa")
y_labels <- c("Tigray", "Afar", "Amhara", "Oromia", "Somale", "Benishangul", "SNN", "Gambela", "Harari", "Addis Abeba", "Dire Dawa")
ggplot(data2, aes(x = Var1, y = Var2, fill = Variations)) +
  geom_tile() +scale_fill_viridis(discrete = FALSE)+geom_tile()+
  labs(title = "",
       x = "Regions",
       y = "Regions")+scale_x_discrete(labels = x_labels)+
  scale_y_discrete(labels = y_labels)
