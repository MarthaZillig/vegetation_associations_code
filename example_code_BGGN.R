

###############################################################################################################

# This is example code for models discussed in Zillig and Fleishman (2022) "TITLE OF PAPER HERE" 

#This code is an example of the final model used for all 19 bird species tested. Here, we show the code for Blue-gray Gnatcatcher. 
#All bird data is available for download on the Forest Service Research Data Archive: https://www.fs.usda.gov/rds/archive/Catalog?authorid=RDS183

################################################################################################################



library(tidyverse)
library(jagsUI)
library(boot)

#Example dataframe for Blue-gray Gnatcatcher (BGGN)

load("output/bggn_example_df.RDS")
load("observers_bggn.RDS")

#Occupancy matrix for the three visits to each site 

bggn_matrix <- bggn_df2 %>% 
  dplyr::select(v1, v2, v3) %>% 
  as.matrix()

#Number of unique survey points 
n_point <- 660

bggn_df2$point <- as.factor(bggn_df2$point)
bggn_df2$point<- as.numeric(bggn_df2$point) #change to a factor for the random effect 

point_matrix <- bggn_df2 %>% dplyr::select(point_year, point) %>% arrange(point_year) %>% dplyr::select(-point_year) %>% as.matrix()

point_vector <- point_matrix %>% as.vector()

n_region <- 5 #number of subregions 

bggn_df2$region <- as.factor(bggn_df2$region)
bggn_df2$region <- as.numeric(bggn_df2$region)

#Observers (random effect)

obs_num <- obs_bggn %>% 
  ungroup() %>% 
  dplyr::select(v1, v2, v3) %>% 
  as.matrix()

n_obs <- 52


#Scaling Variaibles (fixed effects)

bggn_df2$year <- scale(bggn_df2$year) %>% as.vector()

#Functional groups 
bggn_low_sh <- scale(bggn_df2$low_shrub) %>% as.vector() 
bggn_arid_c <- scale(bggn_df2$arid_conifer) %>% as.vector()
bggn_r_tree <- scale(bggn_df2$riparian_trees) %>% as.vector()
bggn_r_sh <- scale(bggn_df2$riparian_shrubs) %>% as.vector()
bggn_con <- scale(bggn_df2$conifers) %>% as.vector() 

#Plant Species 
bggn_sage <- scale(bggn_df2$sagebrush) %>% as.vector() 
bggn_rowo <- scale(bggn_df2$ROWO) %>% as.vector()
bggn_putr <- scale(bggn_df2$PUTR) %>% as.vector()
bggn_pran <- scale(bggn_df2$PRAN) %>% as.vector()
bggn_rabbit <- scale(bggn_df2$rabbitbrush) %>% as.vector()
bggn_potr <- scale(bggn_df2$POTR) %>% as.vector()
bggn_poan <- scale(bggn_df2$POAN) %>% as.vector()
bggn_pimo <- scale(bggn_df2$PIMO) %>% as.vector()
bggn_juniper <- scale(bggn_df2$juniper) %>% as.vector()
bggn_cele <- scale(bggn_df2$CELE) %>% as.vector()
bggn_salix <- scale(bggn_df2$Salix) %>% as.vector()
bggn_pifl <- scale(bggn_df2$PIFL) %>% as.vector()
bggn_fir <- scale(bggn_df2$fir) %>% as.vector() 


#Setting Up Indicator variables 

Q <- 18 #number of indicator variables 

cov_bggn <- cbind(bggn_low_sh, bggn_arid_c, bggn_r_tree, bggn_r_sh, bggn_con, bggn_sage, bggn_rowo, bggn_putr, bggn_pran, bggn_rabbit, bggn_potr, bggn_poan, bggn_pimo, bggn_juniper, bggn_cele, bggn_salix, bggn_pifl, bggn_fir)

cov_bggn2 <- as.data.frame(cov_bggn)

#replace all NAs in veg data with 0

cov_bggn2 <- cov_bggn2 %>% mutate_if(is.numeric, ~replace(., is.na(.), 0))

cov_bggn <- cov_bggn2 %>% as.matrix()


##MODEL

a.vec <- rep(1, n_obs)

C <- bggn_matrix

data <- list(C = C, M =nrow(C), J = ncol(C), facRegion = bggn_df2$region, X = cov_bggn, Q = Q, year = bggn_df2$year, n_point = n_point, point = point_vector, n_obs = n_obs, obs_num = obs_num, a.vec = a.vec, n_region = n_region)

Zst <- apply(C, 1, max, na.rm =TRUE)

inits <- function()list(Z = Zst, alpha0 = rnorm(n_region, 0 ,2), beta0 = dnorm(1))

params <- c("beta0",  "alpha1", "alpha0", "betaT", "theta1", "theta2", "bp.value", "bp.value.n")

#JAGS MODEL

cat("

model { 

  for(i in 1:n_point){
    theta2[i] ~ dnorm(0, sigma.c)
    }
  
  for(i in 1:n_region){
    alpha0[i] ~ dnorm(mu.alpha, tau.alpha)
  }
  
  for(i in 1:n_obs){
    theta1[i] ~ dnorm(0, tau.o)
  }
  
  for(i in 1:M){ #missing observer values 
    for (a in 1:3){
      obs_num[i,a] ~ dcat(b[1:n_obs])  
    }
  }

  
  #Indicator variable priors
    
  for (q in 1:Q) { #Q = number of variables (18)
    for(k in 1:5){
    betaT[k,q]~dnorm(0,tauT)
    
  }
}

  tauT ~ dgamma(0.1,0.1)
  
  for(i in 1:M) {
 
    Z[i] ~ dbern(lambda[i])  #occupancy submodel 
    
    logit(lambda[i]) <- alpha0[facRegion[i]] + inprod(X[i,], betaT[facRegion[i],]) + alpha1*year[i] + theta2[point[i]]
    
    n.new[i]~ dbern(lambda[i]) #calculation of bayes p value
    res.n[i] <- (sqrt(Z[i]) - sqrt(lambda[i]))^2
    res.n.new[i] <- (sqrt(n.new[i]) - sqrt(lambda[i]))^2
    
  
 for (j in 1:J) {
 
      C[i,j] ~ dbin(p[i,j], Z[i]) #detection submodel 
      
      logit(p[i,j]) <- beta0 + theta1[obs_num[i,j]] 
      
      y.new[i,j]~ dbin(p[i,j], Z[i]) #calculation of bayes p values 
      res[i,j] <- (sqrt(C[i,j])-sqrt(p[i,j]*Z[i]))^2
      res.new[i,j] <- (sqrt(y.new[i,j])-sqrt(p[i,j]*Z[i]))^2
      
   
    }
  }
  
  total.res <- sum(res[,])
  total.resnew <- sum(res.new[,])
  bp.value <- total.res>total.resnew
  
  
  total.res.n <- sum(res.n[])
  total.resnew.n <- sum(res.n.new[])
  bp.value.n <- total.res.n>total.resnew.n
  
  
  #Priors 
  
  beta0 ~ dlogis(0,1)
  alpha1 ~ dnorm(0,1)
  sigma.c ~ dunif(0,1)
  tau.o ~ dunif(0,5)
  mu.alpha <- logit(mean.theta)
  mean.theta ~ dunif(0,1)
  tau.alpha <- pow(sd.alpha, -2)
  sd.alpha ~ dunif(0,10)
  b[1:n_obs] ~ ddirch(a.vec)
  

}
", file = "model_veg.txt")

library(jagsUI)


#This will take about 4 hours to run, can shorten the iterations to see example output

out_bggn <- jags(data = data, model.file = "model_veg.txt", inits = inits, parameters.to.save = params, n.chains = 3, n.iter = 50000, n.burnin = 10000, n.adapt = 10000, n.thin = 10)








