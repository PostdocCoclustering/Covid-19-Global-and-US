
library(rjags)

### Step function of pi(t)
change_time <- c("01/23/2020", "02/04/2020", "02/08/2020")
pi0<- c(1.0, 0.9, 0.5, 0.1)


### Needed functions
lognorm.parm <- function(mu0,var0){
  var <- log(var0 / mu0^2 + 1)
  mu <- log(mu0) - var / 2
  list(mu = mu, var = var)
}


gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

### Some default values 
#pi0 = NULL 
change_time = NULL 
exponential = FALSE 
lambda0 = NULL
begin_str = "01/13/2020"
T_fin = 200
nchain = 4
nadapt = 10000
M = 500
thn = 10
nburnin = 200
dic = FALSE
death_in_R = 0.02
beta0 = 0.2586
gamma0 = 0.0821
R0 = beta0/gamma0
gamma0_sd = 0.1
R0_sd = 1
file_add = character(0)
add_death = FALSEsave_files = FALSE 
save_mcmc = FALSE
save_plot_data = FALSE
eps = 1e-10
T_prime <- length(Y)


model1.string <- paste0("
  model{
     for(t in 2:(T_prime+1)){
       Km[t-1,1] <- -beta*pi[t-1]*theta[t-1,1]*theta[t-1,2]
       Km[t-1,9] <- gamma*theta[t-1,2]
       Km[t-1,5] <- -Km[t-1,1]-Km[t-1,9]
       Km[t-1,2] <- -beta*pi[t-1]*(theta[t-1,1]+0.5*Km[t-1,1])*(theta[t-1,2]+0.5*Km[t-1,5])
       Km[t-1,10] <- gamma*(theta[t-1,2]+0.5*Km[t-1,5])
       Km[t-1,6] <- -Km[t-1,2]-Km[t-1,10]
       Km[t-1,3] <- -beta*pi[t-1]*(theta[t-1,1]+0.5*Km[t-1,2])*(theta[t-1,2]+0.5*Km[t-1,6])
       Km[t-1,11] <- gamma*(theta[t-1,2]+0.5*Km[t-1,6])
       Km[t-1,7] <- -Km[t-1,3]-Km[t-1,11]
       Km[t-1,4] <- -beta*pi[t-1]*(theta[t-1,1]+Km[t-1,3])*(theta[t-1,2]+Km[t-1,7])
       Km[t-1,12] <- gamma*(theta[t-1,2]+Km[t-1,7])
       Km[t-1,8] <- -Km[t-1,4]-Km[t-1,12]
       
       alpha[t-1,1] <- theta[t-1,1]+(Km[t-1,1]+2*Km[t-1,2]+2*Km[t-1,3]+Km[t-1,4])/6
       alpha[t-1,2] <- theta[t-1,2]+(Km[t-1,5]+2*Km[t-1,6]+2*Km[t-1,7]+Km[t-1,8])/6
       alpha[t-1,3] <- theta[t-1,3]+(Km[t-1,9]+2*Km[t-1,10]+2*Km[t-1,11]+Km[t-1,12])/6
       
       theta[t,1:3] <- muOfClust[ clust[t], 1:3 ]
       clust[t] ~ dcat(pclust[1:Nclust]) 
       
       R[t-1] ~ dbeta(lambdaR*theta[t,3],lambdaR*(1-theta[t,3]))
       Y[t-1] ~ dbeta(lambdaY*theta[t,2],lambdaY*(1-theta[t,2]))
       
     }
      
      # Priors
      
      for(clustIdx in 1:Nclust){
      muOfClust[clustIdx, 1:3] ~ ddirch(theta0[1:3])
      }
      
      pclust[1:Nclust] ~ ddirch(onesRepNclust)

    theta0[1:3]<-c(",1-Y[1]-R[1],",",Y[1],",", R[1],")
    theta[1,1:3] ~ ddirch(theta0[1:3])
    gamma ~  dlnorm(", lognorm_gamma_parm$mu, ",", 1 / lognorm_gamma_parm$var,")
    R0 ~ dlnorm(", lognorm_R0_parm$mu, ",", 1 / lognorm_R0_parm$var,")
    beta <- R0*gamma
    k ~  dgamma(2,0.0001)
    lambdaY ~ dgamma(2,0.0001)
    lambdaR ~ dgamma(2,0.0001)
  
  }
")


### Data 

set.seed(20192020)

# Hubei province data Jan13 -> Feb 11
# cumulative number of infected

# Number of cases (I)
NI_complete <- c( 41, 41, 41, 45, 62, 131, 200, 270, 375, 444, 549, 729,
                  1052, 1423, 2714, 3554, 4903, 5806, 7153, 9074, 11177,
                  13522,16678,19665,22112,24953,27100,29631,31728,33366) 

# Number of removed
RI_complete <- c(1, 1, 7, 10, 14, 20, 25, 31, 34, 45, 55, 71, 94, 
                 121, 152, 213, 252, 345, 417, 561, 650, 811, 
                 1017, 1261, 1485, 1917, 2260, 2725,3284,3754)
N <- 58.5e6
R <- RI_complete / N
Y <- NI_complete / N - R #Jan13->Feb 11

### Y and R are observed proportions of infected and removed 

Nclust = 2
clust = rep(NA,N) 
clust[which.min(Y)]=1 # smallest value assigned to cluster 1
clust[which.max(Y)]=2 

model.spec <- textConnection(model1.string)

posterior <- jags.model(
  model.spec,
  data = list(
    'Y' = Y,
    'R' = R,
    'T_prime' = T_prime,
    'pi' = pi,
    'Nclust' = Nclust,
    'clust' = clust,
    onesRepNclust = rep(1, 2)
  ),
  n.chains = nchain,
  n.adapt = nadapt
)


update(posterior, nburnin) #burn-in

jags_sample <- jags.samples(
  posterior,
  c('theta','gamma','R0','beta','Y','lambdaY','lambdaR','k'),
  n.iter = M,
  thin = thn
)


