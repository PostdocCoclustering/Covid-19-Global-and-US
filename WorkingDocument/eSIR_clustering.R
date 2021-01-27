
# This is a piece from the eSIR function that implements the state-space 
# model. Specifically it defines the whole model and feeds it to JAGS.

# This implements the Runge-Kutta approximation (Appendix A in the eSIR paper)
model1.string <- paste0("\n  
model{\n     
for(t in 2:(T_prime+1)){\n      

       Km[t-1,1] <- -beta*pi[t-1]*theta[t-1,1]*theta[t-1,2]\n      
       Km[t-1,9] <- gamma*theta[t-1,2]\n       
       Km[t-1,5] <- -Km[t-1,1]-Km[t-1,9]\n\n      
       Km[t-1,2] <- -beta*pi[t-1]*(theta[t-1,1]+0.5*Km[t-1,1])*(theta[t-1,2]+0.5*Km[t-1,5])\n      
       Km[t-1,10] <- gamma*(theta[t-1,2]+0.5*Km[t-1,5])\n       
       Km[t-1,6] <- -Km[t-1,2]-Km[t-1,10]\n\n       
       Km[t-1,3] <- -beta*pi[t-1]*(theta[t-1,1]+0.5*Km[t-1,2])*(theta[t-1,2]+0.5*Km[t-1,6])\n       
       Km[t-1,11] <- gamma*(theta[t-1,2]+0.5*Km[t-1,6])\n       
       Km[t-1,7] <- -Km[t-1,3]-Km[t-1,11]\n\n       
       Km[t-1,4] <- -beta*pi[t-1]*(theta[t-1,1]+Km[t-1,3])*(theta[t-1,2]+Km[t-1,7])\n       
       Km[t-1,12] <- gamma*(theta[t-1,2]+Km[t-1,7])\n       
       Km[t-1,8] <- -Km[t-1,4]-Km[t-1,12]\n\n       
       
       alpha[t-1,1] <- theta[t-1,1]+(Km[t-1,1]+2*Km[t-1,2]+2*Km[t-1,3]+Km[t-1,4])/6\n       
       alpha[t-1,2] <- theta[t-1,2]+(Km[t-1,5]+2*Km[t-1,6]+2*Km[t-1,7]+Km[t-1,8])/6\n       
       alpha[t-1,3] <- theta[t-1,3]+(Km[t-1,9]+2*Km[t-1,10]+2*Km[t-1,11]+Km[t-1,12])/6\n\n       
       
       theta[t,1:3] ~ ddirch(k*alpha[t-1,1:3])\n # Cluster must be indicated      
       
       Y[t-1] ~ dbeta(lambdaY*theta[t,2],lambdaY*(1-theta[t,2]))\n       
       R[t-1] ~ dbeta(lambdaR*theta[t,3],lambdaR*(1-theta[t,3]))\n     
       
       }\n 
              
       theta0[1:3] <- c(", 1 - Y[1] - R[1],",
                        ", Y[1],",
                        ", R[1],")\n    
       
       theta[1,1:3] ~ ddirch(theta0[1:3])\n
       
       gamma ~  dlnorm(", lognorm_gamma_parm$mu,",
                        ", 1/lognorm_gamma_parm$var,")\n
                        
       R0 ~ dlnorm(", lognorm_R0_parm$mu, ",",1/lognorm_R0_parm$var, ")\n
                        
       beta <- R0*gamma\n    
       k ~  dgamma(2,0.0001)\n    
       lambdaY ~ dgamma(2,0.0001)\n    
       lambdaR ~ dgamma(2,0.0001)\n  
       
       }\n")

model.spec <- textConnection(model1.string)
nchain = 4
nadapt = 10000
beta0 = 0.2586
gamma0 = 0.0821
R0 = beta0/gamma0

posterior <- jags.model(model.spec, data = 
                          list(Y = Y, R = R, T_prime = T_prime),
                        n.chains = nchain, n.adapt = nadapt)

update(posterior, nburnin)

jags_sample <- jags.samples(posterior, c("theta", "gamma", "R0", "beta",
                                         "Y", "lambdaY", "lambdaR", "k"),
                            n.iter = M, thin = thn)

  