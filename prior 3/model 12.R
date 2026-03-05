
beta.ch.shape <- 1
beta.ch.scale <- 1.5
beta.cc.shape <- 1
beta.cc.scale <- 1.5
beta.ic.shape <- 1
beta.ic.scale <- 1.5

sigma.shape <- 1
sigma.scale <- 1.5
alpha.shape <- 1
alpha.scale <- 1.5

model.code = nimbleCode({
    
    S[1] <- S0
    Ch[1] <- Ch.0
    Ih[1] <- Ih.0
    Cc[1] <- Cc.0
    Ic[1] <- Ic.0
    A[1] <- A0
    N[1] <- S[1] + Ch[1] + Ih[1] + Cc[1] + Ic[1] + A[1] 
    
    # likelihood
    
    for(t in 1: tau) {
        
        Ch.star.1[t] ~ dbin(1-exp(-(beta.ch)*Ch[t]/N[t]-(beta.cc)*Cc[t]/N[t]-(beta.ic)*Ic[t]/N[t] - sigma), S[t]) 
        Ih.star.1[t] ~ dbin(1-exp(-alpha), Ch[t])
        
        S[t+1] <- S[t] + Admit[t] - Dis[t] - Ch.star.1[t]
        Ch[t+1] <- Ch[t] + Ch.star.1[t] - Ih.star.1[t] - rbinom(1, Ch[t], 1-exp(-rho.1)) 
        Ih[t+1] <- Ih[t] + Ih.star.1[t] - rbinom(1, Ih[t], 1-exp(-rho.2)) 
        Cc[t+1] <- Cc[t] + Cc.star.1[t] - rbinom(1, Cc[t], 1-exp(-rho.3)) 
        Ic[t+1] <- Ic[t] + Ic.star.1[t] - rbinom(1, Ic[t], 1-exp(-rho.4)) 
        A[t+1] <- A[t] + rbinom(1, Ch[t], 1-exp(-rho.1)) + rbinom(1, Ih[t], 1-exp(-rho.2)) +  rbinom(1, Cc[t], 1-exp(-rho.3)) + rbinom(1, Ic[t], 1-exp(-rho.4)) 
        N[t+1] <- S[t+1] + Ch[t+1] + Ih[t+1] + Cc[t+1] + Ic[t+1] + A[t+1]
        
    }
    
    # priors
    beta.ch.shape <- 1
    beta.ch.scale <- 1.5
    beta.cc.shape <- 1
    beta.cc.scale <- 1.5
    beta.ic.shape <- 1
    beta.ic.scale <- 1.5
    
    sigma.shape <- 1
    sigma.scale <- 1.5
    alpha.shape <- 1
    alpha.scale <- 1.5
    
    beta.ch ~ dgamma(shape = beta.ch.shape, scale = beta.ch.scale)
    beta.cc ~ dgamma(shape = beta.cc.shape, scale = beta.cc.scale)
    beta.ic ~ dgamma(shape = beta.ic.shape, scale = beta.ic.scale)
    sigma ~ dgamma(shape = sigma.shape, scale = sigma.scale)
    alpha ~ dgamma(shape = alpha.shape, scale = alpha.scale)
})


# Model fitting

tau = length(MRSA$hospital.infections)

data.list = list(Ch.star.1 = MRSA$hospital.colonization,
                 Ih.star.1 = MRSA$hospital.infections,
                 Cc.star.1 = MRSA$community.colonization,
                 Ic.star.1 = MRSA$community.infections,
                 Admit = MRSA$admissions,
                 Dis = MRSA$discharges
)

tau = length(data.list$Ih.star.1)

constants.list = list(S0 = 3048,  
                      Ch.0 = MRSA$hospital.colonization[1], 
                      Ih.0 = MRSA$hospital.infections[1],  
                      Cc.0 = MRSA$community.colonization[1],  
                      Ic.0 = MRSA$community.infections[1], 
                      A0 = 0, 
                      tau = tau)

inits.list = list(
    beta.ch = rgamma(1, shape = beta.ch.shape, scale = beta.ch.scale),
    beta.cc = rgamma(1, shape = beta.cc.shape, scale = beta.cc.scale),
    beta.ic = rgamma(1, shape = beta.ic.shape, scale = beta.ic.scale),
    sigma = rgamma(1, shape = sigma.shape, scale = sigma.scale),
    alpha = rgamma(1, shape = alpha.shape, scale = alpha.scale),
    rho.1 = 1.3,
    rho.2 = 1.3,
    rho.3 = 10,
    rho.4 = 10
)

mrsa.Model <- nimbleModel(model.code, 
                          constants = constants.list,
                          data = data.list,
                          inits = inits.list)

myConfig <- configureMCMC(mrsa.Model, enableWAIC = TRUE)
myConfig
myMCMC <- buildMCMC(myConfig)
compiled <- compileNimble(mrsa.Model, myMCMC)
output <- runMCMC(compiled$myMCMC, WAIC = TRUE, niter = 60000, nburnin = 10000, setSeed = 1)
samples = output$samples
output$WAIC

# Trace plots for model parameters

plot(samples[,"beta.ch"], type = 'l')
abline(h = mean(samples[,"beta.ch"]), col="red", lwd = 2)
hist(samples[,"beta.ch"], main = "Posterior (beta.ch)")
curve(dgamma(x, shape = beta.ch.shape, scale = beta.ch.scale, log = FALSE), 
      from = 0, to = 10, main = "Prior of beta.ch")

plot(samples[,"beta.cc"], type = 'l')
abline(h = mean(samples[,"beta.cc"]), col="red", lwd = 2)
hist(samples[,"beta.cc"], main = "Posterior (beta.cc)")
curve(dgamma(x, shape = beta.cc.shape, scale = beta.cc.scale, log = FALSE), 
      from = 0, to = 10, main = "Prior of beta.cc")

plot(samples[,"beta.ic"], type = 'l')
abline(h = mean(samples[,"beta.ic"]), col="red", lwd = 2)
hist(samples[,"beta.ic"], main = "Posterior (beta.ic)")
curve(dgamma(x, shape = beta.ic.shape, scale = beta.ic.scale, log = FALSE), 
      from = 0, to = 10, main = "Prior of beta.ic")

plot(samples[,"sigma"], type = 'l')
abline(h = mean(samples[,"sigma"]), col="red", lwd = 2)
hist(samples[,"sigma"], main = "Posterior (sigma)")
curve(dgamma(x, shape = sigma.shape, scale = sigma.scale, log = FALSE), 
      from = 0, to = 10, main = "Prior of sigma")

plot(samples[,"alpha"], type = 'l')
abline(h = mean(samples[,"alpha"]), col="red", lwd = 2)
hist(samples[,"alpha"], main = "Posterior (alpha)")
curve(dgamma(x, shape = alpha.shape, scale = alpha.scale, log = FALSE), 
      from = 0, to = 10, main = "Prior of alpha")
