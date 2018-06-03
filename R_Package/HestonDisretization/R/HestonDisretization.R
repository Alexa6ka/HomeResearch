install.packages("NMOF") # Numerical Methods and Optimization in Finance
library(NMOF)
library(dplyr)

# FULL TRUNCATION EULER DISCRETISATION ####


# TRUE CALCULATIONS
calcHMED.1 <- function(N_MC, N_AVI, S0, K, V0, R, TIME, RHO, KAPPA, THETA, SIGMA){
   dt <- TIME/N_AVI
   payoffs <- c()
   for(j in 1:N_MC){
      V <- V0
      S <- S0
      for(i in 1:N_AVI){
         Z1 <- rnorm(1, mean = 0, sd = 1)
         Z2 <- rnorm(1, mean = 0, sd = 1)
         #full truncation Euler discretisation
         S <- S + R * S * dt + sqrt(max(V, 0)) * S * (RHO * Z1 + sqrt(1 - RHO^2) * Z2) * sqrt(dt)
         V <- V + KAPPA * (THETA - max(V, 0)) * dt + SIGMA * sqrt(max(V, 0)) * Z1 * sqrt(dt)
      }
      payoffs[j] <- max(0, (S - K))
   }
   OptionPrice <- mean(payoffs) * exp(-R * TIME) 
   #cat("dt:", dt, "N_MC:", N_MC, "N_AVI:", N_AVI, "\n")
   return(OptionPrice)
}

calcHMED.2 <- function(N_MC, N_AVI, S0, K, V0, R, TIME, RHO, KAPPA, THETA, SIGMA){
   dt <- TIME/N_AVI
   S_array <- c()
   payoffs <- c()
   for(j in 1:N_MC){
      V <- V0
      S <- S0
      for(i in 1:N_AVI){
         Z1 <- rnorm(1, mean = 0, sd = 1)
         Z2 <- rnorm(1, mean = 0, sd = 1)
         #full truncation
         S <- S + R * S * dt + sqrt(max(V, 0)) * S * (RHO * Z1 + sqrt(1 - RHO^2) * Z2) * sqrt(dt)
         V <- V + KAPPA * (THETA - max(V, 0)) * dt + SIGMA * sqrt(max(V, 0)) * Z1 * sqrt(dt)
      }
      S_array[j] <- S
      #payoffs[j] <- max( (S_array[j] - K), 0)
      payoffs[j] <- max((S - K), 0)
   }
   #cat("N_MC:", N_MC,       "payoffs len:", length(payoffs),  '\n')
   #df <- data.frame(S0, S_array, K, payoffs)
   #print(df)
   #print(mean(payoffs) * exp(-R * TIME))
   
   OptionPrice <- mean(payoffs) * exp(-R * TIME)
   #print(OptionPrice)
   #cat("dt:", dt, "N_MC:", N_MC, "N_AVI:", N_AVI, "\n")
   return(data.frame(S_array, payoffs))
}


N_MC = 2000
N_AVI = 32*6
S0 = 100
K = 70
V0 = 0.0225
R = 0.04
TIME = 6
RHO = -0.5
KAPPA = 2
THETA = 0.04
SIGMA = 0.3





calcHMED.2(N_MC = 2000, N_AVI = 32*6,
           S0 = 100,
           K = 70,
           V0 = 0.0225,
           R = 0.04,
           TIME = 6,
           RHO = -0.5,
           KAPPA = 2,
           THETA = 0.04,
           SIGMA = 0.3)

calcHMED.1(N_MC = 2000, N_AVI = 32*6,
           S0 = 100,
           K = 70,
           V0 = 0.0225,
           R = 0.04,
           TIME = 6,
           RHO = -0.5,
           KAPPA = 2,
           THETA = 0.04,
           SIGMA = 0.3)

callHestoncf(S = 100,
             X = 70,
             v0 = 0.0225,
             r = 0.04,
             tau = 6,
             rho = -0.5,
             k = 2, 
             vT = 0.04,
             sigma = 0.3,
             q = 0)

#########################################



# calculate Heston Model Milstein discretisation TRUE
calcHMMD.1 <- function(N_MC, N_AVI, S0, K, V0, R, TIME, RHO, KAPPA, THETA, SIGMA){
   dt <- TIME/N_AVI
   payoffs <- c()
   for(j in 1:N_MC){
      V <- V0
      S <- S0
      for(i in 1:N_AVI){
         Zv <- rnorm(1, mean = 0, sd = 1)
         Zs <- rnorm(1, mean = 0, sd = 1)
         Zs <- RHO * Zv + sqrt(1 - RHO^2) * Zs
         
         S <- S * exp( (R - 0.5 * max(V, 0)) * dt + sqrt(max(V, 0) * dt) * Zs)
         V <- (sqrt(max(V, 0)) + 0.5 * SIGMA * sqrt(dt) * Zv)^2 + KAPPA * (THETA - max(V, 0)) * dt - 0.25 * SIGMA^2 * dt
      }
      payoffs[j] <- max(0, (S - K))
   }
   OptionPrice <- mean(payoffs) * exp(-R * TIME) 
   #cat("dt:", dt, "N_MC:", N_MC, "N_AVI:", N_AVI, "\n")
   return(OptionPrice)
}

calcHMMD.2 <- function(N_MC, N_AVI, S0, K, V0, R, TIME, RHO, KAPPA, THETA, SIGMA){
   dt <- TIME/N_AVI
   S_array <- c()
   S2_array <- c()
   V_array <- c()
   V2_array <- c()
   payoffs <- c()
   payoffs2 <- c()
   for(j in 1:N_MC){
      V <- V0
      S <- S0
      
      V2 <- V0
      S2 <- S0
      for(i in 1:N_AVI){
         Zv <- rnorm(1, mean = 0, sd = 1)
         Z2 <- rnorm(1, mean = 0, sd = 1)
         Zs <- RHO * Zv + sqrt(1 - RHO^2) * Z2
         
         S <- S + R * S * dt + sqrt(max(V, 0)) * sqrt(dt) * S * Zs + 0.25 * S^2 * dt * (Zs^2 - 1)
         V <- V + KAPPA * (THETA - max(V, 0)) * dt + SIGMA * sqrt(max(V,0)) * sqrt(dt) * Zv + 0.25 * SIGMA^2 * dt * (Zv^2 - 1)
         
         S2 <- S2 * exp( (R - 0.5 * max(V2, 0)) * dt + sqrt(max(V2, 0) * dt) * Zs)
         V2 <- (sqrt(max(V2, 0)) + 0.5 * SIGMA * sqrt(dt) * Zv)^2 + KAPPA * (THETA - max(V2, 0)) * dt - 0.25 * SIGMA^2 * dt
      }
      V_array[j] <- V
      S_array[j] <- S
      V2_array[j] <- V2
      S2_array[j] <- S2
      
      #payoffs[j] <- max( (S_array[j] - K), 0)
      payoffs[j] <- max((S - K), 0)
      payoffs2[j] <- max((S2 - K), 0)
      
   }
   cat("N_MC:", N_MC,       "payoffs len:", length(payoffs),  '\n')
   df <- data.frame(S0, S_array, S2_array, K, payoffs, payoffs2, V_array, V2_array)
   print(df)
   print(mean(payoffs) * exp(-R * TIME))
   
   OptionPrice <- mean(payoffs) * exp(-R * TIME)
   OptionPrice2 <- mean(payoffs2) * exp(-R * TIME)
   
   cat("dt:", dt, "N_MC:", N_MC, "N_AVI:", N_AVI, "\n")
   return(c(OptionPrice, OptionPrice2))
}


calcHMMD.1(N_MC = 2000, N_AVI = 6*32,
           S0 = 100,
           K = 70,
           V0 = 0.0225,
           R = 0.04,
           TIME = 6,
           RHO = -0.5,
           KAPPA = 2,
           THETA = 0.04,
           SIGMA = 0.3)

calcHMMD.2(N_MC = 2000, N_AVI = 6*32,
           S0 = 100,
           K = 70,
           V0 = 0.0225,
           R = 0.04,
           TIME = 6,
           RHO = -0.5,
           KAPPA = 2,
           THETA = 0.04,
           SIGMA = 0.3)

callHestoncf(S = 100,
             X = 70,
             v0 = 0.1,
             r = 0.1,
             tau = 1,
             rho = -0.5,
             k = 1, 
             vT = 0.1,
             sigma = 0.1,
             q = 0)

