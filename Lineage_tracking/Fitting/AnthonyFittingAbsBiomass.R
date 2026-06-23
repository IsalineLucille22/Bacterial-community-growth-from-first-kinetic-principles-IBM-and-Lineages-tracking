library(ggplot2)
library(ggpubr)
library(readxl)
library(dplyr)
library(grid)
library(FME)
library(deSolve)





## =======================================================================
## Model growth simulations OD
## =======================================================================
rm(list = ls())

setwd("/Users/iguex/Library/CloudStorage/OneDrive-UniversitédeLausanne/Images_analysis/Scripts/Lineage_tracking/Fitting")

Data = read_excel("GrowthKinetic_3CBA_triplicate.xlsx", sheet = 3, range = "A1:B143")


Num_iter = length(Data[1,]) - 1
TT = data.matrix(Data[,1])
Mu_max_vect = c()
Mu_max_vect_exp = c()
Mu_max_vect_Monod = c()
Lag_time_vect = c()
Lag_time_vect_Monod = c()
nb_obs = length(data.matrix(Data[,1]))


for(i in 1:Num_iter){#1:Num_iter
  Species_number = i + 1
  Data_temp = Data
  
  nb_data_pts = length(data.matrix(Data_temp[,1]))
  TT = data.matrix(Data_temp[,1])
  Time_Eval = TT
  
  Data_Cells <- data.matrix(Data_temp[,Species_number])*1e15 #Abundance in fg
  
  x_0 = Data_Cells[1]
  
  state = c(x = mean(Data_Cells[1,]))
    
  R = min(3*max(Data_Cells), 3*mean(max(Data_Cells), max(Data_Cells[nb_data_pts],0)))
  state_Monod = c(x = mean(Data_Cells[1,]), R = R)
    
  Data_Cells <- data.frame(
    time = TT,
    x = c(Data_Cells)
  )
  colnames(Data_Cells) <- c('time','x')
    
  #Function for the logistic estimation
  logist <- function(t, state, parms) {
    with(as.list(c(state, parms)), {
      dx <- 1/(1 + (LT/t)^40)*mu_max*x*(1 - x/Ks) #With lag time
      # dx <- mu_max*x*(1 - x/Ks) #Without lag time
      list(dx)
    })
  }
    
  Monod <- function(t, state, parms){
    with(as.list(c(state, parms)), {
      alpha = 1
      R_conc = max(alpha*R, 0)
      dx = 1/(1 + (LT/t)^40)*x*mu_max*R_conc/(R_conc + Ks) #1/(1 + (LT/t)^40)*x*mu_max*R_conc/(R_conc + Ks)
      dR = -1/(1 + (LT/t)^40)*x*(mu_max/yield)*R_conc/(R_conc + Ks) #1/(1 + (LT/t)^40)*x*(mu_max/yield)*R_conc/(R_conc + Ks)
      list(c(dx, dR))
    })
  }
    
  ##===================================
  ## Fitted with logistic model #
  ##===================================
  ## numeric solution 
  ## ODEs system
  parms_init <- c(mu_max = 0.5, Ks = 9e5, LT = 0.1)
  parms_init_Monod <- c(mu_max = 0.5, Ks = 9e2, LT = 0.1, yield = 0.3)
  Times <- TT
  
  # model cost,
  ModelCost2 <- function(P) {
    out <- ode(y = state, func = logist, parms = P, times = TT)
    model = out
    return(modCost(out, Data_Cells)) # object of class modCost
  }
  
  ModelCostMonod <- function(P) {
    out <- ode(y = state_Monod, func = Monod, parms = P, times = TT, atol = 1e-11, rtol = 1e-10)
    model = out
    model = model[,1:2]
    return(modCost(model, Data_Cells)) # object of class modCost
  }
  
  
  Fit <- modFit(f = ModelCost2, p = parms_init, lower = c(0, 0, 0),
                upper = c(2.5, 1e14, 5))

  out <- ode(y = state, func = logist, parms = Fit$par,
             times = Time_Eval)
  
  
  FitMonod <- modFit(f = ModelCostMonod, p = parms_init_Monod, lower = c(0, 0, 0, 0),
                     upper = c(2.5, 1e14, 4, 1))
  
  outMonod <- ode(y = state_Monod, func = Monod, parms = FitMonod$par,
                  times = Time_Eval, atol = 1e-11, rtol = 1e-10)
  
  Param_est = Fit$par
  Mu_max_vect[i] = Param_est[1]
  Lag_time_vect[i] = Param_est[3]
  
  #Linear model, exponential growth
  y = log(Data_Cells$x)
  y = y[1:25]
  x = Data_Cells$time
  x = x[1:25]
  model <- lm(y ~ x) #Linear model
  p = model$coefficients
  Mu_max_vect_exp[i] = p[2]
  plot(x, y, pch = 21, bg = alpha("green", 0.4), col = alpha("green", 0.4))
  abline(model, col = "blue", lwd = 2)
  
  Param_est_Monod = FitMonod$par
  Mu_max_vect_Monod[i] = Param_est_Monod[1]*R[1]/(R[1] + Param_est_Monod[2])
  print(Param_est_Monod[2])
  Lag_time_vect_Monod[i] = Param_est_Monod[3]
  
  # pdf(file = paste("/Users/iguex/Library/CloudStorage/OneDrive-UniversitédeLausanne/Images analysis/Scripts/Lineage tracking/Figures/Fitting/Replicate_", i, ".pdf") ,   # The directory you want to save the file in
  #     width = 4, # The width of the plot in inches
  #     height = 4) # The height of the plot in inches
  
  plot(Data_Cells, xlim = c(0, max(TT)), pch = 21, bg = alpha("green", 0.4), col = alpha("green", 0.4))
  lines(out, col = "red", lty = 2)
  lines(outMonod, col = "blue", lty = 2)
  title(i)
  # dev.off()
  summary(FitMonod)
}
