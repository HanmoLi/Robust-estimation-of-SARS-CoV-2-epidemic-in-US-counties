##########################
# Code for figure 2 in the Extended data figures section
# 
# Author: Hanmo Li & Mengyang Gu
#
# Email: mengyang at pstat.ucsb.edu
##########################

library(deSolve)
library(EpiEstim)

source(file = "./Reproducing results in the paper/Codes/functions_8th_version.R")

noise_index = F

gamma = 0.2
theta = 0.1
delta = 0.0066

n = 100

N = 10^7

I_0 = 1000

R_0 = 1000

D_0 = 0

C_0 = 0

set.seed(3)

if(noise_index){
  beta_t_simulation = exp(-0.7 * seq(1,10,length.out = n))+ rnorm(n, 0, 0.2)
  beta_t_simulation[beta_t_simulation<0] = 0
} else{
  beta_t_simulation = exp(-0.7 * seq(1,10,length.out = n))
}

plot(beta_t_simulation)

betafun = stepfun(2:(n),beta_t_simulation)
plot(betafun)

parameters = list(1:5)
parameters[[1]] = gamma
parameters[[2]] = theta
parameters[[3]] = delta
parameters[[4]] = N
parameters[[5]] = betafun

init = c(
  S = N - (I_0 + R_0 + D_0 + C_0) ,
  I = I_0,
  R = R_0,
  D = D_0,
  C = C_0
)

out_ode = ode(
  y = init,
  times = 1:n,
  func = SIRDC,
  parms = parameters,

)

death_selected = out_ode[,5]

S_t_seq = out_ode[,2]

confirm_selected_smoothed = N - S_t_seq

plot(death_selected)

out_ode_step_size_1 = ode(
  y = init,
  times = 1:n,
  func = SIRDC,
  parms = parameters,
  method = "rk4",
  hini=1,
)

death_selected_rk4_step_1 = out_ode_step_size_1[,5]
infective_rk4_step_1 = out_ode_step_size_1[,3]
susceptive_rk4_step_1 =  out_ode_step_size_1[,2]

plot(death_selected_rk4_step_1)
lines(death_selected)

out_ode_step_size_0.1 = ode(
  y = init,
  times = 1:n,
  func = SIRDC,
  parms = parameters,
  method = "rk4",
  hini=0.1,
)

death_selected_rk4_step_0.1 = out_ode_step_size_0.1[,5]
infective_rk4_step_0.1 = out_ode_step_size_0.1[,3]
susceptive_rk4_step_0.1 =  out_ode_step_size_0.1[,2]

plot(death_selected_rk4_step_0.1)
lines(death_selected)

#########################
# Part 1: Parameter Estimation via our approximation algorithm
#########################

I_0_approx = I_0
R_0_approx = R_0

ratio = confirm_selected_smoothed[1]/(I_0_approx+ R_0_approx+ death_selected[1])

# use ratio to adjust smoothed confirmed cases and get susceptible cases S_t
estimated_confirm = confirm_selected_smoothed/ratio

S_t_seq = N - estimated_confirm

init_for_beta = c(S_t_seq[1], I_0_approx, R_0_approx, death_selected[1], 0)

param_record_approx_for_beta = matrix(0, 5, n) # 5 rows: S_t, I_t, R_t, D_t, C_t
param_record_approx_for_beta[,1] = init_for_beta
param_record_approx_for_beta[1,] = S_t_seq

approx_beta_seq = rep(0, n-1) # record the value of transmission rate

# system.time({
for (i in 1:(n-1)){
  S_t_1 = param_record_approx_for_beta[1,i]
  S_t_2 = param_record_approx_for_beta[1,i+1]
  I_t_1 = param_record_approx_for_beta[2,i]
  R_t_1 = param_record_approx_for_beta[3,i]
  D_t_1 = param_record_approx_for_beta[4,i]
  C_t_1 = param_record_approx_for_beta[5,i]
  
  if(I_t_1<1){
    I_t_1 = 1
  }
  
  beta_t_1_2 = uniroot(find_root_beta, c(0, 10^6), tol = 0.0001, param = c(S_t_1, S_t_2, I_t_1), N = N, gamma=gamma)
  
  I_t_2 = I_t_1 * exp(beta_t_1_2$root*(S_t_1 + S_t_2)/(2*N) - gamma)
  R_t_2 = (2-theta)/(2+theta)*R_t_1 + gamma/(2+theta)*(I_t_1+I_t_2)
  D_t_2 = D_t_1 + delta*theta*(R_t_1+R_t_2)/2
  C_t_2 = C_t_1 + (1-delta)*theta*(R_t_1+R_t_2)/2
  
  param_record_approx_for_beta[2:5, i+1] = c(I_t_2, R_t_2, D_t_2, C_t_2)
  approx_beta_seq[i] = beta_t_1_2$root
}

death_fit = param_record_approx_for_beta[4, ]
infectious_fit = param_record_approx_for_beta[2,]

beta_t_fit = approx_beta_seq

plot(beta_t_fit)
lines(beta_t_simulation)

plot(death_fit)
lines(death_selected)

#######################################
#Part 2: Euler method
#######################################
beta_t_from_euler = calculate_Beta_t_directly(death_selected, N)

beta_t_from_euler[beta_t_from_euler<0] = 0
beta_t_from_euler[beta_t_from_euler>10] = 10

betafun_euler = approxfun(beta_t_from_euler)

plot(beta_t_from_euler)

SIRDC_multi_beta <- function(time, state, parameters) {
  par <- as.list(c(state))
  with(par, {
    
    beta = betafun_euler(time)
    
    dS <- -beta/N * I * S
    dI <- beta/N * I * S - gamma * I
    dR <- gamma * I - theta * R
    dD <- delta * theta * R
    dC <- (1-delta) * theta * R
    list(c(dS, dI, dR, dD, dC))
  })
}

t_euler <- (1):(length(beta_t_from_euler))

I_0_euler = I_0 

R_0_euler = R_0

D_0 = 0
C_0 = 0

init_fitted_euler <- c(S = N-I_0_euler-R_0_euler-D_0-C_0 , I = I_0_euler, R = R_0_euler, D =0, C = 0)

fit_euler <- as.data.frame(ode(y = init_fitted_euler, times = t_euler, func = SIRDC_multi_beta, parms = NULL))


y_limit_death = c(min(death_selected, fit_euler$D, death_fit, na.rm=T), max(death_selected, fit_euler$D, death_fit , na.rm=T))

par(mfrow=c(1,1))

if(noise_index){
  noise_text = "with noise"
}else{
  noise_text = "without noise"
}

#######################################
#Part 2: Draw and save figures
#######################################

file_path = "./Reproducing results in the paper/Results/Extended_fig_2/"

cairo_ps( file = paste0(file_path, "Simulation_results_compare_with_Euler_noise_free.eps"), onefile = FALSE, fallback_resolution = 600, width =3.6 * 5.5, height = 3.6)

par(mfrow=c(1,3))
par( mai=c(0.55,0.55,0.2,0.2),mgp = c(1.9,0.65,0), cex.lab=2, cex.axis=1.8, cex.main=2, cex.sub = 2)  # , mgp = c(1.5,0.2,0)

approx_R_eff = beta_t_fit/0.2 * S_t_seq[1: (length(S_t_seq)-1)]/N

euler_R_eff = beta_t_from_euler[1:(length(beta_t_from_euler)-1)]/0.2 * fit_euler$S[1:(length(fit_euler$S)-1)] / N

lsoda_R_eff = beta_t_simulation/0.2 * out_ode[,2] / N

RK_4_step_1_R_eff = beta_t_simulation/0.2 * susceptive_rk4_step_1/N
RK_4_step_0.1_R_eff = beta_t_simulation/0.2 * susceptive_rk4_step_0.1/N

approx_R_eff_smoothed = data_seven_day_smoothing(approx_R_eff)
real_R_eff_smoothed = data_seven_day_smoothing(lsoda_R_eff)
euler_R_eff_smoothed = data_seven_day_smoothing(euler_R_eff)
RK_4_step_1_R_eff_smoothed = data_seven_day_smoothing(RK_4_step_1_R_eff)
RK_4_step_0.1_R_eff_smoothed = data_seven_day_smoothing(RK_4_step_0.1_R_eff)


y_limit_R_t = range(approx_R_eff_smoothed,real_R_eff_smoothed ,euler_R_eff_smoothed, na.rm=T)

plot(approx_R_eff_smoothed, type="l",ylim = y_limit_R_t, ylab = "Effective reproduction number", xlab = "Days", col="blue") #, main = paste0(county_names[each_index], ", population=", round(N/10^6,2),"M", ", Ratio = ", round(ratio_real,3)))
lines(real_R_eff_smoothed, col = "black", type="p")
lines(approx_R_eff_smoothed, col = "blue", lwd = 1)
lines( euler_R_eff_smoothed, col = "red", lty=2, lwd = 2)
lines(RK_4_step_1_R_eff_smoothed, col = "green", lwd=2)
lines(RK_4_step_0.1_R_eff_smoothed, col = "green", lty=2, lwd=2)
abline(h=1, lty=2,col="black")

y_limit_infectious = range(fit_euler$I/10^5, infectious_fit/10^5, na.rm=T)

plot(infectious_fit/10^5,ylim = y_limit_infectious, type = "l",lwd = 2, col = "blue",xlab = "Days", ylab = expression(paste("Active infectious individuals ", (10^5)))) #,main = paste0(county_name_selected, ", population = ", round(N/10^6,2), "M"))
lines(fit_euler$I/10^5~t_euler, col = "red", lty=2, lwd = 2)
lines(out_ode[,3]/10^5, type = "p")
lines(infective_rk4_step_1/10^5, col ="green", lwd=2)
lines(infective_rk4_step_0.1/10^5, col ="green", lwd=2, lty=2)

plot(death_selected/10^3,ylim = y_limit_death/10^3, type = "p",pch = 2, cex = 1.2, col = "black",xlab = "Days", ylab = expression(paste("Cumulative death toll ", (10^3)))) #,main = paste0(county_name_selected, ", population = ", round(N/10^6,2), "M"))
lines(fit_euler$D/10^3~t_euler, col = "red", lty=2, lwd = 2)
lines(death_fit/10^3, col = "blue", lwd = 4)
lines(death_selected_rk4_step_1/10^3, col = "green", lwd=2)
lines(death_selected_rk4_step_0.1/10^3, col = "green", lty = 2, lwd = 2)
# legend(
#   "bottomright",
#   legend = c("Config 1: Robust estimation", "Config 2: F&J", "Config 3: RK4, step size=1", "Config 4: RK4, step size=0.1", "Config 5: lsoda"),
#   # type = c("l", "l", "p"),
#   pch = c(NA, NA, NA, NA,  2),
#   lty = c(1, 2, 1,2, NA),
#   col = c("blue", "red", "green","green", "black"),
#   lwd = c(4,2,2,2,1.2),
#   cex = 1.5
# )

dev.off()
