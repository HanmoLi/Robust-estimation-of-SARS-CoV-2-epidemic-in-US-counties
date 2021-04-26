##########################
# Code for figure 1 in the Extended data figures section
# 
# Author: Hanmo Li & Mengyang Gu
#
# Email: mengyang at pstat.ucsb.edu
##########################

library(deSolve)
library(stringr)
library(zoo)
library(ggplot2)
library(mFilter)
library(RobustGaSP)
library(matrixStats)
library(dplyr)

source(file = "./Reproducing results in the paper/Codes/functions_8th_version.R")

state_name = "California"
state_name_short = "CA"

county_name_selected = "Imperial"

base = 'https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/'
death = 'time_series_covid19_deaths_'

set.seed(1)

gamma = 0.2
theta = 0.1
delta = 0.0066
# delta = 0.005

# download data from JHU
us_death = read.csv(paste0(base,death,"US.csv"))

death_rate = function(x){x/us_death$Population*10^6}
us_death_rate = us_death %>%
  mutate_at(vars(starts_with("X")), death_rate)

col_names = colnames(us_death)
all_dates = as.Date(paste0(str_sub(col_names[13:dim(us_death)[2]], 2, -1), "20"), tryFormats = c("%m.%d.%Y"))

RData_path = "./Reproducing results in the paper/Data/RData as of Sep 20 2020"

load(paste0(RData_path, "/", state_name_short, "_results_in_GP.RData"))

county_results = results_overall_list[[which(county_names == county_name_selected)]]

start_date = county_results$date[1]
end_date = county_results$date[which(county_results$prediction_indicator==1)[1]-1]

n = sum(county_results$prediction_indicator==0)

N = as.numeric( us_death%>%
  filter(Admin2 == county_name_selected, Province_State == state_name) %>%
  select(Population))

death_selected_raw = us_death%>%
  filter(Admin2 == county_name_selected, Province_State == state_name) %>%
  select(starts_with("x"))
death_selected_raw = as.numeric(death_selected_raw[all_dates %in% seq.Date(start_date, end_date, by=1)])

death_fit = county_results$death_f[county_results$prediction_indicator==0]

death_selected = county_results$death_real[county_results$prediction_indicator==0]

approx_beta_seq = county_results$Beta_estimated[1:(n-1)]

S_t_seq = county_results$suspective_f[1:n]

estimated_confirm = N - S_t_seq

I_0 = county_results$infective_f[1]
R_0 = county_results$resolved_f[1]

#########################
# Part 1: GP for data sampling
#########################

# sample the log daily confirmed cases to get the CI for beta_t and death

input_confirmed_cases= as.matrix(seq(1,n-1 ,1))

output_confirmed_cases=as.matrix(log(diff(estimated_confirm)+1))

m_gp_confirmed=rgasp(input_confirmed_cases, (output_confirmed_cases),
                     # trend = X,
                     # nugget.est=T,
                     num_initial_values=5 ,
                     # kernel_type = "matern_5_2",
                     kernel_type = "pow_exp", alpha=1.9,
                     range.par=1,
                     nugget=.25
)

log_diff_confirm_pred = predict(m_gp_confirmed, input_confirmed_cases,interval_data=F)

log_diff_confirm_pred$upper95 - log_diff_confirm_pred$lower95

plot(log_diff_confirm_pred$mean, type="l", col = "black")
lines(output_confirmed_cases, col="red")
legend("topright", legend = c("Log diff confirmed", "mean of the GP prediction"), lty = c(1,1), col = c("black", "red"))

num_samples=500 # sample 100 samples from adjusted confirmed cases

set.seed(1)

sample_log_diff_confirmed_cases=simulate(m_gp_confirmed,
                                         input_confirmed_cases,
                                         num_sample=num_samples, sample_data=T) #, testing_trend = X)

sample_confirmed_cases = rbind(rep(estimated_confirm[1], num_samples),
                               apply(exp(sample_log_diff_confirmed_cases)-1, 2, cumsum)+estimated_confirm[1])


# eliminate decreasing trend in sampled confirmed cases
for(i in 1:num_samples){
  for(j in 1:(n-1)){
    if (sample_confirmed_cases[j+1, i] < sample_confirmed_cases[j, i]){
      sample_confirmed_cases[j+1, i] = sample_confirmed_cases[j, i]
    }
  }
}

beta_record_matrix = matrix(NA, num_samples, n-1)

Susceptive_5_record_matrix = matrix(NA, num_samples, n)
Death_5_record_matrix = matrix(NA, num_samples, n)
Infectious_5_record_matrix = matrix(NA, num_samples, n)

# Calculate transmission rates for 100 sampled confirmed and death cases
for(sim_index in 1:num_samples){
  
  sample_death_each = death_selected
  sample_confirm_each = sample_confirmed_cases[, sim_index]
  
  daily_sample_confirm_each = diff(sample_confirm_each)
  daily_sample_confirm_each_smoothed = data_seven_day_smoothing(daily_sample_confirm_each)
  confirm_selected_smoothed_sampled = c(sample_confirm_each[1], cumsum(daily_sample_confirm_each_smoothed)+sample_confirm_each[1])
  
  for(j in 1:(n-1)){
    if (confirm_selected_smoothed_sampled[j+1] < confirm_selected_smoothed_sampled[j]){
      confirm_selected_smoothed_sampled[j+1] = confirm_selected_smoothed_sampled[j]
    }
  }
  
  estimated_confirm_sample = confirm_selected_smoothed_sampled #/ratio_sample
  
  S_t_seq_sample = N - estimated_confirm_sample
  
  init_for_beta_sample = c(S_t_seq_sample[1], I_0, R_0, sample_death_each[1], 0)
  
  param_record_approx_all = matrix(0, 5, n)
  param_record_approx_all[,1] = init_for_beta_sample
  param_record_approx_all[1,] = S_t_seq_sample
  
  approx_beta_seq_all = rep(0, n-1)
  
  for (i in 1:(n-1)){
    S_t_1 = param_record_approx_all[1,i]
    S_t_2 = param_record_approx_all[1,i+1]
    I_t_1 = param_record_approx_all[2,i]
    R_t_1 = param_record_approx_all[3,i]
    D_t_1 = param_record_approx_all[4,i]
    C_t_1 = param_record_approx_all[5,i]
    
    if(I_t_1<1){
      I_t_1 = 1
    }
    
    beta_t_1_2 = uniroot(find_root_beta, c(0, 10^6), tol = 0.0001, param = c(S_t_1, S_t_2, I_t_1), N = N, gamma=gamma)
    
    I_t_2 = I_t_1 * exp(beta_t_1_2$root*(S_t_1 + S_t_2)/(2*N) - gamma)
    R_t_2 = (2-theta)/(2+theta)*R_t_1 + gamma/(2+theta)*(I_t_1+I_t_2)
    D_t_2 = D_t_1 + delta*theta*(R_t_1+R_t_2)/2
    C_t_2 = C_t_1 + (1-delta)*theta*(R_t_1+R_t_2)/2
    
    param_record_approx_all[2:5, i+1] = c(I_t_2, R_t_2, D_t_2, C_t_2)
    approx_beta_seq_all[i] = beta_t_1_2$root
  }
  
  beta_record_matrix[sim_index, ] = approx_beta_seq_all
  
  Susceptive_5_record_matrix[sim_index, (1):n] = param_record_approx_all[1,]
  Infectious_5_record_matrix[sim_index, (1):n] = param_record_approx_all[2,]
  Death_5_record_matrix[sim_index, (1):n] = param_record_approx_all[4,]
  
}

death_sample_upper = apply(Death_5_record_matrix, 2, quantile, 0.975)
death_sample_lower = apply(Death_5_record_matrix, 2, quantile, 0.025)

death_sample_mean = apply(Death_5_record_matrix, 2, mean)

I_t_sample_upper = apply(Infectious_5_record_matrix, 2, quantile, 0.975)
I_t_sample_lower = apply(Infectious_5_record_matrix, 2, quantile, 0.025)

I_t_sample_mean = apply(Infectious_5_record_matrix, 2, mean)


########### compare with the ODE method
set.seed(1)

beta_t_from_euler = calculate_Beta_t_directly(death_selected, N)

beta_t_from_euler[beta_t_from_euler<(0)] = 0
beta_t_from_euler[beta_t_from_euler>10] = 10
betafun_euler <- approxfun(1:length(beta_t_from_euler), beta_t_from_euler[1:length(beta_t_from_euler)])

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

RSS <- function(parameters) {
  I_0 = parameters[1]
  R_0 = parameters[2]
  D_0 = death_selected[1]
  C_0 = 0
  init <- c(S = N-(I_0+R_0+D_0+C_0) , I = I_0, R = R_0, D =D_0, C = C_0)
  
  out <- ode(y = init, times = 1:length(beta_t_from_euler), func = SIRDC_multi_beta, parms = NULL)
  index_selected_fitted= 1:length(beta_t_from_euler) #(length(death_selected)-fitted_days_beta+1-5):(length(death_selected)-5)

  fit_D = out[,5]

  loss_cum_death = sum(((death_selected[index_selected_fitted] - fit_D))^2)# * train_Day^2) # weighted loss function
  loss_daily_death = sum(((diff(death_selected[index_selected_fitted])) - (diff(fit_D)))^2 ) #* (2:length(beta_seq))^2 ) # weighted loss function

  return(loss_cum_death + loss_daily_death)
}

if(county_name_selected == "Santa Barbara"){
  I_0_ini =400
  R_0_ini = 400
} else if (county_name_selected == "Imperial"){
  I_0_ini = 800
  R_0_ini = 800
}

set.seed(1)

ui = matrix(c(1,0,-1,0,1,-1), 3,2)
ci = matrix(c(0,0, -(N)  ))

Opt_euler = constrOptim(c(I_0_ini,R_0_ini),
                        RSS,
                        NULL,
                        ui = ui,
                        ci = ci
                        )

Opt_euler_par <- setNames(Opt_euler$par, c('I0','R0'))

t_euler <- 1:length(beta_t_from_euler) # prediction in next 100 days

I_0_euler = as.numeric(Opt_euler_par[1])

R_0_euler = as.numeric(Opt_euler_par[2])

D_0 = death_selected[1]
C_0 = 0

init_fitted_euler <- c(S = N-I_0_euler-R_0_euler-D_0-C_0 , I = I_0_euler, R = R_0_euler, D =death_selected[1], C = 0)

fit_euler <- as.data.frame(ode(y = init_fitted_euler, times = t_euler, func = SIRDC_multi_beta, parms = NULL))

date_death_all = seq.Date(start_date, end_date, by=1)
date_death_euler =  date_death_all[t_euler]
date_death_ode = seq.Date(start_date, end_date-2, by=1)

date_seq_beta = seq.Date(start_date, end_date-1, by=1)
date_seq_beta_smoothed = seq.Date(start_date+3, end_date-1-3, by=1)
date_seq_beta_euler = date_seq_beta[t_euler]

y_limit_death = c(min(death_selected_raw, fit_euler$D, death_fit,death_sample_upper), max(death_selected_raw, fit_euler$D, death_fit,death_sample_upper))

file_path = "./Reproducing results in the paper/Results/Extended_fig_1/"

cairo_ps(file = paste0(file_path, "Real_results_compare_with_Euler_Imperial.eps"), onefile = FALSE, fallback_resolution = 600, width = 16, height = 3)

par(mfrow=c(1,3))
par(cex=0.8, mai=c(0.55,0.55,0.2,0.2),mgp = c(2.2,0.7,0))

R_eff_Euler = beta_t_from_euler[t_euler]/0.2 * fit_euler$S/N

R_eff_Euler_smooth = data_seven_day_smoothing(R_eff_Euler)

R_eff_record_matrix =  beta_record_matrix/0.2 * Susceptive_5_record_matrix[, 2:n] / N 

R_eff_record_matrix_smoothed = t(apply(R_eff_record_matrix, 1, data_seven_day_smoothing))

R_eff_our_fit_sample_upper = apply(R_eff_record_matrix_smoothed, 2, quantile, 0.975)
R_eff_our_fit_sample_lower = apply(R_eff_record_matrix_smoothed, 2, quantile, 0.025)

R_eff_our_fit_sample_mean = apply(R_eff_record_matrix_smoothed, 2, mean)

y_limit_R_t = c(min(R_eff_Euler_smooth, R_eff_our_fit_sample_lower, na.rm = T), max(R_eff_Euler_smooth, R_eff_our_fit_sample_upper, na.rm = T))

plot(R_eff_our_fit_sample_mean~date_seq_beta, type="l",ylim = y_limit_R_t, ylab = "Effective reproduction number", xlab = "Date", col="blue", cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5) #, main = paste0(county_names[each_index], ", population=", round(N/10^6,2),"M", ", Ratio = ", round(ratio_real,3)))
polygon(c(rev(date_seq_beta), date_seq_beta),
        c(rev(R_eff_our_fit_sample_upper), R_eff_our_fit_sample_lower),
        col = 'grey80',
        border = NA)
lines(R_eff_our_fit_sample_mean ~ date_seq_beta, col = "blue")
lines(R_eff_Euler_smooth ~ date_seq_beta_euler, col = "red", lty=2, lwd = 2)
abline(h=1, lty=2,col="black")

y_limit_infectious = range(fit_euler$I, I_t_sample_mean, na.rm=T)

plot(I_t_sample_mean~date_death_all,ylim = y_limit_infectious, type = "l", col = "black",xlab = "Date", ylab = "Active infectious individuals", cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5) #,main = paste0(county_name_selected, ", population = ", round(N/10^6,2), "M"))
polygon(c(rev(date_death_all), date_death_all),
        c(rev(I_t_sample_upper), I_t_sample_lower),
        col = 'grey80',
        border = NA)
lines(fit_euler$I~date_death_euler, col = "red", lty=2, lwd = 2)
lines(I_t_sample_mean~date_death_all, col = "blue")

plot(death_selected_raw~date_death_all,ylim = y_limit_death, type = "l", col = "black",xlab = "Date", ylab = "Cumulative death toll", cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5) #,main = paste0(county_name_selected, ", population = ", round(N/10^6,2), "M"))
polygon(c(rev(date_death_all), date_death_all),
        c(rev(death_sample_upper), death_sample_lower),
        col = 'grey80',
        border = NA)
lines(death_selected_raw~date_death_all, col = "black")
lines(fit_euler$D~date_death_euler, col = "red", lty=2, lwd = 2)
lines(death_sample_mean~date_death_all, col = "blue")
legend("topleft", legend = c("Robust estimation", "F&J", "Observed death toll"), lty = c(1,2,1), col = c("blue", "red", "black"), lwd = c(1,2,1),cex = 1.2)

dev.off()

