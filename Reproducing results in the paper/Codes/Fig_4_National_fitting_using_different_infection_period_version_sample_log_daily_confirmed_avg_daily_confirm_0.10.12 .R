##########################
# Code for figure 4 in the main manuscript
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
library(tidyr)
library(dplyr)
library(foreach)
library(doParallel)
library(grDevices)

source(file = "./Reproducing results in the paper/Codes/functions_8th_version.R")

##for forecast
fitted_days_beta=184 ##how many days from the end day you use to estimate beta_t
predicted_days_train=184  

prediction_length = 90
beta_t_length=184

gamma = 0.2
theta = 0.1
delta = 0.0066

set.seed(1)

base = 'https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/'
death = 'time_series_covid19_deaths_'
confirm = 'time_series_covid19_confirmed_'

covid19_project_url = "https://api.covidtracking.com/v1/states/daily.csv"

us_test = read.csv(covid19_project_url)

# download death and confirmed data from JHU
us_death = read.csv(paste0(base,death,"US.csv"))
us_confirm = read.csv(paste0(base,confirm,"US.csv"))

# clean the death data
us_death = clean_us_death(us_death)

col_names = colnames(us_death)
all_dates = as.Date(paste0(str_sub(col_names[13:dim(us_death)[2]], 2, -1), "20"), tryFormats = c("%m.%d.%Y"))

file_pow_param_path = "./Reproducing results in the paper/Data/"

us_pow_param = read.csv(paste0(file_pow_param_path,"Name_list_of_US_states_updated.csv"))

us_pow_param$state_name = as.character(us_pow_param$state_name)
state_names_all = us_pow_param$state_name

start_date = as.Date("2020-03-21")

end_date = as.Date("2020-09-20")

standard_time_seq = seq.Date(start_date, end_date, by=1)

n_standard  =length(standard_time_seq)

# get the national fitting

N = 321418820 # the population in 50 US states

US_confirmed_overall = rep(0, n_standard)
US_death_overall = rep(0, n_standard)

# us death from 2020-01-22 to today
US_death_overall_all_length = rep(0, length(all_dates))

# sum of state level confirmed and death cases
for(state_index in 1:dim(us_pow_param)[1]){
  
  state_name = us_pow_param$state_name[state_index]
  state_name_short = us_pow_param$state_name_short[state_index]
  
  state_confirmed = us_confirm %>%
    filter(Province_State == state_name)%>%
    dplyr::select(starts_with("x"))
  
  state_confirmed_selected = state_confirmed[all_dates %in% standard_time_seq]
  
  US_confirmed_overall = US_confirmed_overall + as.numeric(apply(state_confirmed_selected, 2, sum))
  
  state_death = us_death %>%
    filter(Province_State == state_name)%>%
    dplyr::select(starts_with("x"))
  
  US_death_overall_all_length = US_death_overall_all_length + as.numeric(apply(state_death,2,sum))
  
  state_death_selected = state_death[all_dates %in% standard_time_seq]
  
  US_death_overall = US_death_overall + as.numeric(apply(state_death_selected,2,sum))
}

confirm_selected = US_confirmed_overall
death_selected = US_death_overall

# eliminate decreasing trend in the county death data

for(i in 1:(n_standard-1)){
  if (death_selected[i+1] < death_selected[i]){
    death_selected[i+1] = death_selected[i]
  }
}

for( i in 1:(n_standard-1)){
  if (confirm_selected[i+1]<confirm_selected[i]){
    confirm_selected[i+1] = confirm_selected[i]
  }
}


# try to calculated the 7-day average using daily cases
daily_confirm_selected = diff(confirm_selected)
daily_confirm_selected_avg = data_seven_day_smoothing(daily_confirm_selected)
unadjusted_confirm_selected_smoothed = c(confirm_selected[1], cumsum(daily_confirm_selected_avg) + confirm_selected[1])

us_test_selected = us_test %>%
  filter(state %in% us_pow_param$state_name_short) %>%
  dplyr::select(date, state, positiveIncrease, totalTestResultsIncrease)

us_test_selected$date = as.Date(as.character(us_test_selected$date), "%Y %m %d")

us_test_selected = us_test_selected %>%
  filter(date %in% standard_time_seq)

# data processing
us_test_selected$totalTestResultsIncrease[us_test_selected$totalTestResultsIncrease<0] = abs(us_test_selected$totalTestResultsIncrease[us_test_selected$totalTestResultsIncrease<0])

us_test_selected$positiveIncrease[us_test_selected$positiveIncrease<0] = abs(us_test_selected$positiveIncrease[us_test_selected$positiveIncrease<0])

us_test_aggregated = us_test_selected %>%
  group_by(date) %>%
  summarise_each(funs(sum), positiveIncrease, totalTestResultsIncrease)

us_test_aggregated$positiveIncrease_7_day_avg = data_seven_day_smoothing(us_test_aggregated$positiveIncrease)

us_test_aggregated$totalTestResultsIncrease_7_day_avg = data_seven_day_smoothing(us_test_aggregated$totalTestResultsIncrease)

us_test_aggregated$positive_rate = us_test_aggregated$positiveIncrease_7_day_avg / us_test_aggregated$totalTestResultsIncrease_7_day_avg

plot(us_test_aggregated$positive_rate~us_test_aggregated$date)


# do not need to reverse the positive rate
daily_positive_rate_selected = us_test_aggregated$positive_rate[us_test_aggregated$date>=(start_date) & us_test_aggregated$date<=end_date]

alpha_ini_seq = seq(0,2,length.out = 1000)

upper_bound_record = rep(NA, length(alpha_ini_seq))

for(alpha_index in 1:length(alpha_ini_seq)){
  c_t_seq_each = daily_positive_rate_selected^(alpha_ini_seq[alpha_index])
  
  daily_confirm_selected_each = diff(confirm_selected)
  
  daily_confirm_selected_each_smoothed = data_seven_day_smoothing(daily_confirm_selected_each * c_t_seq_each[2:n_standard])
  
  confirm_selected_smoothed_each = c( c_t_seq_each[1] * confirm_selected[1] , cumsum(daily_confirm_selected_each_smoothed)+c_t_seq_each[1]*confirm_selected[1])
  
  upper_bound_each = (confirm_selected_smoothed_each[1]*N)/confirm_selected_smoothed_each[n_standard] - death_selected[1]
  
  upper_bound_record[alpha_index] = upper_bound_each
}
upper_final = min(upper_bound_record)


ui = matrix(c(c(1,0,-1,0,1, 0,0),c( 0, 1, 0,-1,-1, 0,0), c(0, 0, 0,0,0, 1, -1)), 7,3)
ci = matrix(c(-10^(10),-10^(10), -log(upper_final),-log(upper_final),0, -10^(10) , -log(2) ))

I_0_ini = 100000
R_0_ini = 100000

if((I_0_ini + R_0_ini)> abs(upper_final)){
  I_0_ini = abs(upper_final)/4
  R_0_ini = abs(upper_final)/4
}

alpha_vector = c(seq(0.1,1,0.1), seq(1.05, 1.95, 0.05))

error_record = rep(NA, length(alpha_vector))
m_approx_record = as.list( rep(NA, length(alpha_vector)))

pb <- txtProgressBar(min = 0, max = length(alpha_vector), style = 3)
for(alpha_index in 1:length(alpha_vector)){
  
  setTxtProgressBar(pb, alpha_index)
  
  m_approx_each = constrOptim(
    c(log(I_0_ini + R_0_ini), log(R_0_ini), log(alpha_vector[alpha_index])),
    loss_approx_beta_three_param_for_constrOptim_daily_avg,
    NULL,
    ui = ui,
    ci = ci,
    death_cases = death_selected,
    confirmed_cases = confirm_selected,
    real_unadjusted_smoothed = unadjusted_confirm_selected_smoothed,
    daily_positive_rate_selected = daily_positive_rate_selected,
    N_population = N,
    training_length = n_standard,
    fixed_global_params = c(gamma, theta, delta),
    penalty = F,
    fitted_days= fitted_days_beta,
    weight_loss=T
  )
  
  error_record[alpha_index] = m_approx_each$value
  m_approx_record[[alpha_index]] = m_approx_each
  
}
m_approx = m_approx_record[[which(error_record==min(error_record))[1]]]

I_0_plus_R_0 = exp(m_approx$par[1])
R_0 = exp(m_approx$par[2])
exp_index = exp(m_approx$par[3])

I_0 = I_0_plus_R_0 - R_0

if (I_0<1){
  I_0 = 1
}
if (R_0<1){
  R_0 = 1
}

#########################
# Part 3: Parameter Estimation: S_t, I_t, R_t, D_t, C_t and beta_t
#########################
c_t_seq = daily_positive_rate_selected^(exp_index)

daily_confirm_selected = diff(confirm_selected)

daily_confirm_selected_smoothed = data_seven_day_smoothing(daily_confirm_selected * c_t_seq[2:n_standard])

confirm_selected_smoothed =  c( c_t_seq[1] * confirm_selected[1] , cumsum(daily_confirm_selected_smoothed)+c_t_seq[1]*confirm_selected[1])

ratio = confirm_selected_smoothed[1]/(I_0+ R_0+ death_selected[1])

ratio_real = confirm_selected[2]/(I_0+ R_0+ death_selected[1])

# use ratio to adjust smoothed confirmed cases and get susceptible cases S_t
estimated_confirm = confirm_selected_smoothed/ratio

for(i in 1:length(estimated_confirm)){
  if (estimated_confirm[i] < unadjusted_confirm_selected_smoothed[i]){
    estimated_confirm[i] = unadjusted_confirm_selected_smoothed[i]
  }
}

S_t_seq = N - estimated_confirm

init_for_beta = c(S_t_seq[1], I_0, R_0, death_selected[1], 0)

param_record_approx_for_beta = matrix(0, 5, n_standard) # 5 rows: S_t, I_t, R_t, D_t, C_t
param_record_approx_for_beta[,1] = init_for_beta
param_record_approx_for_beta[1,] = S_t_seq

approx_beta_seq = rep(0, n_standard-1) # record the value of transmission rate

for (i in 1:(n_standard-1)){
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
  # I_t_2 = uniroot(find_root_I_t_2, c(0, N), tol = 0.0001, param = c(S_t_1, beta_t_1_2$root, I_t_1), N = N, gamma=gamma)
  I_t_2 = I_t_1 * exp(beta_t_1_2$root*(S_t_1 + S_t_2)/(2*N) - gamma)
  R_t_2 = (2-theta)/(2+theta)*R_t_1 + gamma/(2+theta)*(I_t_1+I_t_2)
  D_t_2 = D_t_1 + delta*theta*(R_t_1+R_t_2)/2
  C_t_2 = C_t_1 + (1-delta)*theta*(R_t_1+R_t_2)/2
  
  param_record_approx_for_beta[2:5, i+1] = c(I_t_2, R_t_2, D_t_2, C_t_2)
  approx_beta_seq[i] = beta_t_1_2$root
}

# calculate the smoothed transmission rate using 7 day average
approx_beta_seq_smoothed = rollapply(approx_beta_seq, width = 7, by = 1, FUN = mean, align = "left")

# calculate beta_t from Euler method
beta_t_from_euler = calculate_Beta_t_directly(death_selected, N)
# 
beta_t_from_euler[beta_t_from_euler<0] = 0
beta_t_from_euler[beta_t_from_euler>2] = 2

date_seq_beta = seq.Date(start_date, end_date-1, by=1)
date_seq_beta_smoothed = seq.Date(start_date+3, end_date-1-3, by=1)
date_seq_beta_euler = seq.Date(start_date+2, end_date-1-1, by=1)

y_limit_R_t = c(min(beta_t_from_euler/0.2,approx_beta_seq/0.2), max(beta_t_from_euler/0.2,approx_beta_seq/0.2))
plot(approx_beta_seq/0.2~date_seq_beta, type="p",ylim = y_limit_R_t, ylab = "R_t", xlab = "Date", main = paste0( "US , population=", round(N/10^6,2),"M", ", Ratio = ", round(ratio_real,3)))
lines(approx_beta_seq_smoothed/0.2~date_seq_beta_smoothed, col="blue")
lines(beta_t_from_euler/0.2 ~ date_seq_beta_euler, col = "green")
abline(h=1, lty=2,col="red")
abline(h=0,lty=1,col="black")
legend("topright", legend = c("Daily Approximated R_t", "7-day Averaged Approximated R_t", "R_t from Euler"), pch = c(1,NA,NA), lty=c(NA,1,1), col=c("black", "blue", "green"))

date_seq = seq.Date(start_date, end_date, by=1)

ylimit_death = c(min(param_record_approx_for_beta[4,], death_selected), max(param_record_approx_for_beta[4,], death_selected))

plot(param_record_approx_for_beta[4,]~date_seq,ylim = ylimit_death, type="l", col="blue", xlab = "Date", ylab = "Death Cases", main = paste0( "US, population=", round(N/10^6,2),"M", ", Ratio = ", round(ratio_real,3)))
lines(death_selected~date_seq, col = "red")
legend("topleft", legend = c("Observed Death Cases", "Estimated Death Cases"), lty = c(1,1), col = c("red", "blue"))

#########################
# Part 4: GP for data sampling
#########################

# sample the log daily confirmed cases to get the CI for beta_t and death

input_confirmed_cases= as.matrix(seq(1,n_standard-1 ,1))

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
  for(j in 1:(n_standard-1)){
    if (sample_confirmed_cases[j+1, i] < sample_confirmed_cases[j, i]){
      sample_confirmed_cases[j+1, i] = sample_confirmed_cases[j, i]
    }
  }
}

beta_record_matrix = matrix(NA, num_samples, n_standard-1)

Susceptive_5_record_matrix = matrix(NA, num_samples, n_standard)
Death_5_record_matrix = matrix(NA, num_samples, n_standard)
Infectious_5_record_matrix = matrix(NA, num_samples, n_standard)

# Calculate transmission rates for 100 sampled confirmed and death cases

for(sim_index in 1:num_samples){
  
  # print(sim_index / num_samples)
  
  sample_death_each = death_selected
  sample_confirm_each = sample_confirmed_cases[, sim_index]
  
  daily_sample_confirm_each = diff(sample_confirm_each)
  daily_sample_confirm_each_smoothed = data_seven_day_smoothing(daily_sample_confirm_each)
  confirm_selected_smoothed_sampled = c(sample_confirm_each[1], cumsum(daily_sample_confirm_each_smoothed)+sample_confirm_each[1])
  
  for(j in 1:(n_standard-1)){
    if (confirm_selected_smoothed_sampled[j+1] < confirm_selected_smoothed_sampled[j]){
      confirm_selected_smoothed_sampled[j+1] = confirm_selected_smoothed_sampled[j]
    }
  }
  
  estimated_confirm_sample = confirm_selected_smoothed_sampled #/ratio_sample
  
  S_t_seq_sample = N - estimated_confirm_sample
  
  init_for_beta_sample = c(S_t_seq_sample[1], I_0, R_0, sample_death_each[1], 0)
  
  param_record_approx_all = matrix(0, 5, n_standard)
  param_record_approx_all[,1] = init_for_beta_sample
  param_record_approx_all[1,] = S_t_seq_sample
  
  approx_beta_seq_all = rep(0, n_standard-1)
  
  for (i in 1:(n_standard-1)){
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
  
  Susceptive_5_record_matrix[sim_index, (1):n_standard] = param_record_approx_all[1,]
  Infectious_5_record_matrix[sim_index, (1):n_standard] = param_record_approx_all[2,]
  Death_5_record_matrix[sim_index, (1):n_standard] = param_record_approx_all[4,]

}

Susceptive_4_75_record_matrix = matrix(NA, num_samples, n_standard)
Infectious_4_75_record_matrix = matrix(NA, num_samples, n_standard)
Death_4_75_record_matrix = matrix(NA, num_samples, n_standard)


gamma_updated_1 = 1/4.75

for(sim_index in 1:num_samples){
  approx_beta_seq_all = beta_record_matrix[sim_index, ]
  
  init_for_beta = c(S_t_seq[1], I_0, R_0, death_selected[1], 0)
  
  param_record_approx_all = matrix(NA, 5, n_standard)
  param_record_approx_all[,1] = init_for_beta
  
  for (i in 1:(n_standard-1)){
    S_t_1 = param_record_approx_all[1,i]
    I_t_1 = param_record_approx_all[2,i]
    R_t_1 = param_record_approx_all[3,i]
    D_t_1 = param_record_approx_all[4,i]
    C_t_1 = param_record_approx_all[5,i]
    
    beta_t_1_2 = approx_beta_seq_all[i]
    
    if(I_t_1<1){
      I_t_1 = 1
    }
    
    S_t_2 = uniroot(find_root_S_t_2, c(0, N), tol = 0.0001, param = c(S_t_1, beta_t_1_2, I_t_1), N = N, gamma=gamma_updated_1)
    I_t_2 = I_t_1 * exp(beta_t_1_2*(S_t_1 + S_t_2$root)/(2*N) - gamma_updated_1)
    R_t_2 = (2-theta)/(2+theta)*R_t_1 + gamma_updated_1/(2+theta)*(I_t_1+I_t_2)
    D_t_2 = D_t_1 + delta*theta*(R_t_1+R_t_2)/2
    C_t_2 = C_t_1 + (1-delta)*theta*(R_t_1+R_t_2)/2
    
    param_record_approx_all[, i+1] = c(S_t_2$root, I_t_2, R_t_2, D_t_2, C_t_2)
    
    
  }
  
  Susceptive_4_75_record_matrix[sim_index, 1:n_standard] = param_record_approx_all[1,]
  Infectious_4_75_record_matrix[sim_index, (1):n_standard] = param_record_approx_all[2,]
  Death_4_75_record_matrix[sim_index, (1):n_standard] = param_record_approx_all[4,]

}


Susceptive_4_5_record_matrix = matrix(NA, num_samples, n_standard)
Infectious_4_5_record_matrix = matrix(NA, num_samples, n_standard)
Death_4_5_record_matrix = matrix(NA, num_samples, n_standard)

gamma_updated_2 = 1/4.5

for(sim_index in 1:num_samples){
  approx_beta_seq_all = beta_record_matrix[sim_index, ]
  
  init_for_beta = c(S_t_seq[1], I_0, R_0, death_selected[1], 0)
  
  param_record_approx_all = matrix(NA, 5, n_standard)
  param_record_approx_all[,1] = init_for_beta
  
  for (i in 1:(n_standard-1)){
    S_t_1 = param_record_approx_all[1,i]
    I_t_1 = param_record_approx_all[2,i]
    R_t_1 = param_record_approx_all[3,i]
    D_t_1 = param_record_approx_all[4,i]
    C_t_1 = param_record_approx_all[5,i]
    
    beta_t_1_2 = approx_beta_seq_all[i]
    
    if(I_t_1<1){
      I_t_1 = 1
    }
    
    S_t_2 = uniroot(find_root_S_t_2, c(0, N), tol = 0.0001, param = c(S_t_1, beta_t_1_2, I_t_1), N = N, gamma=gamma_updated_2)
    I_t_2 = I_t_1 * exp(beta_t_1_2*(S_t_1 + S_t_2$root)/(2*N) - gamma_updated_2)
    R_t_2 = (2-theta)/(2+theta)*R_t_1 + gamma_updated_2/(2+theta)*(I_t_1+I_t_2)
    D_t_2 = D_t_1 + delta*theta*(R_t_1+R_t_2)/2
    C_t_2 = C_t_1 + (1-delta)*theta*(R_t_1+R_t_2)/2
    
    param_record_approx_all[, i+1] = c(S_t_2$root, I_t_2, R_t_2, D_t_2, C_t_2)
    
    
  }
  Susceptive_4_5_record_matrix[sim_index, 1:n_standard] = param_record_approx_all[1,]
  Infectious_4_5_record_matrix[sim_index, (1):n_standard] = param_record_approx_all[2,]
  Death_4_5_record_matrix[sim_index, (1):n_standard] = param_record_approx_all[4,]
}

R_eff_7_day_avg_5_infection_period = matrix(NA, num_samples, n_standard-1)
R_eff_7_day_avg_4_75_infection_period = matrix(NA, num_samples, n_standard-1)
R_eff_7_day_avg_4_5_infection_period = matrix(NA, num_samples, n_standard-1)

for(sim_index in 1:num_samples){
  R_eff_5 = beta_record_matrix[sim_index, ] / gamma * Susceptive_5_record_matrix[sim_index, 2:n_standard] / N
  R_eff_5_avg = data_seven_day_smoothing(R_eff_5)
  
  R_eff_4_75 = beta_record_matrix[sim_index, ] / gamma_updated_1 * Susceptive_4_75_record_matrix[sim_index, 2:n_standard] / N
  R_eff_4_75_avg = data_seven_day_smoothing(R_eff_4_75)
  
  R_eff_4_5 = beta_record_matrix[sim_index, ] / gamma_updated_2 * Susceptive_4_5_record_matrix[sim_index, 2:n_standard] / N
  R_eff_4_5_avg = data_seven_day_smoothing(R_eff_4_5)
  
  R_eff_7_day_avg_5_infection_period[sim_index, ] = R_eff_5_avg
  R_eff_7_day_avg_4_75_infection_period[sim_index, ] = R_eff_4_75_avg
  R_eff_7_day_avg_4_5_infection_period[sim_index, ] = R_eff_4_5_avg
}



R_eff_samples_5_overall_df  = data.frame(
  date = standard_time_seq[2:(n_standard)],
  # infection_period = rep("5 days"),
  infection_period = 5,
  mean = apply(R_eff_7_day_avg_5_infection_period, 2, mean),
  upper = apply(R_eff_7_day_avg_5_infection_period, 2, quantile, 0.975),
  lower =  apply(R_eff_7_day_avg_5_infection_period, 2, quantile, 0.025)
)

R_eff_samples_4_75_overall_df  = data.frame(
  date = standard_time_seq[2:(n_standard)],
  # infection_period = rep("5 days"),
  infection_period = 4.75,
  mean = apply(R_eff_7_day_avg_4_75_infection_period, 2, mean),
  upper = apply(R_eff_7_day_avg_4_75_infection_period, 2, quantile, 0.975),
  lower =  apply(R_eff_7_day_avg_4_75_infection_period, 2, quantile, 0.025)
)

R_eff_samples_4_5_overall_df  = data.frame(
  date = standard_time_seq[2:(n_standard)],
  # infection_period = rep("5 days"),
  infection_period = 4.5,
  mean = apply(R_eff_7_day_avg_4_5_infection_period, 2, mean),
  upper = apply(R_eff_7_day_avg_4_5_infection_period, 2, quantile, 0.975),
  lower =  apply(R_eff_7_day_avg_4_5_infection_period, 2, quantile, 0.025)
)




death_real_overall_df = data.frame(
  date = standard_time_seq,
  type = rep("Real Data"),
  # infection_period = "5 days",
  infection_period = 0,
  mean = US_death_overall,
  upper = NA,
  lower =  NA
)


death_samples_5_overall_df = data.frame(
  date = standard_time_seq,
  type = "Estimation",
  # infection_period = rep("5 days"),
  infection_period = 5,
  mean = param_record_approx_for_beta[4,], #apply(Death_5_record_matrix, 2, mean),
  upper = apply(Death_5_record_matrix, 2, quantile, 0.975),
  lower =  apply(Death_5_record_matrix, 2, quantile, 0.025)
)

death_samples_4_75_overall_df = data.frame(
  date = standard_time_seq,
  type = "Estimation",
  # infection_period = rep("4.75 days"),
  infection_period = 4.75,
  mean = apply(Death_4_75_record_matrix, 2, mean),
  upper = apply(Death_4_75_record_matrix, 2, quantile, 0.975),
  lower =  apply(Death_4_75_record_matrix, 2, quantile, 0.025)
)

death_samples_4_5_overall_df = data.frame(
  date = standard_time_seq,
  type = "Estimation",
  # infection_period = rep("4.5 days"),
  infection_period = 4.5,
  mean = apply(Death_4_5_record_matrix, 2, mean),
  upper = apply(Death_4_5_record_matrix, 2, quantile, 0.975),
  lower =  apply(Death_4_5_record_matrix, 2, quantile, 0.025)
)

infectious_samples_5_overall_df = data.frame(
  date = standard_time_seq,
  # infection_period = rep("5 days"),
  infection_period = 5,
  mean = apply(Infectious_5_record_matrix, 2, mean),
  upper = apply(Infectious_5_record_matrix, 2, quantile, 0.975),
  lower =  apply(Infectious_5_record_matrix, 2, quantile, 0.025)
)

infectious_samples_4_75_overall_df = data.frame(
  date = standard_time_seq,
  # infection_period = rep("4.75 days"),
  infection_period = 4.75,
  mean = apply(Infectious_4_75_record_matrix, 2, mean),
  upper = apply(Infectious_4_75_record_matrix, 2, quantile, 0.975),
  lower =  apply(Infectious_4_75_record_matrix, 2, quantile, 0.025)
)

infectious_samples_4_5_overall_df = data.frame(
  date = standard_time_seq,
  # infection_period = rep("4.5 days"),
  infection_period = 4.5,
  mean = apply(Infectious_4_5_record_matrix, 2, mean),
  upper = apply(Infectious_4_5_record_matrix, 2, quantile, 0.975),
  lower =  apply(Infectious_4_5_record_matrix, 2, quantile, 0.025)
)

R_eff_sample_overall_df = rbind(R_eff_samples_5_overall_df, R_eff_samples_4_75_overall_df, R_eff_samples_4_5_overall_df)

death_sample_overall_df = rbind(death_real_overall_df, death_samples_5_overall_df, death_samples_4_75_overall_df, death_samples_4_5_overall_df)

infectious_overall_df = rbind(infectious_samples_5_overall_df, infectious_samples_4_75_overall_df, infectious_samples_4_5_overall_df)

R_eff_sample_overall_df$infection_period = as.factor(R_eff_sample_overall_df$infection_period)

death_sample_overall_df$type = as.factor(death_sample_overall_df$type)
death_sample_overall_df$infection_period = as.factor(death_sample_overall_df$infection_period)

infectious_overall_df$infection_period = as.factor(infectious_overall_df$infection_period)

# figure for R_eff

nation_R_eff_plot = ggplot(R_eff_sample_overall_df, aes(x = date, y = mean))+
  geom_line(aes(group = infection_period, col = infection_period, linetype=infection_period, size = infection_period),show.legend = T) +
  scale_linetype_manual(values=c(rep("solid", 3))) +
  geom_ribbon(aes(ymin=lower,ymax=upper, fill = (infection_period)),alpha=0.25, show.legend = F)+
  scale_size_manual(values=c(0.7, 0.7, 0.7))+ 
  scale_colour_manual(values = c("red", "#009900", "blue"))+
  scale_fill_manual(values = c("#F8766D","#5cd65c","steelblue2"))+
  # scale_fill_manual(values=c("#F8766D",NA,"#00BFC4"))+
  xlab("Date")+
  ylab("Effective reproduction number") + 
  labs(color='Infectious period')+
  guides(colour = FALSE,#guide_legend(reverse=T),
         linetype = FALSE,
         size = FALSE,
         fill = FALSE)+
  geom_hline(yintercept = 1, linetype = "dashed") + 
  theme(text = element_text(size = 20),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.key.width=unit(1,"cm"),
        axis.text.y = element_text(angle=90, hjust=1)) 

nation_infect_plot = ggplot(infectious_overall_df, aes(x = date, y = mean / 10^6 ))+
  geom_line(aes(group = infection_period, col = infection_period, linetype=infection_period, size = infection_period),show.legend = T) +
  scale_linetype_manual(values=c(rep("solid", 3))) +
  geom_ribbon(aes(ymin=lower / 10^6,ymax=upper / 10^6, fill = (infection_period)),alpha=0.25, show.legend = F)+
  scale_size_manual(values=c(0.7, 0.7, 0.7))+ 
  scale_colour_manual(values = c("red", "#009900", "blue"))+
  scale_fill_manual(values = c("#F8766D","#5cd65c","steelblue2"))+
  # scale_fill_manual(values=c("#F8766D",NA,"#00BFC4"))+
  xlab("Date")+
  ylab(expression(paste('Active infectious individuals  ', (10^6)))) + 
  labs(color='Infectious period')+
  guides(colour = FALSE,#guide_legend(reverse=T),
         linetype = FALSE,
         size = FALSE,
         fill = FALSE) + 
  theme(text = element_text(size = 20),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.key.width=unit(1,"cm"),
        axis.text.y = element_text(angle=90, hjust=1)) 

nation_death_with_legend_plot = ggplot(death_sample_overall_df, aes(x = date, y = mean/ 10^3))+
  geom_line(aes( col = infection_period, linetype=type, size = type),show.legend = T) +
  # scale_linetype_manual(values=c(rep("solid", 3), rep(c("dotted"), 1))) +
  geom_ribbon(aes(ymin=lower/ 10^3,ymax=upper/ 10^3, fill = (infection_period)),alpha=0.25, show.legend = T)+
  scale_colour_manual(name = "Infectious period",
                      labels = (c( "4.5 days", "4.75 days", "5 days", "Real data")),
                      breaks = c(4.5,4.75,5,10 ),
                      values = c("red", "#009900", "blue", "black")) +   
  scale_linetype_manual(name = "Type",
                        labels = (c("Estimation", "Observation")),
                        values = (c("solid", "dashed")))+
  scale_size_manual(name = "Type",
                    labels = (c("Estimation", "Observation")),
                    values = (c(0.7, 1)))+
  scale_fill_manual(name = "Infectious period",
                    labels = (c( "4.5 days", "4.75 days", "5 days", "Real data")),
                    breaks = c(4.5,4.75,5,10 ),
                    values=c("#F8766D","#5cd65c","steelblue2", "#555555"))+
  xlab("Date")+
  ylab(expression(paste("Cumulative death toll  ", (10^3)))) + 
  labs(color='Infectious period')+
  guides(colour = guide_legend(reverse=T),
         linetype = guide_legend(reverse=T, override.aes=list(fill=NA)),
         size = guide_legend(reverse=T, override.aes=list(fill=NA)),
         fill = guide_legend(reverse=T)) +
  geom_point(data = death_real_overall_df, aes(x = date, y = mean/ 10^3, group = type),size = 0.1, colours = "#4d2600")+
  theme(text = element_text(size = 20),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.key.width=unit(1,"cm"),
        axis.text.y = element_text(angle=90, hjust=1))

nation_death_plot = ggplot(death_sample_overall_df, aes(x = date, y = mean / 10^3))+
  geom_line(aes( col = infection_period, linetype=type, size = type),show.legend = T) +
  # scale_linetype_manual(values=c(rep("solid", 3), rep(c("dotted"), 1))) +
  geom_ribbon(aes(ymin=lower / 10^3,ymax=upper / 10^3, fill = (infection_period)),alpha=0.25, show.legend = F)+
  scale_colour_manual(name = "Infectious period",
                      labels = (c( "4.5 days", "4.75 days", "5 days", "Real data")),
                      breaks = c(4.5,4.75,5,10 ),
                      values = c("red", "#009900", "blue", "black")) +   
  scale_linetype_manual(name = "Type",
                        labels = (c("Estimation", "Observation")),
                        values = (c("solid", "dashed")))+
  scale_size_manual(name = "Type",
                    labels = (c("Estimation", "Observation")),
                    values = (c(0.7, 1)))+
  scale_fill_manual(name = "Infectious period",
                    labels = (c( "4.5 days", "4.75 days", "5 days", "Real data")),
                    breaks = c(4.5,4.75,5,10 ),
                    values=c("#F8766D","#5cd65c","steelblue2", "#555555"))+
  xlab("Date")+
  ylab(expression(paste("Cumulative death toll  ", (10^3)))  ) + 
  labs(color='Infectious period in the national level')+
  guides(colour = FALSE,#guide_legend(reverse=T),
         linetype = FALSE, #guide_legend(reverse=T),
         size = FALSE,
         fill = FALSE)+ #guide_legend(reverse=T)) + 
  geom_line(data = death_real_overall_df, aes(x = date, y = mean/ 10^3, group = type),linetype = "dashed", size = 1.5) +
  theme(text = element_text(size = 20),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.key.width=unit(1,"cm"),
        axis.text.y = element_text(angle=90, hjust=1)) 


g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

###############################
# Arrange figures and save it as EPS file
###############################

nation_legend<-g_legend(nation_death_with_legend_plot)

save_file_path = "./Reproducing results in the paper/Results/Fig_4/"

ggarrange(nation_R_eff_plot, nation_infect_plot, nation_death_plot, nation_legend, nrow = 1, ncol = 4, widths = c(5,5,5,2), labels = c("a", "b", "c", NA)) %>%
ggsave(file=paste0(save_file_path, "US_different_infection_period.eps"), device=cairo_ps, fallback_resolution = 600, width = 40, height = 12, units = "cm")
