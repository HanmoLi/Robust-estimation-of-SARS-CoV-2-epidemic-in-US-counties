##########################
# Codes for estimating power parameters in US stats.
# The estimation is from March 21, 2020 to Sep 20, 2020.
# 
# Author: Hanmo Li & Mengyang Gu
#
# Contact info: mengyang at pstat.ucsb.edu
#
# Paper: Robust estimation of SARS-CoV-2 epidemic in US counties (https://arxiv.org/pdf/2010.11514.pdf)
##########################

library(deSolve)
library(dplyr)
library(stringr)
library(zoo)
library(mFilter)
library(RobustGaSP)
library(foreach)
library(doParallel)
library(rjson)

source(file = "./Forecast in US counties/Codes/functions_8th_version.R")


start_date_ini = as.Date("2020-3-21")

end_date = as.Date("2020-9-20")

fitted_days_beta=184

set.seed(1)

gamma = 0.2
theta = 0.1
delta = 0.0066

c_t_type = "test positive rate" # test number or test positive rate

# The location where we save the power parameter file
save_pow_param_file_path = "./Forecast in US counties/Data/"

# read the file for optim power parameter in US states
us_pow_param = read.csv(paste0("./Forecast in US counties/Data/Name_list_of_US_states_updated.csv"))
us_pow_param$state_name = as.character(us_pow_param$state_name)
us_pow_param$state_name_short = as.character(us_pow_param$state_name_short)

base = 'https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/'
death = 'time_series_covid19_deaths_'
confirm = 'time_series_covid19_confirmed_'

# download data from JHU
us_death = read.csv(paste0(base,death,"US.csv"))
us_confirm = read.csv(paste0(base,confirm,"US.csv"))

JHU_testing_url = "https://raw.githubusercontent.com/govex/COVID-19/master/data_tables/testing_data/time_series_covid19_US.csv"

JHU_testing <- read.csv(file = JHU_testing_url)

extracted_JHU_date = JHU_testing$date
Year_str_JHU = str_extract_all(extracted_JHU_date, "\\d+")
Year_str_JHU_adjusted = c()

for(i in 1:length(Year_str_JHU)){
  Year_str_JHU_each = Year_str_JHU[[i]]
  Year_str_JHU_each[3] = Year_str_JHU_each[3]
  Year_str_JHU_adjusted = c(Year_str_JHU_adjusted, paste(Year_str_JHU_each,collapse = '-'))
}

JHU_testing$date = as.Date(Year_str_JHU_adjusted,  tryFormats = c("%m-%d-%Y"))

JHU_testing = JHU_testing %>%
  group_by(state) %>%
  mutate(positiveIncrease = c(NA, diff(cases_conf_probable)), totalTestResultsIncrease = c(NA, diff(tests_combined_total ))) %>%
  ungroup() %>%
  filter(date > as.Date("2020-03-10"))

us_death = clean_us_death(us_death)

col_names = colnames(us_death)

extracted_str_date = str_sub(col_names[13:dim(us_death)[2]], 2, -1)
Year_str = str_extract_all(extracted_str_date, "\\d+")
Year_str_adjusted = c()

for(i in 1:length(Year_str)){
  Year_str_each = Year_str[[i]]
  Year_str_each[3] = paste0("20", Year_str_each[3])
  Year_str_adjusted = c(Year_str_adjusted, paste(Year_str_each,collapse = '-'))
}

all_dates = as.Date(Year_str_adjusted,  tryFormats = c("%m-%d-%Y"))

rmse_record_list = as.list(1:dim(us_pow_param)[1])

pow_param_record = rep(NA, dim(us_pow_param)[1])

state_pow_param_df = data.frame(state_name = us_pow_param$state_name, state_name_short = us_pow_param$state_name_short, pow_param = rep(NA, dim(us_pow_param)[1]))

for(state_index in 1:dim(us_pow_param)[1]){
  
  state_name = us_pow_param$state_name[state_index]
  state_name_short = us_pow_param$state_name_short[state_index]
  
  print(paste0("Start calculating power parameter for ", state_name, ": ", state_index, " out of ", dim(us_pow_param)[1] ))
  
  state_test = JHU_testing %>%
    dplyr::filter(state == state_name_short) %>%
    dplyr::select(date, totalTestResultsIncrease, positiveIncrease)
  
  state_test$totalTestResultsIncrease[is.na(state_test$totalTestResultsIncrease)] = 1
  state_test$positiveIncrease[is.na(state_test$positiveIncrease)] = 1
  
  state_test$totalTestResultsIncrease[state_test$totalTestResultsIncrease<0] = abs(state_test$totalTestResultsIncrease[state_test$totalTestResultsIncrease<0])
  
  state_test$positiveIncrease[state_test$positiveIncrease<0] = abs(state_test$positiveIncrease[state_test$positiveIncrease<0])
  
  if (state_name == "Washington"){
    state_test$totalTestResultsIncrease[state_test$date>= as.Date("2020-11-23") & state_test$date<= as.Date("2020-12-05")] = state_test$totalTestResultsIncrease[state_test$date == as.Date("2020-12-05")] / 13
  }
  
  # impute the outlier at March-3-2021 with the 6-day average around it.
  date_seq_temp = c(seq(as.Date("2021-02-28"), as.Date("2021-03-2"), by=1), seq(as.Date("2021-03-4"), as.Date("2021-03-6"), by=1))
  state_test$positiveIncrease[state_test$date==as.Date("2021-03-03")] = mean(state_test$positiveIncrease[state_test$date %in% date_seq_temp])
  
  state_test$daily_test_avg = rep(0, dim(state_test)[1])
  state_test$daily_positive_avg = rep(0, dim(state_test)[1])

  state_test$daily_test_avg = data_seven_day_smoothing(state_test$totalTestResultsIncrease)
  
  state_test$daily_positive_avg = data_seven_day_smoothing(state_test$positiveIncrease)
  
  
  state_test$PositiveRate = state_test$daily_positive_avg / state_test$daily_test_avg
  
  for (i in 2:(dim(state_test)[1])) {
    if(is.na(state_test$PositiveRate[i]) | is.infinite(state_test$PositiveRate[i])){
      state_test$PositiveRate[i] = state_test$PositiveRate[i-1]
    }
  }
  
  state_test$PositiveRate[is.na(state_test$PositiveRate)] = state_test$PositiveRate[which(!is.na(state_test$PositiveRate))[1]]
  
  state_test$PositiveRate = abs(state_test$PositiveRate)
  # change positive rate larger tha 1 to 1
  state_test$PositiveRate[state_test$PositiveRate>1] = 1
  
  state_death = us_death %>%
    dplyr::filter(Province_State == state_name)
  
  state_confirmed = us_confirm %>%
    dplyr::filter(Province_State == state_name)
  
  state_death_sum = apply(state_death[,12:dim(state_death)[2]], 2, sum)
  
  state_confirm_sum = apply(state_confirmed[,12:dim(state_confirmed)[2]], 2, sum)
  
  start_date = all_dates[which(as.numeric(state_confirm_sum)>=5)[1]]

  if (start_date < start_date_ini){
    start_date = start_date_ini
  }
  
  n = as.numeric(end_date - start_date) + 1
  
  state_death_selected = state_death_sum[1 + which(all_dates %in% seq.Date(start_date, end_date, by=1))]
  
  state_confirmed_selected = state_confirm_sum[which(all_dates %in% seq.Date(start_date, end_date, by=1))]
  
  N = as.numeric(state_death_sum[1])
  
  death_selected = as.numeric(state_death_selected)
  confirm_selected = as.numeric(state_confirmed_selected)
  
  # eliminate decreasing trend in the county death data
  for( i in 1:(n-1)){
    if (death_selected[i+1]<death_selected[i]){
      death_selected[i+1] = death_selected[i]
    }
  }
  
  for( i in 1:(n-1)){
    if (confirm_selected[i+1]<confirm_selected[i]){
      confirm_selected[i+1] = confirm_selected[i]
    }
  }
  
  daily_confirm_selected = diff(confirm_selected)
  daily_confirm_selected_avg = data_seven_day_smoothing(daily_confirm_selected)
  unadjusted_confirm_selected_smoothed = c(confirm_selected[1], cumsum(daily_confirm_selected_avg) + confirm_selected[1])

  daily_positive_rate_selected = (state_test$PositiveRate[state_test$date>=(start_date) & state_test$date<=end_date])

  alpha_upper_bound = 2
  
  alpha_ini_seq = seq(0,alpha_upper_bound,length.out = 1000)

  upper_bound_record = rep(NA, length(alpha_ini_seq))
  
  for(alpha_index in 1:length(alpha_ini_seq)){
    c_t_seq_each = daily_positive_rate_selected^(alpha_ini_seq[alpha_index])
    
    daily_confirm_selected_each = diff(confirm_selected)
    
    daily_confirm_selected_each_smoothed = data_seven_day_smoothing(daily_confirm_selected_each * c_t_seq_each[2:n])
    
    confirm_selected_smoothed_each = c( c_t_seq_each[1] * confirm_selected[1] , cumsum(daily_confirm_selected_each_smoothed)+c_t_seq_each[1]*confirm_selected[1])
    
    upper_bound_each = ((confirm_selected_smoothed_each[1]*N)/confirm_selected_smoothed_each[n] - death_selected[1])
    
    upper_bound_record[alpha_index] = upper_bound_each
  }
  upper_final = min(upper_bound_record)
  
  if (upper_final<10){
    upper_final = sort(upper_bound_record[upper_bound_record>10])[1]
    alpha_upper_bound = alpha_ini_seq[upper_bound_record == upper_final]
  }

  ui = matrix(c(c(1,0,-1,0,1, 0,0),c( 0, 1, 0,-1,-1, 0,0), c(0, 0, 0,0,0, 1, -1)), 7,3)
  ci = matrix(c(-10^(10),-10^(10), -log(upper_final),-log(upper_final),0, -10^(10) , -log(2) ))
  
  I_0_ini = 100
  R_0_ini = 100
  
  if((I_0_ini + R_0_ini)> abs(upper_final)){
    I_0_ini = abs(upper_final)/4
    R_0_ini = abs(upper_final)/4
  }
  
  alpha_vector = seq(0.1,alpha_upper_bound-0.1,0.1)#c(0.1,0.3, 0.5, 0.7,0.9)
  
  error_record = rep(NA, length(alpha_vector))
  m_approx_record = as.list( rep(NA, length(alpha_vector)))
  
  pb <- txtProgressBar(min = 0, max = length(alpha_vector), style = 3)
  for(alpha_index in 1:length(alpha_vector)){
    
    setTxtProgressBar(pb, alpha_index)
    
    skip_to_next <- FALSE
    tryCatch({
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
      training_length = n,
      fixed_global_params = c(gamma, theta, delta),
      penalty = F,
      fitted_days= fitted_days_beta,
      weight_loss=T
    )}
    , error = function(e) { skip_to_next <<- TRUE})
    if(skip_to_next) {
      next 
    }else{
      error_record[alpha_index] = m_approx_each$value
      m_approx_record[[alpha_index]] = m_approx_each
    }
  }
  m_approx = m_approx_record[[which(error_record==min(error_record, na.rm = T))[1]]]
  
  I_0_plus_R_0 = exp(m_approx$par[1])
  R_0 = exp(m_approx$par[2])
  exp_index = exp(m_approx$par[3])
  
  I_0 = I_0_plus_R_0 - R_0

  state_pow_param_df$pow_param[state_index] = exp_index
  
  c_t_seq = daily_positive_rate_selected^(exp_index)
  
  daily_confirm_selected = diff(confirm_selected)
  
  daily_confirm_selected_smoothed = data_seven_day_smoothing(daily_confirm_selected * c_t_seq[2:n])
  
  confirm_selected_smoothed =  c( c_t_seq[1] * confirm_selected[1] , cumsum(daily_confirm_selected_smoothed)+c_t_seq[1]*confirm_selected[1])

  if (I_0<1){
    I_0 = 1
  }
  if (R_0<1){
    R_0 = 1
  }
  
  ratio = confirm_selected_smoothed[1]/(I_0+ R_0+ death_selected[1])
  
  ratio_real = confirm_selected[2]/(I_0+ R_0+ death_selected[1])
  
  estimated_confirm = confirm_selected_smoothed/ratio
  
  for(i in 1:length(estimated_confirm)){
    if (estimated_confirm[i] < unadjusted_confirm_selected_smoothed[i]){
      estimated_confirm[i] = unadjusted_confirm_selected_smoothed[i]
    }
  }
  
  S_t_seq = N - estimated_confirm
  
  init_for_beta = c(S_t_seq[1], I_0, R_0, death_selected[1], 0)
  
  param_record_approx_for_beta = matrix(0, 5, n)
  param_record_approx_for_beta[,1] = init_for_beta
  param_record_approx_for_beta[1,] = S_t_seq
  
  approx_beta_seq = rep(0, n-1)
  
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
    if(R_t_1<1){
      R_t_1 = 1
    }
    
    beta_t_1_2 = uniroot(find_root_beta, c(0, 10^6), tol = 0.0001, param = c(S_t_1, S_t_2, I_t_1), N = N, gamma=gamma)

    I_t_2 = I_t_1 * exp(beta_t_1_2$root*(S_t_1 + S_t_2)/(2*N) - gamma)
    R_t_2 = (2-theta)/(2+theta)*R_t_1 + gamma/(2+theta)*(I_t_1+I_t_2)
    D_t_2 = D_t_1 + delta*theta*(R_t_1+R_t_2)/2
    C_t_2 = C_t_1 + (1-delta)*theta*(R_t_1+R_t_2)/2
    
    param_record_approx_for_beta[2:5, i+1] = c(I_t_2, R_t_2, D_t_2, C_t_2)
    approx_beta_seq[i] = beta_t_1_2$root
  }
  
  approx_beta_seq_smoothed = data_seven_day_smoothing(approx_beta_seq)

  date_seq = seq.Date(start_date, end_date, by=1)
  
  ylimit_death = c(min(param_record_approx_for_beta[4,], death_selected), max(param_record_approx_for_beta[4,], death_selected))
  plot(param_record_approx_for_beta[4,]~date_seq,ylim = ylimit_death, type="l", col="blue", xlab = "Date", ylab = "Death Cases", main = paste0(state_name_short, ", population=", round(N/10^6,2),"M", ", pow_param = ", round(exp_index,3)))
  lines(death_selected~date_seq, col = "red")
  legend("topleft", legend = c("Observed Death Cases", "Estimated Death Cases"), lty = c(1,1), col = c("red", "blue"))

  date_seq_beta = seq.Date(start_date, end_date-1, by=1)

  y_limit_R_t = c(min(approx_beta_seq/0.2), max(approx_beta_seq/0.2))
  plot(approx_beta_seq/0.2~date_seq_beta, type="p",ylim = y_limit_R_t, ylab = "R_t", xlab = "Date", main = paste0(state_name, ", population=", round(N/10^6,2),"M", ", Ratio = ", round(ratio_real,3)))
  lines(approx_beta_seq_smoothed/0.2~date_seq_beta, col="blue")
  abline(h=1, lty=2,col="red")
  abline(h=0,lty=1,col="black")
  legend("topright", legend = c("Daily Approximated R_t", "7-day Averaged Approximated R_t"), pch = c(1,NA), lty=c(NA,1), col=c("black", "blue"))
}

state_pow_param_df$pow_param = round(state_pow_param_df$pow_param,3)

write.csv(state_pow_param_df, file = paste0(save_pow_param_file_path, "power_parameters_in_US_states_as_of_Sep_20_2020.csv"), row.names = F)

















