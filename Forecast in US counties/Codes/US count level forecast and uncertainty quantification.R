##########################
# Parallel codes for computing forecasts and CIs in US counties. 
# The forecasts start from Seq 20, 2020.
# 
# Author: Hanmo Li & Mengyang Gu
#
# Contact info: mengyang at pstat.ucsb.edu
#
# Paper: Robust estimation of SARS-CoV-2 epidemic in US counties (https://arxiv.org/pdf/2010.11514.pdf)
##########################

library(deSolve)
library(stringr)
library(zoo)
library(ggplot2)
library(mFilter)
library(RobustGaSP)
library(matrixStats)
library(parallel)
library(foreach)
library(doParallel)
library(dplyr)

source(file = "./Forecast in US counties/Codes/functions_8th_version.R")

start_date_ini = as.Date("2020-03-21")
end_date = as.Date("2020-09-20")

set.seed(1)

# parallel settings
numCores <- detectCores()

gamma = 0.2
theta = 0.1
delta = 0.0066

# for forecast
fitted_days_beta=184 # how many days from the end day you use to estimate beta_t
predicted_days_train=184  

prediction_length = 90
beta_t_length=184

base = 'https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/'
death = 'time_series_covid19_deaths_'
confirm = 'time_series_covid19_confirmed_'

JHU_testing_url = "https://raw.githubusercontent.com/govex/COVID-19/master/data_tables/testing_data/time_series_covid19_US.csv"

file_pow_param_path = "./Forecast in US counties/Data/"
file_sudden_death_path = "./Forecast in US counties/Data/"
file_path_cur= "./Forecast in US counties/Results/"

# counties that have sudden death increase
county_sudden_death = read.csv(paste0(file_sudden_death_path,"counties_with_sudden_death_increase.csv" ) )
county_sudden_death$state = as.character(county_sudden_death$state)

# read the file for optim power parameter in US states
us_pow_param = read.csv(paste0(file_pow_param_path, "power_parameters_in_US_states_as_of_Sep_20_2020.csv"))
us_pow_param$state_name = as.character(us_pow_param$state_name)

# download death and confirmed data from JHU
us_death = read.csv(paste0(base,death,"US.csv"))
us_confirm = read.csv(paste0(base,confirm,"US.csv"))

# clean the death data
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

death_rate = function(x){x/us_death$Population*10^6}
us_death_rate = us_death %>%
  dplyr::mutate_at(vars(starts_with("X")), death_rate)

n_ini = as.numeric(end_date - start_date_ini) + 1

JHU_testing <- read.csv(file = JHU_testing_url)

extracted_JHU_date = JHU_testing$date
Year_str_JHU = str_extract_all(extracted_JHU_date, "\\d+")
Year_str_JHU_adjusted = c()

for(i in 1:length(Year_str_JHU)){
  Year_str_JHU_each = Year_str_JHU[[i]]
  Year_str_JHU_adjusted = c(Year_str_JHU_adjusted, paste(Year_str_JHU_each,collapse = '-'))
}

JHU_testing$date = as.Date(Year_str_JHU_adjusted,  tryFormats = c("%m-%d-%Y"))

JHU_testing = JHU_testing %>%
  group_by(state) %>%
  mutate(positiveIncrease = c(NA, diff(cases_conf_probable)), totalTestResultsIncrease = c(NA, diff(tests_combined_total ))) %>%
  ungroup() %>%
  filter(date > as.Date("2020-03-10"))

avaliable_state_name = c()
avaliable_state_name_short = c()
n_county_death_over_2 = 0
for(state_index in 1:dim(us_pow_param)[1]){
  state_name = as.character(us_pow_param$state_name[state_index])
  state_name_short = as.character(us_pow_param$state_name_short[state_index])
  
  skip_to_next <- FALSE
  tryCatch({    death_and_county_names = get_output_same_time_zone(
    data_type = "death",
    state_name = state_name,
    state_name_short = state_name_short,
    start_date = start_date_ini,
    duration = n_ini,
    training_length = n_ini,
    criterion_death = 2,
    smoothness = F
  )}
  , error = function(e) { skip_to_next <<- TRUE})
  if(skip_to_next) {
    print(paste0(state_name, " No DATA"))
    next 
  }else{
    n_county_death_over_2 = n_county_death_over_2 + length(death_and_county_names[[2]])
    print(paste0(state_name, ": ", length(death_and_county_names[[2]]), " counties"))
    avaliable_state_name = c(avaliable_state_name, state_name)
    avaliable_state_name_short = c(avaliable_state_name_short, state_name_short)
  }
}
n_county_death_over_2

registerDoParallel(numCores)  # use multicore, set to the number of our cores

system.time({
  
  state_results_par = foreach(state_index = c(11, 51), .packages = c("dplyr", "stats", "RobustGaSP", "zoo", "matrixStats", "stringr")) %dopar%{
    
    state_name = avaliable_state_name[state_index]
    state_name_short = avaliable_state_name_short[state_index]
    
    power_parameter = as.numeric(us_pow_param$pow_param[us_pow_param$state_name == state_name])
    
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
    
    death_and_county_names = get_output_same_time_zone(
      data_type = "death",
      state_name = state_name,
      state_name_short = state_name_short,
      start_date = start_date_ini,
      duration = n_ini,
      training_length = n_ini,
      criterion_death = 2,
      smoothness = F
    )
    
    county_names = as.character(death_and_county_names[[2]])
    
    # check which county has death equals to 2
    county_names[which(death_and_county_names[[1]][, dim(death_and_county_names[[1]])[2]] == 2)]
    
    county_index_test = which(death_and_county_names[[1]][, dim(death_and_county_names[[1]])[2]] < 5)
    
    # get populations of selected counties
    SIRDC_county_population = us_death %>%
      filter(Admin2 %in% county_names, Province_State == state_name) %>%
      dplyr::select(Population)
    
    results_overall_list = as.list(rep(NA, length(county_names)))
    
    RMSE_pre_county = rep(NA, length(county_names))
    CI_length_pre_county = rep(NA, length(county_names))
    CI_coverage_pre_county = rep(NA, length(county_names))
    
    time = system.time({
      for (county_index in 1:length(county_names)){
        
        print(paste0("The ", county_index, "th county out of ", length(county_names),": ", county_names[county_index], " in ", state_name))

        county_name_selected = county_names[county_index]
        
        #########################
        # Part 1: Data Cleaning
        #########################
        
        county_death_raw = state_death %>%
          dplyr::filter(Admin2 == county_name_selected)
        
        county_confirm_raw = state_confirmed %>%
          dplyr::filter(Admin2 == county_name_selected)
      
        #cumulative confirmed cases>10
        start_date=all_dates[which(as.numeric(county_confirm_raw[12:length(county_confirm_raw)])>=5)[1]]
        
        if (start_date < start_date_ini){
          start_date = start_date_ini
        }
        
        n = as.numeric(end_date - start_date) + 1
        
        # get county data from real start_date to end_date
        death_with_county_names = get_output_same_time_zone(
          data_type = "death",
          state_name = state_name,
          start_date = start_date,
          state_name_short = state_name_short,
          duration = n ,
          training_length = n,
          criterion_death = 2,
          smoothness = F
        )
        
        county_names = death_with_county_names[[2]]
        
        each_index = which(county_names == county_name_selected)
        
        SIRDC_county_population = us_death %>%
          filter(Admin2 %in% county_names, Province_State == state_name) %>%
          dplyr::select(Population)
        
        state_confirmed_selected = state_confirmed %>%
          dplyr::filter(Admin2 %in% county_names)
        
        state_confirmed_selected = state_confirmed_selected[,12:dim(state_confirmed_selected)[2]]
        
        state_confirmed_selected = state_confirmed_selected[, which(all_dates %in% seq.Date(start_date, end_date, by=1))]
        
        N = SIRDC_county_population$Population[each_index]
        
        death_selected = death_with_county_names[[1]][each_index, ]
        confirm_selected = as.numeric(state_confirmed_selected[each_index,])
        
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
        
        if(paste0(state_name_short, county_name_selected) %in%   paste0(county_sudden_death$state, county_sudden_death$county)){
          daily_death = diff(death_selected)
          sudden_increase_index = which(daily_death==max(daily_death))[1]
          daily_death[(sudden_increase_index-13):(sudden_increase_index-1)] = daily_death[(sudden_increase_index-13):(sudden_increase_index-1)] + daily_death[sudden_increase_index]/14
          daily_death[sudden_increase_index] = daily_death[sudden_increase_index]/14
          death_selected = c(death_selected[1], cumsum(daily_death) + death_selected[1])
        }
        
        daily_confirm_selected = diff(confirm_selected)
        daily_confirm_selected_avg = data_seven_day_smoothing(daily_confirm_selected)
        unadjusted_confirm_selected_smoothed = c(confirm_selected[1], cumsum(daily_confirm_selected_avg) + confirm_selected[1])
        
        daily_positive_rate_selected = rep(0,n)
        
        daily_positive_rate_selected = (state_test$PositiveRate[state_test$date>=(start_date) & state_test$date<=end_date])
        
        # calculate the weights
        c_t_seq =  daily_positive_rate_selected^(power_parameter)
        
        daily_confirm_selected = diff(confirm_selected)
        
        daily_adjusted_confirm_smoothed = data_seven_day_smoothing(daily_confirm_selected * c_t_seq[2:n] )
        
        # adjust the confirmed cases using positive_rate^{power_parameter}
        confirm_selected_smoothed = c(c_t_seq[1] * confirm_selected[1], cumsum(daily_adjusted_confirm_smoothed) + c_t_seq[1] * confirm_selected[1])
        
        # eliminate decreasing trend in the smoothed county confirmed cases
        for(i in 2:n){
          if (confirm_selected_smoothed[i]<confirm_selected_smoothed[i-1]){
            confirm_selected_smoothed[i] = confirm_selected_smoothed[i-1]
          }
        }
        
        #########################
        # Part 2: Optimization
        #########################
        
        # try constrained optimization
        ui = matrix(c(1,0,-1,0,1,-1), 3,2)
        ci = matrix(c(0,0, -((confirm_selected_smoothed[1]*N)/confirm_selected_smoothed[n] -  death_selected[1]-0.1)  ))
        
        I_0_ini = 1000
        R_0_ini = 1000
        
        if((I_0_ini + R_0_ini)> abs(ci[3])){
          I_0_ini = abs(ci[3])/4
          R_0_ini = abs(ci[3])/4
        }
        
        # do optimization
        m_approx = constrOptim(
          c(I_0_ini, R_0_ini),
          loss_approx_beta,
          NULL,
          ui = ui,
          ci = ci,
          death_cases = death_selected,
          confirmed_cases = confirm_selected_smoothed,
          unadjusted_confirm_selected_smoothed = unadjusted_confirm_selected_smoothed,
          N_population = N,
          trianing_length = n,
          fixed_global_params = c(gamma, theta, delta),
          penalty = F,
          fitted_days= fitted_days_beta,
          weight_loss=T
        )
        
        I_0 = m_approx$par[1]
        R_0 = m_approx$par[2]

        #########################
        # Part 3: Parameter Estimation: S_t, I_t, R_t, D_t, C_t and beta_t
        #########################
        
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
        
        param_record_approx_for_beta = matrix(0, 5, n) # 5 rows: S_t, I_t, R_t, D_t, C_t
        param_record_approx_for_beta[,1] = init_for_beta
        param_record_approx_for_beta[1,] = S_t_seq
        
        approx_beta_seq = rep(0, n-1) # record the value of transmission rate
        
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

        approx_beta_seq[approx_beta_seq<0] = 0
        
        # Visualize the point estimations from the approximation algorithm
        
        # calculate the smoothed transmission rate using 7 day average
        # approx_beta_seq_smoothed = rollapply(approx_beta_seq, width = 7, by = 1, FUN = mean, align = "left")
        # 
        # date_seq_beta = seq.Date(start_date, end_date-1, by=1)
        # date_seq_beta_smoothed = seq.Date(start_date+3, end_date-1-3, by=1)
        # 
        # y_limit_R_t = c(min(beta_t_from_euler/0.2,approx_beta_seq/0.2), max(beta_t_from_euler/0.2,approx_beta_seq/0.2))
        # y_limit_R_t = c(min(approx_beta_seq/0.2), max(approx_beta_seq/0.2))
        #
        # plot(approx_beta_seq/0.2~date_seq_beta, type="p",ylim = y_limit_R_t, ylab = "R_t", xlab = "Date", main = paste0(county_names[each_index], ", population=", round(N/10^6,2),"M", ", Ratio = ", round(ratio_real,3)))
        # lines(approx_beta_seq_smoothed/0.2~date_seq_beta_smoothed, col="blue")
        # abline(h=1, lty=2,col="red")
        # abline(h=0,lty=1,col="black")
        # legend("topright", legend = c("Daily Approximated R_t", "7-day Averaged Approximated R_t"), pch = c(1,NA), lty=c(NA,1), col=c("black", "blue"))

        
        # date_seq = seq.Date(start_date, end_date, by=1)
        # date_seq_confirm = seq.Date(start_date-1, end_date,by=1)
        # 
        # ylimit_confirmed = c(min(log(N - param_record_approx_for_beta[1,]), log(confirm_selected[2:length(confirm_selected)])), max(log(N - param_record_approx_for_beta[1,]), log(confirm_selected[2:length(confirm_selected)])))
        # 
        # plot(log(N - param_record_approx_for_beta[1,])~date_seq,ylim = ylimit_confirmed, type="l", col="blue", xlab = "Date", ylab = "log(Confirmed Cases)", main = paste0(county_names[each_index], ", population=", round(N/10^6,2),"M", ", Ratio = ", round(ratio_real,3)))
        # lines(log(confirm_selected[2:length(confirm_selected)])~date_seq, col = "red")
        # legend("bottomright", legend = c("log Observed Confirmed Cases", "log Estimated Confirmed Cases"), lty = c(1,1), col = c("red", "blue"))

        
        # ylimit_IR = c(min(param_record_approx_for_beta[2,], param_record_approx_for_beta[3,]), max(param_record_approx_for_beta[2,], param_record_approx_for_beta[3,]))
        # 
        # plot(param_record_approx_for_beta[2,]~date_seq,ylim = ylimit_IR, type="l", col="blue", xlab = "Date", ylab = "Infective Cases", main = paste0(county_names[each_index], ", population=", round(N/10^6,2),"M", ", Ratio = ", round(ratio_real,3)))
        # lines(param_record_approx_for_beta[3,]~date_seq, col = "red")
        # legend("bottomright", legend = c("Estimated Resolved Cases", "Estimated Infective Cases"), lty = c(1,1), col = c("red", "blue"))

        
        # ylimit_death = c(min(param_record_approx_for_beta[4,], death_selected), max(param_record_approx_for_beta[4,], death_selected))
        # 
        # plot(param_record_approx_for_beta[4,]~date_seq,ylim = ylimit_death, type="l", col="blue", xlab = "Date", ylab = "Death Cases", main = paste0(county_names[each_index], ", population=", round(N/10^6,2),"M", ", Ratio = ", round(ratio_real,3)))
        # lines(death_selected~date_seq, col = "red")
        # legend("topleft", legend = c("Observed Death Cases", "Estimated Death Cases"), lty = c(1,1), col = c("red", "blue"))
        
        #######################################
        #Part 4: Death Prediction in 7 and 21 days
        ##########################################
        
        input_confirmed_cases= as.matrix(seq(1,n-1 ,1))
        
        output_confirmed_cases=as.matrix(log(diff(estimated_confirm)+1))
        
        m_gp_confirmed=rgasp(input_confirmed_cases, (output_confirmed_cases),
                             num_initial_values=5 ,
                             kernel_type = "pow_exp", alpha=1.9,
                             range.par=1,
                             nugget=.25
        )
        
        log_diff_confirm_pred = predict(m_gp_confirmed, input_confirmed_cases,interval_data=F)

        num_samples=500 # sample 100 samples from adjusted confirmed cases
        
        set.seed(each_index)
        
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
        
        sample_confirmed_cases_smoothed_record = matrix(NA, n, num_samples)
        for(sim_index in 1:num_samples){
          sample_confirm_each = sample_confirmed_cases[, sim_index]
          
          daily_sample_confirm_each = diff(sample_confirm_each)
          daily_sample_confirm_each_smoothed = data_seven_day_smoothing(daily_sample_confirm_each)
          confirm_selected_smoothed_sampled = c(sample_confirm_each[1], cumsum(daily_sample_confirm_each_smoothed)+sample_confirm_each[1])
          
          for(j in 1:(n-1)){
            if (confirm_selected_smoothed_sampled[j+1] < confirm_selected_smoothed_sampled[j]){
              confirm_selected_smoothed_sampled[j+1] = confirm_selected_smoothed_sampled[j]
            }
          }
          sample_confirmed_cases_smoothed_record[,sim_index] = confirm_selected_smoothed_sampled
        }
        
        crash_index = which(sample_confirmed_cases_smoothed_record[n,]>N)
        
        if(length(crash_index)>0){
          sample_confirmed_cases_smoothed_record = sample_confirmed_cases_smoothed_record[, -crash_index]
          num_samples = num_samples - length(crash_index)
        }
        
        beta_record_matrix = matrix(0, num_samples, n-1)
        
        Susceptive_record_matrix = matrix(0, num_samples, n)
        Infective_record_matrix = matrix(0, num_samples, n)
        Resolved_record_matrix = matrix(0, num_samples, n)
        Death_record_matrix = matrix(0, num_samples, n)
        ReCovered_record_matrix = matrix(0, num_samples, n)
        
        
        # Calculate transmission rates for 100 sampled confirmed and death cases
        for(sim_index in 1:num_samples){
          
          sample_death_each = death_selected
          confirm_selected_smoothed_sampled = sample_confirmed_cases_smoothed_record[, sim_index]
          
          estimated_confirm_sample = confirm_selected_smoothed_sampled #/ratio_sample
          
          S_t_seq_sample = N - estimated_confirm_sample
          
          Susceptive_record_matrix[sim_index, ] = S_t_seq_sample
          
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
          
          # Susceptive_record_sampled_beta[beta_index, ] = param_record_approx_all[1,]
          Infective_record_matrix[sim_index, ] = param_record_approx_all[2, ]
          Resolved_record_matrix[sim_index, ] = param_record_approx_all[3, ]
          Death_record_matrix[sim_index,] = param_record_approx_all[4, ]
          ReCovered_record_matrix[sim_index,] = param_record_approx_all[5, ]
        }

        seven_day_averaged_beta_record_matrix = matrix(NA, num_samples, n-1)
        for(i in 1:num_samples){
          seven_day_averaged_beta_record_matrix[i, ] = data_seven_day_smoothing(beta_record_matrix[i, ])
        }

        seven_day_averaged_beta_df = data.frame(date = seq.Date(start_date, end_date-1, by=1), mean = apply(seven_day_averaged_beta_record_matrix, 2, mean, na.rm=T), lower = apply(seven_day_averaged_beta_record_matrix, 2, quantile, 0.025, na.rm=T), upper = apply(seven_day_averaged_beta_record_matrix, 2, quantile, 0.975, na.rm=T))
        
        seven_day_avg_beta_limit = c(0, max(seven_day_averaged_beta_df$upper/0.2))
        
        infective_record_df = data.frame(date = seq.Date(start_date, end_date, by=1), mean = apply(Infective_record_matrix, 2, mean, na.rm=T), upper = apply(Infective_record_matrix, 2, quantile, 0.975, na.rm=T), 
                                         lower = apply(Infective_record_matrix, 2, quantile, 0.025, na.rm=T))
        
        infective_limit = c(0, max(infective_record_df$upper))

        prevalence_record_df = data.frame(date = seq.Date(start_date, end_date, by=1), mean = (N - apply(Susceptive_record_matrix, 2, mean, na.rm=T))/N * 100, lower = (N - apply(Susceptive_record_matrix, 2, quantile, 0.975, na.rm=T))/N * 100,
                                          upper = (N - apply(Susceptive_record_matrix, 2, quantile, 0.025, na.rm=T))/N * 100)
        
        
        prevalence_limit = c(0, max(prevalence_record_df$upper))

        ######## extrapolate beta_t in matrix form ################
        
        start_index_beta = which(approx_beta_seq<0.4)[1]
        
        input_beta = as.matrix(seq(start_index_beta,n-1 ,1))
        output_beta =as.matrix(rep(0, n - start_index_beta ))
        
        # here the rgasp function is only used to get Cholesky decomposition with fixed range and nugget parameter
        m_gp_beta_t = rgasp(input_beta, (output_beta), kernel_type = "pow_exp", alpha=1.9, range.par=1, nugget=0.25)
        
        L_beta = m_gp_beta_t@L
        L_beta_inverse = solve(L_beta)
        
        K_beta_inverse = L_beta_inverse %*% t(L_beta_inverse)
        
        one_matrix = matrix(1, n - start_index_beta, 1)
        
        mu_hat_beta = (t(one_matrix) %*% K_beta_inverse %*% t(beta_record_matrix[,start_index_beta:(n-1) ])) /  sum(K_beta_inverse)
        
        simga_2_hat_beta = diag(1/(n-start_index_beta) * (beta_record_matrix[,start_index_beta:(n-1) ] - t(one_matrix %*% mu_hat_beta)) %*% K_beta_inverse %*% t(beta_record_matrix[,start_index_beta:(n-1) ] - t(one_matrix %*% mu_hat_beta)))
        
        input_beta_pred = as.matrix(seq(start_index_beta,n+prediction_length-1 ,1))
        output_beta_pred =as.matrix(rep(0, n+prediction_length - start_index_beta))
        
        # here the rgasp function is only used to get Cholesky decomposition with fixed range and nugget parameter
        m_gp_beta_pred = rgasp(input_beta_pred, (output_beta_pred), kernel_type = "pow_exp", alpha=1.9, range.par=1, nugget=0.25)
        
        L_beta_pred = m_gp_beta_pred@L
        
        K_pred_samples = (L_beta_pred %*% t(L_beta_pred))[((n - start_index_beta + 1):(n+prediction_length-start_index_beta)), (1:(n-start_index_beta))]
        
        mean_beta_pred_given_samples = t(matrix(mu_hat_beta, num_samples, prediction_length)) - K_pred_samples %*% K_beta_inverse %*% t(beta_record_matrix[,start_index_beta:(n-1) ] - t(one_matrix %*% mu_hat_beta))
        
        beta_record_matrix_all = cbind(beta_record_matrix, t(mean_beta_pred_given_samples))

        ######## use extrapolated beta to calculate death prediction ################
        
        
        Susceptive_record_matrix_all = cbind(Susceptive_record_matrix, matrix(NA, num_samples, prediction_length))
        Infective_record_matrix_all = cbind(Infective_record_matrix, matrix(NA, num_samples, prediction_length))
        Resolved_record_matrix_all = cbind(Resolved_record_matrix, matrix(NA, num_samples, prediction_length))
        Death_record_matrix_all = cbind(Death_record_matrix, matrix(NA, num_samples, prediction_length))
        ReCovered_record_matrix_all = cbind(ReCovered_record_matrix, matrix(NA, num_samples, prediction_length))
        
        
        # Calculate transmission rates for 100 sampled confirmed and death cases
        for(sim_index in 1:num_samples){
          
          param_record_approx_all = rbind(Susceptive_record_matrix_all[sim_index, ],
                                          Infective_record_matrix_all[sim_index, ],
                                          Resolved_record_matrix_all[sim_index, ],
                                          Death_record_matrix_all[sim_index, ],
                                          ReCovered_record_matrix_all[sim_index, ])
          
          approx_beta_seq_all = beta_record_matrix_all[sim_index, ]
          
          for (i in n:(n+prediction_length-1)){
            S_t_1 = param_record_approx_all[1,i]
            I_t_1 = param_record_approx_all[2,i]
            R_t_1 = param_record_approx_all[3,i]
            D_t_1 = param_record_approx_all[4,i]
            C_t_1 = param_record_approx_all[5,i]
            
            beta_t_1_2 = approx_beta_seq_all[i]
            
            if(I_t_1<1){
              I_t_1 = 1
            }
            
            S_t_2 = uniroot(find_root_S_t_2, c(0, N), tol = 0.0001, param = c(S_t_1, beta_t_1_2, I_t_1), N = N, gamma=gamma)
            I_t_2 = I_t_1 * exp(beta_t_1_2*(S_t_1 + S_t_2$root)/(2*N) - gamma)
            R_t_2 = (2-theta)/(2+theta)*R_t_1 + gamma/(2+theta)*(I_t_1+I_t_2)
            D_t_2 = D_t_1 + delta*theta*(R_t_1+R_t_2)/2
            C_t_2 = C_t_1 + (1-delta)*theta*(R_t_1+R_t_2)/2
            
            param_record_approx_all[, i+1] = c(S_t_2$root, I_t_2, R_t_2, D_t_2, C_t_2)
            
          }
          
          Susceptive_record_matrix_all[sim_index, (n+1):(n+prediction_length)] = param_record_approx_all[1, (n+1):(n+prediction_length)]
          Infective_record_matrix_all[sim_index, (n+1):(n+prediction_length)] = param_record_approx_all[2, (n+1):(n+prediction_length)]
          Resolved_record_matrix_all[sim_index, (n+1):(n+prediction_length)] = param_record_approx_all[3, (n+1):(n+prediction_length)]
          Death_record_matrix_all[sim_index, (n+1):(n+prediction_length)] = param_record_approx_all[4, (n+1):(n+prediction_length)]
          ReCovered_record_matrix_all[sim_index, (n+1):(n+prediction_length)] = param_record_approx_all[5, (n+1):(n+prediction_length)]
        }

        ######## use GP to fix the death prediction ################
        
        # calculate median and CI for death
        death_record_sampled_beta_df = data.frame(date = seq.Date(start_date, end_date+prediction_length, by=1), mean = apply(Death_record_matrix_all, 2, mean, na.rm=T), upper = apply(Death_record_matrix_all, 2, quantile, 0.975, na.rm=T), 
                                                  lower = apply(Death_record_matrix_all, 2, quantile, 0.025, na.rm=T))

        death_CI_length = death_record_sampled_beta_df$upper - death_record_sampled_beta_df$lower
        
        # calculate sampled variance of death using Law of Large Numbers
        death_sample_variance = (1/(2*1.96) * death_CI_length)^2
        
        input_diff_death = as.matrix(seq(1,n ,1))
        testing_input_diff_death =  as.matrix(seq(1,n+prediction_length ,1))
        output_diff_death = as.matrix(death_selected - death_record_sampled_beta_df$mean[1:n])
        
        # fit GP on (Y-F) with pow_exp and alpha = 1.9
        sink(tempfile())     
        
        index_start=length(input_diff_death)-predicted_days_train+1
        index_end=length(input_diff_death)
        if(index_start<1){
          index_start=1
        }
        m_diff_death = rgasp(
          input_diff_death[index_start:index_end],
          (output_diff_death[index_start:index_end]),
          nugget.est = T,
          num_initial_values = 10,
          kernel_type = "pow_exp",
          zero.mean = "Yes"
        )
        if((1/m_diff_death@beta_hat)<50){
          m_diff_death = rgasp(
            input_diff_death,
            (output_diff_death),
            num_initial_values = 3,
            kernel_type = "pow_exp",
            zero.mean = "Yes",
            range.par=50,
            nugget=0.002
          )
        }
        if((1/m_diff_death@beta_hat)>500){
          m_diff_death = rgasp(
            input_diff_death,
            (output_diff_death),
            num_initial_values = 3,
            kernel_type = "pow_exp",
            zero.mean = "Yes",
            range.par=500,
            nugget=0.002
          )
        }
        sink()
       
        ## interval_data is a new feature implemented in 0.6.0
        m_diff_death_pred=predict(m_diff_death,testing_input_diff_death,interval_data=F,
                                  testing_trend=matrix(0,dim(testing_input_diff_death)[1]))

        # # calculate the adjusted mean ###death_record_sampled_beta_df$mean may be biased as we delete many zero beta, a mean/median one before truncation may be better
        death_record_sampled_beta_df$adjusted_mean = (death_record_sampled_beta_df$mean + m_diff_death_pred$mean)

        # # eliminate the decreasing trend in adjusted mean of death
        for (i in 2:(n+prediction_length)){
          if(death_record_sampled_beta_df$adjusted_mean[i]<death_record_sampled_beta_df$adjusted_mean[i-1]){
            death_record_sampled_beta_df$adjusted_mean[i] = death_record_sampled_beta_df$adjusted_mean[i-1]
          }
        }
        
        date_seq_pred = seq.Date(start_date, end_date+prediction_length, by=1)
        
        death_selected_all = county_death_raw[13:length(county_death_raw)]
        
        # this part of code is to choose and clean real death data with proper length
        if ((which(all_dates==start_date)+n+prediction_length)>length(death_selected_all)){
          death_selected_all_pred = as.numeric(death_selected_all[which(all_dates==start_date):length(death_selected_all)])
          for(i in 2:(length(death_selected_all_pred))){
            if(death_selected_all_pred[i] < death_selected_all_pred[i-1]){
              death_selected_all_pred[i] = death_selected_all_pred[i-1]
            }
          }
          
          if(paste0(state_name_short, county_name_selected) %in%   paste0(county_sudden_death$state, county_sudden_death$county)){
            daily_death_all = diff(death_selected_all_pred)
            sudden_increase_index_all = which(daily_death_all==max(daily_death_all))
            daily_death_all[(sudden_increase_index_all-13):(sudden_increase_index_all-1)] = daily_death_all[(sudden_increase_index_all-13):(sudden_increase_index_all-1)] + daily_death_all[sudden_increase_index_all]/14
            daily_death_all[sudden_increase_index_all] = daily_death_all[sudden_increase_index_all]/14
            death_selected_all_pred = c(death_selected_all_pred[1], cumsum(daily_death_all) + death_selected_all_pred[1])
          }
          
          death_selected_all_pred_aligned = c(death_selected_all_pred, rep(NA, (which(all_dates==start_date)+n+prediction_length) - length(death_selected_all) - 1))
          data_seq_for_real_death = seq.Date(start_date, all_dates[length(all_dates)], by=1)
        } else{
          death_selected_all_pred = as.numeric(death_selected_all[which(all_dates==start_date):(which(all_dates==start_date)+n+prediction_length-1)])
          for(i in 2:(length(death_selected_all_pred))){
            if(death_selected_all_pred[i] < death_selected_all_pred[i-1]){
              death_selected_all_pred[i] = death_selected_all_pred[i-1]
            }
          }
          
          if(paste0(state_name_short, county_name_selected) %in%   paste0(county_sudden_death$state, county_sudden_death$county)){
            daily_death_all = diff(death_selected_all_pred)
            sudden_increase_index_all = which(daily_death_all==max(daily_death_all))
            daily_death_all[(sudden_increase_index_all-13):(sudden_increase_index_all-1)] = daily_death_all[(sudden_increase_index_all-13):(sudden_increase_index_all-1)] + daily_death_all[sudden_increase_index_all]/14
            daily_death_all[sudden_increase_index_all] = daily_death_all[sudden_increase_index_all]/14
            death_selected_all_pred = c(death_selected_all_pred[1], cumsum(daily_death_all) + death_selected_all_pred[1])
          }
          
          death_selected_all_pred_aligned = death_selected_all_pred
          data_seq_for_real_death = seq.Date(start_date, end_date+prediction_length, by=1)
        }
        
        death_record_sampled_beta_df$real_death = death_selected_all_pred_aligned
        
        death_record_sampled_beta_df$prediction_indicator = c(rep(0, n), rep(1, prediction_length))
        
        ## Sample the death
        set.seed(each_index)
        
        m_diff_death_simulate=simulate(m_diff_death,testing_input_diff_death,num_sample=num_samples,
                                       testing_trend= length(testing_input_diff_death))

        sample_death_forecast=(m_diff_death_simulate+t(Death_record_matrix_all))[(n+1):(n+prediction_length),]

        sample_death_forecast_lower_upper=rowQuantiles(sample_death_forecast,probs=c(0.025,0.975))
        
        ## truncate
        sample_death_forecast_lower_upper[which(sample_death_forecast_lower_upper<death_selected_all_pred[n])]=death_selected_all_pred[n]
        
        ## monotone
        for(i in 1:2){
          for(j in 2:dim(sample_death_forecast_lower_upper)[1]){
            if (sample_death_forecast_lower_upper[j, i] < sample_death_forecast_lower_upper[j-1, i]){
              sample_death_forecast_lower_upper[j, i] = sample_death_forecast_lower_upper[j-1, i]
            }
          }
        }
        
        death_record_sampled_beta_df$adjusted_mean[1:n]=death_selected_all_pred[1:n]
        death_record_sampled_beta_df$lower_final=c(death_selected_all_pred[1:n],sample_death_forecast_lower_upper[,1])
        death_record_sampled_beta_df$upper_final=c(death_selected_all_pred[1:n],sample_death_forecast_lower_upper[,2])
        
        death_record_sampled_beta_df$adjusted_mean[(n+1):(n+prediction_length)][which(death_record_sampled_beta_df$adjusted_mean[(n+1):(n+prediction_length)]<death_selected_all_pred[n])]=death_selected_all_pred[n]
        
        ylimit_death = c(min(death_record_sampled_beta_df$lower_final,death_selected_all_pred ), max(death_record_sampled_beta_df$upper_final,death_selected_all_pred))
        
        date_seq_pred = seq.Date(start_date, end_date+prediction_length, by=1)
        date_seq = seq.Date(end_date + 1, end_date+prediction_length, by=1)
        
        # Visualize the death forecast and CIs
        
        plot(death_record_sampled_beta_df$adjusted_mean~date_seq_pred,ylim = ylimit_death, type="l", col="blue", xlab = "Date", ylab = "Death Cases",
        main = paste0(county_names[each_index], ", population=", round(N/10^6,2),"M", ", Ratio = ", round(ratio_real,3)))
        polygon(c(rev(date_seq), date_seq),
                c(rev(death_record_sampled_beta_df$upper_final[(n+1):(n+prediction_length)]), death_record_sampled_beta_df$lower_final[(n+1):(n+prediction_length)]),
                col = 'grey80',
                border = NA)
        lines(death_record_sampled_beta_df$adjusted_mean~date_seq_pred, col = "blue")
        lines(death_selected_all_pred~data_seq_for_real_death, col = "red")
        abline(v = end_date+1, lty=2)
        legend("topleft", legend = c("Observed Death", "Estimated Death"), lty = c(1,1), col = c("red", "blue"))

        #########################
        # Part 6: Record data for plotting figures
        #########################
        
        results_df_overall = data.frame(matrix(0, n+prediction_length, 25))
        colnames(results_df_overall) = c("date", "prediction_indicator","confirmed_real","death_real", "death_f_gp",  "death_lower", "death_upper", "Beta_estimated", "Beta_sampled","Beta_sampled_lower", "Beta_sampled_upper", "Beta_7_day_avg",
                                         "Beta_7_day_avg_lower", "Beta_7_day_avg_upper", "suspective_f", "infective_f", "resolved_f", "death_f", "recovered_f", "Susceptive_sampled", "Susceptive_lower", "Susceptive_upper", "Infective_sampled", "Infective_lower", "Infective_upper")
        
        results_df_overall$date = seq.Date(start_date, end_date+prediction_length, by=1)
        results_df_overall$prediction_indicator = c(rep(0, n), rep(1, prediction_length))
        results_df_overall$confirmed_real = c(confirm_selected[1:length(confirm_selected)], rep(NA, prediction_length))
        
        results_df_overall$death_f_gp = death_record_sampled_beta_df$adjusted_mean
        results_df_overall$death_real = death_record_sampled_beta_df$real_death
        results_df_overall$death_lower = death_record_sampled_beta_df$lower_final
        results_df_overall$death_upper = death_record_sampled_beta_df$upper_final
        
        results_df_overall$Beta_estimated = c(approx_beta_seq, rep(NA, prediction_length+1))
        results_df_overall$Beta_sampled = c(apply(beta_record_matrix, 2, mean, na.rm=T), rep(NA, prediction_length+1))
        results_df_overall$Beta_sampled_lower = c(apply(beta_record_matrix, 2, quantile, 0.025, na.rm=T), rep(NA, prediction_length+1))
        results_df_overall$Beta_sampled_upper = c(apply(beta_record_matrix, 2, quantile, 0.975, na.rm=T), rep(NA, prediction_length+1))
        
        results_df_overall$Beta_7_day_avg = c(apply(seven_day_averaged_beta_record_matrix, 2, mean, na.rm=T) , rep(NA, prediction_length+1))
        results_df_overall$Beta_7_day_avg_lower = c( apply(seven_day_averaged_beta_record_matrix, 2, quantile, 0.025, na.rm=T), rep(NA, prediction_length+1))
        results_df_overall$Beta_7_day_avg_upper = c( apply(seven_day_averaged_beta_record_matrix, 2, quantile, 0.975, na.rm=T), rep(NA, prediction_length+1))

        results_df_overall$suspective_f = c(param_record_approx_for_beta[1,], rep(NA, prediction_length))
        results_df_overall$infective_f = c(param_record_approx_for_beta[2,], rep(NA, prediction_length))
        results_df_overall$resolved_f = c(param_record_approx_for_beta[3,], rep(NA, prediction_length))
        results_df_overall$death_f = c(param_record_approx_for_beta[4,], rep(NA, prediction_length))
        results_df_overall$recovered_f = c(param_record_approx_for_beta[5,], rep(NA, prediction_length))
        
        results_df_overall$Susceptive_sampled =  c(apply(Susceptive_record_matrix[, 1:n], 2, mean, na.rm=T), rep(NA, prediction_length))
        results_df_overall$Susceptive_upper =  c(apply(Susceptive_record_matrix[, 1:n], 2, quantile, 0.975, na.rm=T), rep(NA, prediction_length))
        results_df_overall$Susceptive_lower =  c(apply(Susceptive_record_matrix[, 1:n], 2, quantile, 0.025, na.rm=T), rep(NA, prediction_length))
        
        results_df_overall$Infective_sampled =  c(apply(Infective_record_matrix[, 1:n], 2, mean, na.rm=T), rep(NA, prediction_length))
        results_df_overall$Infective_upper =  c(apply(Infective_record_matrix[, 1:n], 2, quantile, 0.975, na.rm=T), rep(NA, prediction_length))
        results_df_overall$Infective_lower =  c(apply(Infective_record_matrix[, 1:n], 2, quantile, 0.025, na.rm=T), rep(NA, prediction_length))
        
        results_overall_list[[each_index]] = results_df_overall
      }
    })
    
    # save some recorded data for further calculation of RMSE, CI length and C I coverage
    save(county_names, results_overall_list, file = paste0(file_path_cur, state_name_short, "_results_in_GP.RData"))
  }
  
  
})

stopImplicitCluster()


