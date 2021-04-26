
eliminate_abnormal_death = function(data){
  abnormal_index = which(diff(data)<0)
  if(length(abnormal_index)==0){
    return(data)
  } else{
    for (i in 1:length(abnormal_index)){
      data[abnormal_index[i]] = data[abnormal_index[i]+1]
    }
  }
  return(data)
}

data_seven_day_smoothing = function(data, default_method = T){

  if(!default_method){
    data_length = length(data)
    
    data_seven_day_average = rep(NA, data_length)
    
    data_seven_day_average[(1+3):(data_length-3)] = rollapply(data, width = 7, by = 1, FUN = mean, align = "right")
    
    data_seven_day_average[1+2] = mean(data[1:5], na.rm=T)
    data_seven_day_average[data_length-2] = mean(data[(data_length-4):(data_length)], na.rm=T)
    
    data_seven_day_average[1+1] = mean(data[1:3], na.rm=T)
    data_seven_day_average[data_length-1] = mean(data[(data_length-2):(data_length)], na.rm=T)
    
    # data_seven_day_average[1] = mean(data[1:4], na.rm = T)
    # data_seven_day_average[data_length] = mean(data[(data_length-3):data_length], na.rm=T)
    
    data_seven_day_average[1] = mean(data[1], na.rm=T)
    data_seven_day_average[data_length] = mean(data[data_length], na.rm=T)
  } else{
    data_length = length(data)
    
    data_seven_day_average = rep(NA, data_length)
    
    data_seven_day_average[(1+3):(data_length-3)] = rollapply(data, width = 7, by = 1, FUN = mean, align = "right")
    
    data_seven_day_average[1+2] = mean(data[1:6], na.rm=T)
    data_seven_day_average[data_length-2] = mean(data[(data_length-5):data_length], na.rm=T)
    
    data_seven_day_average[1+1] = mean(data[1:5], na.rm=T)
    data_seven_day_average[data_length-1] = mean(data[(data_length-4):data_length], na.rm=T)
    
    # data_seven_day_average[1] = mean(data[1:4], na.rm = T)
    # data_seven_day_average[data_length] = mean(data[(data_length-3):data_length], na.rm=T)
    
    data_seven_day_average[1] = mean(data[1:4], na.rm=T)
    data_seven_day_average[data_length] = mean(data[(data_length-3):data_length], na.rm=T)
  }
  
  return(data_seven_day_average)
}


data_seven_day_smoothing_updated = function(data){
  data_length = length(data)
  
  data_seven_day_average = rep(NA, data_length)
  
  data_seven_day_average[(1+3):(data_length-3)] = rollapply(data, width = 7, by = 1, FUN = mean, align = "right")
  
  data_seven_day_average[1+2] = mean(data[1:5], na.rm=T)
  data_seven_day_average[data_length-2] = mean(data[(data_length-4):(data_length)], na.rm=T)
  
  data_seven_day_average[1+1] = mean(data[1:3], na.rm=T)
  data_seven_day_average[data_length-1] = mean(data[(data_length-2):(data_length)], na.rm=T)
  
  # data_seven_day_average[1] = mean(data[1:4], na.rm = T)
  # data_seven_day_average[data_length] = mean(data[(data_length-3):data_length], na.rm=T)
  
  data_seven_day_average[1] = mean(data[1], na.rm=T)
  data_seven_day_average[data_length] = mean(data[data_length], na.rm=T)
  
  return(data_seven_day_average)
}


clean_us_death = function(data){
  
  
  # data cleaning
  names(data)[names(data) == "Long_"] = "Long"
  names(data)[names(data) == "FIPS"] = "county_fips"
  data$county_fips = as.character(data$county_fips)
  # there are some county fips missing, just redistribute or delete them
  missing_fips = which(is.na(data$county_fips))
  # data$county_fips[3148] = "25007"
  data = data[-missing_fips,]
  
  k = dim(data)[1]
  n = dim(data)[2]
  
  # delete counties with no population
  # no_population_idx = c()
  # for (i in 1:k){
  #   county_fips_iter = data$county_fips[i]
  #   county_fips_length = str_length(county_fips_iter)
  #   diff = 5-county_fips_length
  #   if (diff<5){
  #     data$county_fips[i] = paste0(do.call(paste0,as.list((rep(0,diff)))), county_fips_iter)
  #   }
  #   if(data$Population[i]==0){
  #     no_population_idx = c(no_population_idx, i)
  #   }
  # }
  # 
  # data = data[-no_population_idx,]
  
  # for(i in 1:k){
  #   death_data = data[i, 13:n ]
  #   death_data_modified = eliminate_abnormal_death(death_data)
  #   data[i, 13:n ] = death_data_modified
  # }
  
  return(data)
}

clean_us_death_for_map = function(data){
  
  
  # data cleaning
  names(data)[names(data) == "Long_"] = "Long"
  names(data)[names(data) == "FIPS"] = "county_fips"
  data$county_fips = as.character(data$county_fips)
  # there are some county fips missing, just redistribute or delete them
  missing_fips = which(is.na(data$county_fips))
  # data$county_fips[3148] = "25007"
  data = data[-missing_fips,]
  
  k = dim(data)[1]
  n = dim(data)[2]
  
  # delete counties with no population
  no_population_idx = c()
  for (i in 1:k){
    county_fips_iter = data$county_fips[i]
    
    if(is.na(county_fips_iter)){
      next
    }
    
    county_fips_length = str_length(county_fips_iter)
    diff = 5-county_fips_length
    if (diff<5){
      data$county_fips[i] = paste0(do.call(paste0,as.list((rep(0,diff)))), county_fips_iter)
    }
    if(data$Population[i]==0){
      no_population_idx = c(no_population_idx, i)
    }
  }
  
  data = data[-no_population_idx,]
  
  # for(i in 1:k){
  #   death_data = data[i, 13:n ]
  #   death_data_modified = eliminate_abnormal_death(death_data)
  #   data[i, 13:n ] = death_data_modified
  # }
  
  return(data)
}

clean_us_confirm_for_map = function(data){
  
  
  # data cleaning
  names(data)[names(data) == "Long_"] = "Long"
  names(data)[names(data) == "FIPS"] = "county_fips"
  data$county_fips = as.character(data$county_fips)
  # there are some county fips missing, just redistribute or delete them
  missing_fips = which(is.na(data$county_fips))
  # data$county_fips[3148] = "25007"
  data = data[-missing_fips,]
  
  k = dim(data)[1]
  n = dim(data)[2]
  
  # delete counties with no population
  # no_population_idx = c()
  for (i in 1:k){
    county_fips_iter = data$county_fips[i]
    
    if(is.na(county_fips_iter)){
      next
    }
    
    county_fips_length = str_length(county_fips_iter)
    diff = 5-county_fips_length
    if (diff<5){
      data$county_fips[i] = paste0(do.call(paste0,as.list((rep(0,diff)))), county_fips_iter)
    }
    # if(data$Population[i]==0){
    #   no_population_idx = c(no_population_idx, i)
    # }
  }
  
  # data = data[-no_population_idx,]
  
  # for(i in 1:k){
  #   death_data = data[i, 13:n ]
  #   death_data_modified = eliminate_abnormal_death(death_data)
  #   data[i, 13:n ] = death_data_modified
  # }
  
  return(data)
}

get_output_same_time_zone = function(data_type = "death_rate", state_name, state_name_short,start_date,training_length, duration = 90, criterion_death = 10, smoothness = TRUE){
  
  state_death = us_death %>%
    filter(Province_State == state_name, Admin2 != "Unassigned", Admin2 != paste0("Out of ",state_name), Admin2 != paste0("Out of ",state_name_short))
  
  state_death_rate = us_death_rate %>%
    filter(Province_State == state_name, Admin2 != "Unassigned", Admin2 != paste0("Out of ",state_name), Admin2 != paste0("Out of ",state_name_short))
  
  n_cols = dim(state_death_rate)[2]
  
  if(duration > (n_cols-12)){
    return("Error: the argument Duration is too large")
  }
  
  start_date_index = which(all_dates==start_date)
  
  end_training_date_index = start_date_index + training_length - 1
  
  end_date_index = start_date_index + duration - 1
  
  select_criterion = which(state_death[, (12 + end_training_date_index)] >=criterion_death)
  
  # select_criterion = which(state_death[, dim(state_death)[2]] >=criterion_death)
  
  state_death_rate_selected = state_death_rate[select_criterion,]
  county_name_all = state_death_rate_selected$Admin2
  
  if (data_type == "death_rate"){
    
    output = matrix(0,length(select_criterion), duration)
    
    for (i in 1:length(select_criterion)){
      output[i,] = as.numeric(state_death_rate_selected[ i , (12 + start_date_index):(12 + end_date_index)])
    }
  } else if(data_type == "death"){
    
    state_death_selected = state_death[select_criterion,]
    county_name_all = state_death_selected$Admin2
    
    output = matrix(0,length(select_criterion), duration)
    
    for (i in 1:length(select_criterion)){
      output[i,] = as.numeric(state_death_selected[ i ,  (12 + start_date_index):(12 + end_date_index)])
    }
  }
  
  output[is.na(output)] = 0
  output[is.infinite(output)] = 0
  
  k = dim(output)[1]
  for(i in 1:k){
    for(j in 1:(dim(output)[2]-1)){
      if(output[i, j+1]<output[i, j]){
        output[i, j+1] = output[i, j]
      }
    }
  }
  
  if (smoothness){
    output_smooth_all = matrix(0,k,duration)
    for(i in 1:k){
      output_each = output[i,]
      output_window = zoo(output_each)
      output_smoothed = c(output_each[1:2], rollapply(output_window, width = 5, by = 1, FUN = mean, align = "left"),output_each[(length(output_each)-1):length(output_each)])
      output_smoothed[length(output_smoothed)-1] = mean(c(output_smoothed[length(output_smoothed)-2], output_smoothed[length(output_smoothed)]))
      output_smoothed[2] = mean(c(output_smoothed[1],output_smoothed[3]))
      output_smooth_all[i,] = output_smoothed
    }
    
    
    results = list(1:2)
    results[[1]] = output_smooth_all
    results[[2]] = county_name_all
    
    return(results)
  } else {
    results = list(1:2)
    results[[1]] = output
    results[[2]] = county_name_all
    
    return(results)
  }
}



get_output = function(data_type = "death_rate", state_name = "California", duration = 90, criterion_death_rate = 20, criterion_death = 10){
  
  state_death = us_death %>%
    filter(Province_State == state_name)
  
  state_death_rate = us_death_rate %>%
    filter(Province_State == state_name)
  
  # select_criterion = which((state_death_rate[, dim(state_death_rate)[2]] >=criterion_death_rate) &
  #                            state_death[, dim(state_death)[2]] >=criterion_death)
  
  select_criterion = which(state_death[, dim(state_death)[2]] >=criterion_death)
  
  
  state_death_rate_selected = state_death_rate[select_criterion,]
  county_name_all = state_death_rate_selected$Admin2
  
  start_index_all = c()
  for(i in 1:dim(state_death_rate_selected)[1]){
    start_index_all = c(start_index_all, which(state_death_rate_selected[i,13:dim(state_death_rate_selected)[2]] > 0.3)[1])
  }
  start_index_all_selected = start_index_all[which((dim(state_death_rate_selected)[2]-12-start_index_all)>duration)]
  
  if (data_type == "death_rate"){
    
    # sum((dim(state_death_rate_selected)[2]-12-start_index_all)>duration)
    
    county_index_long_duration = which((dim(state_death_rate_selected)[2]-12-start_index_all)>duration)
    county_name_long_duration = county_name_all[county_index_long_duration]
    output = matrix(0,length(county_index_long_duration), duration)
    
    for (i in 1:length(county_index_long_duration)){
      output[i,] = as.numeric(state_death_rate_selected[county_index_long_duration[i], (0:(duration - 1) + 12 + start_index_all[county_index_long_duration[i]])])
    }
  } else if(data_type == "death"){
    
    state_death_selected = state_death[select_criterion,]
    county_name_all = state_death_selected$Admin2
    
    # sum((dim(state_death_selected)[2]-12-start_index_all)>duration)
    
    county_index_long_duration = which((dim(state_death_selected)[2]-12-start_index_all)>duration)
    county_name_long_duration = county_name_all[county_index_long_duration]
    output = matrix(0,length(county_index_long_duration), duration)
    
    for (i in 1:length(county_index_long_duration)){
      output[i,] = as.numeric(state_death_selected[county_index_long_duration[i], (0:(duration-1) + 12 + start_index_all[county_index_long_duration[i]])])
    }
  }
  k = dim(output)[1]
  while(sum(apply(output,1,diff)<0)>0){
    for (i in 1:k){
      if (sum(diff(output[i,])<0)>0){
        # print(paste0("County ", county_names_death[i],' ',i, " has outliers at "))
        # print(which(diff(output[i,])<0))
        output[i,] = eliminate_abnormal_death(output[i,])
      }
    }
  }
  
  output_smooth_all = matrix(0,k,duration)
  for(i in 1:k){
    output_each = output[i,]
    output_window = zoo(output_each)
    output_smoothed = c(output_each[1:2], rollapply(output_window, width = 5, by = 1, FUN = mean, align = "left"),output_each[(length(output_each)-1):length(output_each)])
    output_smoothed[length(output_smoothed)-1] = mean(c(output_smoothed[length(output_smoothed)-2], output_smoothed[length(output_smoothed)]))
    output_smoothed[2] = mean(c(output_smoothed[1],output_smoothed[3]))
    output_smooth_all[i,] = output_smoothed
  }
  
  
  results = list(1:3)
  results[[1]] = output_smooth_all
  results[[2]] = county_name_long_duration
  results[[3]] = start_index_all_selected
  
  return(results)
}

calculate_Beta_t = function(S_t,death_data,N){
  
  if (length(death_data) != 3){
    print("Error Message: Incorrect Death Data")
    return(NA)
  }
  gamma = 0.2
  theta = 0.1
  delta = 0.01
  
  d_t_plus_1 = death_data[1]
  delta_d_t_plus_2 = death_data[2] -death_data[1]
  delta_delta_d_plus_3 = death_data[3]+death_data[1] - 2*death_data[2]
  
  beta_t = N/S_t*(gamma + (1/theta * delta_delta_d_plus_3 + delta_d_t_plus_2)/(1/theta * delta_d_t_plus_2 + d_t_plus_1))
  return(beta_t)
}

calculate_S_t_plus_1 = function(S_t, beta_t, death_data,N){
  if (length(death_data) != 3){
    print("Error Message: Incorrect Death Data")
    return(NA)
  }
  gamma = 0.2
  theta = 0.1
  delta = 0.01
  
  d_t_plus_1 = death_data[1]
  delta_d_t_plus_2 = death_data[2] -death_data[1]
  # beta_t = calculate_Beta_t(S_t, death_data)
  
  S_t_plus_1 = S_t * (1 - beta_t * 1/(gamma * delta * N) * (1/theta * delta_d_t_plus_2 + d_t_plus_1))
  return(S_t_plus_1)
}

SIRDC <- function(time, state, parms) {
  gamma = parms[[1]]
  theta = parms[[2]]
  delta = parms[[3]]
  N = parms[[4]]
  betafun = parms[[5]]
  
  par <- as.list(c(state))
  with(par, {
    beta <- betafun(time)
    dS <- -beta/N * I * S
    dI <- beta/N * I * S - gamma * I
    dR <- gamma * I - theta * R
    dD <- delta * theta * R
    dC <- (1-delta) * theta * R
    list(c(dS, dI, dR, dD, dC))
  })
}

get_beta_t_same_time_zone = function(county_name_vector,state_name = "California", duration = 90){
  
  d = length(county_name_vector)
  
  beta_t_matrix = matrix(0, d, duration)
  
  for (county_idx in 1:d){
    
    county_name = county_name_vector[county_idx]
    # extract data from the selected place
    data_death = us_death %>% 
      dplyr::filter(Admin2== county_name & Province_State == state_name) %>% 
      dplyr::select(-c(UID:Combined_Key))
    
    N_each = data_death$Population # number of population
    
    death_selected = as.numeric(data_death[1,  which(all_dates==(start_date-8)):dim(data_death)[2]])
    
    ##################
    #data manipulation
    ##################
    # 1. filter the smoothed data
    death_selected_window = zoo(death_selected)
    death_smoothed = c(death_selected[1:2], rollapply(death_selected_window, width = 5, by = 1, FUN = mean, align = "left"),death_selected[(length(death_selected)-1):length(death_selected)])
    death_smoothed[length(death_smoothed)-1] = mean(c(death_smoothed[length(death_smoothed)-2], death_smoothed[length(death_smoothed)]))
    death_smoothed[2] = mean(c(death_smoothed[1],death_smoothed[3]))
    
    
    daily_death_smmothed = diff(death_smoothed)
    
    hpfilter_daily_death = hpfilter(daily_death_smmothed, freq=200)
    filtered_daily_death = hpfilter_daily_death$trend
    
    param_record = data.frame(S = rep(N_each, (length(filtered_daily_death)-3)), beta = rep(0, (length(filtered_daily_death)-3)))
    param_record$beta[1] = calculate_Beta_t(param_record$S[1], filtered_daily_death[2:4],N_each)
    
    # calculate Beta_t from death data
    for (day_index in 2:(length(filtered_daily_death)-3)){
      death_pre_window_1 = filtered_daily_death[day_index + 1:3]
      death_pre_window_2 = filtered_daily_death[day_index + 0:2]
      
      beta_t = param_record$beta[day_index - 1]
      S_t = param_record$S[day_index - 1]
      
      S_t_plus_1 = calculate_S_t_plus_1(S_t, beta_t, death_pre_window_2,N_each)
      beta_t_plus_1 = calculate_Beta_t(S_t_plus_1, death_pre_window_1,N_each)
      
      param_record[day_index, ] = c(S_t_plus_1, beta_t_plus_1)
    }
    
    beta_seq = param_record$beta[1:(dim(param_record)[1]-5)]
    
    if (length(beta_seq) < duration){
      return(print(paste0("The county ",county_name, " has no enough data.")))
    }
    
    beta_t_matrix[county_idx, ] = beta_seq[1:duration]
  }
  
  return(beta_t_matrix)
}

get_beta_t = function(county_name_vector,state_name = "California", duration = 90){
  
  d = length(county_name_vector)
  
  beta_t_matrix = matrix(0, d, duration)
  
  start_date_record = list(1:d)
  
  for (county_idx in 1:d){
    
    county_name = county_name_vector[county_idx]
    # extract data from the selected place
    data_death = us_death %>% 
      dplyr::filter(Admin2== county_name & Province_State == state_name) %>% 
      dplyr::select(-c(UID:Combined_Key))
    
    N = data_death$Population # number of population
    death_rate = as.numeric(data_death[1,2:length(data_death)]/N*10^6)
    
    start_date = all_dates[which(death_rate>0.3)[1]]
    real_start_date = start_date - 1
    
    start_date_record[[county_idx]] = start_date 
    
    death_selected = as.numeric(data_death[1, (1 + which(all_dates==real_start_date)):dim(data_death)[2]])
    
    ##################
    #data manipulation
    ##################
    # 1. filter the smoothed data
    death_selected_window = zoo(death_selected)
    death_smoothed = c(death_selected[1:2], rollapply(death_selected_window, width = 5, by = 1, FUN = mean, align = "left"),death_selected[(length(death_selected)-1):length(death_selected)])
    death_smoothed[length(death_smoothed)-1] = mean(c(death_smoothed[length(death_smoothed)-2], death_smoothed[length(death_smoothed)]))
    death_smoothed[2] = mean(c(death_smoothed[1],death_smoothed[3]))
    
    daily_death_smmothed = diff(death_smoothed)
    
    hpfilter_daily_death = hpfilter(daily_death_smmothed, freq=200)
    filtered_daily_death = hpfilter_daily_death$trend
    
    death_smoothed = death_smoothed[2:length(death_smoothed)]
    
    # 3. hold-out samples
    validation_size = 0
    Day = 1:length(death_smoothed)
    train_Day = 1:(length(death_smoothed)-validation_size)
    # vali_Day = (length(death_smoothed)-validation_size+1):length(death_smoothed)
    
    # param_record = data.frame(S = rep(N - confirm_selected[1]*10 - death_smoothed[1] - death_smoothed[1]/0.01, length(death_smoothed)-3), beta = rep(0, length(death_smoothed)-3))
    param_record = data.frame(S = rep(N, length(death_smoothed)-3), beta = rep(0, length(death_smoothed)-3))
    param_record$beta[1] = calculate_Beta_t(param_record$S[1], filtered_daily_death[2:4],N)
    
    # calculate Beta_t from death data
    for (day_index in 2:(length(death_smoothed)-3)){
      death_pre_window_1 = filtered_daily_death[day_index + 1:3]
      death_pre_window_2 = filtered_daily_death[day_index + 0:2]
      
      beta_t = param_record$beta[day_index - 1]
      S_t = param_record$S[day_index - 1]
      
      S_t_plus_1 = calculate_S_t_plus_1(S_t, beta_t, death_pre_window_2,N)
      beta_t_plus_1 = calculate_Beta_t(S_t_plus_1, death_pre_window_1,N)
      
      param_record[day_index, ] = c(S_t_plus_1, beta_t_plus_1)
    }
    
    beta_seq = param_record$beta[1:(dim(param_record)[1]-5)]
    
    if (length(beta_seq) < duration){
      return(print(paste0("The county ",county_name, " has no enough data.")))
    }
    
    beta_t_matrix[county_idx, ] = beta_seq[1:duration]
  }
  
  results = list(1:d)
  results[[1]] = beta_t_matrix
  results[[2]] = start_date_record
  
  return(results)
}


get_chol <- function(x, beta, kernel_type) {
  R0_00 = abs(outer(x, x, '-'))
  if(kernel_type == "matern_5_2"){
    R = matern_5_2_funct(R0_00, beta)
  } else if(kernel_type == "exp"){
    R = pow_exp_funct(R0_00, beta, alpha_i=1)
  }
  
  rcppeigen_get_chol(R)
}


calculate_Beta_t_directly = function(death_data, population_N, output_S_t = F){
  
  # calculate the beta_t directly
  if (length(death_data)<5){
    print("Error Message: Incorrect death data")
    return(NA)
  }
  
  death_selected_window = zoo(death_data)
  # death_smoothed =  rollapply(death_data, width = 7, by = 1, FUN = mean, align = "left")
  death_smoothed = c(death_data[1:2], rollapply(death_selected_window, width = 5, by = 1, FUN = mean, align = "left"),death_data[(length(death_data)-1):length(death_data)])
  #plot(death_smoothed)
  death_smoothed[length(death_smoothed)-1] = mean(c(death_smoothed[length(death_smoothed)-2], death_smoothed[length(death_smoothed)]))
  death_smoothed[2] = mean(c(death_smoothed[1],death_smoothed[3]))
  #plot(death_smoothed)
  
  daily_death_smmothed = diff(death_smoothed)
  
  hpfilter_daily_death = hpfilter(daily_death_smmothed, freq=200)
  filtered_daily_death = hpfilter_daily_death$trend
  
  filtered_daily_death[filtered_daily_death<0] = 0
  
  death_length = length(filtered_daily_death)
  
  t = length(filtered_daily_death) - 3
  
  nominator = rep(0, death_length-2)
  denominator = rep(0, death_length-2)
  
  for(i in 1:(death_length-2)){
    death_each = filtered_daily_death[i:(i+2)]
    
    d_t_plus_1 = death_each[1]
    delta_d_t_plus_2 = death_each[2] - death_each[1]
    delta_delta_d_plus_3 = death_each[3]+death_each[1] - 2*death_each[2]
    
    nominator[i] = 1/theta * delta_delta_d_plus_3 + delta_d_t_plus_2
    denominator[i] = 1/theta * delta_d_t_plus_2 + d_t_plus_1 + 10^(-5)
  }
  beta_0 = nominator[1] / denominator[1] + gamma
  
  if (length(filtered_daily_death) == 3){
    return(beta_0)
  }
  
  S_part_1 = 1/delta * denominator[1:t] + 1/(delta * gamma) * nominator[1:t]
  S_1_to_t = population_N - cumsum(S_part_1)
  
  
  
  
  beta_1_to_t = (nominator[2:(t+1)] / denominator[2:(t+1)] + gamma) * population_N / S_1_to_t
  
  # beta_0_to_t = c(beta_0, beta_1_to_t)
  
  if(output_S_t){
    results = list(1:2)
    results[[1]] = beta_1_to_t
    results[[2]] = S_1_to_t
    return(results)
  } else{
    return(beta_1_to_t)
  }
}

A_ini_loss_2d = function(param, A_ini, S_1, S_2){
  theta_rotation = param
  
  rotation_matrix = matrix(c(cos(theta_rotation), sin(theta_rotation), -sin(theta_rotation), cos(theta_rotation)),2,2)
  
  A_ini_R = A_ini %*% rotation_matrix
  
  return(-tr(t(A_ini_R) %*% S_1 %*% A_ini_R - 2 * t(A_ini_R) %*% S_2))
}

minmax_normalize <- function(x)
{
  return((x- min(x)) /(max(x)-min(x)))
}

five_day_smoothed_data = function(data){
  data_window = zoo(data)
  data_smoothed = c(data[1:2], rollapply(data_window, width = 5, by = 1, FUN = mean, align = "left"),data[(length(data)-1):length(data)])
  data_smoothed[length(data_smoothed)-1] = mean(c(data_smoothed[length(data_smoothed)-2], data_smoothed[length(data_smoothed)]))
  data_smoothed[2] = mean(c(data_smoothed[1],data_smoothed[3]))
  
  return(data_smoothed)
}

clean_and_five_day_smoothed_data = function(data){
  data_clean = data
  for(i in 1:(length(data)-1)){
    if (data_clean[i+1]<data_clean[i]){
      data_clean[i+1] = data_clean[i]
    }
  }
  
  data_window = zoo(data_clean)
  data_smoothed = c(data_clean[1:2], rollapply(data_window, width = 5, by = 1, FUN = mean, align = "left"),data_clean[(length(data_clean)-1):length(data_clean)])
  data_smoothed[length(data_smoothed)-1] = mean(c(data_smoothed[length(data_smoothed)-2], data_smoothed[length(data_smoothed)]))
  data_smoothed[2] = mean(c(data_smoothed[1],data_smoothed[3]))
  
  return(data_smoothed)
}

clean_data = function(data){
  data_clean = data
  for(i in 1:(length(data)-1)){
    if (data_clean[i+1]<data_clean[i]){
      data_clean[i+1] = data_clean[i]
    }
  }
  return(data_clean)
}

##########################
# functions for estimating beta_t in approximation method
##########################

find_root_S_t_2_given_I_t_2 = function(S_t_2, param, N, gamma){
  S_t_1 = param[1]
  I_t_1 = param[2]
  I_t_2 = param[3]
  
  beta_t_1_2 = (log(I_t_2/I_t_1) + gamma) * (2*N)/(S_t_1 + S_t_2)
  
  results = S_t_2/S_t_1 - exp((beta_t_1_2/N* -(I_t_1+I_t_2)/2))
  return(results)
}

find_root_S_t_2 = function(S_t_2, param, N, gamma){
  S_t_1 = param[1]
  beta_t_1 = param[2]
  I_t_1 = param[3]
  
  I_t_2 = I_t_1 * exp(beta_t_1 * (S_t_1 + S_t_2) / (2*N) - gamma)
  
  results = S_t_2/S_t_1 - exp((beta_t_1/N* -(I_t_1+I_t_2)/2))
  return(results)
}

find_root_I_t_2_second_version = function(I_t_2, param,N, gamma){
  S_t_1 = param[1]
  S_t_2 = param[2]
  beta_t_1 = param[3]
  I_t_1 = param[4]
  
  
  # S_t_2 = S_t_1 * exp((beta_t_1/N * - (I_t_1+I_t_2)/2))
  
  results = I_t_2/I_t_1 - exp(beta_t_1 * (S_t_1 + S_t_2) / (2*N) - gamma)
  return(results)
}

find_root_beta = function(beta_t_1, param,N, gamma){
  S_t_1 = param[1]
  S_t_2 = param[2]
  I_t_1 = param[3]
  
  I_t_2 = I_t_1 * exp(beta_t_1 * (S_t_1 + S_t_2) / (2*N) - gamma)
  
  results = S_t_2/S_t_1 - exp((beta_t_1/N* -(I_t_1+I_t_2)/2))
  return(results)
}

find_root_I_t_2 = function(I_t_2, param,N, gamma){
  S_t_1 = param[1]
  beta_t_1 = param[2]
  I_t_1 = param[3]
  
  S_t_2 = S_t_1 * exp((beta_t_1/N * - (I_t_1+I_t_2)/2))
  
  results = I_t_2/I_t_1 - exp(beta_t_1 * (S_t_1 + S_t_2) / (2*N) - gamma)
  return(results)
}

loss_approx_beta = function(param, death_cases, confirmed_cases,unadjusted_confirm_selected_smoothed, N_population, trianing_length, fixed_global_params, penalty = F,
                            fitted_days=length(death_cases),weight_loss=F){
  # ratio = exp(param[1])
  I_0 = param[1]
  R_0 = param[2]
  
  gamma = fixed_global_params[1]
  theta = fixed_global_params[2]
  delta = fixed_global_params[3]
  
  ratio = confirmed_cases[1]/(I_0+ R_0+ death_cases[1])
  
  estimated_confirm = confirmed_cases/ratio
  
  for(i in 1:length(estimated_confirm)){
    if (estimated_confirm[i] < unadjusted_confirm_selected_smoothed[i]){
      estimated_confirm[i] = unadjusted_confirm_selected_smoothed[i]
    }
  }
  
  S_t_seq = N_population - estimated_confirm
  
  # S_t_seq = N_population - confirmed_cases/ratio
  
  init_for_beta = c(S_t_seq[1], I_0, R_0, death_cases[1], 0)
  
  param_record_approx_for_beta = matrix(0, 5, trianing_length)
  param_record_approx_for_beta[,1] = init_for_beta
  param_record_approx_for_beta[1,] = S_t_seq
  
  approx_beta_seq = rep(0, trianing_length-1)
  
  for (i in 1:(trianing_length-1)){
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
    
    beta_t_1_2 = uniroot(find_root_beta, c(0, 10^6), tol = 0.0001, param = c(S_t_1, S_t_2, I_t_1), N = N_population, gamma=gamma)
    # I_t_2 = uniroot(find_root_I_t_2, c(0, N_population), tol = 0.0001, param = c(S_t_1, beta_t_1_2$root, I_t_1), N = N_population, gamma=gamma)
    
    I_t_2 = I_t_1 * exp(beta_t_1_2$root*(S_t_1 + S_t_2)/(2*N_population) - gamma)
    R_t_2 = (2-theta)/(2+theta)*R_t_1 + gamma/(2+theta)*(I_t_1+I_t_2)
    D_t_2 = D_t_1 + delta*theta*(R_t_1+R_t_2)/2
    C_t_2 = C_t_1 + (1-delta)*theta*(R_t_1+R_t_2)/2
    
    param_record_approx_for_beta[2:5, i+1] = c(I_t_2, R_t_2, D_t_2, C_t_2)
    approx_beta_seq[i] = beta_t_1_2$root
  }
  
  if(penalty){
    negative_index = which(param_record_approx_for_beta[4,] <= death_cases)
    positive_index=  which(param_record_approx_for_beta[4,] > death_cases)
    
    loss = sqrt( 10 * sum(((param_record_approx_for_beta[4,] - death_cases)[negative_index])^2) + 1 * sum(((param_record_approx_for_beta[4,] - death_cases)[positive_index])^2))
    return(loss)
    
  } else{
    #return(sqrt(sum((param_record_approx_for_beta[4,] - death_cases)^2)))
    index_selected_fitted=(length(death_cases)-fitted_days+1):length(death_cases)
    if(fitted_days>length(death_cases)){
      index_selected_fitted=1:length(death_cases)
    }
    
    if(!weight_loss){
      return(sqrt(sum((param_record_approx_for_beta[4, index_selected_fitted] - death_cases[index_selected_fitted])^2)))
    }else{
      weights=seq(length(index_selected_fitted),1,-1)
      return(sqrt(sum( ((param_record_approx_for_beta[4, index_selected_fitted] - death_cases[index_selected_fitted])^2)/weights^2)))
    }
  }
  
}

loss_approx_beta_log_input = function(param, death_cases, confirmed_cases,unadjusted_confirm_selected_smoothed, N_population, trianing_length, fixed_global_params, penalty = F,
                                      fitted_days=length(death_cases),weight_loss=F){
  # ratio = exp(param[1])
  I_0_plus_R_0 = exp(param[1])
  R_0 = exp(param[2])
  
  I_0 = I_0_plus_R_0 - R_0
  
  gamma = fixed_global_params[1]
  theta = fixed_global_params[2]
  delta = fixed_global_params[3]
  
  ratio = confirmed_cases[1]/(I_0+ R_0+ death_cases[1])
  
  estimated_confirm = confirmed_cases/ratio
  
  for(i in 1:length(estimated_confirm)){
    if (estimated_confirm[i] < unadjusted_confirm_selected_smoothed[i]){
      estimated_confirm[i] = unadjusted_confirm_selected_smoothed[i]
    }
  }
  
  S_t_seq = N_population - estimated_confirm
  
  # S_t_seq = N_population - confirmed_cases/ratio
  
  init_for_beta = c(S_t_seq[1], I_0, R_0, death_cases[1], 0)
  
  param_record_approx_for_beta = matrix(0, 5, trianing_length)
  param_record_approx_for_beta[,1] = init_for_beta
  param_record_approx_for_beta[1,] = S_t_seq
  
  approx_beta_seq = rep(0, trianing_length-1)
  
  for (i in 1:(trianing_length-1)){
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
    
    beta_t_1_2 = uniroot(find_root_beta, c(0, 10^6), tol = 0.0001, param = c(S_t_1, S_t_2, I_t_1), N = N_population, gamma=gamma)
    # I_t_2 = uniroot(find_root_I_t_2, c(0, N_population), tol = 0.0001, param = c(S_t_1, beta_t_1_2$root, I_t_1), N = N_population, gamma=gamma)
    
    I_t_2 = I_t_1 * exp(beta_t_1_2$root*(S_t_1 + S_t_2)/(2*N) - gamma)
    R_t_2 = (2-theta)/(2+theta)*R_t_1 + gamma/(2+theta)*(I_t_1+I_t_2)
    D_t_2 = D_t_1 + delta*theta*(R_t_1+R_t_2)/2
    C_t_2 = C_t_1 + (1-delta)*theta*(R_t_1+R_t_2)/2
    
    param_record_approx_for_beta[2:5, i+1] = c(I_t_2, R_t_2, D_t_2, C_t_2)
    approx_beta_seq[i] = beta_t_1_2$root
  }
  
  if(penalty){
    negative_index = which(param_record_approx_for_beta[4,] <= death_cases)
    positive_index=  which(param_record_approx_for_beta[4,] > death_cases)
    
    loss = sqrt( 10 * sum(((param_record_approx_for_beta[4,] - death_cases)[negative_index])^2) + 1 * sum(((param_record_approx_for_beta[4,] - death_cases)[positive_index])^2))
    return(loss)
    
  } else{
    #return(sqrt(sum((param_record_approx_for_beta[4,] - death_cases)^2)))
    index_selected_fitted=(length(death_cases)-fitted_days+1):length(death_cases)
    if(fitted_days>length(death_cases)){
      index_selected_fitted=1:length(death_cases)
    }
    
    if(!weight_loss){
      return(sqrt(sum((param_record_approx_for_beta[4, index_selected_fitted] - death_cases[index_selected_fitted])^2)))
    }else{
      weights=seq(length(index_selected_fitted),1,-1)
      return(sqrt(sum( ((param_record_approx_for_beta[4, index_selected_fitted] - death_cases[index_selected_fitted])^2)/weights^2)))
    }
  }
  
}



loss_approx_beta_three_param = function(param,
                                        death_cases,
                                        confirmed_cases,
                                        real_unadjusted_smoothed,
                                        N_population,
                                        training_length,
                                        daily_positive_rate_selected,
                                        fitted_days,
                                        fixed_global_params,
                                        weight_loss = F,
                                        penalty = F) {
  # ratio = exp(param[1])
  I_0 = exp(param[1])
  R_0 = exp(param[2])
  exp_index=  exp(param[3])
  
  gamma = fixed_global_params[1]
  theta = fixed_global_params[2]
  delta = fixed_global_params[3]
  
  c_t_seq = daily_positive_rate_selected^(exp_index)
  
  daily_confirm_selected = diff(confirmed_cases)
  
  confirm_selected_adjusted = c( c_t_seq[1] * confirmed_cases[1] , cumsum(daily_confirm_selected*c_t_seq[2:training_length]  )+c_t_seq[1]*confirmed_cases[1])
  
  # confirm_selected_adjusted = cumsum(daily_confirm_selected*c_t_seq  )+confirmed_cases[1]
  
  confirmed_cases_smoothed = data_seven_day_smoothing(confirm_selected_adjusted)
  
  ratio = confirmed_cases_smoothed[1]/(I_0 + R_0 + death_cases[1])
  
  estimated_confirm = confirmed_cases_smoothed/ratio
  
  # return(estimated_confirm)
  for(est_confirm_index in 1:length(estimated_confirm)){
    if (estimated_confirm[est_confirm_index] < real_unadjusted_smoothed[est_confirm_index]){
      estimated_confirm[est_confirm_index] = real_unadjusted_smoothed[est_confirm_index]
    }
  }
  
  S_t_seq = N_population - estimated_confirm
  
  init_for_beta = c(S_t_seq[1], I_0, R_0, death_cases[1], 0)
  
  param_record_approx_for_beta = matrix(0, 5, training_length)
  param_record_approx_for_beta[,1] = init_for_beta
  param_record_approx_for_beta[1,] = S_t_seq
  
  approx_beta_seq = rep(0, training_length-1)
  
  for (i in 1:(training_length-1)){
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
    
    beta_t_1_2 = uniroot(find_root_beta, c(0, 10^6), tol = 0.0001, param = c(S_t_1, S_t_2, I_t_1), N = N_population, gamma=gamma)
    # I_t_2 = uniroot(find_root_I_t_2, c(0, N_population), tol = 0.0001, param = c(S_t_1, beta_t_1_2$root, I_t_1), N = N_population, gamma=gamma)
    
    I_t_2 = I_t_1 * exp(beta_t_1_2$root*(S_t_1 + S_t_2)/(2*N_population) - gamma)
    R_t_2 = (2-theta)/(2+theta)*R_t_1 + gamma/(2+theta)*(I_t_1+I_t_2)
    D_t_2 = D_t_1 + delta*theta*(R_t_1+R_t_2)/2
    C_t_2 = C_t_1 + (1-delta)*theta*(R_t_1+R_t_2)/2
    
    param_record_approx_for_beta[2:5, i+1] = c(I_t_2, R_t_2, D_t_2, C_t_2)
    approx_beta_seq[i] = beta_t_1_2$root
  }
  
  if(penalty){
    negative_index = which(param_record_approx_for_beta[4,] <= death_cases)
    positive_index=  which(param_record_approx_for_beta[4,] > death_cases)
    
    loss = sqrt( 10 * sum(((param_record_approx_for_beta[4,] - death_cases)[negative_index])^2) + 1 * sum(((param_record_approx_for_beta[4,] - death_cases)[positive_index])^2))
    return(loss)
    
  } else{
    #return(sqrt(sum((param_record_approx_for_beta[4,] - death_cases)^2)))
    index_selected_fitted=(length(death_cases)-fitted_days+1):length(death_cases)
    if(fitted_days>length(death_cases)){
      index_selected_fitted=1:length(death_cases)
    }
    
    if(!weight_loss){
      return(sqrt(sum((param_record_approx_for_beta[4, index_selected_fitted] - death_cases[index_selected_fitted])^2)))
    }else{
      weights=seq(length(index_selected_fitted),1,-1)
      return(sqrt(sum( ((param_record_approx_for_beta[4, index_selected_fitted] - death_cases[index_selected_fitted])^2)/weights^2)))
    }
  }
  
  
}

loss_approx_beta_three_param_for_constrOptim = function(param,
                                                        death_cases,
                                                        confirmed_cases,
                                                        real_unadjusted_smoothed,
                                                        N_population,
                                                        training_length,
                                                        daily_positive_rate_selected,
                                                        fitted_days,
                                                        fixed_global_params,
                                                        weight_loss = F,
                                                        penalty = F) {
  # ratio = exp(param[1])
  I_0_plus_R_0 = exp(param[1])
  R_0 = exp(param[2])
  exp_index=  exp(param[3])
  
  I_0 = I_0_plus_R_0 - R_0
  
  gamma = fixed_global_params[1]
  theta = fixed_global_params[2]
  delta = fixed_global_params[3]
  
  c_t_seq = daily_positive_rate_selected^(exp_index)
  
  daily_confirm_selected = diff(confirmed_cases)
  
  confirm_selected_adjusted = c( c_t_seq[1] * confirmed_cases[1] , cumsum(daily_confirm_selected*c_t_seq[2:training_length]  )+c_t_seq[1]*confirmed_cases[1])
  
  # confirm_selected_adjusted = cumsum(daily_confirm_selected*c_t_seq  )+confirmed_cases[1]
  
  confirmed_cases_smoothed = data_seven_day_smoothing(confirm_selected_adjusted)
  
  ratio = confirmed_cases_smoothed[1]/(I_0 + R_0 + death_cases[1])
  
  estimated_confirm = confirmed_cases_smoothed/ratio
  
  # return(estimated_confirm)
  for(est_confirm_index in 1:length(estimated_confirm)){
    if (estimated_confirm[est_confirm_index] < real_unadjusted_smoothed[est_confirm_index]){
      estimated_confirm[est_confirm_index] = real_unadjusted_smoothed[est_confirm_index]
    }
  }
  
  S_t_seq = N_population - estimated_confirm
  
  init_for_beta = c(S_t_seq[1], I_0, R_0, death_cases[1], 0)
  
  param_record_approx_for_beta = matrix(0, 5, training_length)
  param_record_approx_for_beta[,1] = init_for_beta
  param_record_approx_for_beta[1,] = S_t_seq
  
  approx_beta_seq = rep(0, training_length-1)
  
  for (i in 1:(training_length-1)){
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
    
    beta_t_1_2 = uniroot(find_root_beta, c(0, 10^6), tol = 0.0001, param = c(S_t_1, S_t_2, I_t_1), N = N_population, gamma=gamma)
    # I_t_2 = uniroot(find_root_I_t_2, c(0, N_population), tol = 0.0001, param = c(S_t_1, beta_t_1_2$root, I_t_1), N = N_population, gamma=gamma)
    
    I_t_2 = I_t_1 * exp(beta_t_1_2$root*(S_t_1 + S_t_2)/(2*N_population) - gamma)
    R_t_2 = (2-theta)/(2+theta)*R_t_1 + gamma/(2+theta)*(I_t_1+I_t_2)
    D_t_2 = D_t_1 + delta*theta*(R_t_1+R_t_2)/2
    C_t_2 = C_t_1 + (1-delta)*theta*(R_t_1+R_t_2)/2
    
    param_record_approx_for_beta[2:5, i+1] = c(I_t_2, R_t_2, D_t_2, C_t_2)
    approx_beta_seq[i] = beta_t_1_2$root
  }
  
  if(penalty){
    negative_index = which(param_record_approx_for_beta[4,] <= death_cases)
    positive_index=  which(param_record_approx_for_beta[4,] > death_cases)
    
    loss = sqrt( 10 * sum(((param_record_approx_for_beta[4,] - death_cases)[negative_index])^2) + 1 * sum(((param_record_approx_for_beta[4,] - death_cases)[positive_index])^2))
    return(loss)
    
  } else{
    #return(sqrt(sum((param_record_approx_for_beta[4,] - death_cases)^2)))
    index_selected_fitted=(length(death_cases)-fitted_days+1):length(death_cases)
    if(fitted_days>length(death_cases)){
      index_selected_fitted=1:length(death_cases)
    }
    
    if(!weight_loss){
      return(sqrt(sum((param_record_approx_for_beta[4, index_selected_fitted] - death_cases[index_selected_fitted])^2)))
    }else{
      weights=seq(length(index_selected_fitted),1,-1)
      return(sqrt(sum( ((param_record_approx_for_beta[4, index_selected_fitted] - death_cases[index_selected_fitted])^2)/weights^2)))
    }
  }
  
  
}


loss_approx_beta_three_param_for_constrOptim_daily_avg = function(param,
                                                        death_cases,
                                                        confirmed_cases,
                                                        real_unadjusted_smoothed,
                                                        N_population,
                                                        training_length,
                                                        daily_positive_rate_selected,
                                                        fitted_days,
                                                        fixed_global_params,
                                                        weight_loss = F,
                                                        penalty = F) {
  # ratio = exp(param[1])
  I_0_plus_R_0 = exp(param[1])
  R_0 = exp(param[2])
  exp_index=  exp(param[3])
  
  I_0 = I_0_plus_R_0 - R_0
  
  gamma = fixed_global_params[1]
  theta = fixed_global_params[2]
  delta = fixed_global_params[3]
  
  c_t_seq = daily_positive_rate_selected^(exp_index)
  
  daily_confirm_selected = diff(confirmed_cases)
  
  daily_confirm_selected_smoothed = data_seven_day_smoothing(daily_confirm_selected*c_t_seq[2:training_length])
  
  confirmed_cases_smoothed = c( c_t_seq[1] * confirmed_cases[1] , cumsum(daily_confirm_selected_smoothed )+c_t_seq[1]*confirmed_cases[1])
  
  # confirm_selected_adjusted = cumsum(daily_confirm_selected*c_t_seq  )+confirmed_cases[1]
  
  # confirmed_cases_smoothed = data_seven_day_smoothing(confirm_selected_adjusted)
  
  ratio = confirmed_cases_smoothed[1]/(I_0 + R_0 + death_cases[1])
  
  estimated_confirm = confirmed_cases_smoothed/ratio
  
  # return(estimated_confirm)
  for(est_confirm_index in 1:length(estimated_confirm)){
    if (estimated_confirm[est_confirm_index] < real_unadjusted_smoothed[est_confirm_index]){
      estimated_confirm[est_confirm_index] = real_unadjusted_smoothed[est_confirm_index]
    }
  }
  
  S_t_seq = N_population - estimated_confirm
  
  init_for_beta = c(S_t_seq[1], I_0, R_0, death_cases[1], 0)
  
  param_record_approx_for_beta = matrix(0, 5, training_length)
  param_record_approx_for_beta[,1] = init_for_beta
  param_record_approx_for_beta[1,] = S_t_seq
  
  approx_beta_seq = rep(0, training_length-1)
  
  for (i in 1:(training_length-1)){
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
    
    beta_t_1_2 = uniroot(find_root_beta, c(0, 10^6), tol = 0.0001, param = c(S_t_1, S_t_2, I_t_1), N = N_population, gamma=gamma)
    # I_t_2 = uniroot(find_root_I_t_2, c(0, N_population), tol = 0.0001, param = c(S_t_1, beta_t_1_2$root, I_t_1), N = N_population, gamma=gamma)
    
    I_t_2 = I_t_1 * exp(beta_t_1_2$root*(S_t_1 + S_t_2)/(2*N_population) - gamma)
    R_t_2 = (2-theta)/(2+theta)*R_t_1 + gamma/(2+theta)*(I_t_1+I_t_2)
    D_t_2 = D_t_1 + delta*theta*(R_t_1+R_t_2)/2
    C_t_2 = C_t_1 + (1-delta)*theta*(R_t_1+R_t_2)/2
    
    param_record_approx_for_beta[2:5, i+1] = c(I_t_2, R_t_2, D_t_2, C_t_2)
    approx_beta_seq[i] = beta_t_1_2$root
  }
  
  if(penalty){
    negative_index = which(param_record_approx_for_beta[4,] <= death_cases)
    positive_index=  which(param_record_approx_for_beta[4,] > death_cases)
    
    loss = sqrt( 10 * sum(((param_record_approx_for_beta[4,] - death_cases)[negative_index])^2) + 1 * sum(((param_record_approx_for_beta[4,] - death_cases)[positive_index])^2))
    return(loss)
    
  } else{
    #return(sqrt(sum((param_record_approx_for_beta[4,] - death_cases)^2)))
    index_selected_fitted=(length(death_cases)-fitted_days+1):length(death_cases)
    if(fitted_days>length(death_cases)){
      index_selected_fitted=1:length(death_cases)
    }
    
    if(!weight_loss){
      return(sqrt(sum((param_record_approx_for_beta[4, index_selected_fitted] - death_cases[index_selected_fitted])^2)))
    }else{
      weights=seq(length(index_selected_fitted),1,-1)
      return(sqrt(sum( ((param_record_approx_for_beta[4, index_selected_fitted] - death_cases[index_selected_fitted])^2)/weights^2)))
    }
  }
  
  
}





loss_approx_beta_CI = function(param, death_cases, confirmed_cases, N_population){
  # ratio = exp(param[1])
  I_0 = param[1]
  R_0 = param[2]
  
  ratio = confirmed_cases[1]/(I_0+ R_0+ death_cases[1]+ 99*death_cases[1])
  
  S_t_seq = N_population - confirmed_cases/ratio
  
  init_for_beta = c(S_t_seq[1], I_0, R_0, death_cases[1], 99*death_cases[1])
  
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
    
    beta_t_1_2 = uniroot(find_root_beta, c(0, 10^6), tol = 0.0001, param = c(S_t_1, S_t_2, I_t_1))
    I_t_2 = uniroot(find_root_I_t_2, c(0, N_population), tol = 0.0001, param = c(S_t_1, beta_t_1_2$root, I_t_1))
    R_t_2 = (2-theta)/(2+theta)*R_t_1 + gamma/(2+theta)*(I_t_1+I_t_2$root)
    D_t_2 = D_t_1 + delta*theta*(R_t_1+R_t_2)/2
    C_t_2 = C_t_1 + (1-delta)*theta*(R_t_1+R_t_2)/2
    
    param_record_approx_for_beta[2:5, i+1] = c(I_t_2$root, R_t_2, D_t_2, C_t_2)
    approx_beta_seq[i] = beta_t_1_2$root
  }
  
  return(sqrt(sum((param_record_approx_for_beta[4,] - death_cases)^2)))
}

calculate_beta_t_approx = function(death_cases, confirmed_cases,N_population, figure = F){
  
  ui = matrix(c(1,0,-1,0,1,-1), 3,2)
  ci = matrix(c(0,0, -((confirmed_cases[1]*N_population)/confirmed_cases[n] - 100 * death_cases[1])  ))
  
  m_approx = constrOptim(
    c(100, 100),
    loss_approx_beta_CI,
    NULL,
    ui = ui,
    ci = ci,
    death_cases = death_cases,
    confirmed_cases = confirmed_cases,
    N_population = N_population
  )
  
  I_0 = m_approx$par[1]
  R_0 = m_approx$par[2]
  
  ratio = confirmed_cases[1]/(I_0+ R_0+ death_cases[1] + 99*death_cases[1])
  
  # ratio_real = confirm_selected[2]/(I_0+ R_0+ death_cases[1]+ 99*death_cases[1])
  
  S_t_seq = N_population - confirmed_cases/ratio
  
  # plot(confirmed_cases/ratio)
  # lines(confirm_selected)
  
  init_for_beta = c(S_t_seq[1], I_0, R_0, death_cases[1], 99*death_cases[1])
  
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
    
    beta_t_1_2 = uniroot(find_root_beta, c(0, 10^6), tol = 0.0001, param = c(S_t_1, S_t_2, I_t_1))
    I_t_2 = uniroot(find_root_I_t_2, c(0, N_population), tol = 0.0001, param = c(S_t_1, beta_t_1_2$root, I_t_1))
    R_t_2 = (2-theta)/(2+theta)*R_t_1 + gamma/(2+theta)*(I_t_1+I_t_2$root)
    D_t_2 = D_t_1 + delta*theta*(R_t_1+R_t_2)/2
    C_t_2 = C_t_1 + (1-delta)*theta*(R_t_1+R_t_2)/2
    
    param_record_approx_for_beta[2:5, i+1] = c(I_t_2$root, R_t_2, D_t_2, C_t_2)
    approx_beta_seq[i] = beta_t_1_2$root
  }
  
  if(figure){
    plot(param_record_approx_for_beta[4,], type = "l", col = "blue")
    lines(death_cases, col = "red")
    legend("topleft", legend = c("Estimated Death", "Real Death"), lty=c(1,1), col = c("blue", "red"))
  }
  
  
  
  results = list(1:3)
  results[[1]] = m_approx
  results[[2]] = approx_beta_seq
  results[[3]] = param_record_approx_for_beta
  return(results)
}

# parallel computing
multiResultClass <- function(beta=NULL,susceptive=NULL, infective=NULL, death = NULL)
{
  me <- list(
    beta = beta,
    susceptive = susceptive,
    infective = infective,
    death = death
  )
  
  ## Set the name for the class
  class(me) <- append(class(me),"multiResultClass")
  return(me)
}

state_overall_result_Class <- function(state_name = NULL, results_overll_list=NULL,county_names=NULL)
{
  me <- list(
    state_name = state_name,
    results_overll_list = results_overll_list,
    county_names = county_names
  )
  
  ## Set the name for the class
  class(me) <- append(class(me),"state_overall_result_Class")
  return(me)
}

pow_param_result_Class <- function(state_name = NULL, power_param=NULL,rmse_record=NULL)
{
  me <- list(
    state_name = state_name,
    power_param = power_param,
    rmse_record = rmse_record
  )
  
  ## Set the name for the class
  class(me) <- append(class(me),"pow_param_result_Class")
  return(me)
}

proportion_death_reduction_Class <- function(state_name = NULL,county_name = NULL,real = NULL, mean=NULL,upper=NULL, lower = NULL)
{
  me <- list(
    state_name = state_name,
    county_name = county_name,
    real = real,
    mean = mean,
    upper = upper,
    lower = lower
  )
  
  ## Set the name for the class
  class(me) <- append(class(me),"proportion_death_reduction_Class")
  return(me)
}

proportion_death_reduction_national_Class <- function(state_name = NULL,county_name = NULL, death_samples_5=NULL,death_samples_4_75=NULL, death_samples_4_5 = NULL, infectious_samples_5=NULL,infectious_samples_4_75=NULL,infectious_samples_4_5 = NULL)
{
  me <- list(
    state_name = state_name,
    county_name = county_name,
    death_samples_5 = death_samples_5,
    death_samples_4_75 = death_samples_4_75,
    death_samples_4_5 = death_samples_4_5,
    infectious_samples_5 = infectious_samples_5,
    infectious_samples_4_75 = infectious_samples_4_75,
    infectious_samples_4_5 = infectious_samples_4_5
  )
  
  ## Set the name for the class
  class(me) <- append(class(me),"proportion_death_reduction_national_Class")
  return(me)
}

bootstrap_matrix = function(data, num_bootstrap){
  boot_strap_order = sample(1:dim(data)[1], num_bootstrap, replace = TRUE)
  results_matrix = data[boot_strap_order,]
  return(results_matrix)
}


death_pred_GP_with_without_linear_Class <-
  function(state_name = NULL,
           county_name = NULL,
           death_pred_7_linear = NULL,
           death_pred_21_linear = NULL,
           death_pred_7_no_linear = NULL,
           death_pred_21_no_linear = NULL,
           death_real_7 = NULL,
           death_real_21 = NULL) {
    me <- list(
      state_name = state_name,
      county_name = county_name,
      death_pred_7_linear = death_pred_7_linear,
      death_pred_21_linear = death_pred_21_linear,
      death_pred_7_no_linear = death_pred_7_no_linear,
      death_pred_21_no_linear = death_pred_21_no_linear,
      death_real_7 = death_real_7,
      death_real_21 = death_real_21
    )
    
    ## Set the name for the class
    class(me) <-
      append(class(me), "death_pred_GP_with_without_linear_Class")
    return(me)
  }

