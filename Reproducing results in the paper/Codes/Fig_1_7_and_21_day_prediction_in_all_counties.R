##########################
# Code for figure 1 in the main manuscript
# 
# Author: Hanmo Li & Mengyang Gu
#
# Email: mengyang at pstat.ucsb.edu
##########################

library(dplyr)
library(stringr)

base = 'https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/'
RData_path= "./Reproducing results in the paper/Data/RData as of Sep 20 2020"
file_pow_param_path = "./Reproducing results in the paper/Data/"
death = 'time_series_covid19_deaths_'

source(file = "./Reproducing results in the paper/Codes/functions_8th_version.R")

# download data from JHU
us_death = read.csv(paste0(base,death,"US.csv"))

death_rate = function(x){x/us_death$Population*10^6}
us_death_rate = us_death %>%
  mutate_at(vars(starts_with("X")), death_rate)

start_date_ini = as.Date("2020-03-21")
end_date = as.Date("2020-09-20")

end_date_char='X9.20.20'

n_ini = as.numeric(end_date - start_date_ini) + 1

col_names = colnames(us_death)
all_dates = as.Date(paste0(str_sub(col_names[13:dim(us_death)[2]], 2, -1), "20"), tryFormats = c("%m.%d.%Y"))

us_pow_param = read.csv(paste0(file_pow_param_path,"Name_list_of_US_states_updated.csv"))

save_pred_death_7_all = NULL
save_pred_death_21_all = NULL
real_pred_death_7_all = NULL
real_pred_death_21_all = NULL
fit_death_all =  NULL
real_death_all = NULL

county_diff_sum = 0

n_US = 0

RMSE_save_all_7=as.list(rep(0,dim(us_pow_param)[1]))
interval_length_save_all_7=as.list(rep(0,dim(us_pow_param)[1]))
coverage_save_all_7=as.list(rep(0,dim(us_pow_param)[1]))

RMSE_save_all_7_death_rate = as.list(rep(0,dim(us_pow_param)[1]))
interval_length_save_all_7_death_rate = as.list(rep(0,dim(us_pow_param)[1]))
coverage_save_all_7_death_rate = as.list(rep(0,dim(us_pow_param)[1]))

RMSE_save_all_21=as.list(rep(0,dim(us_pow_param)[1]))
interval_length_save_all_21=as.list(rep(0,dim(us_pow_param)[1]))
coverage_save_all_21=as.list(rep(0,dim(us_pow_param)[1]))

RMSE_save_all_21_death_rate = as.list(rep(0,dim(us_pow_param)[1]))
interval_length_save_all_21_death_rate = as.list(rep(0,dim(us_pow_param)[1]))
coverage_save_all_21_death_rate = as.list(rep(0,dim(us_pow_param)[1]))

cor_all_7 = c()
cor_all_21 = c()

R_square_all_7 = c()
R_square_all_21 = c()


for(state_index in 1:dim(us_pow_param)[1]){
  state_name = as.character(us_pow_param$state_name[state_index]) 
  state_name_short = as.character(us_pow_param$state_name_short[state_index])
  
  skip_to_next <- FALSE
  
  tryCatch({  load( paste0(RData_path,"/", state_name_short,"_results_in_GP.RData")) }
           , error = function(e) { skip_to_next <<- TRUE})
  if(skip_to_next) { 
    print(paste0(state_name, " No DATA"))
    next }
  
  SIRDC_county_population = us_death %>%
    filter(Admin2 %in% county_names, Province_State == state_name) %>%
    dplyr::select(Population)
  
  RMSE_save_all_7[[state_index]]=rep(0,length(county_names))
  interval_length_save_all_7[[state_index]]=rep(0,length(county_names))
  coverage_save_all_7[[state_index]]=rep(0,length(county_names))
  
  RMSE_save_all_21[[state_index]]=rep(0,length(county_names))
  interval_length_save_all_21[[state_index]]=rep(0,length(county_names))
  coverage_save_all_21[[state_index]]=rep(0,length(county_names))
  
  RMSE_save_all_7_death_rate[[state_index]]=rep(0,length(county_names))
  interval_length_save_all_7_death_rate[[state_index]]=rep(0,length(county_names))
  coverage_save_all_7_death_rate[[state_index]]=rep(0,length(county_names))

  RMSE_save_all_21_death_rate[[state_index]]=rep(0,length(county_names))
  interval_length_save_all_21_death_rate[[state_index]]=rep(0,length(county_names))
  coverage_save_all_21_death_rate[[state_index]]=rep(0,length(county_names))

  save_pred_death_7=NULL
  save_pred_death_21=NULL
  real_pred_death_7=NULL
  real_pred_death_21=NULL
  
  fit_death = NULL
  real_death = NULL
  
  n_state = 0
  
  for(county_index in 1:length(county_names)){
    
    
    if((is.na(results_overall_list)[county_index])){
      print( paste0("Skip county ", county_names[county_index], " in state ", state_name))
      next
    }
    
    
    
    overall_results = results_overall_list[[county_index]]
    
    N = SIRDC_county_population$Population[county_index]
    n = sum(overall_results$prediction_indicator==0)
    
    n_state = n_state + 1
    
    prediction_length =  sum(overall_results$prediction_indicator==1)
    
    prediction_start_index = which(overall_results$prediction_indicator==1)[1]
    
    death_all = us_death %>%
      filter(Admin2==county_names[county_index], Province_State == state_name) %>%
      dplyr::select(starts_with("X"))
    
    index_start_pred=which(names(death_all)==end_date_char)+1
    
    death_selected_compare=as.numeric(death_all[index_start_pred:length(death_all)])
    
    gp_pred_mean=overall_results$death_f_gp[(n+1):length(overall_results$death_f_gp)]
    gp_pred_lower=overall_results$death_lower[(n+1):length(overall_results$death_lower)]
    gp_pred_upper=overall_results$death_upper[(n+1):length(overall_results$death_upper)]
    
    # Death toll
    
    RMSE_7=sqrt(mean((gp_pred_mean[1:7]-death_selected_compare[1:7])^2))
    CI_length_7 = mean(gp_pred_upper[1:7] - gp_pred_lower[1:7])
    CI_Coverage_7 = sum(death_selected_compare[1:7]>=gp_pred_lower[1:7] & death_selected_compare[1:7]<=gp_pred_upper[1:7] ) / 7
    
    RMSE_21=sqrt(mean((gp_pred_mean[1:21]-death_selected_compare[1:21])^2))
    CI_length_21 = mean(gp_pred_upper[1:21] - gp_pred_lower[1:21])
    CI_Coverage_21 = sum(death_selected_compare[1:21]>=gp_pred_lower[1:21] & death_selected_compare[1:21]<=gp_pred_upper[1:21] ) / 21
    
    RMSE_save_all_7[[state_index]][county_index]=RMSE_7
    interval_length_save_all_7[[state_index]][county_index]=CI_length_7
    coverage_save_all_7[[state_index]][county_index]=CI_Coverage_7
    
    RMSE_save_all_21[[state_index]][county_index]=RMSE_21
    interval_length_save_all_21[[state_index]][county_index]=CI_length_21
    coverage_save_all_21[[state_index]][county_index]=CI_Coverage_21
    
    save_pred_death_21=c(save_pred_death_21,gp_pred_mean[1:21])
    save_pred_death_7=c(save_pred_death_7,gp_pred_mean[1:7])
    
    real_pred_death_21=c(real_pred_death_21,death_selected_compare[1:21])
    real_pred_death_7=c(real_pred_death_7,death_selected_compare[1:7])
    
    # Death rate
    
    RMSE_7_death_rate = sqrt(mean(((gp_pred_mean[1:7]-death_selected_compare[1:7])/N*10^6)^2,na.rm=T))
    CI_length_7_death_rate = mean(gp_pred_upper[1:7] - gp_pred_lower[1:7],na.rm=T)/N*10^6
    CI_Coverage_7_death_rate = sum(death_selected_compare[1:7]>=gp_pred_lower[1:7] & death_selected_compare[1:7]<=gp_pred_upper[1:7] ) / 7

    RMSE_21_death_rate = sqrt(mean(((gp_pred_mean[1:21]-death_selected_compare[1:21])/N*10^6)^2,na.rm=T))
    CI_length_21_death_rate = mean(gp_pred_upper[1:21] - gp_pred_lower[1:21],na.rm=T)/N*10^6
    CI_Coverage_21_death_rate = sum(death_selected_compare[1:21]>=gp_pred_lower[1:21] & death_selected_compare[1:21]<=gp_pred_upper[1:21] ) / 21

    RMSE_save_all_7_death_rate[[state_index]][county_index]=RMSE_7_death_rate
    interval_length_save_all_7_death_rate[[state_index]][county_index]=CI_length_7_death_rate
    coverage_save_all_7_death_rate[[state_index]][county_index]=CI_Coverage_7_death_rate

    RMSE_save_all_21_death_rate[[state_index]][county_index]=RMSE_21_death_rate
    interval_length_save_all_21_death_rate[[state_index]][county_index]=CI_length_21_death_rate
    coverage_save_all_21_death_rate[[state_index]][county_index]=CI_Coverage_21_death_rate

    index_start_fit=which(!is.na(overall_results$death_f_gp))[1]
    
    fit_death=c(fit_death,overall_results$death_f[index_start_fit:n])

    real_death=c(real_death,overall_results$death_real[index_start_fit:n])

  }
  
  n_US = n_US + n_state
  
  save_pred_death_21_all = c(save_pred_death_21_all, save_pred_death_21)
  save_pred_death_7_all = c(save_pred_death_7_all, save_pred_death_7)
  
  real_pred_death_21_all = c(real_pred_death_21_all, real_pred_death_21)
  real_pred_death_7_all = c(real_pred_death_7_all, real_pred_death_7)
  
  fit_death_all = c(fit_death_all, fit_death)
  real_death_all = c(real_death_all, real_death)
  
  cor_all_7 = c(cor_all_7, cor(save_pred_death_7, real_pred_death_7))
  cor_all_21 = c(cor_all_21, cor(save_pred_death_21, real_pred_death_21))
  
  seven_day_lm_each = lm(real_pred_death_7~save_pred_death_7)
  twenty_one_day_lm_each = lm(real_pred_death_21~save_pred_death_21)
  
  R_square_all_7 = c(R_square_all_7, summary(seven_day_lm_each)$r.squared)
  R_square_all_21 = c(R_square_all_21, summary(twenty_one_day_lm_each)$r.squared)

}

mean(cor_all_7)
mean(cor_all_21)

mean(R_square_all_7)
mean(R_square_all_21)

plot(save_pred_death_7_all,real_pred_death_7_all)

sqrt(mean((real_pred_death_7_all-save_pred_death_7_all)^2, na.rm = T))
sqrt(mean((real_pred_death_21_all-save_pred_death_21_all)^2, na.rm = T))

seven_day_lm = lm(real_pred_death_7_all~save_pred_death_7_all)

cor(save_pred_death_7_all,real_pred_death_7_all)

plot(save_pred_death_7_all,real_pred_death_7_all, pch=20,
     cex=1,col="blue", xlab = "7 Day Predicted Death Toll", ylab = "7 Day Observed Death Toll", main = paste0("48 US states, R^2=", round(summary(seven_day_lm)$r.squared,3), ", correlation=", round(cor(save_pred_death_7_all,real_pred_death_7_all),3)))
abline(0,1, col="red")

twenty_one_day_lm = lm(real_pred_death_21_all~save_pred_death_21_all)

summary(twenty_one_day_lm)$r.squared


plot(save_pred_death_21_all,real_pred_death_21_all, pch=20, col="blue", xlab = "20 Day Predicted Death", ylab = "20 Day Observed Death", main = paste0("48 US states, R^2=", round(summary(twenty_one_day_lm)$r.squared,3), ", correlation=", round(cor(save_pred_death_21_all,real_pred_death_21_all),3)))
abline(0,1, col="red")

round(cor(save_pred_death_21_all,real_pred_death_21_all),4)
round(summary(twenty_one_day_lm)$r.squared,4)
round(cor(save_pred_death_7_all,real_pred_death_7_all),4)
round(summary(seven_day_lm)$r.squared,4)

fit_lm = lm(real_death_all~fit_death_all)

round(cor(real_death_all,fit_death_all),4)
round(summary(fit_lm)$r.squared,4)


###############################################################
save_file_path = "./Reproducing results in the paper/Results/Fig_1/"

cairo_ps(file = paste0(save_file_path, "pred_7_21_day_as_of_Seq_20.eps"), onefile = FALSE, fallback_resolution = 600, width = 10, height = 10 / 12.5 * 3.3)



state_index=1


state_name = as.character(us_pow_param$state_name[state_index]) 
state_name_short = as.character(us_pow_param$state_name_short[state_index])

skip_to_next <- FALSE

tryCatch({  load( paste0(RData_path,"/", state_name_short,"_results_in_GP.RData")) }
         , error = function(e) { skip_to_next <<- TRUE})
if(skip_to_next) { 
  print(paste0(state_name, " No DATA"))
  next }

SIRDC_county_population = us_death %>%
  filter(Admin2 %in% county_names, Province_State == state_name) %>%
  dplyr::select(Population)

save_pred_death_7=NULL
save_pred_death_21=NULL
real_pred_death_7=NULL
real_pred_death_21=NULL

for(county_index in 1:length(county_names)){
  
  overall_results = results_overall_list[[county_index]]
  
  N = SIRDC_county_population$Population[county_index]
  n = sum(overall_results$prediction_indicator==0)
  prediction_length =  sum(overall_results$prediction_indicator==1)
  
  prediction_start_index = which(overall_results$prediction_indicator==1)[1]
  
  death_all = us_death %>%
    filter(Admin2==county_names[county_index], Province_State == state_name) %>%
    dplyr::select(starts_with("X"))
  
  index_start_pred=which(names(death_all)==end_date_char)+1
  
  death_selected_compare=as.numeric(death_all[index_start_pred:length(death_all)])
  
  gp_pred_mean=overall_results$death_f_gp[(n+1):length(overall_results$death_f_gp)]

  save_pred_death_21=c(save_pred_death_21,gp_pred_mean[1:21])
  save_pred_death_7=c(save_pred_death_7,gp_pred_mean[1:7])
  
  real_pred_death_21=c(real_pred_death_21,death_selected_compare[1:21])
  real_pred_death_7=c(real_pred_death_7,death_selected_compare[1:7])
  
}
par(mfrow = c(1,2))
# par(mgp = c(2.2,0.7,0))
par(cex.lab=1, cex.axis=1, cex.main=1.3, cex.sub=1.5, mai=c(0.55,0.55,0.3,0.3),mgp = c(1.45,0.5,0)) # , mgp = c(1.5,0.2,0)

plot(save_pred_death_7/10^3,real_pred_death_7/10^3,pch=20,
     col=rainbow(50)[state_index],xlab = expression(paste("7-day predicted death toll ", (10^3))) , 
     ylab =expression(paste("7-day held-out death toll ", (10^3))) ,cex=.25,xlim=c(0,7350/10^3),ylim=c(0,7350/10^3),
     main = expression(paste(rho%~~%1,', ',rho[county]%~~%0.9132))
)
abline(0,1, col="black")


for(state_index in 2:dim(us_pow_param)[1]){
  state_name = as.character(us_pow_param$state_name[state_index]) 
  state_name_short = as.character(us_pow_param$state_name_short[state_index])
  
  skip_to_next <- FALSE
  
  tryCatch({  load( paste0(RData_path,"/", state_name_short,"_results_in_GP.RData")) }
           , error = function(e) { skip_to_next <<- TRUE})
  if(skip_to_next) { 
    print(paste0(state_name, " No DATA"))
    next }

  SIRDC_county_population = us_death %>%
    filter(Admin2 %in% county_names, Province_State == state_name) %>%
    dplyr::select(Population)
  
  save_pred_death_7=NULL
  save_pred_death_21=NULL
  real_pred_death_7=NULL
  real_pred_death_21=NULL
  
  for(county_index in 1:length(county_names)){
    
    overall_results = results_overall_list[[county_index]]
    
    N = SIRDC_county_population$Population[county_index]
    n = sum(overall_results$prediction_indicator==0)
    prediction_length =  sum(overall_results$prediction_indicator==1)
    
    prediction_start_index = which(overall_results$prediction_indicator==1)[1]
    
    death_all = us_death %>%
      filter(Admin2==county_names[county_index], Province_State == state_name) %>%
      dplyr::select(starts_with("X"))
    
    index_start_pred=which(names(death_all)==end_date_char)+1
    
    death_selected_compare=as.numeric(death_all[index_start_pred:length(death_all)])

    gp_pred_mean=overall_results$death_f_gp[(n+1):length(overall_results$death_f_gp)]

    save_pred_death_21=c(save_pred_death_21,gp_pred_mean[1:21])
    save_pred_death_7=c(save_pred_death_7,gp_pred_mean[1:7])
    
    real_pred_death_21=c(real_pred_death_21,death_selected_compare[1:21])
    real_pred_death_7=c(real_pred_death_7,death_selected_compare[1:7])
    
  }

  lines(save_pred_death_7/10^3,real_pred_death_7/10^3,pch=20,
        col=rainbow(50)[state_index],cex=0.75,type='p')
}

###################

state_index=1

state_name = as.character(us_pow_param$state_name[state_index]) 
state_name_short = as.character(us_pow_param$state_name_short[state_index])

skip_to_next <- FALSE

tryCatch({  load( paste0(RData_path,"/", state_name_short,"_results_in_GP.RData")) }
         , error = function(e) { skip_to_next <<- TRUE})
if(skip_to_next) { 
  print(paste0(state_name, " No DATA"))
  next }

SIRDC_county_population = us_death %>%
  filter(Admin2 %in% county_names, Province_State == state_name) %>%
  dplyr::select(Population)

save_pred_death_7=NULL
save_pred_death_21=NULL
real_pred_death_7=NULL
real_pred_death_21=NULL

for(county_index in 1:length(county_names)){
  
  overall_results = results_overall_list[[county_index]]
  
  N = SIRDC_county_population$Population[county_index]
  n = sum(overall_results$prediction_indicator==0)
  prediction_length =  sum(overall_results$prediction_indicator==1)
  
  prediction_start_index = which(overall_results$prediction_indicator==1)[1]
  
  death_all = us_death %>%
    filter(Admin2==county_names[county_index], Province_State == state_name) %>%
    dplyr::select(starts_with("X"))
  
  index_start_pred=which(names(death_all)==end_date_char)+1
  
  death_selected_compare=as.numeric(death_all[index_start_pred:length(death_all)])

  gp_pred_mean=overall_results$death_f_gp[(n+1):length(overall_results$death_f_gp)]

  save_pred_death_21=c(save_pred_death_21,gp_pred_mean[1:21])
  save_pred_death_7=c(save_pred_death_7,gp_pred_mean[1:7])
  
  real_pred_death_21=c(real_pred_death_21,death_selected_compare[1:21])
  real_pred_death_7=c(real_pred_death_7,death_selected_compare[1:7])
  
}

plot(save_pred_death_21/10^3,real_pred_death_21/10^3,pch=20,
     col=rainbow(50)[state_index],xlab = expression(paste("21-day predicted death toll ", (10^3))) , 
     ylab =  expression(paste("21-day held-out death toll ", (10^3))),cex=.25,xlim=c(0,7350/10^3),ylim=c(0,7350/10^3),
     main = expression(paste(rho%~~%1,', ',rho[county]%~~%0.9287)))
abline(0,1, col="black")



for(state_index in 2:dim(us_pow_param)[1]){

  state_name = as.character(us_pow_param$state_name[state_index]) 
  state_name_short = as.character(us_pow_param$state_name_short[state_index])
  
  skip_to_next <- FALSE
  
  tryCatch({  load( paste0(RData_path,"/", state_name_short,"_results_in_GP.RData")) }
           , error = function(e) { skip_to_next <<- TRUE})
  if(skip_to_next) { 
    print(paste0(state_name, " No DATA"))
    next }else{
      print(paste0(state_name))
    }
  
  SIRDC_county_population = us_death %>%
    filter(Admin2 %in% county_names, Province_State == state_name) %>%
    dplyr::select(Population)
  
  save_pred_death_7=NULL
  save_pred_death_21=NULL
  real_pred_death_7=NULL
  real_pred_death_21=NULL
  
  for(county_index in 1:length(county_names)){
    
    overall_results = results_overall_list[[county_index]]
    
    N = SIRDC_county_population$Population[county_index]
    n = sum(overall_results$prediction_indicator==0)
    prediction_length =  sum(overall_results$prediction_indicator==1)
    
    prediction_start_index = which(overall_results$prediction_indicator==1)[1]
    
    death_all = us_death %>%
      filter(Admin2==county_names[county_index], Province_State == state_name) %>%
      dplyr::select(starts_with("X"))
    
    index_start_pred=which(names(death_all)==end_date_char)+1
    
    death_selected_compare=as.numeric(death_all[index_start_pred:length(death_all)])

    gp_pred_mean=overall_results$death_f_gp[(n+1):length(overall_results$death_f_gp)]

    save_pred_death_21=c(save_pred_death_21,gp_pred_mean[1:21])
    save_pred_death_7=c(save_pred_death_7,gp_pred_mean[1:7])
    
    real_pred_death_21=c(real_pred_death_21,death_selected_compare[1:21])
    real_pred_death_7=c(real_pred_death_7,death_selected_compare[1:7])
    
  }

  
  lines(save_pred_death_21/10^3,real_pred_death_21/10^3,pch=20,
        col=rainbow(50)[state_index],cex=0.75,type='p')
}
dev.off()