##########################
# Code for figure 6 & 7 in the Extended data figures section
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

gamma = 0.2
theta = 0.1
delta = 0.0066

set.seed(1)

base = 'https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/'
death = 'time_series_covid19_deaths_'
confirm = 'time_series_covid19_confirmed_'

covid19_project_url = "https://api.covidtracking.com/v1/states/daily.csv"

us_test = read.csv(covid19_project_url)

file_path_cur= "./Reproducing results in the paper/Data/RData as of Sep 20 2020"

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

start_date_ini = as.Date("2020-03-21")

end_date = as.Date("2020-09-20")

state_name  = "California"
state_name_short = "CA"

state_population = us_death%>%
  filter( Province_State == state_name)%>%
  select(Admin2, Population)

file_path_rdata = paste0(file_path_cur, "/", state_name_short ,"_results_in_GP.RData")

load(file_path_rdata)

training_data_length = 60

standard_time_seq = seq.Date(end_date - training_data_length + 1, end_date, by=1) #start_date_ini

n_standard = length(standard_time_seq)

prediction_length  = 21

date_seq_all = seq.Date(end_date - training_data_length + 1, end_date+prediction_length, by=1)

save_file_path = "./Reproducing results in the paper/Results/Extended_fig_6_7/"

# pdf(file = paste0(file_path, "FL_predictions_2.pdf"))


cairo_ps(file = paste0(save_file_path, "CA_predictions.eps"), onefile = FALSE, fallback_resolution = 600)


par(mar=c(1.3,1.3,1.3,1.3),mgp=c(3,0.5,0), las=0)

par(mfrow = c(10,5))

for(county_index in 1:50){
  county_name_selected = county_names[county_index]
  
  county_result_selected = results_overall_list[[which(county_names == county_name_selected)]]
  
  N = state_population$Population[state_population$Admin2 == county_name_selected]
  
  n = sum(county_result_selected$prediction_indicator==0)

  projection = county_result_selected$death_f_gp[(n+1):(n+prediction_length)]
  projection_upper = county_result_selected$death_upper[(n+1):(n+prediction_length)]
  projection_lower = county_result_selected$death_lower[(n+1):(n+prediction_length)]
  
  real_death = us_death %>%
    filter(Admin2 == county_name_selected, Province_State == state_name) %>%
    dplyr::select(starts_with("x"))
  
  death_selected = as.numeric(real_death[(all_dates >= (end_date-training_data_length+1)) & (all_dates <= (end_date+prediction_length))])
  
  for(i in 2:length(death_selected)){
    if(death_selected[i] < death_selected[i-1]){
      death_selected[i] = death_selected[i-1]
    }
  }
  
  date_seq_real_death = seq.Date(end_date-training_data_length+1, end_date+prediction_length, by=1)
  
  date_seq_pred = county_result_selected$date[(n+1):(n+prediction_length)]
  
  ylimit_death = range(death_selected, projection_upper)
  
  plot(c(rep(NA,n_standard),projection)~date_seq_all,ylim = ylimit_death, type="l", col="blue",xlab = NULL, ylab = NULL,
       main = paste0(county_name_selected, ", population=", round(N/10^6,2),"M"),cex.lab=0.2, cex.axis=0.5, cex.main=0.8, cex.sub=0.5)
  
  polygon(c(rev(date_seq_pred), date_seq_pred),
          c(rev(projection_upper),projection_lower),
          col = 'grey80',
          border = NA)
  lines(projection~date_seq_pred, col = "blue")
  lines(death_selected~date_seq_real_death, col = "red")
  abline(v = end_date+1, lty=2)

}

dev.off()
