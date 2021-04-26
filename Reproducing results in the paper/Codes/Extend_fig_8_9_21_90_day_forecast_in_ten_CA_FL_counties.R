##########################
# Code for figure 8 & 9 in the Extended data figures section
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

file_path_cur= "./Reproducing results in the paper/Data/RData as of Nov 30 2020"

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

file_pow_param_path = "./Reproducing results in the paper/Data/"

us_pow_param = read.csv(paste0(file_pow_param_path,"Name_list_of_US_states_updated.csv"))
us_pow_param$state_name = as.character(us_pow_param$state_name)
state_names_all = us_pow_param$state_name

start_date_ini = as.Date("2020-03-21")

end_date = as.Date("2020-11-30")

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

save_file_path = "./Reproducing results in the paper/Results/Extended_fig_8_9/"

cairo_ps(file = paste0(save_file_path, "CA_top_ten_21_day_predictions_from_Dec.eps"), onefile = FALSE, fallback_resolution = 600, width = 16, height = 4)


par(mar=c(1.5,1.5,1.5,1.5),mgp=c(3,0.5,0), las=0)

par(mfrow = c(2,5))

# if (state_name_short == "FL"){
#   county_seq = c(43, 6, 28,15,  5, 1, 36, 45, 14, 63)
# }

ten_largest_counties = state_population$Admin2[order(state_population$Population, decreasing = T)][1:10]



for(county_index in ten_largest_counties){
  county_name_selected = county_index
  
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
       main = paste0(county_name_selected, " (", round(N/10^6,2),"M)") ,cex.lab=1.4, cex.axis=1.4, cex.main=1.8, cex.sub=1) # ,cex.lab=1.2, cex.axis=1.2, cex.main=1.6, cex.sub=1 #cex.lab=0.8, cex.axis=0.8, cex.main=1.2, cex.sub=1
  
  polygon(c(rev(date_seq_pred), date_seq_pred),
          c(rev(projection_upper),projection_lower),
          col = 'grey80',
          border = NA)
  lines(projection~date_seq_pred, col = "blue")
  lines(death_selected~date_seq_real_death, col = "red")
  # lines(out_ode_all[,5]~date_seq_pred, col = "green")
  abline(v = end_date+1, lty=2)
  # axis(c(1,2),tck=0)
  # legend("topleft", legend = c("Observed Death", "Estimated Death"), lty = c(1,1), col = c("red", "blue"))
  
}

dev.off()
