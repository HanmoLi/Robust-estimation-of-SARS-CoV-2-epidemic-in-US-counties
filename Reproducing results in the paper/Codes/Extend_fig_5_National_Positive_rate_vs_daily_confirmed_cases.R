##########################
# Code for figure 5 in the Extended data figures section
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

source(file = "./Reproducing results in the paper/Codes/functions_8th_version.R")

##for forecast
fitted_days_beta=184 ##how many days from the end day you use to estimate beta_t
#predicted_days_train=120  
predicted_days_train=184  

prediction_length = 90
# prediction_length=90
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

# parallel settings
numCores <- detectCores()

standard_time_seq = seq.Date(start_date, end_date, by=1)

n_standard  =length(standard_time_seq)

standard_time_seq_extended = seq.Date(start_date-1, end_date, by=1)

n_standard_extened  =length(standard_time_seq_extended)


# get the national fittings

N = 321418820 # the us population

US_confirmed_overall = rep(0, n_standard_extened)
US_death_overall = rep(0, n_standard_extened)

# us death from 2020-01-22 to today
US_death_overall_all_length = rep(0, length(all_dates))

# sum of state level confirmed and death cases
for(state_index in 1:dim(us_pow_param)[1]){
  
  state_name = us_pow_param$state_name[state_index]
  state_name_short = us_pow_param$state_name_short[state_index]
  
  state_confirmed = us_confirm %>%
    filter(Province_State == state_name)%>%
    select(starts_with("x"))
  
  state_confirmed_selected = state_confirmed[all_dates %in% standard_time_seq_extended]
  
  US_confirmed_overall = US_confirmed_overall + as.numeric(apply(state_confirmed_selected, 2, sum))
  
  state_death = us_death %>%
    filter(Province_State == state_name)%>%
    select(starts_with("x"))
  
  US_death_overall_all_length = US_death_overall_all_length + as.numeric(apply(state_death,2,sum))
  
  state_death_selected = state_death[all_dates %in% standard_time_seq_extended]
  
  US_death_overall = US_death_overall + as.numeric(apply(state_death_selected,2,sum))
}

confirm_selected = US_confirmed_overall
death_selected = US_death_overall

# eliminate decreasing trend in the county death data

for(i in 1:(n_standard_extened-1)){
  if (death_selected[i+1] < death_selected[i]){
    death_selected[i+1] = death_selected[i]
  }
}

for( i in 1:(n_standard_extened-1)){
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
  select(date, state, positiveIncrease, totalTestResultsIncrease)

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

save_file_path = "./Reproducing results in the paper/Results/Extended_fig_5/"

p1 = us_test_aggregated %>%
  ggplot(aes(x=date, y=positive_rate)) +
  geom_line(col = "Red", size = 1.1) +
  ylab("7-day averaged test positive rate in the US")+
  xlab("Date")+ 
  theme(text = element_text(size = 20),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.key.width=unit(1,"cm"),
        axis.text.y = element_text(angle=90, hjust=1))

daily_confirm_smoothed_df = data.frame(date = standard_time_seq, value = daily_confirm_selected_avg)

p2 = daily_confirm_smoothed_df %>%
  ggplot(aes(x=date, y=value/10^4)) +
  geom_bar(stat = 'identity', color="white",  fill="#ff8540", width = 0.757) +
  ylab(expression(paste("7-day averaged daily confirmed cases in the US  ", (10^4))))+
  xlab("Date")+ 
  theme(text = element_text(size = 20),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.key.width=unit(1,"cm"),
        axis.text.y = element_text(angle=90, hjust=1))

ggarrange(p2, p1, ncol = 2, heights = 15, widths = 50) #%>%
  ggsave(file=paste0(save_file_path, "US_daily_confirmed_cases_vs_positive_rate.eps"), width = 40, height = 12, units = "cm")

