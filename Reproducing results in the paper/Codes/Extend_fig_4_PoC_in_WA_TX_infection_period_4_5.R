##########################
# Code for figure 4 in the Extended data figures section
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

if (!require("urbnmapr", character.only = TRUE)){
  #make sure you installed the urbnmapr package
  devtools::install_github("UrbanInstitute/urbnmapr")
} else(
  library(urbnmapr)
)

source(file = "./Reproducing results in the paper/Codes/functions_8th_version.R")

gamma = 0.2
theta = 0.1
delta = 0.0066

# state_name = "California"
# state_name_short = "CA"
# county_name_selected = "Santa Barbara"

set.seed(1)

base = 'https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/'
death = 'time_series_covid19_deaths_'
confirm = 'time_series_covid19_confirmed_'

covid19_project_url = "https://api.covidtracking.com/v1/states/daily.csv"

file_pow_param_path = "./Reproducing results in the paper/Data/"
file_path_cur= "./Reproducing results in the paper/Data/RData as of Sep 20 2020"

# download death and confirmed data from JHU
us_death = read.csv(paste0(base,death,"US.csv"))

# clean the death data
us_death = clean_us_death_for_map(us_death)

col_names = colnames(us_death)
all_dates = as.Date(paste0(str_sub(col_names[13:dim(us_death)[2]], 2, -1), "20"), tryFormats = c("%m.%d.%Y"))

start_date_ini = as.Date("2020-03-21")

end_date = as.Date("2020-09-20")

us_pow_param = read.csv(paste0(file_pow_param_path,"Name_list_of_US_states_updated.csv"))
us_pow_param$state_name = as.character(us_pow_param$state_name)
state_names_all = us_pow_param$state_name

dates_extract_prob = c(as.Date("2020-09-20"))

##########################
# Map in WA
##########################

state_prob_df_for_map = data.frame(state = character(),
                                   county = character(),
                                   county_fips = character(),
                                   prob_date_1 = double(),
                                   stringsAsFactors=FALSE)

state_name = "Washington"
state_name_short = "WA"

print(state_name)

state_population = us_death%>%
  filter( Province_State == state_name)%>%
  select(Admin2, Population)

file_path_rdata = paste0(file_path_cur, "/", state_name_short ,"_results_in_GP.RData")

skip_to_next <- FALSE

load(file_path_rdata)

sum_idx = 1

for(county_index in 1:length(county_names)){
  
  county_name_selected = county_names[county_index]
  
  county_fip_each = as.character(us_death %>%
                                   filter( Province_State == state_name, Admin2 == county_name_selected)%>%
                                   select(county_fips))
  
  results_in_selected_county = results_overall_list[[county_index]]
  
  n = sum(results_in_selected_county$prediction_indicator==0)
  
  N = as.numeric(state_population$Population[which(state_population$Admin2==county_name_selected)])
  
  start_date = results_in_selected_county$date[1]
  
  beta_t_daily = results_in_selected_county$Beta_estimated[1:(n-1)]
  
  I_0 = results_in_selected_county$infective_f[1]
  R_0 = results_in_selected_county$resolved_f[1]
  
  death_selected = results_in_selected_county$death_real[1:n]
  
  S_t_seq = results_in_selected_county$suspective_f[1:n]
  
  gamma_updated = 1/4.5
  
  approx_beta_seq_all = beta_t_daily
  
  init_for_beta = c(S_t_seq[1], I_0, R_0, death_selected[1], 0)
  
  param_record_approx_all = matrix(NA, 5, n)
  param_record_approx_all[,1] = init_for_beta
  
  for (i in 1:(n-1)){
    S_t_1 = param_record_approx_all[1,i]
    I_t_1 = param_record_approx_all[2,i]
    R_t_1 = param_record_approx_all[3,i]
    D_t_1 = param_record_approx_all[4,i]
    C_t_1 = param_record_approx_all[5,i]
    
    beta_t_1_2 = approx_beta_seq_all[i]
    
    if(I_t_1<1){
      I_t_1 = 1
    }
    
    S_t_2 = uniroot(find_root_S_t_2, c(0, N), tol = 0.0001, param = c(S_t_1, beta_t_1_2, I_t_1), N = N, gamma=gamma_updated)
    I_t_2 = I_t_1 * exp(beta_t_1_2*(S_t_1 + S_t_2$root)/(2*N) - gamma_updated)
    R_t_2 = (2-theta)/(2+theta)*R_t_1 + gamma_updated/(2+theta)*(I_t_1+I_t_2)
    D_t_2 = D_t_1 + delta*theta*(R_t_1+R_t_2)/2
    C_t_2 = C_t_1 + (1-delta)*theta*(R_t_1+R_t_2)/2
    
    param_record_approx_all[, i+1] = c(S_t_2$root, I_t_2, R_t_2, D_t_2, C_t_2)
    
    
  }
  
  
  R_eff = beta_t_daily/gamma_updated * param_record_approx_all[1, 2:(n)] / N
  
  possibility_infected = param_record_approx_all[2, 2:(n)] * R_eff * gamma_updated / (param_record_approx_all[1, 2:(n)])
  
  possibility_infected_7_day_averaged = data_seven_day_smoothing(possibility_infected) #rollapply(possibility_infected, width = 7, by = 1, FUN = mean, align = "right")
  
  corresponding_date = results_in_selected_county$date[(2):(n)]
  
  probs_vector = c()
  
  for(i in 1:length(dates_extract_prob)){
    if(dates_extract_prob[i] < corresponding_date[1]){
      probs_vector = c(probs_vector, NA)
    } else{
      probs_vector = c(probs_vector, possibility_infected_7_day_averaged[which(corresponding_date==dates_extract_prob[i])])
    }
  }
  
  results_each = c(state_name, county_name_selected, county_fip_each, probs_vector)

  state_prob_df_for_map[sum_idx, ] = results_each
  
  sum_idx = sum_idx + 1
  
}

state_prob_df_for_map$prob_date_1 = as.numeric(state_prob_df_for_map$prob_date_1) * 100

state_prob_df_for_map$prob_date_1[state_prob_df_for_map$prob_date_1<10^(-4)] = 10^(-4)

state_prob_df_for_map$prob_date_1_scaled = log(state_prob_df_for_map$prob_date_1,2 )

range(state_prob_df_for_map$prob_date_1_scaled, na.rm=T)

state_counties = counties %>%
  filter(state_abbv==state_name_short)

paleta <- c("#fefec8", "#fdd88a", "#fdae76", "#f47562", "#d04c67")

state_prob_df_for_map$prob_date_1_f <- cut(state_prob_df_for_map$prob_date_1, breaks = c(0, 0.001, 0.01, 0.1, 1, 100))
niveles <- levels(state_prob_df_for_map$prob_date_1_f)

legend_text = c(expression(""<=10^{-3}),
                expression(paste("(",10^{-3}, ",", 10^{-2}, "]", sep="")),
                expression(paste("(",10^{-2}, ",", 10^{-1}, "]", sep="")),
                expression(paste("(",10^{-1}, ",", 1, "]", sep="")),
                expression("">1))

state_county_infect_prob_joint_WA <- left_join(state_counties,state_prob_df_for_map,  by = "county_fips")




# png(filename = paste0(file_path, "Map-infeciton prob in WA with infection period 4.75.png"), width = 500, height = 450)

WA_PoC_map = state_county_infect_prob_joint_WA %>%
  ggplot(aes(long, lat, group = group)) +
  geom_polygon(aes(fill = (prob_date_1_f))) +
  geom_polygon(color = "black", fill = NA) +
  scale_fill_manual(name = "PoC (%)", values = paleta, breaks = niveles, labels = legend_text)+
  coord_map() +
  labs(fill = expression("Odds (%)")) +
  xlab("Longitude") + ylab("Latitude") + # for the y axis label
  theme(
    legend.position = "right",
    text = element_text(size = 15),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 15),
    legend.key.size = unit(1, "cm"),
    axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    axis.title.y=element_blank(),
    axis.text.y =element_blank(),
    axis.ticks.y=element_blank()
  ) +
  guides(shape = guide_legend(override.aes = list(size = 1)),
         color = guide_legend(override.aes = list(size = 1))) 
# dev.off()

#################################
# the rest two figures in WA
#################################

state_name  = "Washington"
state_name_short = "WA"

state_population = us_death%>%
  filter( Province_State == state_name)%>%
  select(Admin2, Population)

file_path_rdata = paste0(file_path_cur, "/", state_name_short ,"_results_in_GP.RData")

load(file_path_rdata)


# display the 10 best and 10 worest counties in WA
possibility_infected_record_matrix = matrix(NA, length(county_names), as.numeric(end_date - start_date_ini)+1)

standard_time_seq = seq.Date(start_date_ini, end_date, by=1)

for(county_index in 1:length(county_names)){
  
  county_name_selected = county_names[county_index]
  
  county_results_each = results_overall_list[[county_index]]
  
  # n = sum(county_results_each$prediction_indicator==0)
  
  N = as.numeric(state_population$Population[which(state_population$Admin2==county_name_selected)])
  
  start_date = county_results_each$date[1]
  
  # end_date = as.Date("2020-09-20")
  
  n = as.numeric(end_date - start_date+1)
  
  beta_t_daily = county_results_each$Beta_estimated[1:(n-1)]
  
  R_eff = beta_t_daily/0.2 * county_results_each$suspective_f[2:n] / N
  
  possibility_infected = county_results_each$infective_f[2:n] * R_eff / (county_results_each$suspective_f[2:n] * 5)
  
  possibility_infected_7_day_averaged = data_seven_day_smoothing(possibility_infected)#rollapply(possibility_infected, width = 7, by = 1, FUN = mean, align = "right")
  
  time_seq_each = county_results_each$date[(2):(n)]
  
  possibility_infected_record_matrix[county_index, standard_time_seq %in% time_seq_each] = possibility_infected_7_day_averaged
  
}

## top 5 counties

# top 10 counties with biggest probabilities getting infected
top_5_worst_county_names = county_names[order(possibility_infected_record_matrix[,dim(possibility_infected_record_matrix)[2]], decreasing = T)][1:5]

# display the 10 best and 10 worest counties in CA
possibility_infected_record_matrix = matrix(NA, length(top_5_worst_county_names), as.numeric(end_date - start_date_ini)+1)

standard_time_seq = seq.Date(start_date_ini, end_date, by=1)

for(county_index in 1:length(top_5_worst_county_names)){
  
  county_name_selected = top_5_worst_county_names[county_index]
  
  county_results_each = results_overall_list[[which(county_names == county_name_selected)]]
  
  # n = sum(county_results_each$prediction_indicator==0)
  
  N = as.numeric(state_population$Population[which(state_population$Admin2==county_name_selected)])
  
  start_date = county_results_each$date[1]
  
  # end_date = as.Date("2020-09-20")
  
  n = as.numeric(end_date - start_date+1)
  
  gamma_updated = 1/4.5
  
  beta_t_daily = county_results_each$Beta_estimated[1:(n-1)]
  
  I_0 = county_results_each$infective_f[1]
  R_0 = county_results_each$resolved_f[1]
  
  death_selected = county_results_each$death_real[1:n]
  
  S_t_seq = county_results_each$suspective_f[1:n]
  
  approx_beta_seq_all = beta_t_daily
  
  init_for_beta = c(S_t_seq[1], I_0, R_0, death_selected[1], 0)
  
  param_record_approx_all = matrix(NA, 5, n)
  param_record_approx_all[,1] = init_for_beta
  
  for (i in 1:(n-1)){
    S_t_1 = param_record_approx_all[1,i]
    I_t_1 = param_record_approx_all[2,i]
    R_t_1 = param_record_approx_all[3,i]
    D_t_1 = param_record_approx_all[4,i]
    C_t_1 = param_record_approx_all[5,i]
    
    beta_t_1_2 = approx_beta_seq_all[i]
    
    if(I_t_1<1){
      I_t_1 = 1
    }
    
    S_t_2 = uniroot(find_root_S_t_2, c(0, N), tol = 0.0001, param = c(S_t_1, beta_t_1_2, I_t_1), N = N, gamma=gamma_updated)
    I_t_2 = I_t_1 * exp(beta_t_1_2*(S_t_1 + S_t_2$root)/(2*N) - gamma_updated)
    R_t_2 = (2-theta)/(2+theta)*R_t_1 + gamma_updated/(2+theta)*(I_t_1+I_t_2)
    D_t_2 = D_t_1 + delta*theta*(R_t_1+R_t_2)/2
    C_t_2 = C_t_1 + (1-delta)*theta*(R_t_1+R_t_2)/2
    
    param_record_approx_all[, i+1] = c(S_t_2$root, I_t_2, R_t_2, D_t_2, C_t_2)
    
    
  }
  
  
  R_eff = beta_t_daily/gamma_updated * param_record_approx_all[1, 2:(n)] / N
  
  possibility_infected = param_record_approx_all[2, 2:(n)] * R_eff * gamma_updated / (param_record_approx_all[1, 2:(n)])
  
  possibility_infected_7_day_averaged = data_seven_day_smoothing(possibility_infected) #rollapply(possibility_infected, width = 7, by = 1, FUN = mean, align = "right")
  
  time_seq_each = county_results_each$date[(2):(n)]
  
  possibility_infected_record_matrix[county_index, standard_time_seq %in% time_seq_each] = possibility_infected_7_day_averaged
  
}


possibility_infected_top_5_worst = possibility_infected_record_matrix


infection_prob_sep_20 = possibility_infected_top_5_worst[, dim(possibility_infected_top_5_worst)[2]] * 100
names(infection_prob_sep_20) = as.character(top_5_worst_county_names)

# infection_prob_sep_1 = possibility_infected_top_5_worst[, dim(possibility_infected_top_5_worst)[2]-19] * 100
# names(infection_prob_sep_1) = as.character(top_5_worst_county_names)

# infection_prob_sep_1
# infection_prob_sep_20

top_5_worst_df = data.frame(date = rep(standard_time_seq, 5), County = rep(top_5_worst_county_names, each=length(standard_time_seq)), possibility_infected = as.vector(t(possibility_infected_top_5_worst)))

range(top_5_worst_df$possibility_infected*100, na.rm=T)

if(state_name  == "Texas"){
  y_upper_lim = 2.655575
} else if(state_name == "Washington"){
  y_upper_lim = 2.519254e-01
}

# this negative number is due to precision accuracy
top_5_worst_df$possibility_infected[top_5_worst_df$possibility_infected<0] = 0

WA_infection_curve = ggplot(top_5_worst_df, aes(x = date, y = possibility_infected*100,group = County)) + 
  geom_line(aes(color = County, linetype = County), lwd=1, show.legend = F)+
  scale_linetype_manual(values=rep("solid", 5)) +
  xlab("Date") +
  ylab("Probability of contracting COVID-19 (%)")+
  ylim(0,  y_upper_lim) + 
  theme(text = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15),
        axis.text.y = element_text(angle=90, hjust=1)) 

n_standard = length(standard_time_seq)

death_record_df = data.frame(matrix(NA, 10, n_standard+2))

colnames(death_record_df)[1:2] = c("County", "Type")
death_record_df$County = rep(as.character(top_5_worst_county_names), each=2)
death_record_df$Type = rep(c("Real", "Fitted"), 5)

population_record_vector = NULL

for(county_index in 1:5){
  county_name_selected = as.character(top_5_worst_county_names)[county_index]
  
  county_results_each = results_overall_list[[which(county_names == county_name_selected)]]
  
  N = as.numeric(state_population$Population[which(state_population$Admin2==county_name_selected)])
  
  population_record_vector = c(population_record_vector, N)
  
  start_date = county_results_each$date[1]
  
  # end_date = as.Date("2020-09-20")
  
  n = as.numeric(end_date - start_date+1)
  
  print(paste0("Population = ", N, " in ", county_name_selected, ", ", state_name))
  
  county_death_raw = us_death %>%
    filter(Admin2 == county_name_selected, Province_State == state_name) %>%
    select(starts_with("x"))
  
  county_death_raw_selected = as.numeric(county_death_raw[all_dates %in% standard_time_seq])
  
  for(i in 2:length(county_death_raw_selected)){
    if (county_death_raw_selected[i] < county_death_raw_selected[i-1]){
      county_death_raw_selected[i] = county_death_raw_selected[i-1]
    }
  }
  
  death_record_df[2*(county_index-1) + 1, 3:dim(death_record_df)[2]] = county_death_raw_selected / N * 10^6
  
  gamma_updated = 1/4.5
  
  beta_t_daily = county_results_each$Beta_estimated[1:(n-1)]
  
  I_0 = county_results_each$infective_f[1]
  R_0 = county_results_each$resolved_f[1]
  
  death_selected = county_results_each$death_real[1:n]
  
  S_t_seq = county_results_each$suspective_f[1:n]
  
  approx_beta_seq_all = beta_t_daily
  
  init_for_beta = c(S_t_seq[1], I_0, R_0, death_selected[1], 0)
  
  param_record_approx_all = matrix(NA, 5, n)
  param_record_approx_all[,1] = init_for_beta
  
  for (i in 1:(n-1)){
    S_t_1 = param_record_approx_all[1,i]
    I_t_1 = param_record_approx_all[2,i]
    R_t_1 = param_record_approx_all[3,i]
    D_t_1 = param_record_approx_all[4,i]
    C_t_1 = param_record_approx_all[5,i]
    
    beta_t_1_2 = approx_beta_seq_all[i]
    
    if(I_t_1<1){
      I_t_1 = 1
    }
    
    S_t_2 = uniroot(find_root_S_t_2, c(0, N), tol = 0.0001, param = c(S_t_1, beta_t_1_2, I_t_1), N = N, gamma=gamma_updated)
    I_t_2 = I_t_1 * exp(beta_t_1_2*(S_t_1 + S_t_2$root)/(2*N) - gamma_updated)
    R_t_2 = (2-theta)/(2+theta)*R_t_1 + gamma_updated/(2+theta)*(I_t_1+I_t_2)
    D_t_2 = D_t_1 + delta*theta*(R_t_1+R_t_2)/2
    C_t_2 = C_t_1 + (1-delta)*theta*(R_t_1+R_t_2)/2
    
    param_record_approx_all[, i+1] = c(S_t_2$root, I_t_2, R_t_2, D_t_2, C_t_2)
    
    
  }
  
  death_fitted = param_record_approx_all[4,]
  
  if (length(death_fitted) == n_standard){
    death_record_df[2*(county_index-1) + 2, 3:dim(death_record_df)[2]] = death_fitted / N * 10^6
  } else{
    diff_length = n_standard - length(death_fitted)
    
    death_record_df[2*(county_index-1) + 2, (3+diff_length):dim(death_record_df)[2]] = death_fitted / N * 10^6
  }

}

county_name_with_population = paste0(top_5_worst_county_names, " (", round(population_record_vector/10^3, 0), "K)")
death_record_df[,1] = rep(county_name_with_population, each = 2)

death_record_df_long = gather(death_record_df, key="Date", value = "Death", -County, -Type)
death_record_df_long$Date = rep(standard_time_seq, each=10)

death_record_df_long$group = paste0(death_record_df_long$County, death_record_df_long$Type)

death_record_df_long$Type = as.factor(death_record_df_long$Type)

range(death_record_df_long$Death, na.rm = T)

if(state_name  == "Texas"){
  y_upper_lim_death = 3165.126
} else if(state_name == "Washington"){
  y_upper_lim_death = 1034.665
}

WA_death_curve = ggplot(death_record_df_long, aes(x = Date, y = Death,group = group)) + 
  geom_line(aes(color = County, linetype = Type, size = Type), show.legend = T)+
  scale_linetype_manual(values=c("solid", "dotted"))+
  scale_size_manual(values = c(1,1.1))+
  xlab("Date") +
  ylab("Cumulative death per million")+
  guides(colour = guide_legend(order = 1,title="County (Population)"), 
         linetype = guide_legend(order = 2),
         size = FALSE)+
  ylim(0,  y_upper_lim_death) + 
  theme(text = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.key.width=unit(1,"cm"),
        axis.text.y = element_text(angle=90, hjust=1))



#########################
# Map for TX
##########################

state_prob_df_for_map = data.frame(state = character(),
                                   county = character(),
                                   county_fips = character(),
                                   prob_date_1 = double(),
                                   # prob_date_2 = double(),
                                   stringsAsFactors=FALSE)

state_name = "Texas"
state_name_short = "TX"

print(state_name)

state_population = us_death%>%
  filter( Province_State == state_name)%>%
  select(Admin2, Population)

file_path_rdata = paste0(file_path_cur, "/", state_name_short ,"_results_in_GP.RData")

skip_to_next <- FALSE

load(file_path_rdata)

# load(file_path_rdata)

sum_idx = 1

for(county_index in 1:length(county_names)){
  
  county_name_selected = county_names[county_index]
  
  county_fip_each = as.character(us_death %>%
                                   filter( Province_State == state_name, Admin2 == county_name_selected)%>%
                                   select(county_fips))
  
  results_in_selected_county = results_overall_list[[county_index]]
  
  n = sum(results_in_selected_county$prediction_indicator==0)
  
  N = as.numeric(state_population$Population[which(state_population$Admin2==county_name_selected)])
  
  start_date = results_in_selected_county$date[1]
  
  beta_t_daily = results_in_selected_county$Beta_estimated[1:(n-1)]
  
  I_0 = results_in_selected_county$infective_f[1]
  R_0 = results_in_selected_county$resolved_f[1]
  
  death_selected = results_in_selected_county$death_real[1:n]
  
  S_t_seq = results_in_selected_county$suspective_f[1:n]
  
  gamma_updated = 1/4.5
  
  approx_beta_seq_all = beta_t_daily
  
  init_for_beta = c(S_t_seq[1], I_0, R_0, death_selected[1], 0)
  
  param_record_approx_all = matrix(NA, 5, n)
  param_record_approx_all[,1] = init_for_beta
  
  for (i in 1:(n-1)){
    S_t_1 = param_record_approx_all[1,i]
    I_t_1 = param_record_approx_all[2,i]
    R_t_1 = param_record_approx_all[3,i]
    D_t_1 = param_record_approx_all[4,i]
    C_t_1 = param_record_approx_all[5,i]
    
    beta_t_1_2 = approx_beta_seq_all[i]
    
    if(I_t_1<1){
      I_t_1 = 1
    }
    
    S_t_2 = uniroot(find_root_S_t_2, c(0, N), tol = 0.0001, param = c(S_t_1, beta_t_1_2, I_t_1), N = N, gamma=gamma_updated)
    I_t_2 = I_t_1 * exp(beta_t_1_2*(S_t_1 + S_t_2$root)/(2*N) - gamma_updated)
    R_t_2 = (2-theta)/(2+theta)*R_t_1 + gamma_updated/(2+theta)*(I_t_1+I_t_2)
    D_t_2 = D_t_1 + delta*theta*(R_t_1+R_t_2)/2
    C_t_2 = C_t_1 + (1-delta)*theta*(R_t_1+R_t_2)/2
    
    param_record_approx_all[, i+1] = c(S_t_2$root, I_t_2, R_t_2, D_t_2, C_t_2)
    
    
  }
  
  
  R_eff = beta_t_daily/gamma_updated * param_record_approx_all[1, 2:(n)] / N
  
  possibility_infected = param_record_approx_all[2, 2:(n)] * R_eff * gamma_updated / (param_record_approx_all[1, 2:(n)])
  
  possibility_infected_7_day_averaged = data_seven_day_smoothing(possibility_infected) #rollapply(possibility_infected, width = 7, by = 1, FUN = mean, align = "right")
  
  corresponding_date = results_in_selected_county$date[(2):(n)]
  
  probs_vector = c()
  
  for(i in 1:length(dates_extract_prob)){
    if(dates_extract_prob[i] < corresponding_date[1]){
      probs_vector = c(probs_vector, NA)
    } else{
      probs_vector = c(probs_vector, possibility_infected_7_day_averaged[which(corresponding_date==dates_extract_prob[i])])
    }
  }
  
  results_each = c(state_name, county_name_selected, county_fip_each, probs_vector)
  
  state_prob_df_for_map[sum_idx, ] = results_each
  
  sum_idx = sum_idx + 1
  
}

state_prob_df_for_map$prob_date_1 = as.numeric(state_prob_df_for_map$prob_date_1) * 100

state_prob_df_for_map$prob_date_1[state_prob_df_for_map$prob_date_1<10^(-4)] = 10^(-4)

state_prob_df_for_map$prob_date_1_scaled = log(state_prob_df_for_map$prob_date_1,2 )

range(state_prob_df_for_map$prob_date_1_scaled, na.rm=T)

state_counties = counties %>%
  filter(state_abbv==state_name_short)

paleta <- c("#fefec8", "#fdd88a", "#fdae76", "#f47562", "#d04c67")

state_prob_df_for_map$prob_date_1_f <- cut(state_prob_df_for_map$prob_date_1, breaks = c(0, 0.001, 0.01, 0.1, 1, 100))
niveles <- levels(state_prob_df_for_map$prob_date_1_f)

legend_text = c(expression(""<=10^{-3}),
                expression(paste("(",10^{-3}, ",", 10^{-2}, "]", sep="")),
                expression(paste("(",10^{-2}, ",", 10^{-1}, "]", sep="")),
                expression(paste("(",10^{-1}, ",", 1, "]", sep="")),
                expression("">1))

state_county_infect_prob_joint_TX <- left_join(state_counties,state_prob_df_for_map,  by = "county_fips")

TX_PoC_map = state_county_infect_prob_joint_TX %>%
  ggplot(aes(long, lat, group = group)) +
  #scale_fill_gradient2(midpoint=mid, low="#fee8c8", mid="#fdbb84",
  #                     high="red", space ="Lab" ) +
  geom_polygon(aes(fill = (prob_date_1_f))) +
  geom_polygon(color = "black", fill = NA) +
  scale_fill_manual(name = "PoC (%)", values = paleta, breaks = niveles, labels = legend_text)+
  coord_map() +
  labs(fill = expression("Odds (%)")) +
  xlab("Longitude") + ylab("Latitude") + # for the y axis label
  theme(
    legend.position = "right",
    text = element_text(size = 15),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 15),
    legend.key.size = unit(1, "cm"),
    axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    axis.title.y=element_blank(),
    axis.text.y =element_blank(),
    axis.ticks.y=element_blank()
  ) +
  guides(shape = guide_legend(override.aes = list(size = 1)),
         color = guide_legend(override.aes = list(size = 1))) 

###############################
# figures in TX
###############################

state_name  = "Texas"
state_name_short = "TX"

state_population = us_death%>%
  filter( Province_State == state_name)%>%
  select(Admin2, Population)

file_path_rdata = paste0(file_path_cur, "/", state_name_short ,"_results_in_GP.RData")

load(file_path_rdata)


# display the 10 best and 10 worest counties in WA
possibility_infected_record_matrix = matrix(NA, length(county_names), as.numeric(end_date - start_date_ini)+1)

standard_time_seq = seq.Date(start_date_ini, end_date, by=1)

for(county_index in 1:length(county_names)){
  
  county_name_selected = county_names[county_index]
  
  county_results_each = results_overall_list[[county_index]]
  
  # n = sum(county_results_each$prediction_indicator==0)
  
  N = as.numeric(state_population$Population[which(state_population$Admin2==county_name_selected)])
  
  start_date = county_results_each$date[1]
  
  # end_date = as.Date("2020-09-20")
  
  n = as.numeric(end_date - start_date+1)
  
  beta_t_daily = county_results_each$Beta_estimated[1:(n-1)]
  
  R_eff = beta_t_daily/0.2 * county_results_each$suspective_f[2:n] / N
  
  possibility_infected = county_results_each$infective_f[2:n] * R_eff / (county_results_each$suspective_f[2:n] * 5)
  
  possibility_infected_7_day_averaged = data_seven_day_smoothing(possibility_infected)#rollapply(possibility_infected, width = 7, by = 1, FUN = mean, align = "right")
  
  time_seq_each = county_results_each$date[(2):(n)]
  
  possibility_infected_record_matrix[county_index, standard_time_seq %in% time_seq_each] = possibility_infected_7_day_averaged
  
}

## top 5 counties

# top 10 counties with biggest probabilities getting infected
top_5_worst_county_names = county_names[order(possibility_infected_record_matrix[,dim(possibility_infected_record_matrix)[2]], decreasing = T)][1:5]

# display the 10 best and 10 worest counties in CA
possibility_infected_record_matrix = matrix(NA, length(top_5_worst_county_names), as.numeric(end_date - start_date_ini)+1)

standard_time_seq = seq.Date(start_date_ini, end_date, by=1)

for(county_index in 1:length(top_5_worst_county_names)){
  
  county_name_selected = top_5_worst_county_names[county_index]
  
  county_results_each = results_overall_list[[which(county_names == county_name_selected)]]
  
  # n = sum(county_results_each$prediction_indicator==0)
  
  N = as.numeric(state_population$Population[which(state_population$Admin2==county_name_selected)])
  
  start_date = county_results_each$date[1]
  
  # end_date = as.Date("2020-09-20")
  
  n = as.numeric(end_date - start_date+1)
  
  gamma_updated = 1/4.5
  
  beta_t_daily = county_results_each$Beta_estimated[1:(n-1)]
  
  I_0 = county_results_each$infective_f[1]
  R_0 = county_results_each$resolved_f[1]
  
  death_selected = county_results_each$death_real[1:n]
  
  S_t_seq = county_results_each$suspective_f[1:n]
  
  approx_beta_seq_all = beta_t_daily
  
  init_for_beta = c(S_t_seq[1], I_0, R_0, death_selected[1], 0)
  
  param_record_approx_all = matrix(NA, 5, n)
  param_record_approx_all[,1] = init_for_beta
  
  for (i in 1:(n-1)){
    S_t_1 = param_record_approx_all[1,i]
    I_t_1 = param_record_approx_all[2,i]
    R_t_1 = param_record_approx_all[3,i]
    D_t_1 = param_record_approx_all[4,i]
    C_t_1 = param_record_approx_all[5,i]
    
    beta_t_1_2 = approx_beta_seq_all[i]
    
    if(I_t_1<1){
      I_t_1 = 1
    }
    
    S_t_2 = uniroot(find_root_S_t_2, c(0, N), tol = 0.0001, param = c(S_t_1, beta_t_1_2, I_t_1), N = N, gamma=gamma_updated)
    I_t_2 = I_t_1 * exp(beta_t_1_2*(S_t_1 + S_t_2$root)/(2*N) - gamma_updated)
    R_t_2 = (2-theta)/(2+theta)*R_t_1 + gamma_updated/(2+theta)*(I_t_1+I_t_2)
    D_t_2 = D_t_1 + delta*theta*(R_t_1+R_t_2)/2
    C_t_2 = C_t_1 + (1-delta)*theta*(R_t_1+R_t_2)/2
    
    param_record_approx_all[, i+1] = c(S_t_2$root, I_t_2, R_t_2, D_t_2, C_t_2)
    
    
  }
  
  
  R_eff = beta_t_daily/gamma_updated * param_record_approx_all[1, 2:(n)] / N
  
  possibility_infected = param_record_approx_all[2, 2:(n)] * R_eff * gamma_updated / (param_record_approx_all[1, 2:(n)])
  
  possibility_infected_7_day_averaged = data_seven_day_smoothing(possibility_infected) #rollapply(possibility_infected, width = 7, by = 1, FUN = mean, align = "right")
  
  time_seq_each = county_results_each$date[(2):(n)]
  
  possibility_infected_record_matrix[county_index, standard_time_seq %in% time_seq_each] = possibility_infected_7_day_averaged
  
}


possibility_infected_top_5_worst = possibility_infected_record_matrix


infection_prob_sep_20 = possibility_infected_top_5_worst[, dim(possibility_infected_top_5_worst)[2]] * 100
names(infection_prob_sep_20) = as.character(top_5_worst_county_names)

top_5_worst_df = data.frame(date = rep(standard_time_seq, 5), County = rep(top_5_worst_county_names, each=length(standard_time_seq)), possibility_infected = as.vector(t(possibility_infected_top_5_worst)))

range(top_5_worst_df$possibility_infected*100, na.rm=T)

if(state_name  == "Texas"){
  y_upper_lim = 2.655575
} else if(state_name == "Washington"){
  y_upper_lim = 2.519254e-01
}

# this negative number is due to precision accuracy
top_5_worst_df$possibility_infected[top_5_worst_df$possibility_infected<0] = 0

TX_infection_curve = ggplot(top_5_worst_df, aes(x = date, y = possibility_infected*100,group = County)) + 
  geom_line(aes(color = County, linetype = County), lwd=1, show.legend = F)+
  scale_linetype_manual(values=rep("solid", 5)) +
  xlab("Date") +
  ylab("Probability of contracting COVID-19 (%)")+
  ylim(0,  y_upper_lim) + 
  theme(text = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15),
        axis.text.y = element_text(angle=90, hjust=1)) 

# plot the observced and fitted death of top 10 counties in CA

n_standard = length(standard_time_seq)

death_record_df = data.frame(matrix(NA, 10, n_standard+2))

colnames(death_record_df)[1:2] = c("County", "Type")
death_record_df$County = rep(as.character(top_5_worst_county_names), each=2)
death_record_df$Type = rep(c("Real", "Fitted"), 5)

population_record_vector = NULL

for(county_index in 1:5){
  county_name_selected = as.character(top_5_worst_county_names)[county_index]
  
  county_results_each = results_overall_list[[which(county_names == county_name_selected)]]
  
  N = as.numeric(state_population$Population[which(state_population$Admin2==county_name_selected)])
  
  population_record_vector = c(population_record_vector, N)
  
  start_date = county_results_each$date[1]
  
  # end_date = as.Date("2020-09-20")
  
  n = as.numeric(end_date - start_date+1)
  
  print(paste0("Population = ", N, " in ", county_name_selected, ", ", state_name))
  
  county_death_raw = us_death %>%
    filter(Admin2 == county_name_selected, Province_State == state_name) %>%
    select(starts_with("x"))
  
  county_death_raw_selected = as.numeric(county_death_raw[all_dates %in% standard_time_seq])
  
  for(i in 2:length(county_death_raw_selected)){
    if (county_death_raw_selected[i] < county_death_raw_selected[i-1]){
      county_death_raw_selected[i] = county_death_raw_selected[i-1]
    }
  }
  
  death_record_df[2*(county_index-1) + 1, 3:dim(death_record_df)[2]] = county_death_raw_selected / N * 10^6
  
  gamma_updated = 1/4.75
  
  beta_t_daily = county_results_each$Beta_estimated[1:(n-1)]
  
  I_0 = county_results_each$infective_f[1]
  R_0 = county_results_each$resolved_f[1]
  
  death_selected = county_results_each$death_real[1:n]
  
  S_t_seq = county_results_each$suspective_f[1:n]
  
  approx_beta_seq_all = beta_t_daily
  
  init_for_beta = c(S_t_seq[1], I_0, R_0, death_selected[1], 0)
  
  param_record_approx_all = matrix(NA, 5, n)
  param_record_approx_all[,1] = init_for_beta
  
  for (i in 1:(n-1)){
    S_t_1 = param_record_approx_all[1,i]
    I_t_1 = param_record_approx_all[2,i]
    R_t_1 = param_record_approx_all[3,i]
    D_t_1 = param_record_approx_all[4,i]
    C_t_1 = param_record_approx_all[5,i]
    
    beta_t_1_2 = approx_beta_seq_all[i]
    
    if(I_t_1<1){
      I_t_1 = 1
    }
    
    S_t_2 = uniroot(find_root_S_t_2, c(0, N), tol = 0.0001, param = c(S_t_1, beta_t_1_2, I_t_1), N = N, gamma=gamma_updated)
    I_t_2 = I_t_1 * exp(beta_t_1_2*(S_t_1 + S_t_2$root)/(2*N) - gamma_updated)
    R_t_2 = (2-theta)/(2+theta)*R_t_1 + gamma_updated/(2+theta)*(I_t_1+I_t_2)
    D_t_2 = D_t_1 + delta*theta*(R_t_1+R_t_2)/2
    C_t_2 = C_t_1 + (1-delta)*theta*(R_t_1+R_t_2)/2
    
    param_record_approx_all[, i+1] = c(S_t_2$root, I_t_2, R_t_2, D_t_2, C_t_2)
    
    
  }
  
  death_fitted = param_record_approx_all[4,]
  
  if (length(death_fitted) == n_standard){
    death_record_df[2*(county_index-1) + 2, 3:dim(death_record_df)[2]] = death_fitted / N * 10^6
  } else{
    diff_length = n_standard - length(death_fitted)
    
    death_record_df[2*(county_index-1) + 2, (3+diff_length):dim(death_record_df)[2]] = death_fitted / N * 10^6
  }

}

county_name_with_population = paste0(top_5_worst_county_names, " (", round(population_record_vector/10^3, 0), "K)")
death_record_df[,1] = rep(county_name_with_population, each = 2)

death_record_df_long = gather(death_record_df, key="Date", value = "Death", -County, -Type)
death_record_df_long$Date = rep(standard_time_seq, each=10)

death_record_df_long$group = paste0(death_record_df_long$County, death_record_df_long$Type)

death_record_df_long$Type = as.factor(death_record_df_long$Type)

range(death_record_df_long$Death, na.rm = T)

if(state_name  == "Texas"){
  y_upper_lim_death = 3165.126
} else if(state_name == "Washington"){
  y_upper_lim_death = 1034.665
}

TX_death_curve = ggplot(death_record_df_long, aes(x = Date, y = Death,group = group)) + 
  geom_line(aes(color = County, linetype = Type, size = Type), show.legend = T)+
  scale_linetype_manual(values=c("solid", "dotted"))+
  scale_size_manual(values = c(1,1.1))+
  xlab("Date") +
  ylab("Cumulative death per million")+
  guides(colour = guide_legend(order = 1,title="County (Population)"), 
         linetype = guide_legend(order = 2),
         size = FALSE)+
  ylim(0,  y_upper_lim_death)  + 
  theme(text = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.key.width=unit(1,"cm"),
        axis.text.y = element_text(angle=90, hjust=1))  

###############################
# Arrange figures and save it as EPS file
###############################

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

TX_legend<-g_legend(TX_death_curve)

WA_legend<-g_legend(WA_death_curve)

col_1 = ggarrange(TX_PoC_map, WA_PoC_map, ncol = 1, nrow = 2, common.legend = T, legend = "right")

row_1 <- grid.arrange(arrangeGrob(WA_infection_curve + theme(legend.position="none"),
                                  WA_death_curve + theme(legend.position="none"),
                                  nrow=1),
                      WA_legend, nrow=1,ncol=2, widths=c(10, 4))

row_2 <- grid.arrange(arrangeGrob(TX_infection_curve + theme(legend.position="none"),
                                  TX_death_curve + theme(legend.position="none"),
                                  nrow=1),
                      TX_legend, nrow=1,ncol=2, widths=c(10, 4))

save_file_path = "./Reproducing results in the paper/Results/Extended_fig_4/"

ggarrange( col_1, ggarrange(row_1, row_2, nrow = 2, ncol=1), nrow = 1, ncol = 2, widths = c(1,2) )%>%
  ggsave(file=paste0(save_file_path, "Results_in_WA_TX_infection_period_4_5.eps"), width = 36, height = 18, units = "cm")
