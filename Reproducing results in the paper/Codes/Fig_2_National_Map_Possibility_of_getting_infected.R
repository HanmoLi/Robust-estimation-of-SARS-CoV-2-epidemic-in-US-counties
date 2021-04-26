##########################
# Code for figure 2 in the main manuscript
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
library(ggpubr)

if (!require("urbnmapr", character.only = TRUE)){
  #make sure you installed the urbnmapr package
  devtools::install_github("UrbanInstitute/urbnmapr")
} else(
  library(urbnmapr)
)

source(file = "./Reproducing results in the paper/Codes/functions_8th_version.R")

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

us_pow_param = read.csv(paste0(file_pow_param_path,"Name_list_of_US_states_updated.csv"))
us_pow_param$state_name = as.character(us_pow_param$state_name)
state_names_all = us_pow_param$state_name


####################################
# figure 1: Map for Apr 20
####################################

dates_extract_prob = c(as.Date("2020-04-20"))

prob_df_for_map_Apr = data.frame(state = character(),
                             county = character(),
                             county_fips = character(),
                             prob_date_1 = double(),
                             # prob_date_2 = double(),
                             stringsAsFactors=FALSE)

sum_idx = 1

for(state_index in 1:dim(us_pow_param)[1]){
  state_name = us_pow_param$state_name[state_index]
  state_name_short = us_pow_param$state_name_short[state_index]
  
  print(state_name)
  
  state_population = us_death%>%
    filter( Province_State == state_name)%>%
    select(Admin2, Population)
  
  file_path_rdata = paste0(file_path_cur, "/", state_name_short ,"_results_in_GP.RData")
  
  skip_to_next <- FALSE
  
  tryCatch({load(file_path_rdata)},
           error = function(e) {skip_to_next <<- TRUE})
  if(skip_to_next) { next}
  
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
    
    R_eff = beta_t_daily/0.2 * results_in_selected_county$suspective_f[2:(n)] / N
    
    possibility_infected = results_in_selected_county$infective_f[2:(n)] * R_eff / (5*results_in_selected_county$suspective_f[2:(n)])
    
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
    
    results_each = c(state_name, as.character(county_name_selected), county_fip_each, probs_vector)

    prob_df_for_map_Apr[sum_idx, ] = results_each
    
    sum_idx = sum_idx + 1
    
  }
  
}

prob_df_for_map_Apr$prob_date_1 = as.numeric(prob_df_for_map_Apr$prob_date_1) * 100

prob_df_for_map_Apr$prob_date_1[prob_df_for_map_Apr$prob_date_1<10^(-4)] = 10^(-4)

paleta <- c("#fefec8", "#fdd88a", "#fdae76", "#f47562", "#d04c67")

prob_df_for_map_Apr$prob_date_1_f <- cut(prob_df_for_map_Apr$prob_date_1, breaks = c(0, 0.001, 0.01, 0.1, 1, 100))
niveles <- levels(prob_df_for_map_Apr$prob_date_1_f)

legend_text = c(expression(""<=10^{-3}),
                expression(paste("(",10^{-3}, ",", 10^{-2}, "]", sep="")),
                expression(paste("(",10^{-2}, ",", 10^{-1}, "]", sep="")),
                expression(paste("(",10^{-1}, ",", 1, "]", sep="")),
                expression("">1))

us_county_infect_prob_joint_Apr <- left_join(counties,prob_df_for_map_Apr,  by = "county_fips")

g_plot_april_with_legend = us_county_infect_prob_joint_Apr %>%
  ggplot(aes(long, lat, group = group)) +
  geom_polygon(aes(fill = (prob_date_1_f)), show.legend = T) +
  geom_polygon(
    data = urbnmapr::states, mapping = aes(x = long, y = lat, group = group),
    fill = NA, color = 'black', size = 0.5
  ) + # for the y axis label
  theme(
    legend.position = "right",
    text = element_text(size = 20),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    legend.key.size = unit(1, "cm"),
    plot.title = element_text(hjust = 0.5,size=20, face="bold"),
    axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    axis.title.y=element_blank(),
    axis.text.y =element_blank(),
    axis.ticks.y=element_blank()
  ) +
  scale_fill_manual(name = "PoC (%)", values = paleta, breaks = niveles, labels = legend_text)+
  coord_map() +
  labs(fill = expression("PoC (%)")) +
  guides(shape = guide_legend(override.aes = list( size =1)),
         color = guide_legend(override.aes = list(size = 1)))+
  ggtitle(format(dates_extract_prob, "%d %B %Y"))



####################################
# figure 2: Map for Sep 20
####################################


dates_extract_prob = c(as.Date("2020-09-20"))

prob_df_for_map_Sep = data.frame(state = character(),
                             county = character(),
                             county_fips = character(),
                             prob_date_1 = double(),
                             # prob_date_2 = double(),
                             stringsAsFactors=FALSE)

sum_idx = 1

for(state_index in 1:dim(us_pow_param)[1]){
  state_name = us_pow_param$state_name[state_index]
  state_name_short = us_pow_param$state_name_short[state_index]
  
  print(state_name)
  
  state_population = us_death%>%
    filter( Province_State == state_name)%>%
    select(Admin2, Population)
  
  file_path_rdata = paste0(file_path_cur, "/", state_name_short ,"_results_in_GP.RData")
  
  skip_to_next <- FALSE
  
  tryCatch({load(file_path_rdata)},
           error = function(e) {skip_to_next <<- TRUE})
  if(skip_to_next) { next}
  
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
    
    R_eff = beta_t_daily/0.2 * results_in_selected_county$suspective_f[2:(n)] / N
    
    possibility_infected = results_in_selected_county$infective_f[2:(n)] * R_eff / (5*results_in_selected_county$suspective_f[2:(n)])
    
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
    
    results_each = c(state_name, as.character(county_name_selected), county_fip_each, probs_vector)

    prob_df_for_map_Sep[sum_idx, ] = results_each
    
    sum_idx = sum_idx + 1
    
  }
  
}

prob_df_for_map_Sep$prob_date_1 = as.numeric(prob_df_for_map_Sep$prob_date_1) * 100

prob_df_for_map_Sep$prob_date_1[prob_df_for_map_Sep$prob_date_1<10^(-4)] = 10^(-4)

paleta <- c("#fefec8", "#fdd88a", "#fdae76", "#f47562", "#d04c67")

prob_df_for_map_Sep$prob_date_1_f <- cut(prob_df_for_map_Sep$prob_date_1, breaks = c(0, 0.001, 0.01, 0.1, 1, 100))
niveles <- levels(prob_df_for_map_Sep$prob_date_1_f)

us_county_infect_prob_joint_Sep <- left_join(counties,prob_df_for_map_Sep,  by = "county_fips")

g_plot_sep = us_county_infect_prob_joint_Sep %>%
  ggplot(aes(long, lat, group = group)) +
  geom_polygon(aes(fill = (prob_date_1_f)), show.legend = T) +
  geom_polygon(
    data = urbnmapr::states, mapping = aes(x = long, y = lat, group = group),
    fill = NA, color = 'black', size = 0.5
  ) + # for the y axis label
  theme(
    legend.position = "None",
    text = element_text(size = 20),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    legend.key.size = unit(1, "cm"),
    plot.title = element_text(hjust = 0.5,size=20, face = "bold"),
    axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    axis.title.y=element_blank(),
    axis.text.y =element_blank(),
    axis.ticks.y=element_blank()
  ) +
  scale_fill_manual(name = "UVI", values = paleta, breaks = niveles)+
  coord_map() +
  labs(fill = expression("Odds (%)")) +
  guides(shape = guide_legend(override.aes = list( size =1)),
         color = guide_legend(override.aes = list(size = 1)))+
  ggtitle(format(dates_extract_prob, "%d %B %Y"))

g_plot_sep_with_legend =us_county_infect_prob_joint_Sep %>%
  ggplot(aes(long, lat, group = group)) +
  geom_polygon(aes(fill = (prob_date_1_f)), show.legend = T) +
  geom_polygon(
    data = urbnmapr::states, mapping = aes(x = long, y = lat, group = group),
    fill = NA, color = 'black', size = 0.5
  ) + # for the y axis label
  theme(
    legend.position = "right",
    text = element_text(size = 20),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 15),
    legend.key.size = unit(1, "cm"),
    plot.title = element_text(hjust = 0.5,size=20, face = "bold"),
    axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    axis.title.y=element_blank(),
    axis.text.y =element_blank(),
    axis.ticks.y=element_blank()
  ) +
  scale_fill_manual(name = "PoC (%)", values = paleta, breaks = niveles, labels = legend_text)+
  coord_map() +
  guides(shape = guide_legend(override.aes = list( size =1)),
         color = guide_legend(override.aes = list(size = 1)))+
  ggtitle(format(dates_extract_prob, "%d %B %Y"))

###############################
# Arrange figures and save it as EPS file
###############################

save_file_path = "./Reproducing results in the paper/Results/Fig_2/"

ggarrange(g_plot_april_with_legend, g_plot_sep_with_legend,
          ncol=2, nrow=1,widths = c(5,5),  common.legend = TRUE, legend="right") %>%
   ggsave(file=paste0(save_file_path, "US_infection_map_in_APR_SEP.eps"), width = 30, height = 12, units = "cm")
