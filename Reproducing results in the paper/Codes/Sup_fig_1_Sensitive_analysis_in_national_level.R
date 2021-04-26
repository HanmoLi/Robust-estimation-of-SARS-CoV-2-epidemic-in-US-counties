##########################
# Code for figure 1 in the supplementary information
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
library(scales)
library(grid)
library(gridExtra)
library(ggpubr)

source(file = "./Reproducing results in the paper/Codes/functions_8th_version.R")

##for forecast
fitted_days_beta=184 ##how many days from the end day you use to estimate beta_t
  
predicted_days_train=184  
prediction_length = 90

beta_t_length=184

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

# get the national fittings

N = 321418820 # the us population

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
    select(starts_with("x"))
  
  state_confirmed_selected = state_confirmed[all_dates %in% standard_time_seq]
  
  US_confirmed_overall = US_confirmed_overall + as.numeric(apply(state_confirmed_selected, 2, sum))
  
  state_death = us_death %>%
    filter(Province_State == state_name)%>%
    select(starts_with("x"))
  
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

I_0_ini = 10000
R_0_ini = 10000

if((I_0_ini + R_0_ini)> abs(upper_final)){
  I_0_ini = abs(upper_final)/4
  R_0_ini = abs(upper_final)/4
}

alpha_vector = c(seq(0.1,0.5,0.05), seq(0.6, 1.9, 0.1))

error_record = rep(NA, length(alpha_vector))
m_approx_record = as.list( rep(NA, length(alpha_vector)))

gamma = 0.2
theta = 0.1
delta = 0.0066

gamma_vector = c(0.2,1/7, 0.2,0.2)
theta_vector = c(0.1, 0.1, 1/15, 0.1)
delta_vector = c(0.0066, 0.0066, 0.0066, 0.0075)

fix_param_df = data.frame(gamma = gamma_vector,
                          theta = theta_vector,
                          delta = delta_vector)

R_eff_matrix = matrix(NA, dim(fix_param_df)[1], n_standard-1)
Infectious_matrix = matrix(NA, dim(fix_param_df)[1], n_standard)
Death_matrix = matrix(NA, dim(fix_param_df)[1], n_standard)
PoC_matrix = matrix(NA, dim(fix_param_df)[1], n_standard-1)

for(param_index in 1:dim(fix_param_df)[1]){
  
  gamma = fix_param_df$gamma[param_index]
  theta = fix_param_df$theta[param_index]
  delta = fix_param_df$delta[param_index]
  
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
  
  R_eff = approx_beta_seq / gamma * param_record_approx_for_beta[1,2:n_standard] / N
  R_eff_smoothed = data_seven_day_smoothing(R_eff)
  
  Infectious_cases = param_record_approx_for_beta[2,]
  
  Death_toll =  param_record_approx_for_beta[4,]
  
  PoC = Infectious_cases[2:n_standard] * R_eff * gamma / param_record_approx_for_beta[1,2:n_standard]
  PoC_smoothed = data_seven_day_smoothing(PoC)
  
  R_eff_matrix[param_index, ] = R_eff_smoothed
  Infectious_matrix[param_index, ] = Infectious_cases
  Death_matrix[param_index, ] = Death_toll
  PoC_matrix[param_index, ] = PoC_smoothed
  
}

R_eff_df_long = data.frame(
  date = rep(standard_time_seq[2:n_standard], 4),
  type = as.factor(rep(
    c("Setting 1 (default)", "Setting 2", "Setting 3", "Setting 4"),
    each = n_standard - 1
  )),
  value = as.vector(t(R_eff_matrix))
)

Infectious_df_long = data.frame(
  date = rep(standard_time_seq, 4),
  type = as.factor(rep(
    c("Setting 1 (default)", "Setting 2", "Setting 3", "Setting 4"),
    each = n_standard
  )),
  value = as.vector(t(Infectious_matrix))
)

Death_df_long = data.frame(
  date = rep(standard_time_seq, 4),
  settings = rep(
    c(1, 2, 3, 4),
    each = n_standard
  ),
  type = "Estimation",
  value = as.vector(t(Death_matrix))
)

PoC_df_long = data.frame(
  date = rep(standard_time_seq[2:n_standard], 4),
  type = as.factor(rep(
    c("Setting 1 (default)", "Setting 2", "Setting 3", "Setting 4"),
    each = n_standard - 1
  )),
  value = as.vector(t(PoC_matrix))
)

Death_raw = data.frame(date = standard_time_seq,settings = 5, type="Observation", value = death_selected)

Death_df_fill = rbind(Death_df_long,Death_raw )

R_eff = R_eff_df_long %>%
  ggplot(aes(x = date, y = value, group = type, col = type)) +
  geom_line(aes(linetype = type), show.legend = F, size = 1)+
  ylab("Effective reproduction number") +
  scale_y_continuous(label=comma)+
  theme(text = element_text(size=15),
        axis.text.y = element_text(angle=90, hjust=1))

Infectious = Infectious_df_long %>%
  ggplot(aes(x = date, y = value / 10^6, group = type, col = type)) +
  geom_line(aes(linetype = type), show.legend = F, size = 1)+
  ylab(expression(paste('Active infectious individuals  ', (10^6)))) +
  scale_y_continuous(label=comma)+
  theme(text = element_text(size=15),
        axis.text.y = element_text(angle=90, hjust=1))

Death_legends = Death_df_long%>%
  ggplot(aes(
    x = date,
    y = value / 10 ^ 3,
    group = as.factor(settings),
    col = as.factor(settings)
  )) +
  geom_line(aes(linetype = as.factor(settings)), show.legend = T, size = 1) +
  ylab(expression(paste("Cumulative death toll  ", (10^3)))) +
  scale_y_continuous(label = comma) +
  guides(
    colour = guide_legend(title = "Configuration", override.aes=list(shape=NA)),
    linetype = guide_legend(title = "Configuration", override.aes=list(shape=NA)),
    size=guide_legend(title = NULL, override.aes=list(shape=NA, size = 0.8, linetype=3))
  )+
  theme(text = element_text(size=15),
        legend.key.width=unit(1,"cm"),
        axis.text.y = element_text(angle=90, hjust=1))+
  geom_line(data = Death_raw, aes(x = date, y = value / 10^3),show.legend = F,linetype = "dashed", size = 1.5, col = "black")

Death = Death_df_long %>%
  ggplot(aes(
    x = date,
    y = value / 10 ^ 3,
    group = as.factor(settings),
    col = as.factor(settings)
  )) +
  geom_line(aes(linetype = as.factor(settings)), show.legend = F, size = 1) +
  ylab(expression(paste("Cumulative death toll  ", (10^3)))) +
  scale_y_continuous(label = comma) +
  guides(
    colour = FALSE,
    linetype = FALSE
  )+
  theme(text = element_text(size=15),
        legend.key.width=unit(1,"cm"),
        axis.text.y = element_text(angle=90, hjust=1))+
  geom_line(data = Death_raw, aes(x = date, y = value / 10^3),show.legend = F, linetype = "dashed", size = 1.5, col = "black")

PoC = PoC_df_long %>%
  ggplot(aes(x = date, y = value * 100, group = type, col = type)) +
  geom_line(aes(linetype = type), show.legend = F, size = 1)+
  ylab("Probability of contracting COVID-19 (%)")+
  theme(text = element_text(size=15),
        axis.text.y = element_text(angle=90, hjust=1))

##########
# arrange the plots
##########

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

legends<-g_legend(Death_legends)

col_1_2 = ggarrange(R_eff, PoC,Infectious,Death,  ncol = 2, nrow = 2)

save_file_path = "./Reproducing results in the paper/Results/Sup_fig_1/"

ggarrange( col_1_2, legends, nrow = 1, ncol = 2, widths = c(10,2) )%>%
  ggsave(file=paste0(save_file_path, "National_level_sensitive_analysis.eps"), width = 36, height = 18, units = "cm")


