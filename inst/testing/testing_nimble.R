library(tidyverse)
devtools::load_all(".")
library(nimble)
library(ggplot2)

# source('~/research/projects/r_packages/agTrendNimble/R/fit_ssl_nimble.R')
# source('~/research/projects/r_packages/agTrendNimble/R/prep_for_nimble_2.R')
# source('~/research/projects/r_packages/agTrendNimble/R/util_funcs.R')
# source('~/research/projects/r_packages/agTrendNimble/R/agg_funcs.R')

load("~/research/projects/sea_lion_analysis/SSL survey power analysis_2019/ssl_counts_v8.RData")
# source("helper_gam.R")

np <- edpsnp %>% #filter(REGION=="SE AK") %>% 
  arrange(SITE, year) %>% left_join(edps_photo,  by = c("SITE", "REGION", "year")) 

# %>% 
# np <-   filter(np, SITE %in% c("BODEGA ROCK"))

x <- prep_for_nimble(np, timeframe=c(1979,2019), debug=F)

tictoc::tic()
smp <- fit_ssl_nimble(x, niter = 70000, nburnin=20000, thin=10, debug=F)
tictoc::toc()




### Plot all the sites
gr <- rep(1:ceiling(nrow(x$site_data)/4), each=4)[1:nrow(x$site_data)]
for(i in 1:max(gr)){
  s <- x$site_data$site[gr==i]
  idx <- (smp$summary$site%in%s) #& (x$fitting_data$year>=1989)
  df <- smp$summary[idx,]
  p <- ggplot(data=df) + 
    geom_ribbon(aes(x=year, ymin=y_pred_ciL, ymax=y_pred_ciU), alpha=0.2) +
    geom_point(aes(x=year, y=counts), alpha=0.2) +
    geom_path(aes(x=year, y=mu_est), alpha=1, color='darkred', lwd=1.2) +
    # geom_ribbon(aes(x=year, ymin=mu_ciL, ymax=mu_ciU), alpha=0.1, fill='darkred') +
    geom_point(aes(x=year, y=y_real_est), data=df %>% filter(!is.na(counts))) +
    facet_wrap(~site, ncol=2, nrow=2, scales="free_y")
  print(p)
}





### Some aggregations
N_region <- ag_abund(smp$abund, "region")
region_trends <- ag_trend(N_region, c(1999, 2019))


N_region_summary <- summary_agg(N_region)
ggplot(data=N_region_summary %>% filter(year>=1989)) +
  geom_path(aes(x=year, y=Estimate), color='darkred') +
  geom_ribbon(aes(x=year, ymin=CI_predict_lower, ymax=CI_predict_upper), fill='darkred', alpha=0.2) +
  geom_pointrange(aes(x=year, y=Estimate_real, ymin=CI_real_lower, ymax=CI_real_upper), data=N_region_summary %>% filter(survey==1, year>=1989)) +
  geom_path(aes(x=year, y=Est), color='blue', lwd=1.5, data=region_trends$fitted %>% filter(type=='predicted')) +
  facet_wrap(~region, ncol=2, scales='free_y')


smp$abund$total <- 'total'
N_total <- ag_abund(smp$abund, "total")
N_total_summary <- summary_agg(N_total)
total_trend <- ag_trend(N_total, c(1999,2019))
ggplot(data=N_total_summary) +
  geom_path(aes(x=year, y=Estimate), color='darkred') +
  geom_ribbon(aes(x=year, ymin=CI_predict_lower, ymax=CI_predict_upper), fill='darkred', alpha=0.2) +
  geom_pointrange(aes(x=year, y=Estimate_real, ymin=CI_real_lower, ymax=CI_real_upper), data=N_total_summary %>% filter(survey==1)) +
  geom_path(aes(x=year, y=Est), color='blue', lwd=1.5, data=total_trend$fitted %>% filter(type=='predicted'))
