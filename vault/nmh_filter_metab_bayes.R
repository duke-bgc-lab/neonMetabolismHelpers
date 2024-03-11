

nmh_filter_metab_bayes <- function(mod_dir = 'data/model_runs/') {

  est_dir <- glue::glue(mod_dir, q_type, '/Bayes/')
  
  est <- readr::read_csv(glue::glue(est_dir, 'NEON_metab_{q_type}_estimates.csv'))
  diags <- readr::read_csv(glue::glue(est_dir, 'diag_neon_metab_{q_type}.csv'))
  
  good_site_years <- diags %>%
    dplyr::filter(n_days > 365*0.6,
                  err_obs_iid_sigma_Rhat < 1.2,
                  err_proc_iid_sigma_Rhat < 1.2,
                  K600_daily_sigma_Rhat < 1.2,
                  K_median < 100,
                  ER_K_r2 < 0.6,
                  neg_GPP < 60,
                  pos_ER < 60) %>%
    dplyr::mutate(site_year = paste(site, year, sep = '-')) %>%
    dplyr::pull(site_year)
  
  good_est <- est %>%
    dplyr::mutate(year = lubridate::year(date),
                  site_year = paste(site, year, sep = '-')) %>%
    dplyr::filter(site_year %in% good_site_years) %>%
    dplyr::select(-site_year)
  
  write_dir <- 'data/filtered_estimates'
  
  if(!dir.exists(write_dir))
    dir.create(write_dir)
  
  readr::write_csv(good_est,
                   glue::glue('data/filtered_estimates/neon_metab_filt_{q_type}.csv'))
  
  return(good_est)
  
} # end function
