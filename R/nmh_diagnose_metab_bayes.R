#' Run diagnostics on Bayesian model run
#' 
#' @param mod_dir directory of the Bayesian model outputs from streamMetabolizer
#' @param log logical. Should the tracker file populate as it completes the function
#' @returns data frame of diagnostic file 
#' @export

nmh_diagnose_neon_metab_bayes <- function(mod_dir = 'data/model_runs/Bayes/',
                                          log = TRUE) {
  
  mod_dir_fits <- glue::glue(mod_dir,'daily/')
  
  n_mods <- list.files(mod_dir_fits,
                       full.names = TRUE)
  
  out <- data.frame(site = character(),
                    wyear = character(),
                    q_type = character(),
                    z_method = character(),
                    sensor_src = character(),
                    n_days = numeric(),
                    f_days = numeric(),
                    resolution = character(),
                    K600_daily_sigma_Rhat = numeric(),
                    err_obs_iid_sigma_Rhat = numeric(),
                    err_proc_iid_sigma_Rhat = numeric(),
                    K_median = numeric(),
                    K_range = numeric(),
                    neg_GPP = numeric(),
                    pos_ER = numeric(),
                    ER_K_r2 = numeric())
  
  for(i in 1:length(n_mods)) {
    
    file <- gsub(mod_dir_fits,'',n_mods[i])
    file_chunks <- unlist(stringr::str_split(file, '_'))
    
    # reconstruct how the data were compiled
    site <- file_chunks[1]
    wyear <- file_chunks[2]
    q_type <- gsub('Q-','', file_chunks[3])
    z_method <- gsub('z-','',file_chunks[4])
    sensor_src <- gsub('.csv','', gsub('TS-','', file_chunks[5]))
    
    # read in data
    d <- readr::read_csv(glue::glue(mod_dir_fits, file))
    
    # pull out dates
    start_date <- dplyr::pull(d[1, 'date'])
    end_date <- dplyr::pull(d[length(d$date), 'date'])

    if('GPP_50pct' %in% names(d)) {
      days <- d %>%
        dplyr::summarise(n = length(GPP_50pct)) %>%
        dplyr::pull()
      
      GPP_neg <- d %>%
        dplyr::filter(GPP_50pct < -0.5) %>%
        dplyr::summarise(n = length(GPP_50pct),
                         GPP_neg_per = (n/days)*100) %>%
        dplyr::pull(GPP_neg_per)
      
      ER_pos <- d %>%
        dplyr::filter(ER_50pct > 0.5) %>%
        dplyr::summarise(n = length(ER_50pct),
                         ER_pos_per = (n/days)*100) %>%
        dplyr::pull(ER_pos_per)
      
      K <- d %>%
        dplyr::summarise(K_median = median(K600_daily_50pct, na.rm = TRUE),
                         K_min = min(K600_daily_50pct, na.rm = TRUE),
                         K_max = max(K600_daily_50pct, na.rm = TRUE)) %>%
        dplyr::mutate(K_range = K_max - K_min)
      
      ER_K_r2 <- summary(lm(data = d,
                            ER_50pct ~ K600_daily_50pct))$adj.r.squared
    } else  # end if statement
      next
    
    dir_DO <- glue::glue(mod_dir,'mod_and_obs_DO/')
    
    DO <- try(readr::read_csv(glue::glue(dir_DO, '/', site, '_', wyear, '_Q-', q_type, '_z-',z_method, '_TS-',sensor_src,'_mod_and_obs_DO.csv')))
    
    if(inherits(DO, 'try-error')){
      next
    } # end if statement
    
    res <- glue::glue(diff(DO$solar.time) %>%
                        dplyr::first(),
                      ' min')
    
    dir_KQ <- glue::glue(mod_dir,'/KQ_overall')
    
    KQ <- try(readr::read_csv(glue::glue(dir_KQ, '/', site, '_', wyear, '_Q-', q_type, '_z-',z_method, '_TS-',sensor_src, '_KQ_overall.csv')))
    
    if(inherits(KQ, 'try-error')){
      next
    } # end if statement
    
    K600_daily_sigma_Rhat <- KQ %>%
      dplyr::pull(K600_daily_sigma_Rhat)
    
    dir_overall <- glue::glue(mod_dir,'/overall')
    
    overall <- try(readr::read_csv(glue::glue(dir_overall, '/', site, '_', wyear, '_Q-', q_type, '_z-',z_method, '_TS-',sensor_src, '_overall.csv')))
    
    if(inherits(overall, 'try-error')){
      next
    } # end if statement
    
    err_proc_iid_sigma_Rhat <- overall %>%
      dplyr::pull(err_proc_iid_sigma_Rhat)
    
    err_obs_iid_sigma_Rhat <- overall %>%
      dplyr::pull(err_obs_iid_sigma_Rhat)
    
    out <- out %>%
      dplyr::add_row(site = site,
                     wyear = wyear,
                     q_type = q_type,
                     z_method = z_method,
                     sensor_src = sensor_src,
                     n_days = days,
                     f_days = (days/365)*100,
                     resolution = res,
                     K600_daily_sigma_Rhat = K600_daily_sigma_Rhat,
                     err_obs_iid_sigma_Rhat = err_obs_iid_sigma_Rhat,
                     err_proc_iid_sigma_Rhat = err_proc_iid_sigma_Rhat,
                     K_median = K$K_median,
                     K_range = K$K_range,
                     neg_GPP = GPP_neg,
                     pos_ER = ER_pos,
                     ER_K_r2 = ER_K_r2)
  } # end for loop
  
  readr::write_csv(out,
                   glue::glue(mod_dir,'/diag_neon_metab_bayes.csv'))
  
  return(out)
  
} # end function
