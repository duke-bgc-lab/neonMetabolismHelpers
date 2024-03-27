#' Determine the best quality estimates of NEON metabolism Bayesian model run
#' 
#' @author Nick Marzolf, \email{nick.marzolf@@jonesctr.org}
#' @author Wes Slaughter, \email{wslaughter@@berkeley.edu}
#' 
#' @param dir directory location of Bayesian model outputs
#' 
#' @details
#' This function reads in the outputs and passes the data through the diagnostic files output from nmh_diagnose_neon_metab_bayes(). Further, summary
#' statistics from the input time-series file are appended to the output data (e.g. discharge, DO, temperature) and added to a new directory#' 
#' 
#' @examples
#' nmh_estimate_metab_bayes(dir = 'data/model_runs/Bayes/')
#'
#' @export

nmh_estimate_metab_bayes <- function(dir = 'data/model_runs/Bayes/') {

  dir_daily <- glue::glue(dir, 'daily/')
  
  mods <- list.files(dir_daily,
                     full.names = TRUE)
  
  n_mods <- length(list.files(dir_daily,
                              full.names = TRUE))
  
  diagnostic <- readr::read_csv('data/model_runs/Bayes/diag_neon_metab_bayes.csv')
  
  good_years <- diagnostic %>% 
    dplyr::filter(ER_K_r2 < 0.6,
                  K_median < 100)
  
  for(i in 1:n_mods) {
    
    file <- mods[i]
    file_chunks <- unlist(stringr::str_split(file, '_'))
    
    # reconstruct how the data were compiled
    site <- file_chunks[2] %>% 
      substr(18,21)
    wyear <- file_chunks[3]
    q_type <- gsub('Q-','', file_chunks[4])
    z_method <- gsub('z-','',file_chunks[5])
    sensor_src <- gsub('.csv','', gsub('TS-','', file_chunks[6]))
    
    good_year <- diagnostic %>% 
      dplyr::filter(site %in% !!site,
                    wyear %in% !!wyear,
                    q_type %in% !!q_type,
                    z_method %in% !!z_method,
                    sensor_src %in% !!sensor_src)
    
    if(nrow(good_year) == 0)
      next
    
    # read in data
    d <- readr::read_csv(file)
    
    # pull out dates
    start_date <- dplyr::pull(d[1, 'date'])
    end_date <- dplyr::pull(d[length(d$date), 'date'])

    if('GPP_50pct' %in% names(d)) {
      
      d_clean <- d %>%
        dplyr::filter(GPP_50pct > -0.5,
                      ER_50pct < 0.5,
                      GPP_Rhat <= 1.2,
                      ER_Rhat <= 1.2,
                      K600_daily_Rhat <= 1.2) %>%
        dplyr::select(date,
                      GPP = GPP_50pct,
                      GPP.lower = GPP_2.5pct,
                      GPP.upper = GPP_97.5pct,
                      GPP.n_eff = GPP_n_eff,
                      GPP.Rhat = GPP_Rhat,
                      ER = ER_50pct,
                      ER.lower = ER_2.5pct,
                      ER.upper = ER_97.5pct,
                      ER.n_eff = ER_n_eff,
                      ER.Rhat = ER_Rhat,
                      K600 = K600_daily_50pct,
                      K600.lower = K600_daily_2.5pct,
                      K600.upper = K600_daily_97.5pct,
                      K600.n_eff = K600_daily_n_eff,
                      K600.Rhat = K600_daily_Rhat
        )
    }
    
    if(nrow(d_clean) == 0){
      next
    }
    
    dir_DO <- glue::glue(dir,'mod_and_obs_DO/')
    
    DO <- try(readr::read_csv(glue::glue(dir_DO, '/', site, '_', wyear, '_Q-', q_type, '_z-',z_method, '_TS-',sensor_src,'_mod_and_obs_DO.csv')))
    
    if(inherits(DO, 'try-error')){
      next
    }
    
    if(!'date' %in% names(DO)){
      DO <- DO %>%
        dplyr::mutate(date = lubridate::date(solar.time))
    }
    
    res <- glue::glue(diff(DO$solar.time) %>%
                        dplyr::first(),
                      ' min')
    
    DO_clean <- DO %>%
      dplyr::group_by(date) %>%
      dplyr::summarise(temp.water = mean(temp.water, na.rm = TRUE),
                       DO.obs = mean(DO.obs, na.rm = TRUE),
                       DO.sat = mean(DO.sat, na.rm = TRUE),
                       DO.min = min(DO.obs, na.rm = TRUE),
                       DO.max = max(DO.obs, na.rm = TRUE),
                       depth = mean(depth, na.rm = TRUE),
                       discharge = mean(discharge, na.rm = TRUE),
                       shortwave = mean(light, na.rm = TRUE)) %>%
      dplyr::mutate(DO.amp = DO.max - DO.min,
                    DO.psat = (DO.obs/DO.sat)*100)
    
    merged <- dplyr::left_join(d_clean,
                               DO_clean,
                               'date') %>%
      dplyr::mutate(resolution = res,
                    site = site,
                    discharge_type = q_type,
                    mean_depth_type = z_method,
                    sensor_source = sensor_src) %>%
      dplyr::select(site, discharge_type, mean_depth_type, sensor_source,
                    resolution, date, 
                    GPP, GPP.lower, GPP.upper, GPP.n_eff, GPP.Rhat,
                    ER, ER.lower, ER.upper, ER.n_eff, ER.Rhat, 
                    K600, K600.lower, K600.upper, K600.n_eff, K600.Rhat, 
                    DO.obs, DO.sat, DO.amp, DO.psat, depth, temp.water, discharge, shortwave)
    
    write_dir <- glue::glue(dir, 'estimates/')
    
    if(!dir.exists(write_dir))
      dir.create(write_dir)
    
    readr::write_csv(merged,
                     glue::glue(write_dir, 'neon_metab_{site}_{wyear}_Q-{q_type}_Z-{z_method}_TS-{sensor_src}.csv'))
    
  } # end for loop
  
  estimates <- purrr::map_dfr(list.files(write_dir,
                                         full.names = TRUE),
                              readr::read_csv)
  
  readr::write_csv(estimates,
                   glue::glue(dir, 'NEON_metab_estimates.csv'))
  
  return(estimates)
  
} # end function
