nmh_estimate_metab_bayes <- function(mod_dir = 'data/model_runs/',
                                     type = c('raw', 'qaqc', 'simulated')) {
  
  if(!type %in% c('raw','qaqc', 'simulated')){
    stop('Error: please select a discharge input from:\n 1) "raw": Raw NEON input\n 2) "source": NEON discharge evaluated by Rhea et al. (accepted), or\n 3) "simulated": NEON discharge simulations by the Macrosheds project')
  }
  
  file_dir <- glue::glue(mod_dir, type, '/Bayes/')
  
  mod_dir <- glue::glue(file_dir, 'daily')
  
  n_mods <- length(list.files(mod_dir,
                              full.names = TRUE))
  
  for(i in 1:n_mods) {
    
    daily <- list.files(mod_dir,
                        full.names = TRUE)[i]
    
    site <- stringr::str_split(daily, c('_','/'))[[2]][6] %>%
      substr(.,1,4)
    
    start_date <- stringr::str_split(daily, '_')[[1]][3]
    
    end_date <- stringr::str_split(daily, '_')[[1]][4]
    
    year <- substr(start_date, 1,4)
    
    d <- readr::read_csv(daily)
    
    if('GPP_50pct' %in% names(d)) {
      
      d_clean <- d %>%
        dplyr::filter(GPP_50pct > -0.5,
                      ER_50pct < 0.5) %>%
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
    
    dir_DO <- glue::glue(file_dir, 'mod_and_obs_DO')
    
    DO <- try(readr::read_csv(glue::glue(dir_DO, '/', site, '_', start_date, '_', end_date, '_mod_and_obs_DO.csv')))
    
    if(inherits(DO, 'try-error')){
      next
    }
    
    if(!'date' %in% names(DO)){
      DO <- DO %>%
        dplyr::mutate(date = lubridate::date(solar.time))
    }
    
    res <- glue::glue(diff(DO$solar.time) %>%
                        dplyr::first(),
                      'min')
    
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
                    site = site) %>%
      dplyr::select(site, resolution, date, GPP, GPP.lower, GPP.upper, GPP.n_eff, GPP.Rhat,
                    ER, ER.lower, ER.upper, ER.n_eff, ER.Rhat, K600, K600.lower, K600.upper,
                    K600.n_eff, K600.Rhat, DO.obs, DO.sat, DO.amp, DO.psat, depth, temp.water,
                    discharge, shortwave)
    
    write_dir <- glue::glue(file_dir, 'estimates/')
    
    if(!dir.exists(write_dir))
      dir.create(write_dir)
    
    readr::write_csv(merged,
                     glue::glue(write_dir, 'neon_metab_{type}_{site}_{year}.csv'))
    
  } # end for loop
  
  estimates <- purrr::map_dfr(list.files(write_dir,
                                         full.names = TRUE),
                              readr::read_csv)
  
  readr::write_csv(estimates,
                   glue::glue(file_dir, 'NEON_metab_{type}_estimates.csv'))
  
  return(estimates)
  
} # end function