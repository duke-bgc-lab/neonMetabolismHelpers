nmh_eval_metab_mle <- function(dir = 'data/model_runs/',
                               q_type = c('raw', 'source', 'simulated')) {

  dir_type <- glue::glue(dir, q_type)
  
  dir_fits <- glue::glue(dir_type, '/MLE/dat_fit')
  
  fits <- list.files(dir_fits,
                     full.names = TRUE)
  
  mle_diagnostic <- data.frame(
    site = character(),
    year = character(),
    GPP_neg = numeric(),
    ER_pos = numeric(),
    K_median = numeric(),
    K_sd = numeric(),
    no_dat_num = numeric(),
    no_dat_dates = character()
  )
  
  for(i in 1:length(fits)) {
    
    fit <- readr::read_csv(fits[i])
    
    site <- stringr::str_split(fits[i], '/')[[1]][6] %>%
      substr(.,1,4)
    year <- stringr::str_split(fits[i], '/')[[1]][6] %>%
      substr(.,6,9)
    
    GPP_neg <- fit %>%
      dplyr::filter(GPP.daily < -0.5) %>%
      dplyr::summarise(n = n()) %>%
      dplyr::pull(n)
    
    ER_pos <- fit %>%
      dplyr::filter(ER.daily > 0.5) %>%
      dplyr::summarise(n = n()) %>%
      dplyr::pull(n)
    
    K_sum <- fit %>%
      dplyr::filter(K600.daily > 0) %>%
      dplyr::summarise(medianK = median(K600.daily, na.rm = TRUE),
                       sdK = sd(K600.daily, na.rm = TRUE))
    
    no_dat_num <- sum(is.na(fit$GPP.daily))
    
    no_dat_dates <- fit %>%
      dplyr::filter(is.na(GPP.daily)) %>%
      dplyr::pull(date)
    
    mle_diagnostic <- mle_diagnostic %>%
      dplyr::add_row(
        site = site,
        year = year,
        GPP_neg = GPP_neg,
        ER_pos = ER_pos,
        K_median = dplyr::pull(K_sum, medianK),
        K_sd = dplyr::pull(K_sum,sdK),
        no_dat_num = no_dat_num,
        no_dat_dates = stringr::str_c(as.character(no_dat_dates), 
                                      collapse = ' ')
      )
  } # end for loop
  
  readr::write_csv(mle_diagnostic,
                   glue::glue(dir_type, '/MLE/MLE_{q_type}_diagnostics.csv'))
  
  return(mle_diagnostic)
  
} # end function
