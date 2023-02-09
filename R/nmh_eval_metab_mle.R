nmh_eval_metab_mle <- function(dir = 'data/model_runs/MLE/dat_fit/',
                               log = TRUE) {

  fits <- list.files(dir,
                     full.names = TRUE)

  logcode = 'nmh_mle_eval:'
  this_time <- Sys.time()
  if(log){
    cat(paste0(logcode,'-- Evaluating ', length(fits),' MLE model runs to determine priors for Bayesian run at ', this_time),
        file = 'nmh_log.txt', append = TRUE, sep = '\n')
  }

  mle_diagnostic <- data.frame(
    site = character(),
    wyear = character(),
    q_type = character(),
    z_method = character(),
    sensor_src = character(),
    GPP_neg = numeric(),
    ER_pos = numeric(),
    K_median = numeric(),
    K_sd = numeric(),
    no_dat_num = numeric(),
    no_dat_dates = character()
  )

  for(i in 1:length(fits)) {

    file <- gsub(dir,'',fits[i])
    file_chunks <- unlist(stringr::str_split(file, '_'))

    # reconstruct how the data were compiled
    site <- file_chunks[1]
    wyear <- file_chunks[2]
    q_type <- gsub('Q-','', file_chunks[3])
    z_method <- gsub('z-','',file_chunks[4])
    sensor_src <- gsub('.csv','', gsub('TS-','', file_chunks[5]))

    fit <- readr::read_csv(fits[i])

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
        wyear = wyear,
        q_type = q_type,
        z_method = z_method,
        sensor_src = sensor_src,
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
                   'data/model_runs/MLE/MLE_diagnostics.csv')

  this_time <- Sys.time()
  if(log){
    cat(paste0(logcode,'-- Finshed evaluating MLE model runs to determine priors for Bayesian run at ', this_time),
        file = 'nmh_log.txt', append = TRUE, sep = '\n')
  }

  return(mle_diagnostic)

} # end function
