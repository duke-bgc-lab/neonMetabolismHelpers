nmh_model_metab_bayes <- function(input_dir = 'data/sm_input/') {
  
  # TODO: move STAN related things to a separate script to prepare for Bayesian modeling?
  # Bayesian modeling requires rstan; install if not present
  tryCatch({
    library('rstan')
    library('StanHeaders')
  },
  error = function(e){
    install.packages("StanHeaders", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
    install.packages("rstan", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
    library('rstan')
    library('StanHeaders')
  }
  )
  
  # source in helper functions
  source('src/helpers.R')
  source('src/neon_helpers.R')
  
  # define lotic sites we are running the model on
  site_code <- get_neon_site_data()[,'site_code']
  
  # create a data frame with all possible site-years for the model to run
  site_years <- data.frame(
    site = rep(site_code, each = length(2016:2021)),
    year = rep(2016:2021, times = length(site_code))
  )
  
  # get model name, which establishes core hyperparams. these would only be
  # changed as a last resort if the model just refuses to converge. i.e., if
  # we have to settle for a simpler model
  bayes_name_new = streamMetabolizer::mm_name(
    type='bayes',
    pool_K600='binned',
    err_obs_iid=TRUE,
    err_proc_iid=TRUE,
    ode_method='trapezoid',
    deficit_src='DO_mod',
    engine='stan')
  
  has_data <- TRUE
  setup_err <- FALSE
  fit_err <- FALSE
  
  neon_bayes_results <- foreach(m = 1:length(site_years$site), .combine = bind_rows) %dopar% {
    
    # create a matrix to track which site-years have a successful model run
    neon_bayes_row <- matrix('',
                             ncol = 6,
                             nrow = 1)
    
    colnames(neon_bayes_row) = c('Site', 'Year', 'Has Data',
                                 'Setup Error', 'Fit Error', 'Fit Time')
    
    site_id <- site_years[m, 'site']
    run_year <- site_years[m, 'year']
    
    # keep track of our results by site and by year
    neon_bayes_row[m, 'Site'] <- site_id
    neon_bayes_row[m, 'Year'] <- run_year
    
    input_dat <- try(readr::read_csv(glue::glue(input_dir,'/{site_id}_smReady.csv')))
    
    input_dat_sub <- input_dat %>%
      dplyr::filter(lubridate::year(solar.time) == run_year)
    
    if(nrow(input_dat_sub) == 0){
      neon_bayes_row[m, 'Has Data'] <- 'No'
      next
    } else neon_bayes_row[m, 'Has Data'] <- 'Yes'
    
    # pull the earliest and latest date in the modeled run for file naming purposes
    start_date <- lubridate::date(dplyr::first(input_dat_sub$solar.time))
    end_date <- lubridate::date(dplyr::last(input_dat_sub$solar.time))
    
    mle_priors <- readr::read_csv('data/model_runs/source/MLE/MLE_diag.csv') # add in filepath here
    
    median_K600_mle <- mle_priors %>%
      dplyr::filter(site == site_id) %>% # edit this
      dplyr::select(site, year, K_median)
    
    Q_by_day <- tapply(log(input_dat_sub$discharge),
                       substr(input_dat_sub$solar.time, 1, 10),
                       mean)
    
    tryCatch({
      K600_lnQ_node_centers <- seq(from = min(Q_by_day, na.rm = TRUE),
                                   to = max(Q_by_day, na.rm = TRUE),
                                   by = 0.2)
      
      n_nodes <- length(K600_lnQ_node_centers)
      
      bayes_specs <- streamMetabolizer::specs(
        model_name = bayes_name_new,
        K600_lnQ_nodes_centers = K600_lnQ_node_centers,
        K600_lnQ_nodediffs_sdlog = 0.1,
        K600_lnQ_nodes_meanlog = rep(log(12), length(K600_lnQ_node_centers)),  # default; see ?specs
        K600_lnQ_nodes_sdlog = rep(0.1, n_nodes),
        #K600_daily_sigma_mean = 0,          # dont think this is even a settable parameter?
        K600_daily_sigma_sigma = (dplyr::filter(median_K600_mle,
                                                year == run_year)
                                  %>% dplyr::pull(K_median))* 0.02,
        burnin_steps = 1000,
        saved_steps = 500)
    },
    error = function(e){
      setup_err <- TRUE
    }
    )
    
    if(setup_err == TRUE) {
      neon_bayes_row[m, 'Setup Error'] <- TRUE
      next
    } else
      neon_bayes_row[m, 'Setup Error'] <- FALSE
    
    # run the model
    tryCatch({
      dat_metab <- streamMetabolizer::metab(specs = bayes_specs,
                                            data = input_dat_sub)
      dat_fit <- streamMetabolizer::get_fit(dat_metab)
      fit_time <- streamMetabolizer::get_fitting_time(dat_metab)
    },
    error = function(e){
      fit_err <- TRUE
    }
    )
    
    if(fit_err == TRUE) {
      neon_bayes_row[m, 'Fit Error'] <- TRUE
      next
    } else
      neon_bayes_row[m, 'Fit Error'] <- FALSE
    
    neon_bayes_row[m, 'Fit Time'] <- fit_time[3]
    
    # save model output
    # where are we saving data to
    write_dir = 'data/model_runs/source/Bayes'
    
    if(!dir.exists(write_dir))
      dir.create(write_dir)
    
    if(!dir.exists(glue::glue(write_dir, '/daily')))
      dir.create(glue::glue(write_dir, '/daily'))
    if(!dir.exists(glue::glue(write_dir, '/overall')))
      dir.create(glue::glue(write_dir, '/overall'))
    if(!dir.exists(glue::glue(write_dir, '/KQ_overall')))
      dir.create(glue::glue(write_dir, '/KQ_overall'))
    if(!dir.exists(glue::glue(write_dir, '/specs')))
      dir.create(glue::glue(write_dir, '/specs'))
    if(!dir.exists(glue::glue(write_dir, '/datadaily')))
      dir.create(glue::glue(write_dir, '/datadaily'))
    if(!dir.exists(glue::glue(write_dir, '/mod_and_obs_DO')))
      dir.create(glue::glue(write_dir, '/mod_and_obs_DO'))
    
    fn_prefix = paste0('/', site_id, '_', start_date, '_', end_date, '_')
    
    specs_out = data.frame(unlist(get_specs(dat_metab)))
    daily_out = streamMetabolizer::get_data_daily(dat_metab)
    data_out = streamMetabolizer::get_data(dat_metab)
    
    readr::write_csv(dat_fit$daily,
                     paste0(write_dir, '/daily', fn_prefix, 'daily.csv'))
    readr::write_csv(dat_fit$overall,
                     paste0(write_dir, '/overall',fn_prefix, 'overall.csv'))
    readr::write_csv(dat_fit$KQ_overall,
                     paste0(write_dir,'/KQ_overall',fn_prefix, 'KQ_overall.csv'))
    readr::write_csv(specs_out,
                     paste0(write_dir,'/specs', fn_prefix, 'specs.csv'))
    readr::write_csv(daily_out,
                     paste0(write_dir,'/datadaily',fn_prefix, 'datadaily.csv'))
    readr::write_csv(data_out,
                     paste0(write_dir,'/mod_and_obs_DO',fn_prefix, 'mod_and_obs_DO.csv'))
    
    return()
    
  } # end for loop
  readr::write_csv(data.frame(neon_bayes_row),
                   glue::glue(write_dir, '/neon_bayes_row.csv')) #save run results if you want to
  
} # end function