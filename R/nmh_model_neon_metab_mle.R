#' @export
# default directory is the default NEON input files, or data that uses only data from NEON
nmh_model_neon_metab_mle <- function(input_dir = 'data/sm_input/',
                                     log = TRUE) {
  
  # start the log
  logcode <- 'nmh_model_mle:'
  if(log){
    cat(paste0(logcode,' Ready to model metabolism using Maximum likelihood estimation' ),
        file = 'nmh_log.txt', append = TRUE, sep = '\n')
  }
  
  # how many time-series are there to model (max = 27)
  input_ts <- list.files(input_dir,
                         full.names = TRUE)
  
  if(log){
    cat(paste0(logcode,'-- there are ',length(input_ts), ' time-series to model'),
        file = 'nmh_log.txt', append = TRUE, sep = '\n')
  }
  
  # convert to a matrix that we will populate with tracking data if the model does or does not fit
  neon_mle_results <- matrix(ncol = 8,
                             nrow = length(input_ts))
  colnames(neon_mle_results) <- c('Site', 'Water Year', 'q_type','z_method','sensor_src','Has Data','Fit Error', 'Fit Time')
  
  
  # define the model
  init_name <- streamMetabolizer::mm_name(
    type = 'mle',
    pool_K600 = 'none',
    err_obs_iid = TRUE,
    err_proc_iid = FALSE,
    ode_method = 'trapezoid',
    deficit_src = 'DO_mod',
    engine = 'nlm',
    GPP_fun = 'linlight',
    ER_fun = 'constant')
  
  # set all specs, including priors
  init_specs = streamMetabolizer::specs(model_name = init_name)
  
  if(log){
    cat(paste0(logcode,'-- model specifications are defined, starting for loop around each file now'),
        file = 'nmh_log.txt', append = TRUE, sep = '\n')
  }
  
  # loop around site-years
  for(j in 1:length(input_ts)){
    
    # index the site and year
    # models are run for each site year, so we will filter and keep track of each site-year
    file <- gsub(input_dir,'',input_ts[j])
    file_chunks <- unlist(stringr::str_split(file, '_'))
    
    # reconstruct how the data were compiled
    site <- file_chunks[1]
    wyear <- file_chunks[2]
    q_type <- gsub('Q-','', file_chunks[3])
    z_method <- gsub('Z-','',file_chunks[4])
    sensor_src <- gsub('.csv','', 
                       gsub('TS-','', file_chunks[5]))
    
    # populate the results matrix
    neon_mle_results[j,'Site'] <- site
    neon_mle_results[j, 'Water Year'] <- wyear
    neon_mle_results[j,'q_type'] <- q_type
    neon_mle_results[j,'z_method'] <- z_method
    neon_mle_results[j,'sensor_src'] <- sensor_src
    
    # read in the formatted input data
    input_dat <- try(readr::read_csv(input_ts[j]) %>%
                       select(-wateryear))
    
    if(inherits(input_dat, 'try-error')) {
      neon_mle_results[j,'Has Data'] <- FALSE
      
      if(log){
        cat(paste0(logcode, '-- No input time-series at ', site),
            append = TRUE, file = 'nmh_log.txt', sep = '\n')
      }
    } else {
      cat(paste0(logcode, '--  input time-series read in for ', site),
          append = TRUE, file = 'nmh_log.txt', sep = '\n')
    }
    
    # check if there is data to model
    if(sum(is.na(input_dat$DO.obs)) == nrow(input_dat)){
      neon_mle_results[j,'Has Data'] = FALSE
      cat(paste0(logcode, '-- ',site,'-' ,year, ' Has no DO data, jumping to next site-year'),
          append = TRUE, file = 'nmh_log.txt', sep = '\n')
      next
    } else{
      neon_mle_results[j, 'Has Data'] = TRUE
      cat(paste0(logcode, '-- ',site,'-' ,wyear, ' Has data'),
          append = TRUE, file = 'nmh_log.txt', sep = '\n')
    }
    
    # pull the earliest and latest date in the modeled run
    start_date <- dplyr::first(lubridate::date(input_dat$solar.time))
    end_date <- dplyr::last(lubridate::date(input_dat$solar.time))
    
    # tracker: assume all models will be fit successfully
    fit_err <- FALSE
    this_time <- Sys.time()
    
    if(log){
      cat(paste0(logcode, '-- Starting the MLE model run at ', this_time),
          append = TRUE, file = 'nmh_log.txt', sep = '\n')
    }
    
    # fault tolerance: catch an errors and dictate the response to the errors
    tryCatch(
      
      # tries to run the model with specs and data defined above
      {dat_metab = streamMetabolizer::metab(init_specs,
                                            data = input_dat)
      
      # extract the outputs from the model
      dat_fit = streamMetabolizer::get_fit(dat_metab)
      },
      
      # evaluate fault tolerance: did the model fit?
      error = function(e) {
        fit_err <- TRUE
      }
    )
    
    this_time <- Sys.time()
    # did the model fit or not?
    if(fit_err){
      neon_mle_results[j,'Fit Error'] <- 'Fit Error'
      cat(paste0(logcode, '-- MLE run failed at ', this_time),
          append = TRUE, file = 'nmh_log.txt', sep = '\n')
      next
    } else {
      neon_mle_results[j,'Fit Error'] <- 'Fit Successful'
      neon_mle_results[j, 'Fit Time'] <- round(streamMetabolizer::get_fitting_time(dat_metab)[3], 4)
      cat(paste0(logcode, '-- MLE run successful at ', this_time),
          append = TRUE, file = 'nmh_log.txt', sep = '\n')
    }
    
    # save model output
    # where are we saving data to
    write_dir = 'data/model_runs/'
    
    if(!dir.exists(write_dir)){
      dir.create(write_dir)
    }
    
    if(!dir.exists(paste0(write_dir,'MLE')))
      dir.create(paste0(write_dir, 'MLE'))
    if(!dir.exists(paste0(write_dir,'MLE/dat_fit')))
      dir.create(paste0(write_dir, 'MLE/dat_fit'))
    if(!dir.exists(paste0(write_dir,'MLE/specs')))
      dir.create(paste0(write_dir, 'MLE/specs'))
    if(!dir.exists(paste0(write_dir,'MLE/mod_obs_DO')))
      dir.create(paste0(write_dir, 'MLE/mod_obs_DO'))
    
    # write the fit data
    fit_prefix <- paste0(write_dir, 'MLE/dat_fit/', site, '_', wyear, '_Q-', q_type, '_z-', z_method, '_TS-', sensor_src, '_')
    readr::write_csv(dat_fit, paste0(fit_prefix, 'dat_fit.csv'))
    
    # write the specs
    specs_out <- data.frame(unlist(streamMetabolizer::get_specs(dat_metab)))
    specs_prefix <- paste0(write_dir, 'MLE/specs/', site, '_', wyear, '_Q-', q_type, '_z-', z_method, '_TS-', sensor_src, '_')
    readr::write_csv(specs_out, paste0(specs_prefix, 'specs.csv'))
    
    # write the output data
    data_out = streamMetabolizer::get_data(dat_metab)
    dat_prefix <- paste0(write_dir, 'MLE/mod_obs_DO/', site, '_', wyear,'_Q-', q_type, '_z-', z_method, '_TS-', sensor_src, '_')
    readr::write_csv(data_out, paste0(dat_prefix, 'mod_and_obs_DO.csv'))
    
    if(log){
      cat(paste0(logcode, '-- MLE data saved to data/model_runs/MLE'),
          append = TRUE, file = 'nmh_log.txt', sep = '\n')
    }
    
  } # end for loop
  
  # return the matrix summarizing the model fits
  readr::write_csv(data.frame(neon_mle_results),
                   glue::glue('data/model_runs/MLE/neon_MLE_results.csv'))
  
  return(neon_mle_results)
  
} # end function

