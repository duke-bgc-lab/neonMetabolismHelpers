nmh_model_metab_bayes <- function(input_dir = 'data/sm_input/',
                                  type = 'raw') {
  
  if(!type %in% c('raw','source', 'simulation')){
    stop('Error: please select a discharge input from:\n 1) "raw": Raw NEON input\n 2) "source": NEON discharge evaluated by Rhea et al. (accepted), or\n 3) "simulation": NEON discharge simulations by the Macrosheds project')
  }
  
  # TODO: source in helper functions
  source('R/nmh_internals.R')
  
  # Define the model name, which establishes core hyperparameters. These would only be
  # changed as a last resort if the model just refuses to converge. i.e., if
  # we have to settle for a simpler model
  bayes_name_new = streamMetabolizer::mm_name(
    type = 'bayes',
    pool_K600 = 'binned',
    err_obs_iid = TRUE,
    err_proc_iid = TRUE,
    ode_method = 'trapezoid',
    deficit_src = 'DO_mod',
    engine = 'stan')
  
  # Define fault tolerances
  has_data <- TRUE
  setup_err <- FALSE
  fit_err <- FALSE
  
  # print into console that we are ready to model!
  writeLines(glue::glue('bayes modelling environment setup complete, ',
                        'starting site-year model runs \n\n\n'))
  
  # run through site-years in parallel
  # this loop contains various steps to prepare modeling
  
  # run through sit years in parallel
  
  neon_bayes_results <- foreach::foreach(m = 1:length(site_years$site),
                                .combine = bind_rows,
                                .packages = c('dplyr',
                                              'glue',
                                              'lubridate')) %dopar% {
          print("_________________________________________________")
          print(type)
          print("_________________________________________________")
          
          # Step 0: define starting conditions and initiate tracker
          
          # Define which site and which year we are modeling on
          site_id <- site_years[m, 'site']
          run_year <- site_years[m, 'year']
          
          # Define where the model outputs will be saved
          write_dir = glue::glue('data/model_runs/{type}/Bayes')
          
          # and create the directory if it doesn't exist
          if(!dir.exists(write_dir))
            dir.create(write_dir)
          
          # define cores and times
          session_id = Sys.getpid()
          session_date = Sys.Date()
          session_time = paste(Sys.time(), Sys.timezone())
          
          # create a tracker text file for each site-year model run and save it to the write_dir
          tracker_fn = glue::glue("{site_id}_{run_year}__{session_id}_tracker.txt")
          tracker_fp = file.path(write_dir, "trackers", tracker_fn)
          
          # and create the directory for the tracker files if it doesn't exist
          if(!dir.exists(glue::glue(write_dir, '/trackers')))
            dir.create(glue::glue(write_dir, '/trackers'))
          
          # create the file and write the meta-data for the model run
          cat(glue::glue("Tracker Log, R Session ID: {session_id}\n    datetime {session_date} {session_time}\n"), file=tracker_fp, sep="\n")
          cat(glue::glue('site: {site_id} \n year: {run_year} \n discharge: {type}'), file=tracker_fp, sep = "\n", append=TRUE)
          
          # Step 1: read in data
          # create a matrix to track which site-years have a successful model run
          neon_bayes_row <- matrix('',
                                   ncol = 7,
                                   nrow = 1)
         
          # we will track each site-year and record if there: ,
          colnames(neon_bayes_row) = c('Site',
                                       'Year',
                                       'Has Data',       # is data (T/F)
                                       'Setup Error',    # a K-Q relationship can be defined (T/F)
                                       'Fit Error',      # there was an error fitting the model (T/F)
                                       'Fit Error Text', # record the text from a model error (text)
                                       'Fit Time')       # and if the model run, how long did it take (numeric)
          
          
          writeLines(glue::glue('\n\n\nrunning model for\n\n    site:{site_id}\n    year: {run_year}\n\n\n'))
          
          # keep track of our results by site and by year
          neon_bayes_row[1, 'Site'] <- site_id
          neon_bayes_row[1, 'Year'] <- run_year
          
          # read in the data
          input_dat <- tryCatch({
            
            # read in site model input data
            readr::read_csv(glue::glue(input_dir,'/{type}/{site_id}_{type}_smReady.csv'))
          },
          # fault tolerance if the file couldn't be read in; not, jump to the next site
          error = function(e) {
            writeLines(glue('ERROR: input data failed to read in for {site_id}'))
            cat('Step_1_Error',
                file=tracker_fp, sep = "\n", append=TRUE)
            next
          }
          )
          
          # if the data exists, record to the tracker file that it was successful
          cat('Step_1_Success',
              file = tracker_fp, sep = "\n", append=TRUE)
          
          # Step 2: filter data to the run year
          input_dat_sub <- input_dat %>%
            dplyr::filter(lubridate::year(solar.time) == run_year)
          
          # if there are no rows in the dataframe, there is no data from that site-year and jump to the next
          if(nrow(input_dat_sub) == 0){
            neon_bayes_row[1, 'Has Data'] <- FALSE
            cat('Step_2_Error',                            # write the success/fail to the tracker file
                file = tracker_fp, sep = "\n", append=TRUE)
            # next
          } else {
            neon_bayes_row[1, 'Has Data'] <- TRUE
            cat('Step_2_Success',                          # write the success/fail to the tracker file
                file = tracker_fp, sep = "\n", append=TRUE)
          }

          # Step 3: check that input data is good (not all NAs)
          min_ratio <- 1 # 1 = 100% NAs in a column of the dataframe
          for(i in 1:ncol(input_dat_sub[-1])){
            
            # define the column names that have data
            vars <- names(input_dat_sub[-1])
            
            # which column to examine first
            var <- vars[i]
            
            # how many observations of that variable
            length <- input_dat_sub[,var] %>%
              nrow()
            
            # how many NAs in that variable
            nas <- input_dat_sub[,var] %>%
              dplyr::filter(is.na(.)) %>%
              nrow()
            
            # what is the percentage of NAs in that column
            ratio <- (nas/length)
            
            # if the ratio is 100% NAs, record that in the tracker file and jump to the next
            if(ratio == min_ratio) {
              cat('Step_3_Error', file = tracker_fp, sep = "\n", append=TRUE)
              cat(glue::glue('---- ', '{site_id} has no {var} data, jumping to next site-year'),
                  file = tracker_fp, sep = "\n", append=TRUE)
              writeLines(glue::glue('{site_id} has no {var} data, jumping to next site-year'))
            } else { # or keep going
              cat(glue::glue('Step_3_{i}_Success'),
                  file = tracker_fp, sep = "\n", append=TRUE)
              cat(glue::glue('---- ', '{site_id} has {round(as.numeric(ratio)*100,2)}% NA of {var} data'),
                  file = tracker_fp, sep = "\n", append=TRUE)
            }
          } # end for loop
          
          # Step 4: define KQ priors and model specifications
          # pull the earliest and latest date in the modeled run for file naming purposes
          tryCatch({
            # define start and end dates from the filtered dataset
            start_date <- lubridate::date(dplyr::first(input_dat_sub$solar.time))
            end_date <- lubridate::date(dplyr::last(input_dat_sub$solar.time))
            
            # read in results from MLE model run
            mle_priors <- readr::read_csv(glue::glue('data/model_runs/{type}/MLE/MLE_{type}_diagnostics.csv'))
            
            # subset to pull gas exchange results
            median_K600_mle <- mle_priors %>%
              dplyr::filter(site == site_id) %>%
              dplyr::select(site, year, K_median)
            
            # create an object that sorts discharge in log space
            Q_by_day <- tapply(log(input_dat_sub$discharge),
                               substr(input_dat_sub$solar.time, 1, 10),
                               mean)
            
            # define each node of log discharge by 0.2 units
            K600_lnQ_node_centers <- seq(from = min(Q_by_day, na.rm = TRUE),
                                         to = max(Q_by_day, na.rm = TRUE),
                                         by = 0.2)
            
            # and how many nodes are there
            n_nodes <- length(K600_lnQ_node_centers)
            
            # define the model specs and other hyperparameters
            bayes_specs <- streamMetabolizer::specs(
              model_name = bayes_name_new,
              K600_lnQ_nodes_centers = K600_lnQ_node_centers,
              K600_lnQ_nodediffs_sdlog = 0.1,
              K600_lnQ_nodes_meanlog = rep(log(12), length(K600_lnQ_node_centers)),  # default; see ?specs
              K600_lnQ_nodes_sdlog = rep(0.1, n_nodes),
              K600_daily_sigma_sigma = (dplyr::filter(median_K600_mle,
                                                      year == run_year) %>%
                                          dplyr::pull(K_median))* 0.02,
              burnin_steps = 1000,
              saved_steps = 500)
          },
          # if we can't define the K-Q relationship or model specs, there model setup fails
          error = function(e){
            setup_err <- TRUE
          }
          )
          
          # record in the tracker if the model setup is successful, and if not, jump to the next run
          if(setup_err == TRUE) {
            neon_bayes_row[1, 'Setup Error'] <- TRUE
            cat('Step_4_Error', file = tracker_fp, sep = "\n", append=TRUE)
            next
          } else {
            neon_bayes_row[1, 'Setup Error'] <- FALSE
            cat('Step_4_Success', file = tracker_fp, sep = "\n", append=TRUE)
          }
          
          # Step 5: run the model
          tryCatch({
            # run the model
            dat_metab <- streamMetabolizer::metab(specs = bayes_specs,
                                                  data = input_dat_sub)
            # if the model returns a warning, fit error is TRUE and we move
            # to the next site-year
            model_error <- dat_metab@fit$errors
            
            # if an error exists, model_error object will contain some character object(s), of some length
            if(length(model_error) > 0) {
              fit_err <- TRUE
            } else {
              fit_err <- FALSE
              model_error <- NA
              
              fit_time <- streamMetabolizer::get_fitting_time(dat_metab)
              neon_bayes_row[1, 'Fit Time'] <- fit_time[3]
              dat_fit <- streamMetabolizer::get_fit(dat_metab)
            }
          },
          # if an error (of literal class 'error') this code will run instead:
          error = function(e){
            fit_err <- TRUE
            model_error <- "streamMetabolizer::metab() function ERROR"
          },
          # this will run at end of tryCatch no matter what
          finally = {
            neon_bayes_row[1, 'Fit Error'] <- fit_err
            neon_bayes_row[1, 'Fit Error Text'] <- model_error
          }
          )
          
          # logging of model run errors
          if(fit_err) {
            writeLines(glue::glue("fit error on\n site: {site_id}\n year: {run_year}"))
            cat('Step_5_fitError', file=tracker_fp, sep = "\n", append=TRUE)
            cat(glue::glue('---- {model_error}'), file=tracker_fp, sep = "\n", append=TRUE)
            next
          } else {
            cat('Step_5_fitSuccess', file=tracker_fp, sep = "\n", append=TRUE)
          }
          
          # Step 6: write model output
          # define fault tolerance
          write_err <- FALSE
          
          # create, if necessary, the directory structure of the output files
          tryCatch({
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
            
            # define file name prefix
            fn_prefix = paste0('/', site_id, '_', start_date, '_', end_date, '_')
            
            # create objects to save from the model outputs
            specs_out = data.frame(unlist(streamMetabolizer::get_specs(dat_metab)))
            daily_out = data.frame(streamMetabolizer::get_data_daily(dat_metab))
            data_out = data.frame(streamMetabolizer::get_data(dat_metab))
            
            # write the model output objects to the directory needed
            readr::write_csv(data.frame(dat_fit$daily),
                             paste0(write_dir, '/daily', fn_prefix, 'daily.csv'))
            readr::write_csv(data.frame(dat_fit$overall),
                             paste0(write_dir, '/overall',fn_prefix, 'overall.csv'))
            readr::write_csv(data.frame(dat_fit$KQ_overall),
                             paste0(write_dir,'/KQ_overall',fn_prefix, 'KQ_overall.csv'))
            readr::write_csv(specs_out,
                             paste0(write_dir,'/specs', fn_prefix, 'specs.csv'))
            readr::write_csv(daily_out,
                             paste0(write_dir,'/datadaily',fn_prefix, 'datadaily.csv'))
            readr::write_csv(data_out,
                             paste0(write_dir,'/mod_and_obs_DO',fn_prefix, 'mod_and_obs_DO.csv'))
          },
          # fault tolerance: if there was an error writing these files, flag the error
          error = function(e) {
            write_err <- TRUE
          }
          )
          
          # and record the error in the tracker file
          if(write_err) {
            cat('Step_6_Error', file=tracker_fp, sep = "\n", append=TRUE)
          } else
            cat('Step_6_Success', file=tracker_fp, sep = "\n", append=TRUE)
          
          # write the site-year tracker information
          return(neon_bayes_row)
                                                
      } # end site-year foreach()
  
  writeLines("---- NEON Metabolism Bayes Model Runs Complete -----")
  readr::write_csv(data.frame(neon_bayes_results),
                   glue::glue(write_dir, '/neon_bayes_row.csv')) #save run results if you want to
  
  # TODO: add a column(s) to site_years to look at general stats on model runs etc.
  
  return(neon_bayes_results)
  
} # end function