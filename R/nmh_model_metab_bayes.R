#' @export
nmh_model_metab_bayes <- function(input_dir = 'data/sm_input/',
                                  log = TRUE) {


  logcode <- 'nmh_model_bayes:'
  if(log){
    cat(paste0(logcode,' Ready to model metabolism using Bayesian model structure'),
        file = 'nmh_log.txt', append = TRUE, sep = '\n')
  }

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
  has_data <- FALSE
  setup_err <- FALSE
  fit_err <- FALSE

  # run through site-years in parallel
  # this loop contains various steps to prepare modeling

  input_ts <- list.files(input_dir,
                         full.names = TRUE)

  if(log){
    cat(paste0(logcode,' There are ', length(input_ts), ' site - wateryears to run'),
        file = 'nmh_log.txt', append = TRUE, sep = '\n')
  }

  # run through sit years in parallel
  neon_bayes_results <- foreach::foreach(m = 1:length(input_ts),
                                .combine = bind_rows,
                                .packages = 'dplyr') %dopar% {

          # Step 0: get information from each input time-series
          file <- gsub(input_dir,'',input_ts[m])
          file_chunks <- unlist(stringr::str_split(file, '_'))

          # reconstruct how the data were compiled
          site <- file_chunks[1]
          wyear <- file_chunks[2]
          q_type <- gsub('Q-','', file_chunks[3])
          z_method <- gsub('Z-','',file_chunks[4])
          sensor_src <- gsub('.csv','', gsub('TS-','', file_chunks[5]))

          # Define where the model outputs will be saved
          write_dir = glue::glue('data/model_runs/Bayes/')

          # and create the directory if it doesn't exist
          if(!dir.exists(write_dir))
            dir.create(write_dir)

          # define cores and times
          session_id = Sys.getpid()
          session_date = Sys.Date()
          session_time = paste(Sys.time(), Sys.timezone())

          # create a tracker text file for each site-year model run and save it to the write_dir
          tracker_fn = glue::glue("{site}_{wyear}_{q_type}_{z_method}_{sensor_src}_{session_id}_tracker.txt")
          tracker_fp = file.path(write_dir, "trackers", tracker_fn)

          # and create the directory for the tracker files if it doesn't exist
          if(!dir.exists(glue::glue(write_dir, '/trackers')))
            dir.create(glue::glue(write_dir, '/trackers'))

          # create the file and write the meta-data for the model run
          cat(glue::glue("Tracker Log, R Session ID: {session_id}\n    datetime {session_date} {session_time}\n"), file=tracker_fp, sep="\n")
          cat(glue::glue('site: {site} \n year: {wyear} \n discharge: {q_type}, depth method: {z_method}, sensor data from" {sensor_src}'),
              file=tracker_fp, sep = "\n", append=TRUE)

          # Step 1: read in data
          # create a matrix to track which site-years have a successful model run
          neon_bayes_row <- matrix(ncol = 10,
                                   nrow = 1)

          # we will track each site-year and record if there: ,
          colnames(neon_bayes_row) = c('Site',
                                       'WYear',
                                       'q_type',
                                       'z_method',
                                       'sensor_src',
                                       'Has Data',       # is data (T/F)
                                       'Setup Error',    # a K-Q relationship can be defined (T/F)
                                       'Fit Error',      # there was an error fitting the model (T/F)
                                       'Fit Error Text', # record the text from a model error (text)
                                       'Fit Time')       # and if the model run, how long did it take (numeric)


          writeLines(glue::glue('\n\n\nrunning model for\n\n    site:{site}\n    year: {wyear}\n\n\n'))

          # keep track of our results by site and by year
          neon_bayes_row[1, 'Site'] <- site
          neon_bayes_row[1, 'WYear'] <- wyear
          neon_bayes_row[1, 'q_type'] <- q_type
          neon_bayes_row[1, 'z_method'] <- z_method
          neon_bayes_row[1, 'sensor_src'] <- sensor_src

          # read in the data
          input_dat <- tryCatch({

            # read in site model input data
            readr::read_csv(input_ts[m]) %>%
              dplyr::select(-wateryear)
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
          if(log){
            cat(paste0(logcode,' Data read in successfully'),
                file = 'nmh_log.txt', append = TRUE, sep = '\n')
            cat('Step_1_Success',
              file = tracker_fp, sep = "\n", append=TRUE)
          }

          # check if there is data to model
          if(sum(is.na(input_dat$DO.obs)) == nrow(input_dat)){
            neon_bayes_row[1, 'Has Data'] <- FALSE
            cat('Step_2_Error',                            # write the success/fail to the tracker file
                file = tracker_fp, sep = "\n", append=TRUE)
            cat(paste0(logcode, '-- ',site,'-' ,wyear, ' Has no DO data'),
                append = TRUE, file = 'nmh_log.txt', sep = '\n')
          } else{
            neon_bayes_row[1, 'Has Data'] <- TRUE
            cat(paste0(logcode, '-- ',site,'-' , wyear, ' Has data'),
                append = TRUE, file = 'nmh_log.txt', sep = '\n')
            cat('Step_2_Success',                          # write the success/fail to the tracker file
                file = tracker_fp, sep = "\n", append=TRUE)
          }

          # Step 3: check that input data is good (not all NAs)
          min_ratio <- 1 # 1 = 100% NAs in a column of the dataframe
          for(i in 1:ncol(input_dat[-1])){

            # define the column names that have data
            vars <- names(input_dat[-1])

            # which column to examine first
            var <- vars[i]

            # how many observations of that variable
            length <- input_dat[,var] %>%
              nrow()

            # how many NAs in that variable
            nas <- input_dat[,var] %>%
              dplyr::filter(is.na(.)) %>%
              nrow()

            # what is the percentage of NAs in that column
            ratio <- (nas/length)

            # if the ratio is 100% NAs, record that in the tracker file and jump to the next
            if(ratio == min_ratio) {
              cat('Step_3_Error', file = tracker_fp, sep = "\n", append=TRUE)
              cat(glue::glue('---- ', '{site} has no {var} data, jumping to next site-year'),
                  file = tracker_fp, sep = "\n", append=TRUE)
              writeLines(glue::glue('{site} has no {var} data, jumping to next site-year'))
            } else { # or keep going
              cat(glue::glue('Step_3_{i}_Success'),
                  file = tracker_fp, sep = "\n", append=TRUE)
              cat(glue::glue('---- ', '{site} has {round(as.numeric(ratio)*100,2)}% NA of {var} data'),
                  file = tracker_fp, sep = "\n", append=TRUE)
            }
          } # end for loop

          # Step 4: define KQ priors and model specifications
          # pull the earliest and latest date in the modeled run for file naming purposes
          tryCatch({
            # define start and end dates from the filtered dataset
            start_date <- lubridate::date(dplyr::first(input_dat$solar.time))
            end_date <- lubridate::date(dplyr::last(input_dat$solar.time))

            # read in results from MLE model run
            mle_priors <- readr::read_csv(glue::glue('data/model_runs/MLE/MLE_diagnostics.csv'))

            # subset to pull gas exchange results
            median_K600_mle <- mle_priors %>%
              dplyr::filter(site %in% !!site,
                            wyear %in% !!wyear,
                            q_type %in% !!q_type,
                            z_method %in% !!z_method,
                            sensor_src %in% !!sensor_src) %>%
              dplyr::pull(K_median)

            minQ <- input_dat %>%
                dplyr::mutate(logQ = log(discharge)) %>%
                dplyr::filter(is.finite(logQ)) %>%
                dplyr::summarise(minQ = min(discharge, na.rm = TRUE)) %>%
                dplyr::pull()

            logQ <- ifelse(is.finite(log(input_dat$discharge)),
                           log(input_dat$discharge),
                           minQ)

            # create an object that sorts discharge in log space
            Q_by_day <- tapply(logQ,
                               substr(input_dat$solar.time, 1, 10),
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
              K600_daily_sigma_sigma = median_K600_mle* 0.02,
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
            cat(paste0(logcode, '-- failed to define K600-Q Bayesian prior, model setup failed'),
                append = TRUE, file = 'nmh_log.txt', sep = '\n')
            next
          } else {
            neon_bayes_row[1, 'Setup Error'] <- FALSE
            cat('Step_4_Success', file = tracker_fp, sep = "\n", append=TRUE)
            cat(paste0(logcode, '-- defined K600-Q Bayesian prior, model setup succeeded'),
                append = TRUE, file = 'nmh_log.txt', sep = '\n')
          }

          # Step 5: run the model
          tryCatch({
            # run the model
            dat_metab <- streamMetabolizer::metab(specs = bayes_specs,
                                                  data = input_dat)
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
              neon_bayes_row[1, 'Fit Time'] <- round(fit_time[3], 4)
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
            writeLines(glue::glue("fit error on\n site: {site}\n year: {wyear}"))
            cat('Step_5_fitError', file=tracker_fp, sep = "\n", append=TRUE)
            cat(glue::glue('---- {model_error}'), file=tracker_fp, sep = "\n", append=TRUE)
            cat(paste0(logcode, '-- Modeling returned an error for ', site.,'-',wyear),
                append = TRUE, file = 'nmh_log.txt', sep = '\n')
            next
          } else {
            cat('Step_5_fitSuccess', file=tracker_fp, sep = "\n", append=TRUE)
            cat(paste0(logcode, '-- Model fit successfully for ', site, '-', wyear),
                append = TRUE, file = 'nmh_log.txt', sep = '\n')
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
            fn_prefix <- paste0('/',site, '_', wyear, '_Q-', q_type, '_z-', z_method, '_TS-', sensor_src, '_')

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

      } # end foreach() loop

  writeLines("---- NEON Metabolism Bayes Model Runs Complete -----")
  readr::write_csv(data.frame(neon_bayes_results),
                   glue::glue(write_dir, '/neon_bayes_row.csv')) #save run results if you want to

  this_time <- Sys.time()
  if(log) {
    cat(paste0(logcode, '-- all model runs completed at ', this_time),
      append = TRUE, file = 'nmh_log.txt', sep = '\n')
    }

  return(neon_bayes_results)

} # end function
