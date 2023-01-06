model_neon_metab_mle <- function(input_dir = 'data/sm_input/',
                                 type = c('raw', 'source', 'simulated')) {

    library(glue)
    library(tidyverse)
    library(lubridate)
    library(streamMetabolizer)
    library(dplyr)

    if(!type %in% c('raw','source', 'simulation')){
        stop('Error: please select a discharge input from:\n 1) "raw": Raw NEON input\n 2) "source": NEON discharge evaluated by Rhea et al. (accepted), or\n 3) "simulation": NEON discharge simulations by the Macrosheds project')
    }

    source('src/helpers.R')
    source('src/neon_helpers.R')

    # define lotic sites we are running the model on
    site_code <- get_neon_site_data()[,'site_code']

    # create a data frame with all possible site-years for the model to run
    site_years <- data.frame(
        site = rep(site_code, each = length(2016:2021)),
        year = rep(2016:2021, times = length(site_code))
    )

    # convert to a matrix that we will populate with tracking data if the model does or does not fit
    neon_mle_results = matrix('',
                              ncol = 5,
                              nrow = nrow(site_years))
    colnames(neon_mle_results) = c('Site', 'Year', 'Has Data','Fit Error', 'Fit Time')

    # define the model
    init_name = streamMetabolizer::mm_name(
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

    # loop around site-years
    for(j in 1:length(site_years$site)){

        # index the site and year
        # models are run for each site year, so we will filter and keep track of each site-year
        site <- site_years[j,'site']
        year <- site_years[j,'year']

        # populate the results matrix
        neon_mle_results[j,'Site'] <- site
        neon_mle_results[j,'Year'] <- year

        # read in the formatted input data
        input_dat <- try(read_csv(glue(input_dir,'{type}/{site}_{type}_smReady.csv')))

        if(inherits(input_dat, 'try-error')) {
            neon_mle_results[j,'Has Data'] = FALSE
            next
        }

        # subset to the year we are modeling from
        input_dat_sub <- input_dat %>%
            filter(lubridate::year(solar.time) == year)

        if(nrow(input_dat_sub) == 0){
            neon_mle_results[j,'Has Data'] = FALSE
            next
        } else{
            neon_mle_results[j, 'Has Data'] = TRUE
        }

        # pull the earliest and latest date in the modeled run
        start_date <- first(lubridate::date(input_dat_sub$solar.time))
        end_date <- last(lubridate::date(input_dat_sub$solar.time))

        # tracker: assume all models will be fit successfully
        fit_err <- FALSE

        # fault tolerance: catch an errors and dictate the response to the errors
        tryCatch(

            # tries to run the model with specs and data defined above
            {dat_metab = streamMetabolizer::metab(init_specs,
                                                  data = input_dat_sub)

            # extract the outputs from the model
            dat_fit = streamMetabolizer::get_fit(dat_metab)
            },

            # evalute fault tolerance: did the model fit?
            error = function(e) {
                fit_err <- TRUE
            }#,

            # warning = function(w) {
            #   fit_err <<- TRUE
            # }
        )


        # did the model fit or not?
        neon_mle_results[j,'Fit Error'] <- ifelse(fit_err == TRUE,
                                        'Fit Error',
                                        'Fit successful')

        # if the model fit, how long did it take?
        neon_mle_results[j,'Fit Time'] <- ifelse(fit_err == TRUE,
                                        'NA',
                                        get_fitting_time(dat_metab)[3])

        # if there was an error fitting a site-year, jump to the next level in the for loop
        if(fit_err)
            next

        # save model output
        # where are we saving data to
        write_dir = glue('data/model_runs/{type}/')

        if(!dir.exists(write_dir)){
            dir.create(write_dir)
        }

        if(!dir.exists(glue(write_dir,'/MLE')))
            dir.create(glue(write_dir, '/MLE'))
        if(!dir.exists(glue(write_dir,'/MLE/dat_fit')))
            dir.create(glue(write_dir, '/MLE/dat_fit'))
        if(!dir.exists(glue(write_dir,'/MLE/specs')))
            dir.create(glue(write_dir, '/MLE/specs'))
        if(!dir.exists(glue(write_dir,'/MLE/mod_obs_DO')))
            dir.create(glue(write_dir, '/MLE/mod_obs_DO'))

        # write the fit data
        fn_prefix <- paste0(write_dir, 'MLE/dat_fit/', site, '_', start_date, '_', end_date, '_')
        write_csv(dat_fit, paste0(fn_prefix, 'dat_fit.csv'))

        # write the specs
        specs_out <- data.frame(unlist(get_specs(dat_metab)))
        fn_prefix <- paste0(write_dir, 'MLE/specs/', site, '_', start_date, '_', end_date, '_')
        write_csv(specs_out, paste0(fn_prefix, 'specs.csv'))

        # write the output data
        data_out = get_data(dat_metab)
        fn_prefix = paste0(write_dir, 'MLE/mod_obs_DO/', site, '_', start_date, '_', end_date, '_')
        write_csv(data_out, paste0(fn_prefix, 'mod_and_obs_DO.csv'))

    } # end for loop

    # return the matrix summarizing the model fits
    write_csv(data.frame(neon_mle_results),
              glue('data/model_runs/{type}/MLE/neon_{type}_MLE_results.csv'))

    return(neon_mle_results)

} # end function
