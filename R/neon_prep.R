prep_daily_discharge <- function(site_data, domain = 'neon') {
    # pull in site data
    site_code <- read_csv('data/site_data/site_data.csv') %>%
    # ensure domain is neon
      filter(domain == 'neon') %>%
    # extract site codes
      pull(site_code)

    # loop through data
    for(i in 1:length(site_code)){

      site <- site_code[i]

      d <- try(read_csv(glue('data/discharge/Q_simulations/simulated_discharge_2022-01-25/Q_by_site/{site}.csv')) %>%
                 filter(year(date) > 2014) %>%
                 mutate(discharge_m3s_sim = discharge_Ls_sim/1000) %>%
                 select(site_code, date, discharge_m3s_sim))

      if(inherits(d, 'try-error')) {
        next
      }

      write_csv(d,
                glue('data/sm_ready_dailyQ/{site}_daily_Q.csv'))
    }
}

prep_neon_metab_inputs <- function(dir = 'data/raw',
                                   type = 'raw') {

    # make sure supplied q type is in accepted options
    if(all(!type %in% c('raw','source', 'simulation'))){
        stop('Error: please select a discharge input from:\n 1) "raw": Raw NEON input\n 2) "source": NEON discharge evaluated by Rhea et al. (accepted), or\n 3) "simulation": NEON discharge simulations by the Macrosheds project')
    }

    {
        library(tidyverse)
        library(dplyr)
        library(lubridate)
        library(feather)
        library(glue)
        # Metabolism functions
        library(StreamPULSE)
        library(streamMetabolizer)
        # Time-series interpolation
        library(imputeTS)
        library(padr)

        # source functions to complete downstream function
        source('src/helpers.R')
        source('src/neon_helpers.R')
        source('src/streampulse_helpers.R')
    }


    if(any(type == 'simulation')) {
        neon_Q_sim <- get_neon_q_simulated()
    }


    # read in site data file
    site_data <- get_neon_site_data(arg = 'deets')

    # create a vector of lotic site codes
    site_id <- site_data$site_code

    # specify directories where the raw data is saved
    sp_dir <- glue(dir, '/streampulse')          # DO and temperature data, from streampulse
    q_dir <- glue(dir, '/Continuous discharge')  # Discharge from NEON
    bp_dir <- glue(dir, '/Barometric pressure')  # Barometric pressure from NEON
    light_dir <- glue(dir, '/Photosynthetically active radiation at water surface/') # PAR from NEON

    # for loop across each site
    # this will manipulate data at the site level within the function
    for(i in 1:length(site_id)) {

        # define site and its geographic coordinates
        site <- site_data$site_code[i]
        lat <- site_data$latitude[i]
        lon <- site_data$longitude[i]

        # read in the temperature and DO sensor data
        sp_data <- read_csv(glue(sp_dir,'/{site}.csv'))

        # manipulate sensor data
        sp_data_clean <- sp_data %>%
            dplyr::filter(is.na(flagtype)) %>%                # remove flagged data; see text for specifics and StreamPULSE portal for definitions
            pivot_wider(names_from = 'variable',              # convert from long to wide format
                        values_from = 'value') %>%
            mutate(difftime = c(NA, diff(DateTime_UTC))) %>%  # determine difference in time between each row
            filter(difftime <= 180,
                   is.finite(DO_mgL),
                   is.finite(WaterTemp_C)) %>%
            mutate(DO_mgL = na_kalman(DO_mgL, maxgap = 12),
                   WaterTemp_C = na_kalman(WaterTemp_C, maxgap = 12)) %>%                # remove gaps that are greater than 3 hours (10800 seconds)
            pad('15 min') %>%                                 # interpolate time stamps in these gaps at 15 min intervals
            dplyr::select(-difftime)


        # define start and end dates for each site- will be part of file naming convention
        start_date <- first(sp_data_clean[1,'DateTime_UTC'])
        end_date <- last(sp_data_clean[nrow(sp_data_clean),'DateTime_UTC'])

        # read in barometric pressure data
        bp <- try(read_feather(glue(bp_dir, '/{site}/BP_1min.feather')))

        # fault tolerance: if NEON data doesn't exist, use streamPULSE function to find nearest pressure data
        if(inherits(bp, 'try-error')) {
            print(paste0(site, ' bp failed to load'))

            bp_15 <- FindandCollect_airpres(lat = lat,                     # define where the site is
                                            lon = lon,
                                            start_datetime = start_date,   # and what dates to look for
                                            end_datetime = end_date) %>%
                select(-air_temp) %>%                                      # remove columns
                rename(BP_15min = air_mb)                                  # standardize naming conventions
        } else {
            bp_15 <- bp %>%
                mutate(startDateTime = lubridate::ymd_hms(startDateTime,
                                                          tz = 'UTC'),
                       DateTime_UTC = lubridate::round_date(startDateTime,      # NEON data read at 1 minute interval, aggregate to 15 minute windows
                                                            '15 minutes')) %>%
                group_by(DateTime_UTC) %>%                                      # group by 15 minute window
                summarise(BP_15_mean = mean(staPresMean,                        # calculate mean bp during each 15 minute window
                                            na.rm = TRUE)) %>%
                mutate(BP_15min = BP_15_mean*10)                                # unit conversion
        }

        # Join sensor data to barometric pressure data
        use <- full_join(sp_data_clean,
                         bp_15,
                         by = 'DateTime_UTC')

        # calculate DO.sat
        use <- use %>%
            mutate(DO.sat = calc_DO_sat(temp.water = WaterTemp_C,   # streamMetabolizer function to estimate DO saturation
                                        pressure.air = BP_15min,    # based on bp, temp, salinity
                                        salinity.water = 0,
                                        model = 'garcia-benson'))

        # read in light data
        par <- try(read_feather(glue(light_dir, '{site}/PARWS_5min.feather')))

        if(inherits(par, 'try-error')) {
            print('No NEON measured PAR at this site, will esimate from streamMetabolizer::calc_light()')
        } else {

            par_15 <- par %>%
                dplyr::filter(horizontalPosition %in% c('102', '103'),
                              PARFinalQF == 0) %>%       # select the PAR sensors to use
                mutate(DateTime_UTC = lubridate::round_date(startDateTime,       # like bp, light is aggregated to 15 min intervals
                                                            '15 minutes')) %>%
                group_by(DateTime_UTC) %>%
                summarise(PAR_15min = sum(PARMean, na.rm = TRUE))
        }

        # join light data to sensor data
        use <- full_join(use,
                         par_15,
                         by = 'DateTime_UTC')

        if(type == 'raw') {
            # TODO: copy from NM repo of old code
            # read in discharge file
            discharge <- try(read_feather(glue(q_dir, '/{site}/csd_continuousDischarge.feather')))

            # fault tolerance: did the data load?
            if(inherits(discharge, 'try-error')) {
                print(paste0(site, ' failed to load discharge'))
                next
            }

            q_final <- discharge %>%
                dplyr::select(endDate,
                              equivalentStage,
                              maxpostDischarge) %>%
                mutate(DateTime_UTC = lubridate::round_date(endDate,
                                                            '15 minutes')) %>%
                group_by(DateTime_UTC) %>%
                summarise(Q_15min = mean(maxpostDischarge, na.rm = TRUE),
                          #z_15min = mean(equivalentStage, na.rm = TRUE)
                )
        }

        if(type == 'source') {
            # append discharge data
            # read in evaluation from Rhea et al. (in review)
            # this function pulls the hydroshare dataset for Rhea et al. (in review)
            q_eval <- get_neon_eval_q()

            # read in discharge file
            discharge <- try(read_feather(glue(q_dir, '/{site}/csd_continuousDischarge.feather')))

            # fault tolerance: did the data load?
            if(inherits(discharge, 'try-error')) {
                print(paste0(site, ' failed to load discharge'))
            }

            # TOMB uses USGS discharge data, we'll call dataRetrival to access those data
            if(site == 'TOMB'){

                q_tomb <- dataRetrieval::readNWISuv('02469761',                            # site code
                                                    c('00060', '00065'),                   # parameter codes
                                                    startDate = '2015-01-01')
                datetime_cor <- data.frame(datetime_15 = seq.POSIXt(min(q_tomb$dateTime),  # date range to access
                                                                    max(q_tomb$dateTime),
                                                                    by = '15 min')) %>%
                    mutate(dateTime = ymd_hms(paste0(date(datetime_15),
                                                     ' ',
                                                     hour(datetime_15), ':00:00')))
                discharge_filled <- q_tomb %>%
                    full_join(.,
                              datetime_cor,
                              by = 'dateTime') %>%
                    arrange(datetime_15)

                # convert USGS discharge file into the same units and naming conventions as NEON
                q_sub <- discharge_filled %>%
                    mutate(usgs_discharge_liter = X_00060_00000*28.3168,    # convert cfs to L/s
                           usgsStage_m = X_00065_00000*0.3048) %>%          # convert ft to m
                    select(endDate = datetime_15,
                           maxpostDischarge = usgs_discharge_liter,
                           equivalentStage = usgsStage_m) %>%
                    mutate(siteID = site)

                discharge <- q_sub %>%
                    select(siteID, endDate, equivalentStage, maxpostDischarge) %>%
                    mutate(year = lubridate::year(endDate),
                           month = lubridate::month(endDate),
                           maxpostDischarge = maxpostDischarge/1000)       # convert to m3/s
            } else {
                # apply Rhea et al. NEON q evaluations (tiers assigned to data quality by month and year) to raw NEON discharge
                # to ensure use of only Q data which meets or exceeds minimum standards
                qaqc_keep = c('Tier1', 'Tier2')
                discharge <- apply_neon_eval_q(q_eval = q_eval,
                                               q_df = discharge,
                                               site = site,
                                               qaqc_keep = qaqc_keep)
            }

            q_final <- discharge %>%
                dplyr::select(endDate,
                              equivalentStage,
                              maxpostDischarge) %>%
                mutate(DateTime_UTC = lubridate::round_date(endDate,
                                                            '15 minutes')) %>%
                group_by(DateTime_UTC) %>%
                summarise(Q_15min = mean(maxpostDischarge, na.rm = TRUE),
                          #z_15min = mean(equivalentStage, na.rm = TRUE)
                )
        }

        if(type == 'simulation') {
            # TODO: re-run the get_neon_q_simulated() to get the sub-daily data and write a pipeline to process it

            q_final <- neon_Q_sim %>%
                rename(site_id = site) %>%
                dplyr::filter(site == site_id) %>%
                dplyr::select(DateTime_UTC = datetime,
                              Q_15min = Q_predicted)
        }

        # compile the final output file for each site
        out <- left_join(use,                   # join sensor data
                         q_final,                    # and discharge file
                         by = 'DateTime_UTC') %>%
            arrange(DateTime_UTC) %>%                # sorts dates from old to new, for some reason merging put things out of order
            pad %>%                                  # fills in missing 15 minute rows with NAs
            dplyr::mutate(solar.time = convert_UTC_to_solartime(DateTime_UTC,        # streamMetabolizer uses solar time, using internal function to convert to solar.time
                                                                longitude = lon,
                                                                time.type = 'mean'),
                          light = ifelse(is.na(PAR_15min),                    # if light data is missing from NEON, use internal functions to gap fill light data
                                         calc_light(solar.time = solar.time,
                                                    latitude = lat,
                                                    longitude = lon),
                                         PAR_15min),
                          #depth = z_15min) %>%
                          mean_depth = calc_depth(Q_15min,     # convert discharge to mean depth in the study reach; this uses hydraulic scaling coefficients from Leopold and Maddock 1953
                                                  c = 0.409,   # coefficients from Raymond et al. 2012
                                                  f = 0.294)   # coefficients from Raymond et al. 2012
            ) %>%
            # format for streamMetabolizer input, including column headers; these are required input file types
            dplyr::select(solar.time,
                          DO.obs = DO_mgL,
                          DO.sat,
                          depth = mean_depth,
                          discharge = Q_15min,
                          temp.water = WaterTemp_C,
                          light) %>%
            full_days()                             # checks that each day has 96 observations (every 15 minutes)


        if(!dir.exists(glue('data/sm_input/{type}'))){
            print('Directory does not exist, creating now')
            dir.create(glue('data/sm_input/{type}'))
        }

        # write a CSV file for each site
        write_csv(out,
                  glue('data/sm_input/{type}/{site}_{type}_smReady.csv'))
    } # end for loop

} # end function
