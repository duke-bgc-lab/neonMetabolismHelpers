nmh_prep_metab_inputs <- function(dir = 'data/raw',
                                  q_type = c('raw','qaqc', 'simulation'),
                                  z_method = c('model', 'meas')) {
  
  if(!q_type %in% c('raw','qaqc', 'simulation')) {
    stop('Error: please select a discharge input from:\n 1) "raw": Raw NEON input\n 2) "qaqc": NEON discharge evaluated by Rhea et al. (accepted), or\n 3) "simulation": NEON discharge simulations by the Macrosheds project')
  }
  
  if(!z_method %in% c('model', 'meas')){
    stop('No mean depth method assigned, please select either \n "1) model: use values from Raymond et al. (2012)" or \n "2) meas: site-specific coefficients"')
  }
  
  # TODO: how does sourcing nmh_internals.R work? Is the following line necessary?
  source('sandbox/nmh_internals.R')
  
  if(q_type == 'qaqc') {
    # read in evaluation from Rhea et al. (accepted)
    # this function pulls the hydroshare dataset for Rhea et al. (in review)
    q_eval <- nmh_get_neon_Q_eval()
  }
  
  # get the simulated Q data from MacroSheds portal
  if(q_type == 'simulation') {
    neon_Q_sim <- get_neon_q_simulated()
  }
  
  # read in site data file
  site_data <- get_neon_site_data(arg = 'deets')
  
  # create a vector of lotic site codes
  site_id <- site_data$site_code
  
  # specify directories where the raw data is saved
  sp_dir <- glue::glue(dir, '/streampulse')          # DO and temperature data, from streampulse
  q_dir <- glue::glue(dir, 'neon', '/Continuous discharge')  # default to raw discharge from NEON
  bp_dir <- glue::glue(dir, '/Barometric pressure')  # Barometric pressure from NEON
  light_dir <- glue::glue(dir, '/Photosynthetically active radiation at water surface/') # PAR from NEON
  
  # q_dir is different based on "q type"
  if(q_type == 'simulated') {
    q_dir <- file.path(dir, 'macrosheds', 'simulated')  # Discharge from NEON Q simulations from MacroSheds scientist Mike Vlah
  } else if(q_type == 'qaqc') {
    q_dir <- file.path('data', 'munged', 'qaqc')  # Discharge filtered via NEON Q evaluation from Rhea et al. 2023
  }
  
  # for loop across each site
  # this will manipulate data at the site level within the function
  for(i in 1:length(site_id)){
    
    # define site and its geographic coordinates
    site <- site_data$site_code[i]
    lat <- site_data$latitude[i]
    lon <- site_data$longitude[i]
    
    # read in the temperature and DO sensor data
    sp_data <- readr::read_csv(glue::glue(sp_dir,'/{site}.csv'))
    
    # manipulate sensor data
    sp_data_clean <- sp_data %>%
      dplyr::filter(is.na(flagtype)) %>%                # remove flagged data; see text for specifics and StreamPULSE portal for definitions
      tidyr::pivot_wider(names_from = 'variable',              # convert from long to wide format
                         values_from = 'value') %>%
      dplyr::mutate(difftime = c(NA, diff(DateTime_UTC))) %>%  # determine difference in time between each row
      dplyr::filter(difftime <= 180,                           # filter time-stamps less than 3 hr in length
                    is.finite(DO_mgL),                         # and non-finite data from numeric strings
                    is.finite(WaterTemp_C)) %>%
      dplyr::mutate(DO_mgL = na_kalman(DO_mgL,                 # gap fill using kalman filter for gaps <12 15 minute periods (= 3 hr)
                                       maxgap = 12),
                    WaterTemp_C = na_kalman(WaterTemp_C,
                                            maxgap = 12)) %>%
      padr::pad('15 min') %>%                                 # interpolate time stamps in these gaps at 15 min intervals
      dplyr::select(-difftime)
    
    
    # define start and end dates for each site- will be part of file naming convention
    start_date <- dplyr::first(sp_data_clean[1,'DateTime_UTC'])
    end_date <- dplyr::last(sp_data_clean[nrow(sp_data_clean),'DateTime_UTC'])
    
    # read in barometric pressure data
    bp <- try(feather::read_feather(glue::glue(bp_dir, '/{site}/BP_1min.feather')))
    
    # fault tolerance: if NEON data doesn't exist, use streamPULSE function to find nearest pressure data
    if(inherits(bp, 'try-error')) {
      print(paste0(site, ' bp failed to load'))
      
      bp_15 <- FindandCollect_airpres(lat = lat,                     # define where the site is
                                      lon = lon,
                                      start_datetime = start_date,   # and what dates to look for
                                      end_datetime = end_date) %>%
        dplyr::select(-air_temp) %>%                                      # remove columns
        dplyr::rename(BP_15min = air_mb)                                  # standardize naming conventions
    } else {
      bp_15 <- bp %>%
        dplyr::mutate(startDateTime = lubridate::ymd_hms(startDateTime,
                                                         tz = 'UTC'),
                      DateTime_UTC = lubridate::round_date(startDateTime,      # NEON data read at 1 minute interval, aggregate to 15 minute windows
                                                           '15 minutes')) %>%
        dplyr::group_by(DateTime_UTC) %>%                                      # group by 15 minute window
        dplyr::summarise(BP_15_mean = mean(staPresMean,                        # calculate mean bp during each 15 minute window
                                           na.rm = TRUE)) %>%
        dplyr::mutate(BP_15min = BP_15_mean*10)                                # unit conversion
    }
    
    # Join sensor data to barometric pressure data
    use <- dplyr::full_join(sp_data_clean,
                            bp_15,
                            by = 'DateTime_UTC')
    
    # calculate DO.sat
    use <- use %>%
      dplyr::mutate(DO.sat = streamMetabolizer::calc_DO_sat(temp.water = WaterTemp_C,   # streamMetabolizer function to estimate DO saturation
                                                            pressure.air = BP_15min,    # based on bp, temp, salinity
                                                            salinity.water = 0,
                                                            model = 'garcia-benson'))
    
    # read in light data
    par <- try(feather::read_feather(glue::glue(light_dir, '{site}/PARWS_5min.feather')))
    
    if(inherits(par, 'try-error')) {
      print('No NEON measured PAR at this site, will esimate from streamMetabolizer::calc_light()')
    } else {
      par_15 <- par %>%
        dplyr::filter(horizontalPosition %in% c('102', '103'),
                      PARFinalQF == 0) %>%       # select the PAR sensors to use
        dplyr::mutate(DateTime_UTC = lubridate::round_date(startDateTime,       # like bp, light is aggregated to 15 min intervals
                                                           '15 minutes')) %>%
        dplyr::group_by(DateTime_UTC) %>%
        dplyr::summarise(PAR_15min = sum(PARMean, na.rm = TRUE))
      
      # join light data to sensor data
      use <- dplyr::full_join(use,
                              par_15,
                              by = 'DateTime_UTC')
    }
    
    # compile discharge data
    if(q_type == 'raw') {
      # read in discharge file
      discharge <- try(feather::read_feather(glue::glue(q_dir, '/{site}/csd_continuousDischarge.feather')))
      
      # fault tolerance: did the data load?
      if(inherits(discharge, 'try-error')) {
        print(paste0(site, ' failed to load raw NEON discharge'))
        next
      }
      
      q_final <- discharge %>%
        dplyr::select(endDate,
                      maxpostDischarge) %>%
        dplyr::mutate(DateTime_UTC = lubridate::round_date(endDate,
                                                           '15 minutes')) %>%
        dplyr::group_by(DateTime_UTC) %>%
        dplyr::summarise(Q_15min = mean(maxpostDischarge/1000, na.rm = TRUE))
    }
    
    if(q_type == 'qaqc') {
      # append discharge data
      # read in discharge file
      raw_Q <- try(feather::read_feather(glue::glue(q_dir, '/{site}/csd_continuousDischarge.feather')))
      
      # fault tolerance: did the data load?
      if(inherits(raw_Q, 'try-error')) {
        print(paste0(site, ' failed to load Rhea QAQC NEON discharge'))
      }
      
      # TOMB uses USGS discharge data, we'll call dataRetrival to access those data
      if(site == 'TOMB'){
        
        discharge <- nmh_get_tomb_q()
        #   q_tomb <- dataRetrieval::readNWISuv('02469761',                            # site code
        #                                       c('00060', '00065'),                   # parameter codes
        #                                       startDate = '2015-01-01')
        #   
        #   datetime_cor <- data.frame(datetime_15 = seq.POSIXt(min(q_tomb$dateTime),  # date range to access
        #                                                       max(q_tomb$dateTime),
        #                                                       by = '15 min')) %>%
        #     dplyr::mutate(dateTime = ymd_hms(paste0(date(datetime_15),
        #                                             ' ',
        #                                             lubridate::hour(datetime_15), ':00:00')))
        #   discharge_filled <- q_tomb %>%
        #     dplyr::full_join(.,
        #                      datetime_cor,
        #                      by = 'dateTime') %>%
        #     dplyr::arrange(datetime_15)
        #   
        #   # convert USGS discharge file into the same units and naming conventions as NEON
        #   q_sub <- discharge_filled %>%
        #     dplyr::mutate(usgs_discharge_liter = X_00060_00000*28.3168,    # convert cfs to L/s
        #                   usgsStage_m = X_00065_00000*0.3048) %>%          # convert ft to m
        #     dplyr::select(endDate = datetime_15,
        #                   maxpostDischarge = usgs_discharge_liter,
        #                   equivalentStage = usgsStage_m) %>%
        #     dplyr::mutate(siteID = site)
        #   
        #   discharge <- q_sub %>%
        #     dplyr::select(siteID, endDate, equivalentStage, maxpostDischarge) %>%
        #     dplyr::mutate(year = lubridate::year(endDate),
        #                   month = lubridate::month(endDate),
        #                   maxpostDischarge = maxpostDischarge/1000)       # convert to m3/s
      } else {
        # apply Rhea et al. NEON q evaluations (tiers assigned to data quality by month and year) to raw NEON discharge
        # to ensure use of only Q data which meets or exceeds minimum standards
        qaqc_keep = c('Tier1', 'Tier2')
        
        # this function passes through Rhea et al.
        # and converts from L/s to m3/s
        discharge <- apply_neon_eval_q(q_eval = q_eval,
                                       q_df = raw_Q,
                                       site = site,
                                       qaqc_keep = qaqc_keep)
      }
      
      q_final <- discharge %>%
        dplyr::select(endDate,
                      equivalentStage,
                      maxpostDischarge) %>%
        dplyr::mutate(DateTime_UTC = lubridate::round_date(endDate,
                                                           '15 minutes')) %>%
        dplyr::group_by(DateTime_UTC) %>%
        dplyr::summarise(Q_15min = mean(maxpostDischarge, na.rm = TRUE))
    }
    
    if(q_type == 'simulation') {
      q_final <- neon_Q_sim %>%
        dplyr::rename(site_id = site) %>%
        dplyr::filter(site == site_id) %>%
        dplyr::mutate(Q_15min = Q_predicted/1000) %>%
        dplyr::select(DateTime_UTC = datetime,
                      Q_15min)
    }
    
    if(z_method == 'model') {
      c <- 0.409
      f <- 0.294
    }  
    
    if(z_method == 'meas') {
      coefs <- readr::read_csv('data/NEON_site_scaling_coefs.csv')
      
      good_fits <- coefs %>% 
        dplyr::filter(r2_depth > 0.1) %>% 
        dplyr::pull(site)
      
      if(site %in% good_fits) {
        
        c <- coefs %>% 
          dplyr::filter(!!site == site) %>% 
          dplyr::select(c) %>% 
          dplyr::pull()
        
        f <- coefs %>% 
          dplyr::filter(site %in% !!site) %>% 
          dplyr::select(f) %>% 
          dplyr::pull()
      } else {
        cat(glue::glue('No hydraulic scaling coefficients available at {site}\n Using default values c = 0.409, f = 0.294'))
        c <- 0.409
        f <- 0.294
      }
    }
    
    # compile the final output dataset
    # if NEON provided light, we will go from here:
    if('PAR_15min' %in% names(use)) {
      out <- dplyr::left_join(use,                        # join sensor data
                              q_final,                    # and discharge file
                              by = 'DateTime_UTC') %>%
        dplyr::arrange(DateTime_UTC) %>%                # sorts dates from old to new, for some reason merging put things out of order
        padr::pad('15 min') %>%                                         # fills in missing 15 minute rows with NAs
        dplyr::mutate(solar.time = streamMetabolizer::convert_UTC_to_solartime(DateTime_UTC,        # streamMetabolizer uses solar time, using internal function to convert to solar.time
                                                                               longitude = lon,
                                                                               time.type = 'mean'),
                      mean_depth = streamMetabolizer::calc_depth(Q_15min,                      # convert discharge to mean depth in the study reach; this uses hydraulic scaling coefficients from Leopold and Maddock 1953
                                                                 c = c, f = f)) %>%    # coefficients from Raymond et al. 2012
        dplyr::select(solar.time,
                      DO.obs = DO_mgL,
                      DO.sat,
                      depth = mean_depth,
                      discharge = Q_15min,
                      temp.water = WaterTemp_C,
                      light = PAR_15min) %>%
        full_days()
    } else { # if NEON does not have light data, we will calculate it here
      out <- dplyr::left_join(use,                   # join sensor data
                              q_final,                    # and discharge file
                              by = 'DateTime_UTC') %>%
        dplyr::arrange(DateTime_UTC) %>%                # sorts dates from old to new, for some reason merging put things out of order
        padr::pad('15 min') %>%
        dplyr::mutate(solar.time = streamMetabolizer::convert_UTC_to_solartime(DateTime_UTC,lon,'mean'),
                      depth = streamMetabolizer::calc_depth(Q_15min,
                                                            c = c, f = f),
                      light = streamMetabolizer::calc_light(solar.time, lat, lon)) %>%
        dplyr::select(solar.time,
                      DO.obs = DO_mgL,
                      DO.sat = DO.sat,
                      depth,
                      discharge = Q_15min,
                      temp.water = WaterTemp_C,
                      light)
    }
    
    # create the save directory if need be
    if(!dir.exists(glue('data/sm_input/{q_type}'))){
      print('Directory does not exist, creating now')
      dir.create(glue('data/sm_input/{q_type}'))
    }
    
    # write a CSV file for each site
    readr::write_csv(out,
                     glue::glue('data/sm_input/{q_type}/{site}_{q_type}_smReady.csv'))
  } # end for loop
  
} # end function
