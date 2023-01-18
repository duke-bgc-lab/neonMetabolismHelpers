nmh_get_tomb_q <- function(site = 'TOMB') {
  
  q_tomb <- dataRetrieval::readNWISuv('02469761',                            # site code
                                      c('00060', '00065'),                   # parameter codes for discharge and depth
                                      startDate = '2015-01-01')
  
  datetime_cor <- data.frame(datetime_15 = seq.POSIXt(min(q_tomb$dateTime),  # date range to access
                                                      max(q_tomb$dateTime),
                                                      by = '15 min')) %>%
    dplyr::mutate(dateTime = lubridate::ymd_hms(paste0(lubridate::date(datetime_15),
                                            ' ',
                                            lubridate::hour(datetime_15), ':00:00')))
  
  discharge_filled <- q_tomb %>%
    dplyr::full_join(.,
                     datetime_cor,
                     by = 'dateTime') %>%
    dplyr::arrange(datetime_15)
  
  # convert USGS discharge file into the same units and naming conventions as NEON
  q_sub <- discharge_filled %>%
    dplyr::mutate(usgs_discharge_liter = X_00060_00000*28.3168,    # convert cfs to L/s
                  usgsStage_m = X_00065_00000*0.3048) %>%          # convert ft to m
    dplyr::select(endDate = datetime_15,
                  maxpostDischarge = usgs_discharge_liter,
                  equivalentStage = usgsStage_m) %>%
    dplyr::mutate(siteID = site)
  
  discharge <- q_sub %>%
    dplyr::select(siteID, endDate, equivalentStage, maxpostDischarge) %>%
    dplyr::mutate(year = lubridate::year(endDate),
                  month = lubridate::month(endDate),
                  maxpostDischarge = maxpostDischarge/1000)       # convert to m3/s
  
  return(discharge)
} # end function
