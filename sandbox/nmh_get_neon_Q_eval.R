get_neon_eval_q <- function(dest_fp = 'data/raw/qaqc/', dest_fn = 'neon_q_eval.csv') {
  destfile <- paste0(dest_fp, dest_fn)
  
  
  if(!dir.exists(dest_fp)){
    # create direcotry if doesnt exists
    print(paste('Directory does not exist, creating dir at:', dest_fp))
    dir.create(dest_fp)
  }
  
  # read in the discharge evaluation from Rhea et al. (in prep)
  download.file(
    url = 'https://www.hydroshare.org/resource/03c52d47d66e40f4854da8397c7d9668/data/contents/neon_q_eval.csv',
    destfile = destfile
  )
  
  q_eval <- read.csv(destfile)
  
  # filter the sites with bad discharge data
  good_q <- q_eval %>%
    filter(!final_qual %in% c('drift', 'regression_flag', 'Tier3'))
  return(good_q)
}

apply_neon_eval_q <- function(q_eval = NULL, q_df = NULL, dir = 'data/raw', site = NULL,
                              qaqc_keep = c('Tier1', 'Tier2')) {
  
  if(is.null(site)) {
    stop(glue('no site name supplied to function'))
  }
  
  if(is.null(q_df)) {
    tryCatch(
      expr = {
        q_dir <- file.path(dir, 'Continuous discharge')  # Discharge from NEON
        writeLines(glue('looking for {site} discharge data at {q_dir}'))
        
        # fault tolerance: did the data load?
        discharge <- read_feather(glue::glue(q_dir, '/{site}/csd_continuousDischarge.feather'))
      },
      error = function(e) {
        stop('must supply NEON raw continous discharge data as df or in dir filepath')
      }
    )
  }
  
  if(is.null(q_eval)) {
    writeLines(glue('retrieving NEON Q evaluation (Rhea et al.) from hydroshare'))
    q_eval <- get_neon_eval_q()
  }
  
  q_eval_bysite <- q_eval %>%
    dplyr::rename(siteID = site) %>%
    dplyr::filter(
      siteID == site,
      final_qual %in% qaqc_keep,
    ) %>%
    dplyr::select(year, month) %>%
    dplyr::mutate(
      month_year = paste0(month, '_', year)
    )
  
  discharge <- q_df %>%
    select(siteID, endDate, equivalentStage, maxpostDischarge) %>%
    mutate(year = lubridate::year(endDate),
           month = lubridate::month(endDate),
           maxpostDischarge = maxpostDischarge/1000, # convert to m3/s
           month_year = paste0(month, '_', year))
  
  q_evaluated <- discharge %>%
    merge(q_eval_bysite)
  
  # checks``
  qaqc_check <- unique(q_evaluated$month_year) %in% unique(q_eval_bysite$month_year)
  
  if(FALSE %in% qaqc_check) {
    stop('QAQC failed somehow, check q_eval and q_df input data')
  }
  
  return(q_evaluated)
}
