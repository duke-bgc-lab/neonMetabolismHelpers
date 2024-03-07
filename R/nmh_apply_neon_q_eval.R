#' Apply discharge evaluation (Rhea et al. 2023) to raw Continuous discharge data product
#' 
#' @param q_eval object containing discharge evaluation details. If not included, will be downloaded using get_neon_eval_q()
#' @param q_df data frame of continuous discharge product
#' @param raw_q_dir directory of raw Continuous discharge data product
#' @param write_dir directory for where to save the evaluated discharge product
#' @param site 4 character NEON site code to specify
#' @param qaqc_keep which tier of discharge evaluation should be kept, Tier 1 or Tier 2?
#' @param q_write logical. Write data to csv file or return the dataframe?
#' @export 
nmh_apply_neon_q_eval <- function(q_eval = NULL, 
                                  q_df = NULL, 
                                  raw_q_dir = 'data/raw/neon', 
                                  write_dir = 'data/munged/qaqc', 
                                  site = NULL,
                                  qaqc_keep = c('Tier1', 'Tier2'), 
                                  q_write = FALSE) {

  if(is.null(site)) {
    stop(glue('no site name supplied to function'))
  }

  if(is.null(q_df)) {
    tryCatch(
      expr = {
        q_dir <- file.path(raw_q_dir)  # Discharge from NEON
        writeLines(glue::glue('looking for {site} discharge data at {q_dir}'))

        # fault tolerance: did the data load?
        q_df <- feather::read_feather(glue::glue(q_dir, '/{site}/csd_continuousDischarge.feather'))
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

  # checks
  qaqc_check <- unique(q_evaluated$month_year) %in% unique(q_eval_bysite$month_year)

  if(FALSE %in% qaqc_check) {
    stop('QAQC failed somehow, check q_eval and q_df input data')
  }

  if(nrow(q_evaluated) == 0) {
    writeLines(paste0('no discharge data returned which passed QAQC evaluation',
                      ', returning empty dataframe. \n\n', 'this could be a correct rersult, if none of the data',
                      'from your selected site and time frame pass the QAQC evaulation'))
    return(q_evaluated)
  }

  if(q_write == FALSE) {
    return(q_evaluated)
  } else {
    q_eval_fp <- file.path(write_dir, '/{site}/csd_continuousDischarge.feather')
    writeLines(paste('writing evaluated discharge file to', q_eval_fp))
    write.csv(q_evaluated, q_eval_fp)
  }
}
