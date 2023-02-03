#' @export
nmh_get_neon_data <- function(product_codes = 'all', q_type = 'raw',
                              dest_fp = NULL, # file path where all data is saved to
                              log = TRUE,
                              site_filter = NULL, startdate = NA, enddate = NA,
                              neon_api_token = NA, stream_only = TRUE, forceParallel = FALSE,
                              check.size = TRUE,
                              quietly = FALSE
                              ) {
                              ## q_evaluation_fp = 'data/raw/macrosheds', # file path where NEON Q evaluation data is saved to and read from
                              ## q_evaluated_fp = 'data/munged' # file path where evaluated Q is saved to and read from

    # require macrosheds package to be loaded, and if user does not have it, option to install form github
    pkg_namespace_require(pkg = 'macrosheds', quietly = quietly)

    # start logging for this function
    logcode = 'nmh_get:'
    if(log) {
        cat("neonMetabolismHelpers Logging File", file="nmh_log.txt", sep="\n")
        cat(paste0(logcode, "now logging nmh_get_neon_data() request"), file="nmh_log.txt", append=TRUE,sep="\n")
        cat(paste0(logcode, "all logging from this request preceded by the following logging code prefix - ", logcode), file="nmh_log.txt", append=TRUE,sep="\n")

        this_time = Sys.time()
        cat(paste0(logcode, "system datetime at beginning,", this_time), file="nmh_log.txt", append=TRUE,sep="\n")
    }

    # decide folder structure
    data_folder <- 'data'
    data_raw_folder <- file.path(data_folder, 'raw')
    data_munged_folder <- file.path(data_folder, 'munged')

    # checking for user defined filepath
    if(is.null(dest_fp)) {
      # if none, save to working directory in "data" folder
        dest_fp <- file.path(getwd(), data_folder)
    }

    if(log) {
        cat(paste0(logcode, "all raw data to be saved in", dest_fp), file="nmh_log.txt", append=TRUE,sep="\n")
    }

    # if 'data' folder does not exist
    if(!dir.exists(dest_fp)){
        # create 'data' folder
        print(paste('Directory does not exist, creating directory here:', dest_fp))
        dir.create(dest_fp, recursive = TRUE)

        if(log) {
            cat(paste0(logcode, "created directory at", dest_fp), file="nmh_log.txt", append=TRUE,sep="\n")
        }
    }

    # default is to retrieve all data necessary for NEON metabolism, this is
    if(product_codes == 'all') {
      # products are: discharge, light, barometric pressure, water temp, & dissolved oxygen
      product_codes <- c('DP4.00130.001', 'DP1.20042.001', 'DP1.00004.001', 'DP1.20053.001', 'DP1.20264.001','DP1.20288.001')
    }

    products_n = length(product_codes)
    writeLines(paste('retrieving NEON data for', products_n, 'data products'))

    if(log) {
        cat(paste0(logcode, "retreiving", products_n, "NEON data products:"), file="nmh_log.txt", append=TRUE,sep="\n")
        cat(paste0(logcode, "\n    ", product_codes), file="nmh_log.txt", append=TRUE,sep="\n")
        cat(paste0(logcode, 'looping through product codes now'))
    }


    for(product_code in product_codes) {

      if(log) {
          cat(paste0(logcode, "-- retrieving ", product_code), file="nmh_log.txt", append=TRUE,sep="\n")
      }

      tryCatch(
        expr = {

          # discharge handling (raw, qaqc, or simulated)
          if(product_code == 'DP4.00130.001') {
            if(!q_type %in% c('qaqc', 'raw', 'simulated')){
              stop("invalid q_type selection, must be 'raw' 'qaqc' or 'simulated'")

              if(log) {
                 cat(paste0(logcode, "-- ERROR invalid q_type selection, must be 'raw' 'qaqc' or 'simulated'"), file="nmh_log.txt", append=TRUE,sep="\n")
              }
            } else {
              if(log) {
                 cat(paste0(logcode, "-- discharge data type:", q_type), file="nmh_log.txt", append=TRUE,sep="\n")
              }
            }
          }


           # if method is QC, we need to download the NEON Q evaluation data
           if(product_code == 'DP4.00130.001' & q_type == 'qaqc'){
             # we should download this to the
             neon_eval_q_dir <- file.path(data_raw_folder, 'macrosheds', 'qaqc')
             neon_eval_q_fn <- 'neon_q_eval.csv'
             neon_eval_q_fp <- file.path(neon_eval_q_dir, neon_eval_q_fn)

             # downloads q_eval and saves as df
             neon_eval_q_df <- nmh_get_neon_q_eval(
               dest_fp = neon_eval_q_dir,
               dest_fn = neon_eval_q_fn,
               download = TRUE
             )

             if(log) {
                cat(paste0(logcode, "-- saving NEON discharge evaluation data at", neon_eval_q_fp), file="nmh_log.txt", append=TRUE,sep="\n")
             }

             # also need to check if "raw" discharge exists, and retrieve if it is not there
             product_name <- neonUtilities::getProductInfo(dpID = product_code, token = neon_api_token)$productName
             raw_q_data_fp <- file.path(data_raw_folder, 'neon', product_name)
             # also set "evaluated" w filepath
             ## munged_q_data_fp <- file.path(data_munged_folder, 'neon', product_name)

             # look for raw NEON Q data, necessary for QAQC
             if(!dir.exists(raw_q_data_fp)) {
                if(log) {
                   cat(paste0(logcode, "-- ERROR raw NEON discharge not found, cant run QAQC method", neon_eval_q_fp), file="nmh_log.txt", append=TRUE,sep="\n")
                }

                stop(paste0('it looks like you do not have a directory at ', raw_q_data_fp,
                           ' \nmake sure you have run this function to retrieve "raw" NEON discharge',
                           ' before running the QAQC version, which filters raw NEON discharge data\n',
                           ' through a rigorous monthly evaluation of the usability of NEON discharge data'))
             }

           } # end evaluate q conditional

           # if the product is discharge and the q_type is simulated
           if(product_code == 'DP4.00130.001' & q_type == 'simulated'){

             if(log) {
                cat(paste0(logcode, "-- retrieving simulated Q data"), file="nmh_log.txt", append=TRUE,sep="\n")
             }

             # get simulated discharge
             qsim <- try(
               nmh_get_neon_Q_sim()
             )

             if(inherits(qsim, 'try-error')) {
                if(log) {
                   cat(paste0(logcode, "-- ERROR simulated Q retrieval failed"), file="nmh_log.txt", append=TRUE,sep="\n")
                }
                 stop(paste('q simulation data download failed, try anohter q type'))
             } else {
                 # Create the file path for where the files will be saved. This can be changed
                 # if you  want to save the files in a different location
                 raw_qsim_dest <- file.path(data_raw_folder, 'macrosheds', 'simulated')
                 raw_qsim_fp <- paste0(raw_qsim_dest, '/q_simulation.csv')

                 if(!dir.exists(raw_qsim_dest)) {
                     # create direcotry if doesnt exists
                     print(paste('Directory does not exist, creating directory here:', raw_qsim_dest))
                     dir.create(raw_qsim_dest, recursive = TRUE)
                 }

                if(log) {
                   cat(paste0(logcode, "-- saved simulated NEON Q data to", raw_qsim_fp), file="nmh_log.txt", append=TRUE,sep="\n")
                }

                 writeLines(paste('download simulated NEON q data to:\n', "    ", raw_qsim_fp))
                 write.csv(qsim, raw_qsim_fp)
             }

           } else {

             if(log) {
                cat(paste0(logcode, "-- retrieving", product_code), file="nmh_log.txt", append=TRUE,sep="\n")
             }

             # get NEON product name, used for filepath
             product_name <- neonUtilities::getProductInfo(dpID = product_code, token = neon_api_token)$productName
             data_fp <- file.path(data_raw_folder, 'neon', product_name)

             writeLines(paste('retrieving NEON product:', product_name,
                          'and saving results at:\n', data_fp))

             if(!dir.exists(data_fp)){
                 # create direcotry if doesnt exists
                 print(paste('Directory does not exist, creating directory here:', data_fp))
                 dir.create(data_fp, recursive = TRUE)
             }

             if(log) {
               cat(paste0(logcode, "-- saved NEON product:"," \n", product_code,
                          " \n", product_name, "\n to ", data_fp), file="nmh_log.txt", append=TRUE,sep="\n")
             }

             writeLines(paste('checking available NEON sites with data for product', product_name))

             # call to NEON for avilable data of code type
             avail_sets <- get_avail_neon_product_sets(product_code)

             # sites form this record
             avail_sites <- unique(avail_sets$site_name)

             # filter, by default, to only stream sites (only NEON sites included in MacroSheds)
             if(stream_only) {
               writeLines('this package is meant for stream and river data only -- filtering available sites\n',
                          'to only stream and river type NEON sites')

               neon_sites <- macrosheds::ms_load_sites() %>%
                 dplyr::filter(network == 'neon') %>%
                 dplyr::pull("site_code")

               avail_sites <- avail_sites[avail_sites %in% neon_sites]

               possible_site_n <- length(unique(neon_sites))

               if(possible_site_n != 27) {
                 warning(paste0('product query found ', possible_site_n,
                                   'possible NEON sites (from macrosheds data). Number should be 27 as of January 2023'))
               }

               ## avail_sites_n <- length(unique(avail_sites))
               ## writeLines(paste0('product query found ', available_site_n, 'sites with data out of',
               ##                   possible_site_n, 'total NEON sites (user supplied filter not yet applied)'))
             }



             if(!is.null(site_filter)) {
                 print('values supplied to site_filter argument, filtering queryable sites to those in site_filter values')
                 avail_sites <- avail_sites[avail_sites %in% site_filter]
             }

             avail_sites_n <- length(unique(avail_sites))
             if(avail_sites_n == 0) {
                 stop('no available sites')
             } else {
               # download product for each site
               writeLines(paste('\nquerying for NEON product code:', product_code,
                                '\n  total sites to query:', avail_sites_n))
             }

             if(log) {
                cat(paste0(logcode, "-- querying ", avail_sites_n, " sites"), file="nmh_log.txt", append=TRUE,sep="\n")
             }

             for(j in 1:length(avail_sites)) {

                 writeLines(paste('  querying site', j, 'of', avail_sites_n))

                 # Defines what site to retrieve
                 site_name <- avail_sites[j]
                 writeLines(paste('  site name:', site_name))

                 if(log) {
                    cat(paste0(logcode, "-- querying ", site_name, "\n"), file="nmh_log.txt", append=TRUE,sep="\n")
                 }

                 # Download data product for the site
                 tryCatch(
                     expr = {
                         data_pile <- neonUtilities::loadByProduct(dpID = product_code,
                                                                   site = site_name,
                                                                   package='basic',
                                                                   startdate = startdate,
                                                                   enddate = enddate,
                                                                   token = neon_api_token,
                                                                   forceParallel = forceParallel,
                                                                   check.size = check.size
                                                                   )
                     },
                     error = function(e) {
                         writeLines(paste('ERROR at', site_name, 'trying again with expanded package'))
                         data_pile <- try(neonUtilities::loadByProduct(dpID = product_code,
                                                                       site = site_name,
                                                                       package='expanded',
                                                                       startdate = startdate,
                                                                       enddate = enddate,
                                                                       token = neon_api_token,
                                                                       forceParallel = forceParallel,
                                                                       check.size = check.size
                                                                       ))

                     }
                 )

                 if(inherits(data_pile, 'try-error')) {
                     writeLines(paste('no data downloaded skipping', site_name))

                     if(log) {
                        cat(paste0(logcode, "-- query failed, no data download, skipping"), file="nmh_log.txt", append=TRUE,sep="\n")
                     }

                     next
                 } else {
                     # Create the file path for where the files will be saved. This can be changed
                     # if you  want to save the files in a different location
                     raw_data_dest <- file.path(data_fp, site_name)

                     # Save the list of all dataframes for the site as individual feather files
                     # in the file path raw_data_dest
                     serialize_list_to_dir(data_pile, raw_data_dest)

                     if(log) {
                       cat(paste0(logcode, "-- saved\n product:", product_code,
                                  "\n site:", site_name, " to ",
                                  raw_data_dest), file="nmh_log.txt", append=TRUE,sep="\n")
                     }
                 }

                 # run this to make a new evaluated Q folder for each site
                 if(product_code == 'DP4.00130.001' & q_type == 'qaqc'){

                   # also set "evaluated" w filepath
                   munged_q_data_fp <- file.path(data_munged_folder, 'neon', product_name, site_name)

                   nmh_apply_neon_q_eval(
                     q_eval = neon_eval_q_df,
                     q_df = NULL,
                     raw_q_dir = raw_q_data_fp,
                     write_dir = munged_q_data_fp,
                     site = site_name,
                     q_write = TRUE
                     )

                    if(log) {
                       cat(paste0(logcode, "-- evaluated Q saved to", munged_q_data_fp), file="nmh_log.txt", append=TRUE,sep="\n")
                    }
                 } # end evaluate q conditional

             } # end available sites loop

           } # end 'if simulated discharge' conditional
        }, # end try-expr section
        error = function(e) {
          print(e)
          writeLines('product failed')
          next
        } # end try-error section
      ) # end try
    } # end product loop
}

## nmh_get_neon_data(product_codes = 'DP4.00130.001', q_type = 'raw', site_filter = 'ARIK', startdate = '2017-01', enddate = '2017-12')
