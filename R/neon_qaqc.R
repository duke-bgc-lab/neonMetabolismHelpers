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
                discharge <- read_feather(glue(q_dir, '/{site}/csd_continuousDischarge.feather'))
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

    discharge <- discharge %>%
        select(siteID, endDate, equivalentStage, maxpostDischarge) %>%
        mutate(year = lubridate::year(endDate),
               month = lubridate::month(endDate),
               maxpostDischarge = maxpostDischarge/1000, # convert to m3/s
               month_year = paste0(month, '_', year)
        )

    q_evaluated <- discharge %>%
        merge(q_eval_bysite)

    # checks``
    qaqc_check <- unique(q_evaluated$month_year) %in% unique(q_eval_bysite$month_year)
    if(FALSE %in% qaqc_check) {
        stop('QAQC failed somehow, check q_eval and q_df input data')
    }

    return(q_evaluated)
}


get_neon_site_data <- function(arg = 'details') {
    # create the site codes
    us_states <- USAboundaries::us_states()

    site_data <- macrosheds::ms_download_site_data() %>%
        filter(domain == 'neon',
               site_type == 'stream_gauge')

    site_geo <- site_data %>%
        sf::st_as_sf(coords = c('longitude', 'latitude'), crs = 4326)

    region_code = sf::st_join(site_geo, us_states) %>%
        dplyr::pull(state_abbr)

    sites <- site_data %>%
        dplyr::pull(site_code)

    # build a data frame with sites, site codes, and dates that will carry out of the loop
    # sp_code is the look-up site name in streampulse
    site_deets <- data.frame(
        site_code = sites,
        sp_code = paste0(region_code, '_', sites, 'temp'),
        start_date = NA,
        end_date = NA
    )
    if(arg == 'details') {
        return(site_deets)
    } else {
        return(site_data)
    }
}

# Mean Channel Depth from Q Helper Funcs
# draft functions
batch_hydraulic_geometry_from_ws_area <- function(data, model) {
    # get model coefficients
    neon_model <- left_join(data, model, by = c('HUC02'))

    neon_hg <- neon_model %>%
        mutate(
            slope = as.numeric(slope),
            intercept = as.numeric(intercept),
            ws_area = as.numeric(ws_area),
            # log10(hg) = log10(ws_area)*slope + intercept
            hg = 10^(log10(ws_area)*slope + intercept),
            log10_hg = log10(ws_area)*slope + intercept
        )

    # calculate metric in (m)
    # NOTE: make sure all area is km^2 and all metrics are in cm
    return(neon_hg)
}

model_depth_from_Q <- function(q, hg_model) {
    model <- hg_model %>%
        mutate(siteID = site_code)

    # get model coefficients
    neon_model <- left_join(q, model, by = c('siteID'))
    neon_q_depth <- neon_model %>%
        filter(
            metric == 'depth',
            ## condition == 'baseflow'
        ) %>%
        select(
            site_code = siteID,
            HUC02,
            datetime = collectDate,
            discharge = sectionFlow,
            slope,
            intercept,
            metric,
            condition
        ) %>%
        mutate(
            discharge = as.numeric(discharge),
            slope = as.numeric(slope),
            intercept = as.numeric(intercept),
            depth = 10^(log10(discharge)*slope + intercept)
        )

    neon_composite_depth <- neon_q_depth %>%
        dplyr::group_by(site_code, datetime, condition) %>%
        dplyr::summarise(
            n = dplyr::n(),
            depth = mean(depth, na.rm = TRUE),
            discharge = mean(discharge, na.rm = TRUE),
            .groups = "drop") %>%
        dplyr::filter(n > 1L) %>%
        pivot_wider(
            id_cols = c('site_code', 'datetime', 'discharge'),
            names_from = condition,
            values_from = depth
        ) %>%
        dplyr::group_by(site_code, datetime) %>%
        summarize(
            comp_depth = (bankfull + baseflow)/2,
            bankfull_depth = mean(bankfull, na.rm = TRUE),
            baseflow_depth = mean(baseflow, na.rm = TRUE),
            discharge = mean(discharge, na.rm = TRUE)
        ) %>%
        ungroup() %>%
        pivot_longer(
            cols = ends_with('_depth'),
            names_to = 'condition',
            values_to = 'depth')


    ## max_q <- max(neon_composite_depth$discharge)
    ## min_q <- min(neon_composite_depth$discharge)
    ## neon_rolling_q <- neon_composite_depth %>%
    ##   mutate(
    ##     roll_depth = (max_q)
    ##   )

    # calculate metric in (m)
    # NOTE: make sure all area is km^2 and all metrics are in cm
    return(neon_composite_depth)
}

# testing
## stream_width_sample <- sample_neon_product(product_codes = stream_width_pc)
## some NEON products are incomplete, and require 'check-and-fill' workflows to
## provide extant data where possible, two such are
##  pressure: provide NOAA gauge pressure where not in NEON
##  Q: provide USGS dataRetrieval Q for TOMB specifically
#

read_all_neon_feathers <- function(file_path, by_site = TRUE){

    if(by_site == TRUE){


        # full paht ot each NEON site directory in data location
        site_dirs <- list.files(file_path, full.names = TRUE)
        # get full path of all files inside all NEON site directories
        neon_files <- list.files(site_dirs, full.names = TRUE)
        # get just file names fpr each of these files
        file_names <- list.files(site_dirs)

        neon_list <- lapply(neon_files, read_feather)
        names(neon_list) <- file_names

        return(neon_list)

    } else{

        neon_files <- list.files(file_path, full.names = TRUE, recursive = TRUE)
        file_names <- str_split_fixed(neon_files, '.feather', n = Inf)[,1]
        file_names <- str_split_fixed(file_names, '[/]', n = Inf)
        inv_file_names <- file_names[,dim(file_names)[2]]
        site_name <- file_names[,dim(file_names)[2]-1]

        inv_file_names <- unique(inv_file_names)
        site_names <- unique(site_name)

        #categoricalCodes
        #readme
        #validation
        #variables

        data_files <-  inv_file_names[!grepl('categoricalCodes|readme|validation|variables', inv_file_names)]

        final_list <- list()
        for(i in 1:length(data_files)){

            all_files <- tibble()
            for(s in 1:length(site_names)){

                path <- glue('{f}/{s}/{i}.feather',
                             f = file_path,
                             i = data_files[i],
                             s = site_names[s])

                one_site <- read_feather(path)

                all_files <- rbind(all_files, one_site)
            }

            all_site_files <- list(all_files)

            names(all_site_files) <- data_files[i]

            final_list <- c(final_list, all_site_files)
        }

        meta_data_files <-  inv_file_names[grepl('categoricalCodes|readme|validation|variables', inv_file_names)]

        for(i in 1:length(meta_data_files)){

            path <- glue('{f}/{s}/{i}.feather',
                         f = file_path,
                         i = meta_data_files[i],
                         s = site_names[1])

            meta_file <- read_feather(path)

            meta_list <- list(meta_file)
            names(meta_list) <- meta_data_files[i]

            final_list <- c(final_list, meta_list)
        }
        return(final_list)
    }

}

read_neon_feathers <- function(file_path, by_site){

    if(by_site == TRUE){

        neon_files <- list.files(file_path, full.names = TRUE)

        file_names <- list.files(file_path)

        neon_list <- lapply(neon_files, read_feather)

        names(neon_list) <- file_names

        return(neon_list)

    } else{

        neon_files <- list.files(file_path, full.names = TRUE, recursive = TRUE)
        file_names <- str_split_fixed(neon_files, '.feather', n = Inf)[,1]
        file_names <- str_split_fixed(file_names, '[/]', n = Inf)
        inv_file_names <- file_names[,dim(file_names)[2]]
        site_name <- file_names[,dim(file_names)[2]-1]

        inv_file_names <- unique(inv_file_names)
        site_names <- unique(site_name)

        #categoricalCodes
        #readme
        #validation
        #variables

        data_files <-  inv_file_names[!grepl('categoricalCodes|readme|validation|variables', inv_file_names)]

        final_list <- list()
        for(i in 1:length(data_files)){

            all_files <- tibble()
            for(s in 1:length(site_names)){

                path <- glue('{f}/{s}/{i}.feather',
                             f = file_path,
                             i = data_files[i],
                             s = site_names[s])

                one_site <- read_feather(path)

                all_files <- rbind(all_files, one_site)
            }

            all_site_files <- list(all_files)

            names(all_site_files) <- data_files[i]

            final_list <- c(final_list, all_site_files)
        }

        meta_data_files <-  inv_file_names[grepl('categoricalCodes|readme|validation|variables', inv_file_names)]

        for(i in 1:length(meta_data_files)){

            path <- glue('{f}/{s}/{i}.feather',
                         f = file_path,
                         i = meta_data_files[i],
                         s = site_names[1])

            meta_file <- read_feather(path)

            meta_list <- list(meta_file)
            names(meta_list) <- meta_data_files[i]

            final_list <- c(final_list, meta_list)
        }
        return(final_list)
    }

}
