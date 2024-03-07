rad2deg <- function(x){
  return(x * 180 / pi)
}

full_days <- function(df) {
  df_sum <- df %>%
    dplyr::group_by(date = lubridate::date(solar.time)) %>%
    dplyr::summarise(obs = dplyr::n(),
                     min = min(solar.time),
                     max = max(solar.time))
  
  bad_days <- df_sum %>%
    dplyr::filter(obs == 97) %>%
    dplyr::pull(max)
  
  df_out <- df %>%
    dplyr::filter(!solar.time %in% bad_days) %>%
    dplyr::mutate(date = lubridate::date(solar.time)) %>%
    dplyr::select(-date)
  return(df_out)
}

serialize_list_to_dir <- function(l, dest){
  
  #l must be a named list
  #dest is the path to a directory that will be created if it doesn't exist
  
  #list element classes currently handled: data.frame, character
  
  elemclasses = lapply(l, class)
  
  handled = lapply(elemclasses,
                   function(x) any(c('character', 'data.frame') %in% x))
  
  if(! all(unlist(handled))){
    stop('Unhandled class encountered')
  }
  
  dir.create(dest, showWarnings=FALSE, recursive=TRUE)
  
  for(i in 1:length(l)){
    
    if('data.frame' %in% elemclasses[[i]]){
      
      fpath = paste0(dest, '/', names(l)[i], '.feather')
      feather::write_feather(l[[i]], fpath)
      
    } else if('character' %in% elemclasses[[i]]){
      
      fpath = paste0(dest, '/', names(l)[i], '.txt')
      readr::write_file(l[[i]], fpath)
    }
  }
}

read_neon_feathers <- function(file_path, by_site){
  
  if(by_site == TRUE){
    
    neon_files <- list.files(file_path, full.names = TRUE)
    
    file_names <- list.files(file_path)
    
    neon_list <- lapply(neon_files, feather::read_feather)
    
    names(neon_list) <- file_names
    
    return(neon_list)
    
  } else{
    
    neon_files <- list.files(file_path, full.names = TRUE, recursive = TRUE)
    file_names <- stringr::str_split_fixed(neon_files, '.feather', n = Inf)[,1]
    file_names <- stringr::str_split_fixed(file_names, '[/]', n = Inf)
    inv_file_names <- file_names[,dim(file_names)[2]]
    site_name <- file_names[,dim(file_names)[2]-1]
    
    inv_file_names <- unique(inv_file_names)
    site_names <- unique(site_name)
    
    data_files <- inv_file_names[!grepl('categoricalCodes|readme|validation|variables', inv_file_names)]
    
    final_list <- list()
    for(i in 1:length(data_files)){
      
      all_files <- tibble()
      for(s in 1:length(site_names)){
        
        path <- glue::glue('{f}/{s}/{i}.feather',
                           f = file_path,
                           i = data_files[i],
                           s = site_names[s])
        
        one_site <- feather::read_feather(path)
        
        all_files <- rbind(all_files, one_site)
      }
      
      all_site_files <- list(all_files)
      
      names(all_site_files) <- data_files[i]
      
      final_list <- c(final_list, all_site_files)
    }
    
    meta_data_files <- inv_file_names[grepl('categoricalCodes|readme|validation|variables', inv_file_names)]
    
    for(i in 1:length(meta_data_files)){
      
      path <- glue::glue('{f}/{s}/{i}.feather',
                         f = file_path,
                         i = meta_data_files[i],
                         s = site_names[1])
      
      meta_file <- feather::read_feather(path)
      
      meta_list <- list(meta_file)
      names(meta_list) <- meta_data_files[i]
      
      final_list <- c(final_list, meta_list)
    }
    return(final_list)
  }
}

get_avail_neon_product_sets <- function(prodcode_full) {
  
  # returns a tibble with url, site_name, component columns
  
  avail_sets = tibble::tibble()
  
  req = httr::GET(paste0("http://data.neonscience.org/api/v0/products/",
                         prodcode_full))
  
  txt = httr::content(req, as="text")
  
  neondata = jsonlite::fromJSON(txt, simplifyDataFrame=TRUE, flatten=TRUE)
  
  urls = unlist(neondata$data$siteCodes$availableDataUrls)
  
  avail_sets = stringr::str_match(urls, '(?:.*)/([A-Z]{4})/([0-9]{4}-[0-9]{2})') %>%
    tibble::as_tibble(.name_repair='unique') %>%
    dplyr::rename(url=`...1`, site_name=`...2`, component=`...3`)
  
  return(avail_sets)
}

FindandCollect_airpres = function(lat, lon, start_datetime, end_datetime) {
  #get df of all available air pressure stations
  tf = tempfile()
  utils::download.file("ftp://ftp.ncdc.noaa.gov/pub/data/noaa/isd-history.txt", tf, mode="wb")
  noaa.sites <- utils::read.fwf(tf, skip = 22, header = F,
                                widths = c(6,-1,5,-45, 8, 9,-8, 9, 8), 
                                comment.char = "",
                                col.names = c("USAF", "WBAN", "LAT", "LON", "BEGIN", "END"),
                                flush = TRUE, 
                                colClasses=c('USAF'='character', 'WBAN'='character'))
  noaa.sites <- na.omit(noaa.sites)
  
  # narrow them down to those within 5 lats/longs
  noaa.sites <- noaa.sites %>%
    dplyr::mutate(LAT = as.numeric(as.character(LAT))) %>%
    dplyr::mutate(LON = as.numeric(as.character(LON))) %>%
    dplyr::filter(LAT < (lat + 5) & LAT > (lat - 5) & LON < (lon + 5) & LON > (lon - 5))
  
  # filter by coverage, order by distance
  pt1 <- cbind(rep(lon, length.out = length(noaa.sites$LAT)),
               rep(lat, length.out = length(noaa.sites$LAT)))
  pt2 <- cbind(noaa.sites$LON, noaa.sites$LAT)
  dist <- diag(geosphere::distm(pt1, pt2, fun=geosphere::distHaversine))/1000
  noaa.sites$dist <- dist
  tmp <- which((as.numeric(substr(noaa.sites$END,1,4)) >=
                  as.numeric(substr(end_datetime, 1, 4))) &
                 as.numeric(substr(noaa.sites$BEGIN,1,4)) <=
                 as.numeric(substr(start_datetime, 1, 4)))
  noaa.sites <- noaa.sites[tmp,]
  noaa.sites <- noaa.sites[with(noaa.sites, order(dist)),]
  yrs <- seq(as.numeric(substr(start_datetime, 1, 4)),
             as.numeric(substr(end_datetime, 1, 4)), 
             by = 1)
  
  for (i in 1:length(noaa.sites$dist)) {
    k <- i
    available <- vector(mode = 'logical', length = length(yrs))
    USAF <- as.character(noaa.sites$USAF[i])
    if(nchar(as.character(noaa.sites$WBAN[i])) == 5){
      WBAN <- as.character(noaa.sites$WBAN[i])
    } else {
      WBAN <- paste0(0,as.character(noaa.sites$WBAN[i]))
    }
    
    y <- as.data.frame(matrix(NA, nrow = 1, ncol = 12))
    
    for(j in 1:length(yrs)){
      tf = tempfile()
      res = tryCatch(suppressWarnings(download.file(paste0("ftp://ftp.ncdc.noaa.gov/pub/data/noaa/isd-lite/",
                                                           yrs[j], "/", USAF, "-", WBAN, "-", yrs[j], ".gz"), tf, mode="wb")),
                     error=function(e){
                       # message('NCDC download failed; trying next closest station')
                       return('download failed')
                     })
      if(exists('res') && res == 'download failed'){
        break #try next station
      }
      x = read.table(tf)
      x[x==-9999] = NA
      if(length(which(!is.na(x$V7))) >= 0.9 * length(x$V7)) {
        available[j] <- TRUE
        y <- rbind(x,y)
      } else {
        break #too many NAs, move to next station
      }
    }
    if(length(yrs) == length(which(available))){
      break #got one
    }
  }
  y <- y[!is.na(y$V1),]
  colnames(y) = c("y","m","d","h","air_temp","dewtemp","air_kPa","winddir","sindspeed","skycover","precip1h","precip6h")
  y$air_kPa = y$air_kPa/100
  y$air_temp = y$air_temp/10
  y$DateTime_UTC = readr::parse_datetime(paste0(y$y,"-",
                                                sprintf("%02d",y$m),"-",sprintf("%02d",y$d)," ",sprintf("%02d",y$h),
                                                ":00:00"), "%F %T")
  y <- y[with(y, order(DateTime_UTC)),]
  y = tibble::as_tibble(y) %>% select(DateTime_UTC,air_temp,air_kPa)
  ss = tibble::tibble(DateTime_UTC=seq(y$DateTime_UTC[1],
                                       y$DateTime_UTC[nrow(y)], by=900))
  xx = dplyr::left_join(ss, y, by = "DateTime_UTC")
  xx = dplyr::mutate(xx, 
                     air_temp=zoo::na.approx(air_temp),
                     air_kPa=zoo::na.approx(air_kPa))
  daterng = c(start_datetime, end_datetime)
  xtmp = xx %>%
    dplyr::filter(DateTime_UTC>=daterng[1] & DateTime_UTC<=daterng[2]) %>%
    dplyr::mutate(air_mb = air_kPa*10) # convert to mbar
  
  return(select(xtmp, DateTime_UTC, air_mb, air_temp))
}

prep_daily_discharge <- function(site_data, domain = 'neon') {
  # pull in site data
  site_code <- readr::read_csv('data/site_data/site_data.csv') %>%
    # ensure domain is neon
    dplyr::filter(domain == 'neon') %>%
    # extract site codes
    dplyr::pull(site_code)
  
  # loop through data
  for(i in 1:length(site_code)){
    
    site <- site_code[i]
    
    d <- try(read_csv(glue('data/discharge/Q_simulations/simulated_discharge_2022-01-25/Q_by_site/{site}.csv')) %>%
               dplyr::filter(year(date) > 2014) %>%
               dplyr::mutate(discharge_m3s_sim = discharge_Ls_sim/1000) %>%
               dplyr::select(site_code, date, discharge_m3s_sim))
    
    if(inherits(d, 'try-error')) {
      next
    }
    
    readr::write_csv(d,
                     glue('data/sm_ready_dailyQ/{site}_daily_Q.csv'))
  }
}

get_neon_product <- function(product_codes = 'DP0.20288.001',
                             dest_fp = NULL, 
                             site_filter = NULL) {
  
  # checking for user defined filepath
  if(is.null(dest_fp)) {
    dest_fp <- file.path(getwd(), 'data', 'raw')
  }
  
  if(!dir.exists(dest_fp)){
    # create direcotry if doesnt exists
    print(paste('Directory does not exist, creating directory here:', dest_fp))
    dir.create(dest_fp)
  }
  
  for(product_code in product_codes) {
    # get NEON product name, used for filepath
    product_name <- neonUtilities::getProductInfo(product_code)$productName
    data_fp <- file.path(dest_fp, product_name)
    products_n = length(product_codes)
    
    if(!dir.exists(data_fp)){
      # create direcotry if doesnt exists
      print(paste('Directory does not exist, creating directory here:', data_fp))
      dir.create(data_fp)
    }
    
    writeLines(paste('retrieving NEON data for', products_n, 'data products',
                     'and saving results at:\n', data_fp))
    
    # call to NEON for avilable data of code type
    avail_sets <- get_avail_neon_product_sets(product_code)
    
    # sites form this record
    avail_sites <- unique(avail_sets$site_name)
    
    if(!is.null(site_filter)) {
      print('values suppleid to site_filter argument, filtering queryable sites to those in site_filter values')
      avail_sites <- avail_sites[avail_sites %in% site_filter]
    }
    
    avail_sites_n <- length(avail_sites)
    if(avail_sites_n == 0) {
      stop('no available sites')
    }
    
    # download product for each site
    writeLines(paste('\nquerying for NEON product code:', product_code,
                     '\n  total sites to query:', avail_sites_n))
    for(j in 1:length(avail_sites)){
      writeLines(paste('  querying site', j, 'of', avail_sites_n))
      
      # Defines what site to retrieve
      site_name <- avail_sites[j]
      writeLines(paste('  site name:', site_name))
      
      # Download data product for the site
      tryCatch(
        expr = {
          data_pile <- neonUtilities::loadByProduct(dpID = product_code,
                                                    site = site_name,
                                                    package='basic',
                                                    check.size=FALSE)
        },
        error = function(e) {
          writeLines(paste('ERROR at', site_name, 'trying again with expanded package'))
          data_pile <- try(neonUtilities::loadByProduct(dpID = product_code,
                                                        site = site_name,
                                                        package='expanded',
                                                        check.size=FALSE))
        }
      )
      
      if(inherits(data_pile, 'try-error')) {
        writeLines(paste('no data downloaded skipping', site_name))
        next
      } else {
        # Create the file path for where the files will be saved. This can be changed
        # if you  want to save the files in a different location
        raw_data_dest <- file.path(data_fp, site_name)
        
        # Save the list of all dataframes for the site as individual feather files
        # in the file path raw_data_dest
        serialize_list_to_dir(data_pile, raw_data_dest)
      }
    }
  }
}

sample_neon_product <- function(product_codes = 'DP0.20288.001', product_name = 'neon_data',
                                dest_fp = NULL, site_filter = NULL) {
  # checking for user defined filepath
  if(is.null(dest_fp)) {
    dest_fp <- file.path(getwd(), 'data', 'raw')
  }
  
  data_fp <- file.path(dest_fp, product_name)
  products_n = length(product_codes)
  
  writeLines(paste('retrieving NEON data for', products_n, 'data products',
                   'and saving results at:\n', data_fp))
  
  for(product_code in first(product_codes)) {
    
    # call to NEON for avilable data of code type
    avail_sets <- get_avail_neon_product_sets(product_code)[1,]
    
    # sites form this record
    avail_sites <- unique(avail_sets$site_name)
    
    if(!is.null(site_filter)) {
      print('values suppleid to site_filter argument, filtering queryable sites to those in site_filter values')
      avail_sites <- avail_sites[avail_sites$site_name %in% site_filter,]
    }
    
    avail_sites_n <- length(avail_sites)
    if(avail_sites_n == 0) {
      stop('no available sites')
    }
    
    # download product for each site
    writeLines(paste('\nquerying for NEON product code:', product_code,
                     '\n  total sites to query:', avail_sites_n))
    for(j in 1:length(avail_sites)){
      writeLines(paste('  querying site', j, 'of', avail_sites_n))
      
      # Defines what site to retrieve
      site_name <- avail_sites[j]
      
      # Download data product for the site
      data_pile <- neonUtilities::loadByProduct(dpID = product_code,
                                                site = site_name,
                                                package='basic',
                                                check.size=FALSE)
      
      # Create the file path for where the files will be saved. This can be changed
      # if you  want to save the files in a different location
      ## raw_data_dest <- file.path(data_fp, site_name)
      # Save the list of all dataframes for the site as individual feather files
      # in the file path raw_data_dest
      ## serialize_list_to_dir(data_pile, raw_data_dest)
      
      return(data_pile)
    }
  }
}


read_all_neon_feathers <- function(file_path, by_site = TRUE){
  
  if(by_site == TRUE){
    
    # full paht ot each NEON site directory in data location
    site_dirs <- list.files(file_path, full.names = TRUE)
    
    # get full path of all files inside all NEON site directories
    neon_files <- list.files(site_dirs, full.names = TRUE)
    
    # get just file names fpr each of these files
    file_names <- list.files(site_dirs)
    
    neon_list <- lapply(neon_files, feather::read_feather)
    names(neon_list) <- file_names
    
    return(neon_list)
    
  } else {
    
    neon_files <- list.files(file_path, full.names = TRUE, recursive = TRUE)
    file_names <- stringr::str_split_fixed(neon_files, '.feather', n = Inf)[,1]
    file_names <- stringr::str_split_fixed(file_names, '[/]', n = Inf)
    inv_file_names <- file_names[,dim(file_names)[2]]
    site_name <- file_names[,dim(file_names)[2]-1]
    
    inv_file_names <- unique(inv_file_names)
    site_names <- unique(site_name)
    
    data_files <-  inv_file_names[!grepl('categoricalCodes|readme|validation|variables', inv_file_names)]
    
    final_list <- list()
    for(i in 1:length(data_files)){
      
      all_files <- tibble()
      for(s in 1:length(site_names)){
        
        path <- glue('{f}/{s}/{i}.feather',
                     f = file_path,
                     i = data_files[i],
                     s = site_names[s])
        
        one_site <- feather::read_feather(path)
        
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
      
      meta_file <- feather::read_feather(path)
      
      meta_list <- list(meta_file)
      names(meta_list) <- meta_data_files[i]
      
      final_list <- c(final_list, meta_list)
    }
    return(final_list)
  }
  
}

get_streampulse_data <- function(site_deets, 
                                 site_data,
                                 dest_fp = NULL) {
  
  # checking for user defined filepath
  if(is.null(dest_fp)) {
    dest_fp <- file.path(getwd(), 'data', 'raw', 'streampulse_2022')
  }
  
  if(!dir.exists(dest_fp)){
    # create direcotry if doesnt exists
    print(paste('Directory does not exist, creating directory here:', dest_fp))
    dir.create(dest_fp)
  }
  
  # loop to download and save data from streampulse ----
  for(i in 1:nrow(site_deets)) {
    
    # define site codes and inputs parameters
    site_code = dplyr::pull(site_deets, site_code)[i]
    sp_code = dplyr::pull(site_deets, sp_code)[i]
    
    # define lat and long, based on NEON site data
    long <- site_data %>%
      dplyr::filter(site_code == !!site_code) %>%
      dplyr::pull(field_longitude)
    
    lat <- site_data %>%
      dplyr::filter(site_code == !!site_code) %>%
      dplyr::pull(field_latitude)
    
    # request data from streampulse
    # each NEON site has temperature and DO data at 15 minute intervals
    streampulse_data <- try({
      StreamPULSE::request_data(sitecode = sp_code)
    })
    
    print(paste('downloading', site_code, 'data from streampulse'))
    file_path <- glue::glue('data/raw/streampulse_2022/{site_code}_2022.csv')
    
    readr::write_csv(streampulse_data$data, file_path)
  } # end for loop
  
} # end function



pkg_namespace_require <- function(pkg = 'macrosheds', pkg_gh = "https://github.com/MacroSHEDS/macrosheds.git", quietly = FALSE) {
  res <- try(requireNamespace(pkg, quietly = quietly))
  if(res) {
    library(macrosheds)
  } else {
    if(!quietly) {
      if(menu(c("Yes", "No"),
              title= paste("Are you sure you want to install package", pkg)) == "1") {
        inst <- try(devtools::install_github(pkg_gh))
        
        if(class('inst') == 'try-error') {
          warning('devtools::install_github did not work, trying remotes::install_github')
          inst <- try(remotes::install_github(pkg_gh))
        }
        
      } else {
        warning("no installation of required package")
      }
    }
  }
}


calc_DO_sat <- function (temp.water, pressure.air, salinity.water = u(0, "PSU"), 
                         model = "garcia-benson", ...) {
  with.units <- any(sapply(list(temp.water, pressure.air), 
                           is.unitted)) || (if (!missing(salinity.water)) 
                             is.unitted(salinity.water)
                             else FALSE)
  if (with.units) {
    verify_units(temp.water, "degC")
    verify_units(pressure.air, "mb")
    verify_units(salinity.water, "PSU")
  }
  temp.water <- v(temp.water)
  pressure.air <- v(pressure.air)
  salinity.water <- v(salinity.water)
  o2.at.sat <- LakeMetabolizer::o2.at.sat.base(temp = temp.water, 
                                               baro = pressure.air, salinity = salinity.water, model = model, 
                                               ...)
  if (with.units) {
    return(u(o2.at.sat, "mgO2 L^-1"))
  }
  else {
    return(o2.at.sat)
  }
}
