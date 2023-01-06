FindandCollect_airpres = function(lat, lon, start_datetime, end_datetime) {
  #get df of all available air pressure stations
  tf = tempfile()
  download.file("ftp://ftp.ncdc.noaa.gov/pub/data/noaa/isd-history.txt", tf, mode="wb")
  noaa.sites <- read.fwf(tf, skip = 22, header = F,
                         # widths = c(6,-1,5,-1,30, 5, 3, 6, 8, 9, 8, 9, 8), comment.char = "",
                         widths = c(6,-1,5,-45, 8, 9,-8, 9, 8), comment.char = "",
                         col.names = c("USAF", "WBAN", "LAT", "LON", "BEGIN", "END"),
                         # col.names = c("USAF", "WBAN", "STATION NAME", "CTRY", "ST", "CALL", "LAT", "LON", "ELEV(M)", "BEGIN", "END"),
                         flush = TRUE, colClasses=c('USAF'='character', 'WBAN'='character'))
  noaa.sites <- na.omit(noaa.sites)
  #narrow them down to those within 5 lats/longs
  noaa.sites <- noaa.sites %>%
    mutate(LAT = as.numeric(as.character(LAT))) %>%
    mutate(LON = as.numeric(as.character(LON))) %>%
    filter(LAT < (lat + 5) & LAT > (lat - 5) & LON < (lon + 5) & LON > (lon - 5))
  #filter by coverage, order by distance
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
  yrs <- seq(as.numeric(substr(start_datetime, 1, 4)),as.numeric(substr(end_datetime, 1, 4)), by = 1)
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
      }else {
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
  xx = left_join(ss, y, by = "DateTime_UTC")
  xx = mutate(xx, air_temp=zoo::na.approx(air_temp),
              air_kPa=zoo::na.approx(air_kPa))
  daterng = c(start_datetime, end_datetime)
  xtmp = xx %>%
    filter(DateTime_UTC>=daterng[1] & DateTime_UTC<=daterng[2]) %>%
    mutate(air_mb = air_kPa*10) # convert to mbar
  # select(xtmp, DateTime_UTC, air_kPa, air_temp)
  # print(noaa.sites[k,])
  return(select(xtmp, DateTime_UTC, air_mb, air_temp))
}

get_streampulse_data <- function(site_deets, site_data,
                                 dest_fp = NULL) {

    # checking for user defined filepath
    if(is.null(dest_fp)) {
        dest_fp <- file.path(getwd(), 'data', 'raw', 'streampulse')
    }

    if(!dir.exists(dest_fp)){
        # create direcotry if doesnt exists
        print(paste('Directory does not exist, creating directory here:', dest_fp))
        dir.create(dest_fp)
    }

    # loop to download and save data from streampulse ----
    for(i in 1:nrow(site_deets)) {

        # define site codes and inputs parameters
        site_code = pull(site_deets, site_code)[i]
        sp_code = pull(site_deets, sp_code)[i]

        # define lat and long, based on NEON site data
        long <- site_data %>%
            dplyr::filter(site_code == !!site_code) %>%
            dplyr::pull(longitude)

        lat <- site_data %>%
            dplyr::filter(site_code == !!site_code) %>%
            dplyr::pull(latitude)

        # request data from streampulse
        # each NEON site has temperature and DO data at 15 minute intervals
        streampulse_data = try({
            StreamPULSE::request_data(sitecode = sp_code)
        })

        print(paste('downloading', site_code, 'data from streampulse'))
        file_path <- glue('data/raw/streampulse/{site_code}.csv')

        # if there is an error, let us know
        # if(inherits(streampulse_data, 'try-error')){
        #   results[i, 2] = 'request error'
        #   next
        # }

        write_csv(streampulse_data$data, file_path)
    }
}
