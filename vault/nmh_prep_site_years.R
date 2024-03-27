#' @export
nmh_prep_site_years <- function(years = c('2016', '2021'),
                                site = 'all') {
  
  if(site[1] == 'all') {
    site_codes <- get_neon_site_data()[,'site_code']
  } else {
    site_codes <- unique(site)
  }
  
  if(length(years) > 1) {
    year <- rep(years[1]:years[2], 
                times = length(site_codes)) 
  } else {
    year <- rep(years[1],
                times = length(site_codes))
  }
  
  site_years <- data.frame(site_code = rep(unique(site_codes), 
                                           each = length(unique(year))),
                           year = year)
  
} # end function
