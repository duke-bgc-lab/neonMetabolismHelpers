
get_neon_site_data <- function(download = TRUE) {
  
  if(download){
    dest_dir <- 'data/'
    dest_fn <- '1_neon_site_data.csv'
    
    destfile <- file.path(dest_dir, dest_fn)
    
    if(!dir.exists(dest_dir)){
      dir.create(dest_dir)
    }
    
    download.file(url = 'https://www.neonscience.org/sites/default/files/NEON_Field_Site_Metadata_20231026.csv',
                  destfile = destfile)
    
    site_data <- readr::read_csv(destfile)
    
    table_1 <- site_data %>% 
      dplyr::filter(field_site_subtype %in% c('Wadeable Stream', 'Non-wadeable River')) %>% 
      dplyr::select(Domain = field_domain_id,
                    `Site Code` = field_site_id,
                    `Site Name` = field_site_name,
                    Latitude = field_latitude,
                    Longitude = field_longitude,
                    Elevation = field_mean_elevation_m,
                    `Watershed Area (km2)` = field_watershed_size_km2) %>% 
      sf::st_as_sf(.,
                   coords = c('Longitude', 'Latitude'),
                   crs = 4326)
    
    readr::write_csv(table_1,
                     'data/1_neon_site_data.csv')
    
  }  
}
