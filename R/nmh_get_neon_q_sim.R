nmh_get_neon_q_sim <- function(
    gdrive = "MacroSheds",
    gdrive_filename = 'modeled_discharge_2022-10-13',
    file_ext = '.zip') {

  # NOTE: this all must be replaced with most up to date public data source, a la
  # Vlah et al. 2023 -- https://doi.org/10.5194/egusphere-2023-1178
  # figshare: https://figshare.com/articles/dataset/Composite_discharge_series_for_each_NEON_stream_river_site_/23206592
  
  # first ----
  # download THE MOST RECENT NEON discharge simulations from Google Drive
  # for me: NCSU GDrive/Shared Drives/macroSHEDS/data/Q_simulations
  # save to data folder within RProject
  
  # connect to MS gdrive and locate desired file
  ms_gdrive = googledrive::drive_find(gdrive)
  
  ## ms_gdrive = googledrive::shared_drive_find(gdrive)
  gdrive_file = paste0(gdrive_filename, file_ext)
  ms_sim_q = googledrive::drive_find(gdrive_file)
  
  # download into temp folder, unzip
  temp <- tempfile(fileext = file_ext)
  googledrive::drive_download(ms_sim_q$id, temp, overwrite = TRUE)
  out <- unzip(temp, exdir = tempdir())
  
  print('ATTENTION: printing README:')
  print(readLines(out[1]))
  
  # read data in ----
  this_dir <- file.path(tempdir(), gdrive_filename, 'q_lm_outdata', 'predictions')
  print(paste('loading in site data from', this_dir))
  
  files <- list.files(this_dir, full.names = TRUE)
  sites <- stringr::str_split_fixed(files, pattern = '/',n=5)[,5] %>%
    stringr::str_remove(.,'.csv')
  names(files) = sites
  
  # the map function allows us to read in multiple files at once
  site_discharge_lm <- purrr::map_dfr(
    files,
    readr::read_csv,
    .id = 'site')
  
  ## # read data in ----
  ## this_dir <- file.path(tempdir(), gdrive_filename, 'lstm_outdata', 'predictions')
  ## print(paste('loading in site data from', this_dir))
  
  ## # the map function allows us to read in multiple files at once
  ## site_discharge_lstm <- purrr::map_dfr(
  ##               list.files(this_dir,
  ##                          full.names = TRUE), readr::read_csv)
  
  site_discharge <- site_discharge_lm
  
  return(site_discharge)
}
