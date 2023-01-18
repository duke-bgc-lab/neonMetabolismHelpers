nmh_get_scaling_coefs <- function(dest_fp = 'data/raw/macrosheds/', 
                                  dest_fn = 'neon_hyd_scaling_coefs.csv') {
  destfile <- paste0(dest_fp, dest_fn)
  
  
  if(!dir.exists(dest_fp)){
    # create direcotry if doesnt exists
    print(paste('Directory does not exist, creating dir at:', dest_fp))
    dir.create(dest_fp)
  }
  
  
  # read in the discharge evaluation from Rhea et al. (in prep)
  download.file(
    url = 'https://raw.githubusercontent.com/nmarzolf91/NEON-reaeration/master/data/derived/scaling_coefs/NEON_site_scaling_coefs.csv',
    destfile = destfile
  )
  
  coefs <- readr::read_csv(destfile)
  
  return(coefs)
}

