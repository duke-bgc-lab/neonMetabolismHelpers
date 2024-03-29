#' Download hydraulic scaling equations determined for NEON wadable sites from my personal github repo
#' 
#' @author Nick Marzolf, \email{nick.marzolf@@jonesctr.org}
#' @author Wes Slaughter, \email{wslaughter@@berkeley.edu}

nmh_get_scaling_coefs <- function(dest_fp = 'data/raw/macrosheds/', 
                                  dest_fn = 'neon_hyd_scaling_coefs.csv') {
  
  destfile <- paste0(dest_fp, dest_fn)
  
  if(!dir.exists(dest_fp)){
    # create direcotry if doesnt exists
    print(paste('Directory does not exist, creating dir at:', dest_fp))
    dir.create(dest_fp)
  }
  
  
  # NOTE: applies to this and other similar function (1) only room in this town for one nmh_get_scaling_coef?
  # and (2) should get this data from best available public src associated w publication, or, hosted on figshare with
  # consent of original authors
  download.file(
    url = 'https://raw.githubusercontent.com/nmarzolf91/NEON-reaeration/master/data/derived/scaling_coefs/NEON_site_scaling_coefs.csv',
    destfile = destfile
  )
  
  coefs <- readr::read_csv(destfile)
  
  return(coefs)
}

