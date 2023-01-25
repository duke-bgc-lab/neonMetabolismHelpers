#' @export
nmh_get_neon_q_eval <- function(dest_file = 'data/raw/qaqc/', dest_fn = 'neon_q_eval.csv', download = FALSE) {
  if(download) {
    destfile <- paste0(dest_fp, dest_fn)


    if(!dir.exists(dest_fp)){
      # create direcotry if doesnt exists
      print(paste('Directory does not exist, creating dir at:', dest_fp))
      dir.create(dest_fp, recursive = TRUE)
    }

    # read in the discharge evaluation from Rhea et al. (in prep)
    download.file(
      url = 'https://www.hydroshare.org/resource/03c52d47d66e40f4854da8397c7d9668/data/contents/neon_q_eval.csv',
      destfile = destfile
    )

    q_eval <- read.csv(destfile)
    return(q_eval)

  } else {
    return(nmh_q_eval)
  }
}
