#' Download simulated NEON discharge data from Vlah et al. 2023 -- https://doi.org/10.5194/egusphere-2023-1178
#'
#' @author Nick Marzolf, \email{@@}
#' @author Wes Slaughter, \email{wslaughter@@berkeley.edu}
#'
#' @param gdrive character.
#' @details
#' @seealso [nmh_get_hydraulic_coef_data()], [nmh_get_neon_data()], [nmh_get_neon_q_eval()], [nmh_get_scaling_coefs()], [nmh_get_tomb_q()]
#' @examples
#' nmh_get_neon_q_sim(
#' )
#' @export
nmh_get_neon_q_sim_figshare <- function(dest_fp = 'data/raw/macrosheds/') {
  
  if(!dir.exists(dest_fp)){
    dir.create(dest_fp)
  }
  
  download.file(url = 'https://figshare.com/ndownloader/files/40935107/composite_series.zip', 
                destfile = glue::glue(dest_fp, 'neon_Q_sims.zip'))
  
  
  sim_dir <- glue::glue(dest_fp,'simulated/')
  
  if(!dir.exists(sim_dir)){
    dir.create(sim_dir)
  }
  
  sim_files <- unzip(zipfile = glue::glue(dest_fp, 'neon_Q_sims.zip'),
                     exdir = sim_dir)
  
  sims <- sim_files %>% 
    rlang::set_names() %>% 
    purrr::map_dfr(.,
                   readr::read_csv,
                   .id = 'site') %>% 
    dplyr::mutate(site = substr(site, 49,52))
  
  return(sims)
  
}
