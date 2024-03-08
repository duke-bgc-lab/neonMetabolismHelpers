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