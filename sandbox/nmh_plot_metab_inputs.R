nmh_plot_metab_inputs <- function(input_dir = 'data/sm_input/') {
  
  q_types <- list.files(input_dir)
  
  list <- list()
  for(i in 1:length(q_types)) {
    
    type_use <- q_types[i]
    
    file_dir <- glue::glue(input_dir, type_use)
    
    files <- list.files(file_dir, full.names = TRUE)
    
    list[[i]] <- map_dfr(files,
                         readr::read_csv,
                         .id = 'type_use') 
  }  # end for loop
  
  d <- unlist(list)

  plot <- ggplot2::ggplot(d %>%
                            tdiyr::pivot_longer(cols = 2:7,
                                                values_to = 'value',
                                                names_to = 'vars'),
                          aes(x = solar.time,
                              y = value,
                              color = type))+
    ggplot2::geom_point(alpha = 0.5)+
    ggplot2::facet_grid(vars ~ type,
                        scales = 'free_y')+
    ggplot2::ggtitle(site)
  
  save_dir <- 'figures/sm_input/'
  if(!dir.exists(save_dir))
    dir.create(save_dir)
  
  # TODO: check site name saving convention
  ggplot2::ggsave(plot,
                  glue::glue(save_dir, '/{site}_sm_input.png'))
  
} # end function