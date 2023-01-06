diagnose_neon_metab_bayes <- function(mod_dir = 'data/model_runs/',
                                      type = c('raw', 'source', 'simulation')) {

    library(tidyverse)
    library(dplyr)
    library(ggplot2)
    library(glue)

    # add if statement to check for 'type'

    mod_dir_type <- glue::glue(mod_dir, type)
    mod_dir_fits <- glue::glue(mod_dir_type,'/Bayes/daily')

    n_mods <- list.files(mod_dir_fits,
                         full.names = TRUE)

    out <- data.frame(site = character(),
                      year = character(),
                      n_days = numeric(),
                      f_days = numeric(),
                      resolution = character(),
                      K600_daily_sigma_Rhat = numeric(),
                      err_obs_iid_sigma_Rhat = numeric(),
                      err_proc_iid_sigma_Rhat = numeric(),
                      K_median = numeric(),
                      K_range = numeric(),
                      neg_GPP = numeric(),
                      pos_ER = numeric(),
                      ER_K_r2 = numeric())

    for(i in 1:length(n_mods)) {

        file <- list.files(mod_dir_fits,
                           full.names = TRUE)[i]

        d <- read_csv(file)

        # TODO: check these
        site <- str_split(file, c('_','/'))[[2]][6] %>%
            substr(.,1,4)

        start_date <- str_split(file, c('_','/'))[[2]][6] %>%
            substr(.,6,15)
        end_date <- str_split(file, c('_','/'))[[2]][6] %>%
            substr(.,17,26)

        year <- substr(start_date, 1,4)

        if('GPP_50pct' %in% names(d)) {
            days <- d %>%
                summarise(n = length(GPP_50pct)) %>%
                pull()

            GPP_neg <- d %>%
                filter(GPP_50pct < -0.5) %>%
                summarise(n = length(GPP_50pct),
                          GPP_neg_per = (n/days)*100) %>%
                pull(GPP_neg_per)

            ER_pos <- d %>%
                filter(ER_50pct > 0.5) %>%
                summarise(n = length(ER_50pct),
                          ER_pos_per = (n/days)*100) %>%
                pull(ER_pos_per)

            K <- d %>%
                summarise(K_median = median(K600_daily_50pct, na.rm = TRUE),
                          K_min = min(K600_daily_50pct, na.rm = TRUE),
                          K_max = max(K600_daily_50pct, na.rm = TRUE)) %>%
                mutate(K_range = K_max - K_min)

            ER_K_r2 <- summary(lm(data = d,
                                  ER_50pct ~ K600_daily_50pct))$adj.r.squared
        } else  # end if statement
            next

        dir_DO <- glue::glue(mod_dir_type,'/Bayes/mod_and_obs_DO')

        DO <- try(read_csv(glue(dir_DO, '/', site, '_', start_date, '_', end_date, '_mod_and_obs_DO.csv')))

        if(inherits(DO, 'try-error')){
            next
        } # end if statement

        res <- glue(diff(DO$solar.time) %>%
                        first(),
                    'min')

        dir_KQ <- glue::glue(mod_dir_type,'/Bayes/KQ_overall')

        KQ <- try(read_csv(glue(dir_KQ, '/', site, '_', start_date, '_', end_date, '_KQ_overall.csv')))

        if(inherits(KQ, 'try-error')){
            next
        } # end if statement

        K600_daily_sigma_Rhat <- KQ %>%
            pull(K600_daily_sigma_Rhat)

        dir_overall <- glue::glue(mod_dir_type,'/Bayes/overall')

        overall <- try(read_csv(glue(dir_overall, '/', site, '_', start_date, '_', end_date, '_overall.csv')))

        if(inherits(overall, 'try-error')){
            next
        } # end if statement

        err_proc_iid_sigma_Rhat <- overall %>%
            pull(err_proc_iid_sigma_Rhat)
        err_obs_iid_sigma_Rhat <- overall %>%
            pull(err_obs_iid_sigma_Rhat)

        out <- out %>%
            add_row(site = site,
                    year = year,
                    n_days = days,
                    f_days = (days/365)*100,
                    resolution = res,
                    K600_daily_sigma_Rhat = K600_daily_sigma_Rhat,
                    err_obs_iid_sigma_Rhat = err_obs_iid_sigma_Rhat,
                    err_proc_iid_sigma_Rhat = err_proc_iid_sigma_Rhat,
                    K_median = K$K_median,
                    K_range = K$K_range,
                    neg_GPP = GPP_neg,
                    pos_ER = ER_pos,
                    ER_K_r2 = ER_K_r2)
    } # end for loop

    write_csv(out,
              glue::glue(mod_dir_type,'/Bayes/diag_neon_metab_{type}.csv'))

    return(out)
}
