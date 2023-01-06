eval_neon_metab_mle <- function(dir = 'data/model_runs/',
                                type = c('raw', 'source', 'simulation')) {

    library(tidyverse)
    library(dplyr)
    library(ggplot2)
    library(glue)
    library(streamMetabolizer)

    dir_type <- glue(dir, type)

    dir_fits <- glue(dir_type, '/MLE/dat_fit')

    fits <- list.files(dir_fits,
                       full.names = TRUE)

    mle_diagnostic <- data.frame(
        site = character(),
        year = character(),
        GPP_neg = numeric(),
        ER_pos = numeric(),
        K_median = numeric(),
        K_sd = numeric(),
        no_dat_num = numeric(),
        no_dat_dates = character()
    )

    for(i in 1:length(fits)) {

        fit <- read_csv(fits[i])

        site <- str_split(fits[i], '/')[[1]][6] %>%
            substr(.,1,4)
        year <- str_split(fits[i], '/')[[1]][6] %>%
            substr(.,6,9)

        GPP_neg <- fit %>%
            filter(GPP.daily < -0.5) %>%
            summarise(n = n()) %>%
            pull(n)

        ER_pos <- fit %>%
            filter(ER.daily > 0.5) %>%
            summarise(n = n()) %>%
            pull(n)

        K_sum <- fit %>%
            filter(K600.daily > 0) %>%
            summarise(medianK = median(K600.daily, na.rm = TRUE),
                      sdK = sd(K600.daily, na.rm = TRUE))

        no_dat_num <- sum(is.na(fit$GPP.daily))

        no_dat_dates <- fit %>%
            filter(is.na(GPP.daily)) %>%
            pull(date)

        mle_diagnostic <- mle_diagnostic %>%
            add_row(
                site = site,
                year = year,
                GPP_neg = GPP_neg,
                ER_pos = ER_pos,
                K_median = pull(K_sum,medianK),
                K_sd = pull(K_sum,sdK),
                no_dat_num = no_dat_num,
                no_dat_dates = str_c(as.character(no_dat_dates), collapse = ' ')
                )
    } # end for loop

    write_csv(mle_diagnostic,
              glue(dir_type, '/MLE/MLE_{type}_diagnostics.csv'))
} # end function
