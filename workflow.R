
remotes::install_github('MacroSHEDS/macrosheds')
library(macrosheds)
library(unitted)
source('R/nmh_internals.R')

# get data ----
source('R/nmh_get_neon_data.R')

nmh_get_neon_data(product_codes = 'all',
                  quietly = TRUE,
                  check.size = FALSE,
                  neon_api_token = 'eyJ0eXAiOiJKV1QiLCJhbGciOiJFUzI1NiJ9.eyJhdWQiOiJodHRwczovL2RhdGEubmVvbnNjaWVuY2Uub3JnL2FwaS92MC8iLCJzdWIiOiJtYWNyb3NoZWRzLnByb2plY3RAZ21haWwuY29tIiwic2NvcGUiOiJyYXRlOnB1YmxpYyIsImlzcyI6Imh0dHBzOi8vZGF0YS5uZW9uc2NpZW5jZS5vcmcvIiwiZXhwIjoxODMyMjY5MjQ2LCJpYXQiOjE2NzQ1ODkyNDYsImVtYWlsIjoibWFjcm9zaGVkcy5wcm9qZWN0QGdtYWlsLmNvbSJ9.aBPsQZFZU8fgzTtcI78GjS5lJN9yt9JDZETQUXS6cVTt2IiTQIn7vUwalSm0yP5acL8wD8tcEyuCNf9qTaVlBA')

# prep data ----
source('R/nmh_prep_metab_inputs.R')
source('R/nmh_get_scaling_coefs.R')
source('R/nmh_get_neon_q_eval.R')
source('R/nmh_apply_neon_q_eval.R')
source('R/nmh_get_neon_q_sim.R')
# nmh_prep_metab_inputs(sensor_src = 'neon',
#                       q_type = 'raw',
#                       z_method = 'model')
# nmh_prep_metab_inputs(sensor_src = 'neon',
#                       q_type = 'raw',
#                       z_method = 'meas')
# nmh_prep_metab_inputs(sensor_src = 'neon',
#                       q_type = 'qaqc',
#                       z_method = 'model')
# nmh_prep_metab_inputs(sensor_src = 'neon',
#                       q_type = 'qaqc',
#                       z_method = 'meas')
# nmh_prep_metab_inputs(sensor_src = 'neon',
#                       q_type = 'simulated',
#                       z_method = 'model')
# nmh_prep_metab_inputs(sensor_src = 'neon',
#                       q_type = 'simulated',
#                       z_method = 'meas')

# nmh_prep_metab_inputs(sensor_src = 'streampulse',
#                       q_type = 'raw',
#                       z_method = 'model')
# nmh_prep_metab_inputs(sensor_src = 'streampulse',
#                       q_type = 'raw',
#                       z_method = 'meas')
# nmh_prep_metab_inputs(sensor_src = 'streampulse',
#                       q_type = 'qaqc',
#                       z_method = 'model')
# nmh_prep_metab_inputs(sensor_src = 'streampulse',
#                       q_type = 'qaqc',
#                       z_method = 'meas')
# nmh_prep_metab_inputs(sensor_src = 'streampulse',
#                       q_type = 'simulated',
#                       z_method = 'model')
# nmh_prep_metab_inputs(sensor_src = 'streampulse',
#                       q_type = 'simulated',
#                       z_method = 'meas')

# model MLE ----
source('R/nmh_model_neon_metab_mle.R')
nmh_model_neon_metab_mle()

# eval mle ----
source('R/nmh_eval_metab_mle.R')
nmh_eval_metab_mle()


# model bayes

source('R/nmh_model_metab_bayes.R')

library(dplyr)
library(foreach)
library(doParallel)

unregister_dopar <- function() {
    env <- foreach:::.foreachGlobals
    rm(list = ls(env), pos = env)
}

unregister_dopar()

# set up parallel processing environment info
n.cores <- parallel::detectCores()

if(.Platform$OS.type == 'windows') {
    cl <- parallel::makeCluster(n.cores, type = "PSOCK")
} else {
    cl <- parallel::makeCluster(n.cores, type = "FORK")
}

doParallel::registerDoParallel(cl)
foreach::getDoParWorkers()

nmh_model_metab_bayes()

stopCluster(cl)
unregister_dopar()
