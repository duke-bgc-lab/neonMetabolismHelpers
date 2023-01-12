nmh_prep_bayes_model <- function() {
  tryCatch({
    library('rstan')
    library('StanHeaders')
    
  },
  error = function(e){
    install.packages("StanHeaders", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
    install.packages("rstan", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
    library('rstan')
    library('StanHeaders')
  }
  )
  
  # load streamMetabolizer 0.12.0
  # version 0.12.0 is required, so easiest path is to install from github
  # also loads dependencies not on CRAN (including lakemetabolizer)
  tryCatch({
    library('unitted')
    library('streamMetabolizer')
  },
  error = function(e){
    remotes::install_github("https://github.com/appling/unitted")
    remotes::install_github("https://github.com/USGS-R/streamMetabolizer")
    library('unitted')
    library('streamMetabolizer')
  }
  )
}