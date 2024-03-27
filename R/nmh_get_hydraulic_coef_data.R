#' Access file with hydraulic scaling coefficients for 24 NEON wadable sites
#' 
#' @author Nick Marzolf, \email{nick.marzolf@@jonesctr.org}
#' @author Wes Slaughter, \email{wslaughter@@berkeley.edu}
#' 
#' @details
#' This function loads Data Citation item 3, which contains scaling coefficients, error estimates, and their power law scaling fits for each site (Leopold and Maddock 1953).
#'
#' @export

nmh_get_hydraulic_coef_data <- function() {
  return(nmh_hydraulic_coef_data)
}
