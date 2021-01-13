#' @importFrom stats coef lm median pnorm sd time 
#' @importFrom utils data tail 
#' 

.onAttach <- function(library, pkgname)
{
  info <-utils::packageDescription(pkgname)
  package <- info$Package
  version <- info$Version
  date <- info$Date
  packageStartupMessage(
      paste(package, version, paste("(",date, ")", sep=""), "\n")
  )
}

if(getRversion() >= "2.15.1")  utils::globalVariables(
  c('.','CI_predict', 'CI_predict_lower', 'CI_predict_upper', 'CI_real', 
    'CI_real_lower', 'CI_real_upper', 'Estimate_predict', 'Estimate_real', 'K', 
    'N', 'N_predict', 'N_real', 'Na', 'Ns', 'Nsurv', 'Nt', 'alpha', 'beta0', 
    'beta1', 'num_non_zero', 'obl', 'p', 'region', 'site', 'site_num',
    'sp_scale', 'survey', 'survey_years', 'year', 'z', 'z_pred','fitted_pred', 
    'fitted_real', 'growth_pred', 'growth_real')
  )