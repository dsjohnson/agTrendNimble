#' @title Create data list for NIMBLE analysis
#' @description Create a named list with data structures necessary to run the MCMC in \code{NIMBLE}
#' @param x Data containing the Steller sea lion counts and sites
#' @param timeframe A 2-vector containing the start and years for the modeling effort
#' @param model.cuts A 2-vector specifying the numbers of positive counts necessary for each model type. 
#' E.g., \code{model.cuts=c(5,10)} implies that a constant model (over years) 
#' is used for sites with <=5 positive counts, a linear model is used for the number of positive counts in 
#' (5, 10], and a low-rank Gaussian process is used for the number of positive counts >10. 
#' @param max.factor Maximum value for imputed counts calculated as \code{max.factor * max(counts)}. 
#' Currently this is not implemented, so, it is ignored
#' @param debug Logical. Should debug mode be entered. Probably just for Devin to use. 
#' @author Devin S. Johnson
#' @import dplyr tidyr mgcv
#' @export
#' 
prep_for_nimble <- function(x, timeframe, model.cuts=c(5,10), max.factor=3, debug=FALSE){
  
  ##############################################################################
  ####                Process/Summarize data                                ####
  ##############################################################################
  if(debug==1) browser()
  colnames(x) <- tolower(colnames(x))
  if(!all(c('year', 'site','counts')%in%colnames(x))) stop("There are no 'year' ('YEAR'), 'site' ('SITE'), or 'counts' ('COUNTS') columns. Please rename the appropriate columns!")
  if(missing(timeframe)) timeframe <- range(x$year)
  
  # Remove counts outside timeframe
  x <- filter(x, year>=timeframe[1], year<=timeframe[2], !is.na(.data$counts)) %>% 
    arrange(site, year)
  if(!'obl'%in%colnames(x)) x$obl <- 0
  
  # Summarize site-level information
  site_data <- x %>% group_by(site, region) %>% 
    summarize(
      first = min(year),
      last = max(year),
      num_survey = n(),
      num_non_zero = sum(.data$counts>0),
      num_zero = sum(.data$counts==0),
      min_count = min(.data$counts, na.rm=TRUE),
      max_count = max(.data$counts, na.rm=TRUE),
      mean_count = mean(.data$counts, na.rm=TRUE),
      sd_count = sd(.data$counts, na.rm=TRUE),
      survey_years = list(unique(year)),
      .groups='drop'
    ) %>% filter(num_non_zero>=2)
  
  model.cuts <- c(0,model.cuts,timeframe[2]-timeframe[1]+1)
  site_data <- mutate(site_data,
                      model = cut(num_non_zero, model.cuts, labels=c('const','lin','tprs'))
                      )
  
  # Add missing values for unsurveyed years and sites
  x <- filter(x, site%in%site_data$site) %>% 
    complete(nesting(site, region), year=timeframe[1]:timeframe[2], fill=list(obl=0)) %>% 
    arrange(site, year) %>% 
    mutate(
      site_num = as.integer(factor(site)),
      time = as.integer(year-min(year)+1),
      idx = 1:n()
    ) 

  ##############################################################################
  ####                Create lists for NIMBLE run                           ####
  ##############################################################################
  if(debug==2) browser()
  # Constants
  
  xred <- x %>% filter(!is.na(.data$counts))
  
  constants <- list(
    N=nrow(x), 
    Nsurv = nrow(xred),
    Ns=nrow(site_data), 
    Na=8, 
    Nt=length(timeframe[1]:timeframe[2]),
    site=x$site_num, # Data
    time=x$time, # Data
    idx = xred$idx
    )
  
  # Pen. spline smoothing matrix
  sm <- mgcv::smoothCon(mgcv::s(time), data.frame(time=timeframe[1]:timeframe[2]))[[1]]
  null.sp <- tail(1:ncol(sm$X), sm$null.space.dim)
  Xf <- as.matrix(sm$X[,null.sp])
  S <- sm$S[[1]][-null.sp,-null.sp]
  L <- mgcv::mroot(solve(S))
  K <- sm$X[,-null.sp] %*% L
  K <- (diag(nrow(Xf))-Xf%*%solve(crossprod(Xf))%*%t(Xf))%*%K
  # rw.order <- 2
  # K <- iar_basis(length(timeframe[1]:timeframe[2]), rw.order)
  
  data <- list(
    c = rep(0.5, constants$Nsurv),
    y = ifelse(xred$counts==0, 0, 1),
    z = ifelse(xred$counts==0, NA, xred$counts),
    obl=xred$obl, # Data
    beta1 = ifelse(site_data$model=='const',0,NA), # parameter vector (Ns)
    sig_beta1 = ifelse(site_data$model=='const', 1, NA),
    K=K, #known matrix (Nt x Na)
    p = ifelse(site_data$model=='const', 0, 0), # parameter vector (Ns)
    alpha=sapply(
      site_data$model, 
      function(x){if(x=='tprs') return(rep(NA,ncol(K))) else return(rep(0,ncol(K)))}
      ), # parameter matrix (Nt x Ns)
    tau=ifelse(site_data$model=='tprs', NA, 0), # parameter vector (Ns)
    lambda_phi = -log(0.05)/site_data$max_count,
    sp_scale=1
  )
  
  beta0 <- rep(NA, nrow(site_data))
  beta1 <- rep(NA, nrow(site_data))
  phi <- rep(NA, nrow(site_data))
  for(i in 1:nrow(site_data)){
    df <- xred[xred$site==site_data$site[i],]
    if(site_data$model[i]=='const'){
      fff <- lm(counts ~ 1, data=df)
    } else{
      fff <- lm(counts ~ I(time-1), data=df)
      beta1[i] <- coef(fff)[2]
    }
    beta0[i] <- coef(fff)[1]
    phi[i] <- suppressWarnings(summary(fff)$sigma)
  }
  
  inits <- list(
    beta0 = beta0,
    beta1 = beta1,
    sig_beta1 = 2*abs(beta1),
    z = ifelse(xred$counts==0, 0, NA),
    alpha = sapply(
      site_data$model, 
      function(x){if(x=='tprs') return(rep(0,ncol(K))) else return(rep(NA,ncol(K)))}
    ),
    tau = ifelse(site_data$model=='tprs', 1, NA),
    phi = ifelse(phi>0, phi, 0.001),
    p = ifelse(site_data$model=='const', NA, NA),
    gamma = -0.03903366
  )
  
  monitors <- sort(c(
    'mu', 'alpha', 'beta0', 'beta1', 'tau',
    'phi','p','gamma', 'y_pred',
    'sig_beta1', 'prob_0', 'y_real'
    ))
  
  
  out <- list(site_data=site_data, fitting_data=x, constants=constants,
              data=data, monitors=monitors, inits=inits)
  if(debug==3) browser()
  return(out)
  
}