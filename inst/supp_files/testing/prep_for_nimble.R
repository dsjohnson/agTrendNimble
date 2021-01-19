#' @title Create data list for NIMBLE analysis
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
  x <- filter(x, year>=timeframe[1], year<=timeframe[2], !is.na(counts)) %>% 
    arrange(site, year)
  if(!'obl'%in%colnames(x)) x$obl <- 0
  
  # Summarize site-level information
  site_data <- x %>% group_by(site, region) %>% 
    summarize(
      first = min(year),
      last = max(year),
      num_survey = n(),
      num_non_zero = sum(counts>0),
      num_zero = sum(counts==0),
      min_count = min(counts, na.rm=TRUE),
      max_count = max(counts, na.rm=TRUE),
      mean_count = mean(counts, na.rm=TRUE),
      sd_count = sd(counts, na.rm=TRUE),
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
      time = as.integer(year-min(year)+1)
    ) 
  
  upper <- select(site_data, site, max_count, mean_count) %>% 
    mutate(
      upper=max.factor*max_count,
      ) %>% select(-max_count)
  
  x <- left_join(x, upper, by='site') %>% 
    mutate(
      zi=ifelse(is.na(counts), mean_count, counts),
      lu1 = ifelse(counts==0, -Inf, counts-0.499999),
      lu2 =  ifelse(is.na(counts), upper, counts+0.499999),
      idx=1:nrow(x)
    ) %>% select(-mean_count)
  
  ##############################################################################
  ####                Create lists for NIMBLE run                           ####
  ##############################################################################
  if(debug==2) browser()
  # Constants
  
  xred <- x %>% filter(!is.na(counts))
  
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
    lu = select(xred, lu1, lu2) %>% as.matrix(),
    y = rep(1,constants$Nsurv),
    obl=xred$obl, # Data
    beta1 = ifelse(site_data$model=='const',0,NA), # parameter vector (Ns)
    sig_beta1 = ifelse(site_data$model=='const', 1, NA),
    K=K, #known matrix (Nt x Na)
    p = ifelse(site_data$model=='const', 0, NA), # parameter vector (Ns)
    alpha=sapply(
      site_data$model, 
      function(x){if(x=='tprs') return(rep(NA,ncol(K))) else return(rep(0,ncol(K)))}
      ), # parameter matrix (Nt x Ns)
    tau=ifelse(site_data$model=='tprs', NA, 0), # parameter vector (Ns)
    sp_scale=1
  )
  
  inits <- list(
    beta0 = site_data$mean_count,
    sig_beta0 = site_data$mean_count,
    beta1 = ifelse(site_data$model=='const', NA, 0),
    sig_beta1 = ifelse(site_data$model=='const', NA, 1),
    z = xred$counts*exp(xred$obl*0.03903366),
    alpha = sapply(
      site_data$model, 
      function(x){if(x=='tprs') return(rep(0,ncol(K))) else return(rep(NA,ncol(K)))}
    ),
    tau = ifelse(site_data$model=='tprs', 1, NA),
    phi = site_data$sd_count^2,
    p = ifelse(site_data$model=='const', NA, 0),
    gamma = 0.03903366
  )
  
  monitors <- sort(c(
    'z', 'mu', 'alpha', 'beta0', 'beta1', 'tau',
    'v', 'phi','p','gamma', 'z_pred', 'y_pred', 'sig_beta0', 
    'sig_beta1', 'prob_0'
    ))
  
  
  out <- list(site_data=site_data, fitting_data=x, constants=constants,
              data=data, monitors=monitors, inits=inits)
  if(debug==3) browser()
  return(out)
  
}