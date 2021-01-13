#' @title Aggregate abundances based on a factor variable
#' @param x A data set containing columns '\code{agg.var}', '\code{y_pred}' 
#' and '\code{y_real}'. The last two are obtained from a calls to 
#' \code{agTrendNimble::fit_ssl_nimble}.
#' @param agg.var Variable used for aggregation
#' @import dplyr purrr
#' @export
#' 
ag_abund <- function(x, agg.var){
  results <- x %>% group_by(.data[[agg.var]]) %>% nest() %>%
    mutate(
      survey_years = map(data, ~{sort(reduce(.x$survey_years, union))}),
      N_predict = map(data, ~{reduce(.x$y_pred, `+`)}),
      N_real = map(data, ~{reduce(.x$y_real, `+`)})
    ) %>% ungroup() %>% select(-data)
  attr(results, "years") <- attr(x, "years")
  results
}


#' @title Summarize aggregation results
#' @param x An aggregated abundance object from \code{agg_abund} 
#' @param ci.prob Confidence (credible) interval value
#' @param smooth Logical. Should smoothing of piecewise credible intervals be used. This is primarily used
#' for plotting. No statistical inference is made from smoothed intervals. 
#' @import tidyr
#' @importFrom coda HPDinterval mcmc
#' @export
#' 
summary_ag <- function(x, ci.prob=0.95, smooth=TRUE){
  df <- data.frame(year=attr(x,"years"))
  x$Estimate_predict <- map(x$N_predict, ~{
    out <- data.frame(Estimate=apply(.x, 2, median))
    bind_cols(df, out)
  })
  x$CI_predict <- map(x$N_predict, ~{
    cidf <- coda::HPDinterval(coda::mcmc(.x), prob=ci.prob) %>% as.data.frame() %>% 
      `colnames<-`(c("CI_predict_lower", "CI_predict_upper"))
    if(smooth){
      time <- 1:nrow(cidf)
      for(i in 1:2){
        fff <- mgcv::gam(cidf[,i]~s(time, k=nrow(cidf)-1), method='REML')
        cidf[,i] <- predict(fff) 
      }
    }
    cidf
  })
  x <- select(x, -N_predict)
  x$Estimate_real <- map(x$N_real, ~{data.frame(Estimate_real=apply(.x, 2, median))})
  x$CI_real <- map(x$N_real, ~{
    coda::HPDinterval(coda::mcmc(.x), prob=ci.prob) %>% as.data.frame() %>% 
      `colnames<-`(c("CI_real_lower", "CI_real_upper"))
  })
  x <- select(x, -N_real)
  surv <- select(x, 1:2) %>% unnest(cols=survey_years) %>% 
    rename(year=survey_years) %>% mutate(survey=1)
  x <- x %>% select(-survey_years) %>% unnest(cols = c(Estimate_predict, CI_predict, Estimate_real, CI_real))
  x <- left_join(x, surv, by=colnames(surv)[1:2]) %>% 
    mutate(
      survey = ifelse(is.na(survey), 0, survey),
      CI_real_upper = ifelse(!survey, CI_predict_upper, CI_real_upper),
      CI_real_lower = ifelse(!survey, CI_predict_lower, CI_real_lower)
    )
  return(x)
}


#' @title Summarize aggregation results
#' @param x An aggregated abundance object from \code{ag_abund} 
#' @param timeframe A 2-vector giving the time frames for the calculated growth trends
#' @param ci.prob Confidence (credible) interval value
#' @importFrom coda HPDinterval mcmc
#' @export
#' 
ag_trend <- function(x, timeframe, ci.prob=0.95){
  if(is.null(dim(timeframe))) timeframe <- t(as.matrix(timeframe))
  nms <- pull(x,1)
  cn <- colnames(x)[1]
  yrs <- attr(x, 'years')
  min_yr <- min(yrs)
  max_yr <- max(yrs)
  zeros_pred <- map_lgl(x$N_predict, ~{any(.x==0)})
  zeros_real <- map_lgl(x$N_real, ~{any(.x==0)})  
  warning(paste0("There are abundnces of 0 animals for: ", nms[zeros_pred | zeros_real], ". Adding '0.5' to allow calculation."))
  ln_Np <- map2(x$N_predict, zeros_pred, ~{log(as.matrix(.x + 0.5*.y))})
  ln_Nr <- map2(x$N_real, zeros_real, ~{log(as.matrix(.x + 0.5*.y))})
  out <- x[,1]
  idx <- yrs>=timeframe[1] & yrs<=timeframe[2]
  H <- cbind(1,1:sum(idx))
  P <-  H %>% {solve(crossprod(.), t(.))}
  out$trend_pred <- map(ln_Np, ~{.x[,idx]}) %>% map(~{t(P%*%t(.x))}) %>% 
    map(coda::mcmc)
  out$fitted_pred <- map(out$trend_pred, ~{exp(t(H%*%t(.x)))}) %>% 
    map(coda::mcmc)
  out$growth_pred <-  map(out$trend_pred, ~{100*(exp(.x[,2])-1)}) %>% 
    map(coda::mcmc)
  out$trend_real <- map(ln_Nr, ~{.x[,idx]}) %>% map(~{t(P%*%t(.x))}) %>% 
    map(coda::mcmc)
  out$fitted_real <- map(out$trend_real, ~{exp(t(H%*%t(.x)))}) %>% 
    map(coda::mcmc)
  out$growth_real <-  map(out$trend_real, ~{100*(exp(.x[,2])-1)}) %>% 
    map(coda::mcmc)
  summ_fitted_pred <- select(out, .data[[cn]], fitted_pred) %>% 
    mutate(
      year = rep(list(yrs[idx]),nrow(x)),
      type = 'predicted',
      Est = map(fitted_pred, ~{apply(.x,2,median)}),
      CI = map(fitted_pred, ~{data.frame(coda::HPDinterval(.x, prob=ci.prob))})
    ) %>% select(-fitted_pred) %>% unnest(cols=c('year', 'Est', 'CI'))
  summ_fitted_real <- select(out, .data[[cn]], fitted_real) %>% 
    mutate(
      year = rep(list(yrs[idx]),nrow(x)),
      type = 'realized',
      Est = map(fitted_real, ~{apply(.x,2,median)}),
      CI = map(fitted_real, ~{data.frame(coda::HPDinterval(.x, prob=ci.prob))})
    ) %>% select(-fitted_real) %>% unnest(cols=c('year', 'Est', 'CI'))
  summary_fitted <- bind_rows(summ_fitted_pred, summ_fitted_real)
  summ_growth_pred <- select(out, .data[[cn]], growth_pred) %>% 
    mutate(
      type='predicted',
      Est = map(growth_pred, median),
      CI = map(growth_pred, ~{data.frame(coda::HPDinterval(.x, prob=ci.prob))})
    ) %>% select(-growth_pred) %>% unnest(cols=c('Est', 'CI'))
  summ_growth_real <- select(out, .data[[cn]], growth_real) %>% 
    mutate(
      type='realized',
      Est = map(growth_real, median),
      CI = map(growth_real, ~{data.frame(coda::HPDinterval(.x, prob=ci.prob))})
    ) %>% select(-growth_real) %>% unnest(cols=c('Est', 'CI'))
  summary_growth <- bind_rows(summ_growth_pred, summ_growth_real)
  
  return(list(growth=summary_growth, fitted=summary_fitted, sample=out))
}
