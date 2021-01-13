#' @title Fit agTrend SSL model using MCMC with NIMBLE
#' 
#' @param x A named list generated from the \code{prep_for_nimble} that prepares the raw data
#' @param niter Number of iterations for the MCMC, this includes the burnin and is prior to thinning.
#' @param nburnin Number of burnin iterations. These will be discarded from \code{niter} before returning results.
#' @param thin Thinning of the MCMC sample for storage. Total number of retuned samples is \code{(niter-nburnin)/thin}.
#' @param debug Jump into the fucntion for debugging.
#' @param ... Additional arguments.
#' @import nimble
#' @export
#' 
#' @author Devin S. Johnson
fit_ssl_nimble <- function(x, niter = 110000, nburnin=10000, thin=10, debug=FALSE, ...){
  if(debug==1) browser()
  ### Model code
  # agt_ssl_code <- nimbleCode({
  #   
  #   for(j in 1:Nsurv){
  #     y[j] ~ dinterval(z[j], lu_adj[j,1:2])
  #     z[j] ~ dnorm(mu[idx[j]], var=v[idx[j]])
  #     lu_adj[j,1:2] <- exp(obl[j]*gamma)*lu[j,1:2]
  #   }
  #   
  #   for(i in 1:N){
  #     mu[i] <- beta0[site[i]] + beta1[site[i]]*(time[i]-1) + eta[site[i],time[i]]
  #     v_base[i] <- log(1 + exp(-sp_scale*abs(mu[i])))/sp_scale + max(0, mu[i])
  #     v[i] <- phi[site[i]]*v_base[i]^p[site[i]]
  #     z_pred[i] ~ dnorm(mu[i], var=v[i])
  #     y_pred[i] <- max(0, round(z_pred[i]))
  #     prob_0[i] <- pnorm(0.5, mu[i], sd=sqrt(v[i]))
  #   }
  #   for(s in 1:Ns){
  #     eta[s,1:Nt] <- K[1:Nt,1:Na]%*%alpha[1:Na,s]
  #     for(k in 1:Na){alpha[k,s] ~ dnorm(0, sd=tau[s])}
  #     tau[s] ~ dexp(1.0E-6)
  #     beta0[s] ~ dnorm(0, sd=sig_beta0[s])
  #     sig_beta0[s] ~ dexp(-log(0.05)/10000)
  #     beta1[s] ~ dnorm(0, sd=sig_beta1[s])
  #     sig_beta1[s] ~ dexp(-log(0.05)/10000)
  #     phi[s] ~ dexp(-log(0.05)/10000)
  #     p[s] ~ dexp(-log(0.05)/2)
  #   }
  #   gamma ~ dnorm(0.03903366,sd=0.01068773)
  #   
  # })
  agt_ssl_code <- nimbleCode({
    
    for(j in 1:Nsurv){
      y[j] ~ dinterval(z[j], c[j])
      z[j] ~ dnorm(exp(obl[j]*gamma)*mu[idx[j]], var=exp(2*obl[j]*gamma)*v[idx[j]])
      y_real[j] <- max(0, exp(-obl[j]*gamma)*z[j])
    }
    for(i in 1:N){
      mu[i] <- beta0[site[i]] + beta1[site[i]]*(time[i]-1) + eta[site[i],time[i]]
      v_base[i] <- log(1 + exp(-sp_scale*abs(mu[i])))/sp_scale + max(0, mu[i])
      v[i] <- phi[site[i]]*phi[site[i]]*v_base[i]^p[site[i]]
      z_pred[i] ~ dnorm(mu[i], var=v[i])
      y_pred[i] <- max(0, z_pred[i])
      prob_0[i] <- pnorm(0.5, mu[i], sd=sqrt(v[i]))
    }
    for(s in 1:Ns){
      eta[s,1:Nt] <- K[1:Nt,1:Na]%*%alpha[1:Na,s]
      for(k in 1:Na){alpha[k,s] ~ dnorm(0, sd=tau[s])}
      tau[s] ~ dexp(1.0E-8)
      beta0[s] ~ dnorm(0, var=1.0E8)
      beta1[s] ~ dnorm(0, sd=sig_beta1[s])
      sig_beta1[s] ~ dexp(1.0E-8)
      phi[s] ~ dexp(lambda_phi[s])
      p[s] ~ dexp(-log(0.05)/2)
    }
    gamma ~ dnorm(-0.03903366,sd=0.01068773)
    
  })
  
  message("building model...")
  suppressMessages(
    agt_ssl_model <- nimbleModel(
      code = agt_ssl_code,
      constants = x$constants,
      data = x$data,
      inits = x$inits
    )
  )
  message("building MCMC sampler...")
  agt_ssl_mcmc <- buildMCMC(agt_ssl_model, monitors=x$monitors, print=FALSE) 
  message("compling code...this will take a couple of minutes...")
  suppressMessages(agt_ssl_model_c <- compileNimble(agt_ssl_model))
  suppressMessages(agt_ssl_mcmc_c <- compileNimble(agt_ssl_mcmc))
  message("running MCMC...this will be a while...go get some coffee...")
  if(debug==2) browser()
  suppressMessages(smp <- runMCMC(agt_ssl_mcmc_c, 
                                  niter = niter, nchains = 1, nburnin=nburnin, thin=thin,
                                  samplesAsCodaMCMC = TRUE))
  message("sampler done...just doing some post-processing...")
  nms_idx <- as.numeric(factor(map_chr(strsplit(colnames(smp), "\\["), ~{.x[[1]]})))
  smp <- lapply(unique(nms_idx), function(x, smp, nms_idx){smp[,nms_idx==x]}, 
                smp=smp, nms_idx=nms_idx)
  names(smp) <- sort(x$monitors)
  smp <- list(
    abund = smp[names(smp)%in%c('y_pred','y_real','prob_0','mu')], 
    param = smp[!names(smp)%in%c('y_pred','y_real','prob_0','mu')]
  )
  y_real <- smp$abund$y_pred
  idx <- x$constants$idx
  for(i  in 1:length(idx)){
    y_real[,idx[i]] <- smp$abund$y_real[,i]
  }
  colnames(y_real)[idx] <- colnames(smp$abund$y_real)
  smp$abund$y_real <- y_real
  # for(i in 1:nrow(x$fitting_data)){
  #   if(is.na(x$fitting_data$counts[i])){
  #     next
  #   } else{
  #     smp$abund$y_real[,i] <- x$fitting_data$counts[i]
  #   }
  # }
  smp$summary <- x$fitting_data %>% select(-site_num,-time,-idx)
  for(i in 1:length(smp$abund)){
    ci <- round(coda::HPDinterval(smp$abund[[i]]),1)
    est <- round(apply(smp$abund[[i]], 2, median),1)
    tmp <- tibble(est, ci[,1], ci[,2]) %>% 
      `colnames<-`(paste0(names(smp$abund)[i], c('_est','_ciL','_ciU')))
    smp$summary <- cbind(smp$summary,tmp)
    idx <- rep(x$site_data$site, each=ncol(smp$abund[[i]])/nrow(x$site_data))
    smp$abund[[i]] <- split_col(smp$abund[[i]], idx)
  }
  smp$abund <- as_tibble(smp$abund) %>% bind_cols(x$site_data,.)
  attr(smp$abund, "years") <- sort(unique(x$fitting_data$year))
  
  if(debug==2) browser()
  
  return(smp)
  
}