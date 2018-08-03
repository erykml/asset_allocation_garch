# Author: Eryk Lewinson

#-------------------------------------------------------------------------------
# Long-term asset allocation strategies based on GARCH models - 
# a simulation exercise
# ------------------------------------------------------------------------------
# This script contains the functions used for obtaining the results in the 
# separate notebook.
#-------------------------------------------------------------------------------

# User defined functions ----

## loading packages
load_packages <- function(pkg) {
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg))
        install.packages(new.pkg, dependencies = TRUE)
    sapply(pkg, require, character.only = TRUE)
}

## estimating GARCH models
garch_ar_expanding_estimation <- function(data, start, length, garch_type,
                                          window_length, distr){
    for (i in 1:length){
    ### garch specification
        garch.spec <- ugarchspec(variance.model = list(model = garch_type,
                                                       garchOrder = c(1,1),
                                                       variance.targeting = FALSE),
                                 mean.model = list(armaOrder = c(1,0)), 
                                 distribution.model = distr)
        
        ### garch estimation on expanding window
        garch.fit <- ugarchfit(spec = garch.spec, 
                               data = data[1:(start + i - 1)], 
                               solver = "hybrid", 
                               solver.control = list(trace = 0))
        
        ### saving results into a data.frame
        coefs_input <- c(garch.fit@fit$coef, 
                         residuals(garch.fit)[length(residuals(garch.fit))],
                         (sigma(garch.fit)[length(sigma(garch.fit))])^2)
        
        if(i == 1)
        {
            garchcoefs <- matrix(0, length, length(coefs_input)) %>% 
                          as.data.frame()
        }
        
        garchcoefs[i, ] <- coefs_input
    }
    colnames(garchcoefs) <- c(names(garch.fit@fit$coef), 'e_T', 'sigma2_T')
    return(garchcoefs)
}

garch_ar_fixed_estimation <- function(data, start, length, window_length,
                                      garch_type, distr){
    ### garch specification
    garch.spec <- ugarchspec(variance.model = list(model = garch_type,
                                                   garchOrder = c(1,1),
                                                   variance.targeting = FALSE),
                             mean.model = list(armaOrder = c(1,0)), 
                             distribution.model = distr)
    
    ### garch estimation on fixed window
    garch.fit <- ugarchfit(spec = garch.spec, 
                           data = data[1:start], 
                           solver = "hybrid", 
                           solver.control = list(trace = 0))
    
    ### saving results into a data.frame
    coefs_input <- c(garch.fit@fit$coef, 
                     residuals(garch.fit)[length(residuals(garch.fit))],
                     (sigma(garch.fit)[length(sigma(garch.fit))])^2)
    
    garchcoefs <- matrix(coefs_input, length, length(coefs_input),
                         byrow = TRUE) %>% 
                  as.data.frame() %>% 
                  'colnames<-'(c(names(garch.fit@fit$coef), 'e_T', 'sigma2_T'))
    
    return(garchcoefs)
    
}

garch_ar_rolling_estimation <- function(data, start, length, window_length, 
                                        garch_type, distr){
    for (i in 1:length){
        
        ### garch specification
        garch.spec <- ugarchspec(variance.model = list(model = garch_type,
                                                       garchOrder = c(1,1),
                                                       variance.targeting = FALSE),
                                 mean.model = list(armaOrder = c(1,0)), 
                                 distribution.model = distr)
        ### garch estimation on rolling window
        garch.fit <- ugarchfit(spec = garch.spec, 
                               data = data[(start - window_length + i - 1):
                                               (start + i - 1)], 
                               solver = "hybrid", 
                               solver.control = list(trace = 0))
        
        ###csaving results into a data.frame
        coefs_input <- c(garch.fit@fit$coef, 
                         residuals(garch.fit)[length(residuals(garch.fit))],
                         (sigma(garch.fit)[length(sigma(garch.fit))])^2)
        
        if(i == 1)
        {
            garchcoefs <- matrix(0, length, length(coefs_input)) %>% 
                as.data.frame()
        }
        
        garchcoefs[i, ] <- coefs_input
    }
    colnames(garchcoefs) <- c(names(garch.fit@fit$coef), 'e_T', 'sigma2_T')
    return(garchcoefs)
}

## forecasting using garch estimation results
garch_forecast <- function(est_results, rtns_t, h_t, resid_t, garch_type){
    if(grepl('gjr', garch_type) == 1){
        ### gjr garch forecasting 
        mu <-  (est_results$mu + est_results$ar1 * rtns_t) / 100
        sigma <- (est_results$omega + est_results$alpha1 * resid_t^2 +
                      est_results$beta1 * h_t + 
                      est_results$gamma1 * resid_t^2 * (resid_t<0)) / 10000
    } else {
        ### garch forecasting    
        mu <-  (est_results$mu + est_results$ar1 * rtns_t) / 100
        sigma <- (est_results$omega + est_results$alpha1 * resid_t^2 +
            est_results$beta1 * h_t) / 10000 
    }
    return(data.frame(cond_mean = mu, cond_var = sigma))
}

## solving the asset allocation problem
asset_allocation <- function(garch_sim, cond_var, n_sims, length, curr_sim, 
                             lambda, mean_rf){
    uf <- rep(1, n_sims) %>% cbind()
    opt_u <- rep(-Inf, n_sims) %>% cbind()
    port_w <- rep(0, n_sims) %>% cbind()
    
    for(i in seq((length - curr_sim), 1, -1))
    {
        for(alpha in seq(0, 1, 0.01))
        {
            U <- (1 / (1 - lambda)) * (((1 + alpha * garch_sim[, (i+1)] + 
                                             mean_rf * (1 - alpha)) ^ (1 - lambda)) * uf)
            
            x <- cbind(rep(1, n_sims), cond_var[, (i+1)], garch_sim[, i])
            coefs <- solve(t(x) %*% x) %*% t(x) %*% U
            cu_yhat <- x %*% coefs
            
            port_w[cu_yhat > opt_u] <- alpha
            opt_u[cu_yhat > opt_u] <- cu_yhat[cu_yhat > opt_u]
        }
        
        uf <- ((1 + port_w * garch_sim[, (i+1)] + (1 - port_w) *  mean_rf) ^ 
                   (1 - lambda)) * uf
        
        opt_u <- rep(-Inf, n_sims) %>% cbind()
        port_w = rep(0, n_sims) %>% cbind()
    }
    
    for(alpha in seq(0, 1, 0.01))
    {
        U <- (1 / (1 - lambda)) * (((1 + alpha * garch_sim[, 1] + 
                                         mean_rf * (1 - alpha)) ^ (1 - lambda)) * uf)
        
        x <- rep(1, n_sims) %>% cbind()
        coefs <- solve(t(x) %*% x) %*% t(x) %*% U
        cu_yhat <- x %*% coefs
        
        port_w[cu_yhat > opt_u] <- alpha
        opt_u[cu_yhat > opt_u] <- cu_yhat[cu_yhat > opt_u]
    }
    return(port_w[1])    
} 

## solving the asset allocation problem for the last time period
asset_allocation_last <- function(garch_sim, cond_var, n_sims, length, curr_sim, 
                                  lambda, mean_rf){
    opt_u <- rep(-Inf, n_sims) %>% cbind()
    port_w <- rep(0, n_sims) %>% cbind()
    
    for(alpha in seq(0, 1, 0.01))
    {
        U <- (1 / (1 - lambda)) * (((1 + alpha * garch_sim[, 1] + 
                                         mean_rf * (1 - alpha)) ^ (1 - lambda)))
        
        x <- rep(1, n_sims) %>% cbind()
        coefs <- solve(t(x) %*% x) %*% t(x) %*% U
        cu_yhat <- x %*% coefs
        
        port_w[cu_yhat > opt_u] <- alpha
        opt_u[cu_yhat > opt_u] <- cu_yhat[cu_yhat > opt_u]
    }
    return(port_w[1])  
} 

## calculating portfolio returns
port_rtns <- function(garch_wghts, rtns, mean_rf, start, length){
    n <- length(rtns)
    port_rtns <- garch_wghts * exp(rtns[(start + 1) : n] / 100) + 
        (1 - garch_wghts) * rep(exp(mean_rf), length) - 1
    return(port_rtns)
}

## garch simulation function
garch_simulations <- function(rtns, mean_rf, n_sims, est_type, garch_type, 
                              distr, start, length, window_length, lambda){
    result_weights <- rep(NA, length)
    
    est_results <- eval(parse(text = paste0('garch_ar_', est_type, '_estimation',
                                     '(data = rtns, start = start,' , 
                                     'length = length,', 
                                     'window_length = window_length,',
                                     'garch_type = garch_type, distr = distr)')))
    for (l in seq(length))
    {
        ### result matrices
        cond_var <- matrix(0, n_sims, (length - l + 1)) 
        garch_sim <- matrix(0, n_sims, (length - l + 1)) 
        resids <- matrix(0, n_sims, (length - l + 1)) 
        
        ### generating the random component based on the selected distribution
        if(distr == 'norm')
        { rand_var <- matrix(rnorm(n_sims * (length - l + 1)),
                             n_sims, 
                             (length - l + 1))
        }
        
        if(distr == 'snorm')
        { rand_var <- matrix(rdist(distribution = "snorm", 
                                   n = n_sims * (length - l + 1), mu = 0, 
                                   sigma = 1, skew = est_results$skew[l]),
                             n_sims, 
                             (length - l + 1))
        }
        
        if(distr == 'std')
        { rand_var <- matrix(rdist(distribution = "std",
                                   n = n_sims * (length - l + 1), 
                                   shape = est_results$shape[l]),
                             n_sims, 
                             (length - l + 1))
        }
        
        if(distr == 'sstd')
        { rand_var <- matrix(rdist(distribution = "sstd",
                                   n = n_sims * (length - l + 1), 
                                   shape = est_results$shape[l],
                                   skew = est_results$skew[l]),
                             n_sims, 
                             (length - l + 1))
        }
        
        if(distr == 'ged')
        { rand_var <- matrix(rdist(distribution = "ged",
                                   n = n_sims * (length - l + 1), 
                                   shape = est_results$shape[l]),
                             n_sims, 
                             (length - l + 1))
        }
        
        if(distr == 'sged')
        { rand_var <- matrix(rdist(distribution = "sged",
                                   n = n_sims * (length - l + 1), 
                                   shape = est_results$shape[l],
                                   skew = est_results$skew[l]),
                             n_sims, 
                             (length - l + 1))
        }
        
        
        ### inputting the values for the first step
        preds <- garch_forecast(est_results[l, ], 
                                rtns_t = rtns[start + l - 1], 
                                h_t = est_results$sigma2_T[l], 
                                resid_t = est_results$e_T[l], 
                                garch_type)
            
        cond_var[, 1] <- preds$cond_var 
        garch_sim[, 1] <- preds$cond_mean + sqrt(preds$cond_var) * rand_var[, 1]
        resids[, 1] <- sqrt(preds$cond_var) * rand_var[, 1]
        
        if(l != length)
        {
            ### generating simulated paths
            for(tt in seq(2, (length - l + 1)))
            {
                ### conditional mean forecast
                cond_mean <- est_results$mu[l] / 100 + 
                             est_results$ar1[l] * garch_sim[, (tt-1)]
                ### conditional variance forecast
                cond_var[, tt] <- est_results$omega[l] / 10000 + 
                                  est_results$alpha1[l] * (resids[, (tt - 1)]^2) +
                                  est_results$beta1[l] * cond_var[, (tt - 1)] +
                                  ifelse(grepl('gjr', garch_type) == 1, 1, 0) * 
                                  ifelse(is.null(est_results$gamma1[l]), 0,
                                         est_results$gamma1[l]) * 
                                  (resids[, (tt - 1)]^2) *
                                  (resids[, (tt - 1)] < 0)
                resids[, tt] <- sqrt(cond_var[, tt]) * rand_var[, tt]
                ### simulations
                garch_sim[, tt] <- cond_mean + resids[, tt]
            }
            
            result_weights[l] <- asset_allocation(garch_sim, cond_var, n_sims,
                                                  length, curr_sim = l, lambda, 
                                                  mean_rf = mean_rf)
        } else {
            result_weights[l] <- asset_allocation_last(garch_sim, cond_var, 
                                                       n_sims, length, 
                                                       curr_sim = l, lambda, 
                                                       mean_rf = mean_rf)
        }
    }
    
    port_rtns <- port_rtns(result_weights, rtns, mean_rf, start, length)
    
    return(port_rtns)
    
}

## modified Sharpe ratio

mSR <- function(rtns_to_compare, rtns_benchmark, rf){
    msr <- (sd(rtns_benchmark) / sd(rtns_to_compare)) * 
           (mean(rtns_to_compare) - mean(rf)) - 
           (mean(rtns_benchmark) - mean(rf))
    return(msr)
}