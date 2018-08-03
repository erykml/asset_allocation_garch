# Author: Eryk Lewinson

#-------------------------------------------------------------------------------
# Long-term asset allocation strategies based on GARCH models - 
# a simulation exercise
# ------------------------------------------------------------------------------
# A short exercise connected to investigating the forecasting performance of 
# different GARCH models in a long-term asset allocation problem. 
# 
# For this project I estimate various univariate AR-(GJR)GARCH models (with 
# different distributions of innovations and estimation types, i.e., fixed/
# /rolling/expanding window). I limit the scope of this exercise to only one
# risky asset (in this case the SP500 index) due to the curse of dimensionality.
# For the risk-free rate I take the Federal Fund Rate. I also adjust the returns 
# for inflation, using CPI data from the FRED database.
# 
# I use the estimated parameters of the GARCH models to simulate the evolution 
# of asset returns using the conditional mean/variance decomposition of asset
# returns as in Wilhelmsson (2013). Having the simulated returns' paths, 
# I calculate the optimal weights by directly maximizing the expected power 
# utility function, as in Barberis (2000), Brandt and van Binsbergen (2007).
#
# To evaluate the performance of different strategies (strategy is in fact the
# same, but based of different model specifications) I employ the modified
# Sharpe ratio. I choose the most basic model (Gaussian GARCH model with the 
# fixed window estimation) as the benchmark and measure the outperformance of
# the more complex models.
# 
# The following script produces the results of the abovementioned exercise.
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

# Loading packages ----

packages <-
    c("tseries", "dplyr", "rugarch", "foreach", "doParallel", 
      "PerformanceAnalytics", "fBasics", "reshape2", "ggplot2")

load_packages(packages)
rm(packages)

# Loading and manipulating data ----

## downloading sp500 stock prices

sp <- get.hist.quote(instrument = "^GSPC", 
                     start = "1987-12-22", 
                     end = "2017-03-01", 
                     quote = "Close", 
                     provider = "yahoo", 
                     compression = "w") %>% 
      as.data.frame()

sp_p <- data.frame(date = row.names(sp), 
                 price = sp$Close)  

## calculating returns

n <- dim(sp_p)[1]
rtns <- sp_p$price[2:n] / sp_p$price[1:(n - 1)] - 1
rtns <- data.frame(date = sp_p$date[2:n], 
                   rtns = rtns)
rtns$date <- as.Date(rtns$date)
rm(sp, sp_p, n)

## loading risk-free rate and inflation

rf_df <- read.csv('FEDFUNDS.csv')
n_rf <- dim(rf_df)[1]
rtns_rf <- rf_df$FEDFUNDS[2:n_rf] / rf_df$FEDFUNDS[1:(n_rf - 1)] - 1
rf <- data.frame(date = rf_df$DATE[2:n_rf],
                 rf = rtns_rf)
rf$date <- as.Date(rf$date)
rm(n_rf, rf_df, rtns_rf)

cpi_df <- read.csv('CPIAUCSL.csv')
n_cpi <- dim(cpi_df)[1]
rtns_cpi <- cpi_df$CPIAUCSL[2:n_cpi] / cpi_df$CPIAUCSL[1:(n_cpi - 1)] - 1
infl <- data.frame(date = cpi_df$DATE[2:n_cpi],
                   infl = rtns_cpi)
infl$date <- as.Date(infl$date)
rm(n_cpi, cpi_df, rtns_cpi)

## selecting sample period

sample_start <- as.Date('1988-01-01')
estimation_end <- as.Date('1999-12-31')
sample_end <- as.Date('2000-12-31')

## selecting only valid dates

rtns %>% dplyr::filter(date >= sample_start & date <= sample_end) -> rtns
rf %>% dplyr::filter(date >= sample_start & date <= sample_end) -> rf
infl %>% dplyr::filter(date >= sample_start & date <= sample_end) -> infl

## expanding rf and infl from monthly to weekly

rf_weekly <- rf %>% 
                mutate(date = substr(date, 1, 7), 
                       rf = rf * 12 / 52)

infl_weekly <- infl %>% 
               mutate(date = substr(date, 1, 7),
                      infl = infl * 12 / 52)

rf <- data.frame(date = substr(rtns$date, 1, 7)) %>% 
      left_join(., rf_weekly, by = 'date') %>% 
      mutate(date = rtns$date)

infl <- data.frame(date = substr(rtns$date, 1, 7)) %>% 
        left_join(., infl_weekly, by = 'date') %>% 
        mutate(date = rtns$date)

rm(rf_weekly, infl_weekly)

## creating expost returns

rf$rf <- (1 + rf$rf) / (1 + infl$infl) - 1
rtns$rtns <- (1 + rtns$rtns) / (1 + infl$infl) - 1

rm(infl)

## global parameters

length <- 52 # no. of forecasting periods
rolling_window <- 52 * 12 # length of rolling window
start <- sum(rtns$date <= estimation_end) # last period of estimation sample
mean_rf <- mean(rf$rf[1:start]) # mean risk-free rate

## summary statistics of sp500 returns
n <- length(rtns$rtns)

summary_stat <- data.frame(period = c('estimation_sample', 'evaluation_sample'),
                           mean = c((((1 + mean(rtns$rtns[1:start])) ^ 52) - 1) * 100,
                                    (((1 + mean(rtns$rtns[(start+1):n])) ^ 52) - 1) * 100),
                           std = c(sd(100 * rtns$rtns[1:start]) * sqrt(52),
                                   sd(100 * rtns$rtns[(start+1):n]) * sqrt(52)),
                           skewness = c(skewness(rtns$rtns[1:start]),
                                        skewness(rtns$rtns[(start+1):n])),
                           kurtosis = c(kurtosis(rtns$rtns[1:start], method="moment"),
                                        kurtosis(rtns$rtns[(start+1):n], method="moment")),
                           jarque_bera = c(jarqueberaTest(rtns$rtns[1:start])@test$p.value,
                                           jarqueberaTest(rtns$rtns[(start+1):n])@test$p.value),
                           acf = c(acf(rtns$rtns[1:start],plot=F,lag.max=1)$acf[2],
                                   acf(rtns$rtns[(start+1):n],plot=F,lag.max=1)$acf[2]),
                           acf_2 = c(acf((rtns$rtns[1:start])^2,plot=F,lag.max=1)$acf[2],
                                     acf((rtns$rtns[(start+1):n])^2,plot=F,lag.max=1)$acf[2])
                           )
rm(n)

write.csv(summary_stat, 'summary_statistics.csv', row.names = FALSE)

# Running  simulations ----

## Preparing the grid of all simulation inputs

model_combinations <- expand.grid(c('sGARCH', 'gjrGARCH'),
                                  c('norm', 'snorm', 'std', 
                                    'sstd', 'ged', 'sged'), 
                                  c('fixed', 'rolling', 'expanding'), 
                                  stringsAsFactors = FALSE) %>% 
                      as.data.frame() %>% 
                      'colnames<-'(c('model', 'distr', 'est_type')) %>% 
                      mutate(window = ifelse(est_type == 'rolling',
                                             rolling_window, 0)) 
## setting up parallel computation
    
no_cores <- detectCores() - 1
registerDoParallel(makeCluster(no_cores))

## running the simulations

sim_results <- foreach(model_indices = seq(dim(model_combinations)[1]),
                       .combine = data.frame,
                       .export = as.vector(lsf.str()), 
                       .packages = c('dplyr', 'tseries', 'rugarch'))  %dopar%  
    garch_simulations(rtns = 100 * rtns$rtns, 
                      mean_rf = mean_rf,
                      n_sims = 100000,
                      est_type = model_combinations$est_type[model_indices], 
                      garch_type = model_combinations$model[model_indices],
                      distr = model_combinations$distr[model_indices],
                      start = start, 
                      length = length,
                      window_length = model_combinations$window[model_indices],
                      lambda = 5)

## renaming simulation results
sim_results %>% 'colnames<-'(paste(model_combinations$est_type,
                                   model_combinations$model,
                                   model_combinations$distr, 
                                   sep = '_'
                                   )) -> sim_results

write.csv(sim_results, 'sim_results.csv')

## calculating summary statistics of realized portfolio returns
portfolio_rtns_stats <- data.frame(model = character(),
                                   mean = double(),
                                   sd = double(), 
                                   skew = double(),
                                   kurt = double(),
                                   SR = double(),
                                   stringsAsFactors = FALSE)

for (i in seq(dim(sim_results)[2]))
{
    portfolio_rtns_stats[i, ] <- c(colnames(sim_results)[i],
                                   (((1 + mean(sim_results[, i])) ^ 52) - 1) * 100,
                                   sd(sim_results[, i]) * sqrt(52) * 100,
                                   skewness(sim_results[, i]),
                                   kurtosis(sim_results[, i], method="moment"),
                                   (mean(sim_results[, i]) - mean_rf) / sd(sim_results[, i])
                                   )
}

write.csv(portfolio_rtns_stats, 'portfolio_rtns_stats.csv')

## evaluating the results

msr_results <- sapply(sim_results[, -1], function(x) mSR(x, sim_results[, 1], 
                                          rf$rf[(start + 1):length(rf$rf)])) %>% 
               data.frame(variant = names(.),
                          mSR = .,
                          row.names = NULL)
write.csv(msr_results, 'msr_results.csv')

# Additional plots and figures (load averything from start) ----

## downloading sp500 stock prices

sp <- get.hist.quote(instrument = "^GSPC", 
                     start = "1987-12-22", 
                     end = "2017-03-01", 
                     quote = "Close", 
                     provider = "yahoo", 
                     compression = "w") %>% 
      as.data.frame()

sp_p <- data.frame(date = row.names(sp), 
                   price = sp$Close)  

## calculating returns

n <- dim(sp_p)[1]
rtns <- sp_p$price[2:n]/sp_p$price[1:(n-1)] - 1
rtns <- data.frame(date = sp_p$date[2:n], 
                   rtns = rtns)
rtns$date <- as.Date(rtns$date)
rm(sp, sp_p, n)

## load simulation results
to_add <- dplyr::filter(rtns, date >= '2000-01-01' &
                              date <= '2001-01-01') %>% 
          'colnames<-'(c('date', 'actual_rtns'))

sim_results <- read.csv('sim_results.csv')[, -1] %>% 
               cbind(., to_add)

colnames(sim_results) <- colnames(sim_results) %>% 
                         gsub('_', ' ', .) %>% 
                         gsub('norm', 'N', .) %>% 
                         gsub('snorm', 'SN', .) %>%
                         gsub('std', 'T', .) %>%
                         gsub('sstd', 'ST', .) %>%
                         gsub('ged', 'GED', .) %>%
                         gsub('sged', 'SGED', .) %>% 
                         gsub('fixed ', '', .) %>%
                         gsub('rolling ', '', .) %>% 
                         gsub('expanding ', '', .) %>% 
                         gsub(' ', '-', .) %>% 
                         gsub('actual-rtns', 'Actual returns', .)

rm(to_add)

sim_results_fixed_long <- melt(sim_results[,c(1:12, 37:38)], id="date")  
sim_results_rolling_long <- melt(sim_results[,c(13:24, 37:38)], id="date")  
sim_results_expanding_long <- melt(sim_results[,c(25:36, 37:38)], id="date")  

Sys.setlocale("LC_TIME", "English")

ggplot(data=sim_results_fixed_long,
       aes(x=date, y=value, colour=variable)) +
       geom_line() +
       ggtitle("Portfolio returns over the evaluation sample (fixed window)") +
       xlab('Date') +
       ylab('Returns') +
       theme(plot.title = element_text(lineheight=.8, face="bold"))

ggplot(data=sim_results_rolling_long,
       aes(x=date, y=value, colour=variable)) +
       geom_line() +
       ggtitle("Portfolio returns over the evaluation sample (rolling window)") +
       xlab('Date') +
       ylab('Returns') +
       theme(plot.title = element_text(lineheight=.8, face="bold"))

ggplot(data=sim_results_expanding_long,
       aes(x=date, y=value, colour=variable)) +
       geom_line() +
       ggtitle("Portfolio returns over the evaluation sample (expanding window)") +
       xlab('Date') +
       ylab('Returns') +
       theme(plot.title = element_text(lineheight=.8, face="bold"))

 