#' This function implenets the posterior sampling strategy described in the article
#' "Bayesian Time-to-event Density Regression with Application to High Dimensional Data".

#' @references Maity, A. K., Pati, D., and Mallick, B. K. (2024) 
#'             "Bayesian Time-to-event Density Regression with Application to High Dimensional Data"
#'             
#'             
#'@param ct survival response, a \eqn{n*2} matrix with first column as response and second column as right censored indicator,

#'  1 is event time and 0 is right censored.

#'@param X Matrix of covariates, dimension \eqn{n*p}.

#'@param nburnin Number of burn-in MCMC samples. Default is 1000.

#'@param nmc Number of posterior draws to be saved. Default is 5000.

#'@param thin Thinning parameter of the chain. Default is 1 (no thinning).

#'@param method.mcmc posterior MCMC sampling methods. There are two options -- 
#' Metropolis-Hastings and Slice sampling. 
#' 
#' @param prior prior specifications on the regression parameters. There are two 
#' options -- Laplace prior (or Bayesian lasso prior) and Horseshoe prior

#'@param Xtest test design matrix.

#'@param cttest test survival response.



#'@return \item{BetaHat}{Posterior mean of \eqn{\beta} for survival model, a \eqn{p} by 1 vector.}

#'\item{BetaMedian}{Posterior median of \eqn{beta} for survival model, a \eqn{p} by 1 vector.}

#'\item{Sigma2Hat}{Posterior mean of variance covariance matrix.}

#'\item{gammaHat}{Posterior mean of \eqn{\gamma}}

#'\item{BetaSamples}{Posterior samples of \eqn{\beta} for survival model.}

#'\item{Sigma2Samples}{Posterior samples of \eqn{\sigma^2}}
#'
#'\item{gamma.initial}{Starting values of \eqn{\gamma}}
#'
#'\item{gamma.samples}{Posterior samples of \eqn{\gamma}}

#'\item{SurvivalHat}{Predictive survival probability.}

#'\item{LogTimeHat}{Predictive log time.}

#'\item{LambdaSamples}{Posterior samples of \eqn{\lambda}.}

#'\item{TauSamples}{Posterior samples of \eqn{\tau}.}

#'\item{DIC}{DIC of the model}

#'\item{WAIC}{WAIC for binary model.}

#'\item{LeftCI}{The left bounds of the credible intervals for BetaHat.}

#'\item{RightCI}{The right bounds of the credible intervals for BetaHat.}


#' @examples
#' 
#' 


gph <- function (ct, X, nburnin = 1000, nmc = 5000, thin = 1, method.mcmc = c("MH", "slice"),
                    prior = c("laplace", "horseshoe"), Xtest = NULL, cttest = NULL)
{
  
  niter   <- nburnin + nmc
  effsamp <- (niter - nburnin)/thin
  
  n <- nrow(X)
  p <- ncol(X)
  
  time         <- ct[, 1]
  status       <- ct[, 2]
  censored.id  <- which(status == 0)
  n.censored   <- length(censored.id)  # number of censored observations
  n.observed   <- n - n.censored
  X.censored   <- X[censored.id, ]
  X.observed   <- X[-censored.id, ]
  y <- logtime <- log(time)   # for coding convenience, since the whole code is written with y
  
  Gamma <- basis(y, knots)
  Hamma <- hasis(y, knots)
  N     <- ncol(Hamma)
  
  
  beta      <- rep(0, p)
  # beta      <- beta.t
  lambda    <- rep(1, p)
  tau       <- 1
  lambda.sq <- 1^2
  tau.sq    <- rep(1, p)^2
  sigma_sq  <- 1
  r         <- 1
  delta     <- 1.78
  # gamma     <- c(a * knots[1]^5 + b, 5 * a * knots^4)
  # gamma     <- sample(c(-4, 4), size = N, replace = TRUE) + gamma
  gamma     <- c(1, rep(1, length(knots)))
  gamma.int <- gamma
  s2.gamma  <- 0.15
  adapt_par <- c(10, 50, 0.5, 0.95, 0.05)
  # The first element determines after which iteration to begin adaptation. The second gives the frequency
  # with which updating occurs. The third gives the proportion of the previous states to include when
  # updating (by default 1/2). Finally, the fourth element indicates when to stop adapting
  # (default after 75% of the iterations). The fifth is the porposrtion of mixing used for new paprameter sampling
  # see ROberts and Rosenthal (2009) JCGS
  
  
  ncluster    <- detectCores()    # number of clusters
  lower.bound <- c(-Inf, rep(0, (N - 1)))
  I_n         <- diag(n)
  l0          <- rep(0, p)
  l1          <- rep(1, n)
  l2          <- rep(1, p)
  Q_star      <- t(X) %*% X
  
  
  if(is.null(Xtest))
  {
    Xtest <- X
    ntest <- n
    cttest<- ct
  } else {
    ntest <- nrow(Xtest)
  }
  
  betaout     <- matrix(0, p, effsamp)
  sigmaSqout  <- rep(0, effsamp)
  gammaout    <- matrix(0, N, effsamp)
  trace       <- array(dim = c(niter, length(gamma)))
  predsurvout <- matrix(0, ntest, effsamp)
  logtimeout  <- matrix(0, ntest, effsamp)
  loglikelihood.out <- rep(0, effsamp)
  likelihood.out    <- matrix(0, n, effsamp)
  Gout        <- matrix(0, n, effsamp)
  yout        <- matrix(0, n, effsamp)
  
  
  
  method.mcmc <- match.arg(method.mcmc)
  prior       <- match.arg(prior)
  
  
  if(method.mcmc == "MH")
  {
    cl <- makeCluster(ncluster)  # number of clusters
  }
  
  
  
  
  for (i in 1:niter)
  {
    
    #     ## Sample sigma square
    #     E_1 = max(t((Gamma %*% gamma) - X %*% beta) %*% ((Gamma %*% gamma) - X %*% beta), (1e-10))
    #     # E_2 = max(t(beta) %*% beta), (1e-10))
    # # E_2 = max(sum(beta^2/((tau * lambda))^2), (1e-10))
    # E_2 = max(sum(beta^2/tau.sq), (1e-10))
    #     sigma_sq = 1/stats::rgamma(1, (n + p)/2, scale = 2/(E_1 + E_2))
    
    
    
    ## Impute the Censored data ##
    mean.impute <- X.censored %*% beta
    sd.impute   <- sqrt(sigma_sq)
    ## update censored data ##
    xi.censored <- rnorm(n.censored, mean = mean.impute, sd = sd.impute)  # response
    # xi.observed <- rnorm(n.observed, mean = X.observed %*% beta, sd = sd.impute)
    ## Impute the Censored data using inverse function ##
    time.censored <- rep(0, n.censored)
    # time.censored <- (1/a)^(1/3) * (xi - b)^(1/3)
    
    # interpolation.fit <- sortedXyData(y[-censored.id], xi.observed)  # sorted according to y
    # for(j in 1:n.censored)
    # {
    #   time.censored[j] <- NLSstClosestX(interpolation.fit, xi.censored[j])
    # }
    
    G <- Gamma %*% gamma
    interpolation.fit <- sortedXyData(y[-censored.id], G[-censored.id])  # sorted according to y
    for(j in 1:n.censored)
    {
      time.censored[j] <- NLSstClosestX(interpolation.fit, xi.censored[j])
    }
    
    y[censored.id] <- time.censored
    Gamma <- basis(y, knots)  # construct basis matrix with the imputed data
    Hamma <- hasis(y, knots)
    
    
    
    ## Sample beta
    # Sigma.beta <- sigma_sq * diag(p)
    # A          <- t(X) %*% X + chol2inv(chol(Sigma.beta))
    # Ainv       <- chol2inv(chol(A))
    # Sigma.beta <- sigma_sq * Ainv
    # mean.beta  <- as.vector(Ainv %*% t(X) %*% (Gamma %*% gamma))
    #
    # beta <- as.vector(rmvnorm(n = 1, mean = mean.beta, sigma = Sigma.beta))
    
    
    if(prior == "horseshoe")
    {
      ## Update beta according to horseshoe
      if (p > n)
      {
        lambda_star = tau * lambda
        U = as.numeric(lambda_star^2) * t(X)
        u = stats::rnorm(l2, l0, lambda_star)
        v = X %*% u + stats::rnorm(n)
        v_star = solve((X %*% U + I_n), (((Gamma %*% gamma)/sqrt(sigma_sq)) - v))
        beta = sqrt(sigma_sq) * (u + U %*% v_star)
      }
      else if (p < n)
      {
        lambda_star = tau * lambda
        L = chol((1/sigma_sq) * (Q_star + diag(1/as.numeric(lambda_star^2), p, p)))
        v = solve(t(L), t(t(Gamma %*% gamma) %*% X)/sigma_sq)
        mu = solve(L, v)
        u = solve(L, stats::rnorm(p))
        beta = mu + u
      }
      # print(as.vector(beta))
      
      
      ## Sample lambda
      eta = 1/(lambda^2)
      upsi = stats::runif(p, 0, 1/(1 + eta))
      tempps = beta^2/(2 * sigma_sq * tau^2)
      ub = (1 - upsi)/upsi
      Fub = 1 - exp(-tempps * ub)
      Fub[Fub < (1e-04)] = 1e-04
      up = stats::runif(p, 0, Fub)
      eta = -log(1 - up)/tempps
      lambda = 1/sqrt(eta)
      
      
      
      ## sample tau
      tempt = sum((beta/lambda)^2)/(2 * sigma_sq)
      et = 1/tau^2
      utau = stats::runif(1, 0, 1/(1 + et))
      ubt = (1 - utau)/utau
      Fubt = stats::pgamma(ubt, (p + 1)/2, scale = 1/tempt)
      Fubt = max(Fubt, 1e-08)
      ut = stats::runif(1, 0, Fubt)
      et = stats::qgamma(ut, (p + 1)/2, scale = 1/tempt)
      tau = 1/sqrt(et)
      
    } else if (prior == "laplace") {
      
      
      ## Update beta according to Bayesian lasso
      if (p > n)
      {
        lambda_star = tau.sq
        U = as.numeric(lambda_star) * t(X)
        u = stats::rnorm(l2, l0, lambda_star)
        v = X %*% u + stats::rnorm(n)
        v_star = solve((X %*% U + I_n), (((Gamma %*% gamma)/sqrt(sigma_sq)) - v))
        beta = sqrt(sigma_sq) * (u + U %*% v_star)
      }
      else if (p < n)
      {
        lambda_star = tau.sq
        L = chol((1/sigma_sq) * (Q_star + diag(1/as.numeric(lambda_star), p, p)))
        v = solve(t(L), t(t(Gamma %*% gamma) %*% X)/sigma_sq)
        mu = solve(L, v)
        u = solve(L, stats::rnorm(p))
        beta = mu + u
      }
      
      
      ## Update tau^2
      tau.sq <- 1/rinvGauss(p, nu = sqrt(lambda.sq * sigma_sq/beta^2), lambda = lambda.sq)
      
      ## Update lambda^2
      lambda.sq <- rgamma(1, shape = r + p, rate = delta + 0.5 * sum(tau.sq))
      
    }
    
    
    
    
    
    
    if(method.mcmc == "MH")
    {
      ######################################################################
      ## Update gamma using Metropolis Hastings; use directional sampling ##
      ######################################################################
      
      prop_sigma <- rep(s2.gamma, N)
      trace[i, ] <- gamma
      
      # adaptive markov chain according to adaptive parameters
      # reference: Roberts and Rosenthal (2009) Examples of Adaptive MCMC, Journal of
      #            Computational and Graphical Statistics, 18(2), 349 -- 367.
      if (i > adapt_par[1] && i%%adapt_par[2] == 0 && i < (adapt_par[4] * niter)) {
        len <- floor(i * adapt_par[3]):i
        for(k in 1:N)
        {
          x  <- trace[len, k]
          Nl <- length(len)
          p_sigma <- 2.38^2 * (Nl - 1) * var(x)/Nl
          p_sigma <- makePositiveDefinite(p_sigma)
          if (!(0 %in% p_sigma))
            prop_sigma[k] <- p_sigma
        }
      }
      
      gamma.current <- gamma
      gamma.proposed <- rep(0, length(gamma))
      
      
      MH <- function(k)
      {
        if(runif(1) < adapt_par[5])
        {
          gamma.temp <- rnorm(n = 1, mean = gamma.current[k], sd = sqrt(s2.gamma))
          
          gamma.proposed[k] <- gamma.temp
          gamma.proposed[-k] <- gamma.current[-k]
          
          
          # log.kernel <- min(0,
          #                   sum(dnorm(Gamma %*% gamma.proposed, mean = X %*% beta, sd = sqrt(sigma_sq),
          #                             log = TRUE)) -
          #                     sum(dnorm(Gamma %*% gamma.current, mean = X %*% beta, sd = sqrt(sigma_sq),
          #                               log = TRUE)) +
          #                     sum(log(abs(Hamma %*% gamma.proposed))) -
          #                     sum(log(abs(Hamma %*% gamma.current))) +
          #                     dmvnorm(gamma.proposed, sigma = CM, log = TRUE) -
          #                     dmvnorm(gamma.current, sigma = CM, log = TRUE) +
          #                     log(1) - sum(log(1 + exp(-2 * gamma.proposed[-1]))) -
          #                     log(1) + sum(log(1 + exp(-2 * gamma.current[-1]))) +
          #                     dnorm(gamma.current[k], mean = gamma.proposed[k], sd = sqrt(s2.gamma), log = TRUE) -
          #                     dnorm(gamma.proposed[k], mean = gamma.current[k], sd = sqrt(s2.gamma), log = TRUE)
          # )
          
          log.kernel <- min(0,
                            sum(status * dnorm(Gamma %*% gamma.proposed, mean = X %*% beta,
                                               log = TRUE)) -
                              sum(status * dnorm(Gamma %*% gamma.current, mean = X %*% beta,
                                                 log = TRUE)) +
                              sum((1 - status) * pnorm(Gamma %*% gamma.proposed, mean = X %*% beta,
                                                       lower.tail = FALSE, log.p = TRUE)) -
                              sum((1 - status) * pnorm(Gamma %*% gamma.current, mean = X %*% beta,
                                                       lower.tail = FALSE, log.p = TRUE)) +
                              sum(log(abs(Hamma %*% gamma.proposed))) -
                              sum(log(abs(Hamma %*% gamma.current))) +
                              dmvnorm(gamma.proposed, sigma = CM, log = TRUE) -
                              dmvnorm(gamma.current, sigma = CM, log = TRUE) +
                              log(1) - sum(log(1 + exp(-20 * gamma.proposed[-1]))) -
                              log(1) + sum(log(1 + exp(-20 * gamma.current[-1]))) +
                              dnorm(gamma.current[k], mean = gamma.proposed[k], sd = sqrt(s2.gamma), log = TRUE) -
                              dnorm(gamma.proposed[k], mean = gamma.current[k], sd = sqrt(s2.gamma), log = TRUE)
          )
          
        } else {
          gamma.temp <- rnorm(n = 1, mean = gamma.current[k], sd = sqrt(prop_sigma[k]))
          
          gamma.proposed[k] <- gamma.temp
          gamma.proposed[-k] <- gamma.current[-k]
          
          
          # log.kernel <- min(0,
          #                   sum(dnorm(Gamma %*% gamma.proposed, mean = X %*% beta, sd = sqrt(sigma_sq),
          #                             log = TRUE)) -
          #                     sum(dnorm(Gamma %*% gamma.current, mean = X %*% beta, sd = sqrt(sigma_sq),
          #                               log = TRUE)) +
          #                     sum(log(abs(Hamma %*% gamma.proposed))) -
          #                     sum(log(abs(Hamma %*% gamma.current))) +
          #                     dmvnorm(gamma.proposed, sigma = CM, log = TRUE) -
          #                     dmvnorm(gamma.current, sigma = CM, log = TRUE) +
          #                     log(1) - sum(log(1 + exp(-2 * gamma.proposed[-1]))) -
          #                     log(1) + sum(log(1 + exp(-2 * gamma.current[-1]))) +
          #                     dnorm(gamma.current[k], mean = gamma.proposed[k], sd = prop_sigma[k], log = TRUE) -
          #                     dnorm(gamma.proposed[k], mean = gamma.current[k], sd = prop_sigma[k], log = TRUE)
          # )
          
          
          log.kernel <- min(0,
                            sum(status * dnorm(Gamma %*% gamma.proposed, mean = X %*% beta,
                                               log = TRUE)) -
                              sum(status * dnorm(Gamma %*% gamma.current, mean = X %*% beta,
                                                 log = TRUE)) +
                              sum((1 - status) * pnorm(Gamma %*% gamma.proposed, mean = X %*% beta,
                                                       lower.tail = FALSE, log.p = TRUE)) -
                              sum((1 - status) * pnorm(Gamma %*% gamma.current, mean = X %*% beta,
                                                       lower.tail = FALSE, log.p = TRUE)) +
                              sum(log(abs(Hamma %*% gamma.proposed))) -
                              sum(log(abs(Hamma %*% gamma.current))) +
                              dmvnorm(gamma.proposed, sigma = CM, log = TRUE) -
                              dmvnorm(gamma.current, sigma = CM, log = TRUE) +
                              log(1) - sum(log(1 + exp(-20 * gamma.proposed[-1]))) -
                              log(1) + sum(log(1 + exp(-20 * gamma.current[-1]))) +
                              dnorm(gamma.current[k], mean = gamma.proposed[k], sd = prop_sigma[k], log = TRUE) -
                              dnorm(gamma.proposed[k], mean = gamma.current[k], sd = prop_sigma[k], log = TRUE)
          )
          
          
          
        }
        
        gamma.out.par <- gamma.current[k]
        
        if(log(runif(1)) < log.kernel)
        {
          gamma.out.par <- gamma.proposed[k]
        }
        
        return(gamma.out.par)
        
      }
      
      
      registerDoParallel(cl)
      gamma <- foreach(k = 1:N, .combine = c,
                       .packages = c("mvtnorm", "msm", "MHadaptive"),
                       .export = c("gamma.proposed", "CM")
      ) %dopar% MH(k)
      
    } else if(method.mcmc == "slice"){
      ##################################################
      ## Update gamma using elliptical slice sampling ##
      ##################################################
      
      loglikelihood <- function(gamma)
      {
        # logl <- sum(dnorm(Gamma %*% gamma, mean = X %*% beta, sd = sqrt(sigma_sq),
        #               log = TRUE)) + sum(log(abs(Hamma %*% gamma))) - sum(log(1 + exp(-2 * gamma[-1])))
        logl <- sum(status * dnorm(Gamma %*% gamma, mean = X %*% beta,
                                   log = TRUE)) +
          sum((1 - status) * pnorm(Gamma %*% gamma, mean = X %*% beta,
                                   lower.tail = FALSE, log.p = TRUE)) + sum(log(abs(Hamma %*% gamma))) -
          sum(log(1 + exp(-2 * gamma[-1])))
        return(logl)
      }
      
      gamma.current <- gamma
      rsample <- mvrnorm(n = 1, mu = rep(0, N), Sigma = CM)
      theta <- runif(n = 1, max = 2 * pi)
      theta.min <- theta - 2 * pi
      theta.max <- theta
      llk.threshold <- loglikelihood(gamma.current) + log(runif(1))
      gamma.proposed <- gamma.current * cos(theta) + rsample * sin(theta)
      while(loglikelihood(gamma.proposed) < llk.threshold)
      {
        if(theta < 0)
        {
          theta.min <- theta
        } else {
          theta.max <- theta
        }
        theta <- runif(n = 1, min = theta.min, max = theta.max)
        gamma.proposed <- gamma.current * cos(theta) + rsample * sin(theta)
      }
      
      gamma <- gamma.proposed
    }
    
    
    
    
    mean.t <- Xtest %*% beta
    predictive.survivor <- pnorm(basis(log(cttest[, 1]), knots) %*% gamma, mean = mean.t, 
                                 lower.tail = FALSE)
    
    # Survival time prediction
    logt <- rep(0, ntest)
    interpolation.fit <- sortedXyData(y, G)  # sorted according to y
    for(j in 1:ntest)
    {
      logt[j] <- NLSstClosestX(interpolation.fit, mean.t[j])
    }
    
    
    loglikelihood <- sum(status * dnorm(Gamma %*% gamma, mean = X %*% beta, log = TRUE)) +
      sum((1 - status) * pnorm(Gamma %*% gamma, mean = X %*% beta, lower.tail = FALSE, log.p = TRUE)) +
      sum(log(abs(Hamma %*% gamma))) 
    likelihood    <- exp(loglikelihood)
    
    
    
    if (i%%100 == 0)
    {
      print(i)
      # print(loglikelihood)
    }
    
    
    if (i > nburnin && i%%thin == 0)
    {
      betaout[, (i - nburnin)/thin]     <- beta
      sigmaSqout[(i - nburnin)/thin]    <- sigma_sq
      gammaout[, (i - nburnin)/thin]    <- gamma
      predsurvout[ ,(i - nburnin)/thin] <- predictive.survivor
      logtimeout[, (i - nburnin)/thin]  <- logt
      loglikelihood.out[(i - nburnin)/thin] <- loglikelihood
      likelihood.out[, (i - nburnin)/thin]  <- likelihood
      Gout[, (i - nburnin)/thin]        <- G 
      yout[, (i - nburnin)/thin]        <- y
    }
  }
  
  # if(method.mcmc == "MH")
  # {
  #   stopCluster(cl)
  # }
  
  pMean    <- apply(betaout, 1, mean)
  pMedian  <- apply(betaout, 1, stats::median)
  pSigma   <- mean(sigmaSqout)
  pgamma   <- apply(gammaout, 1, mean)
  pPS      <- apply(predsurvout, 1, mean)
  pLogtime <- apply(logtimeout, 1, mean)
  pLoglikelihood <- mean(loglikelihood.out)
  plikelihood    <- apply(likelihood.out, 1, mean)
  pG       <- apply(Gout, 1, mean)
  py       <- apply(yout, 1, mean)
  
  alpha        <- 0.05
  left         <- floor(alpha * effsamp/2)
  right        <- ceiling((1 - alpha/2) * effsamp)
  gammaSort    <- apply(gammaout, 1, sort, decreasing = F)
  left.points  <- gammaSort[left, ]
  right.points <- gammaSort[right, ]
  
  
  Gamma <- basis(py, knots)  
  Hamma <- hasis(py, knots)
  
  
  loglikelihood.posterior <- sum(status * dnorm(Gamma %*% pgamma, mean = X %*% pMean, log = TRUE)) +
    sum((1 - status) * pnorm(Gamma %*% pgamma, mean = X %*% pMean, lower.tail = FALSE, log.p = TRUE)) +
    sum(log(abs(Hamma %*% pgamma)))
  
  
  DIC  <- -4 * pLoglikelihood + 2 * loglikelihood.posterior
  lppd <- sum(log(plikelihood))
  WAIC <- -2 * (lppd - 2 * (loglikelihood.posterior - pLoglikelihood))
  
  
  
  result = list(BetaHat = pMean, BetaMedian = pMedian, Sigma2Hat = pSigma, gammaHat = pgamma,
                BetaSamples = betaout, Sigma2Samples = sigmaSqout, gamma.initial = gamma.int,
                gamma.samples = gammaout, SurvivalHat = pPS, LogTimeHat = pLogtime, DIC = DIC,
                WAIC = WAIC, LeftCI = left.points, RightCI = right.points)
  return(result)
}
