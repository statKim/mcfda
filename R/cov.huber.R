cov.huber <- function(Lt,Ly,newt=NULL,
                      domain=NULL,
                      weig=NULL,
                      corf=NULL, # correlation function(theta,x,y)
                      mu=NULL,
                      sig2e=NULL,
                      sig2x=NULL,
                      pfunc=NULL,
                      theta0=NULL,
                      lb=NULL,
                      ub=NULL,
                      D=NULL,
                      kernel = "gauss")
{
    
    if(is.null(corf)){
        corf <- function(x,y,theta) matern(x,y,nu=theta)
        D <- 1
    }
    else
    {
        if(is.null(theta0) && is.null(D)) stop('The dimension D must be specified')
    }
    
    if(is.null(mu))
    {
        mu <- meanfunc(Lt, Ly, kernel = kernel, method = 'HUBER')   # Huber option
    }
    
    if(is.function(mu)) mu.hat <- lapply(Lt,mu)
    else mu.hat <- predict(mu, unlist(Lt))       # Huber option
    
    if(is.null(sig2e)) sig2e <- sigma2(Lt,Ly)
    
    if(is.null(sig2x))
    {
        # sig2x <- varfunc(Lt,Ly,mu=mu,sig2=sig2e)
        sig2x <- varfunc(Lt, Ly, mu = mu, sig2 = sig2e, kernel = kernel, method = "HUBER")   # Huber option
    }
    
    if(is.function(sig2x))
    {
        var.hat <- lapply(Lt, sig2x)
    }
    else
    {
        var.hat <- predict(sig2x,Lt)
    }
    
    if(is.null(domain))
    {
        t.vec <- unlist(Lt)
        domain <- c(min(t.vec),max(t.vec))
    }
    
    th.est <- estimate.theta(Lt,Ly,
                              D=D,
                              var.hat=var.hat,
                              mu.hat=mu.hat,
                              method='LS',
                              rho=corf,
                              weig=weig,
                              pfunc=pfunc,
                              theta.lb=lb,
                              theta.ub=ub,
                              theta0=theta0,
                              domain=domain)$LS
    
    rslt <- list(sig2e=sig2e,
                 theta=th.est,
                 mu.hat=mu.hat,
                 domain=domain,
                 mu=mu,
                 sig2x=sig2x,
                 rho=function(x,y) corf(x,y,th.est),
                 method='HUBER')
    class(rslt) <- 'covfunc'
    
    if(!is.null(newt))
        rslt$fitted <- predict(rslt,newt)
    
    return(rslt)
}




### Local polynomial kernel smoothing with huber loss (mean estimator)
# Lt : a list of vectors or a vector containing time points for all curves
# Ly : a list of vectors or a vector containing observations for all curves
# newt : a vector containing time points to estimate
# bw : bandwidth
# kernel : a kernel function for kernel smoothing ("epan", "gauss" are supported.)
# loss : a loss function for kernel smoothing("L2" is squared loss, "Huber" is huber loss.)
#   For loss = "Huber", it uses `rlm()` in `MASS` and fits the robust regression with Huber loss. 
#   So additional parameters of `rlm()` can be applied. (k2, maxit, ...)
local_kern_smooth <- function(Lt, Ly, newt = NULL, bw = NULL, kernel = "epanechnikov", loss = "L2", ...) {
  if (is.list(Lt) | is.list(Ly)) {
    Lt <- unlist(Lt)
    Ly <- unlist(Ly)
  }
  
  if (is.null(newt)) {
    newt <- Lt
  }
  if (is.list(newt)) {
    newt <- unlist(newt)
  }
  
  # # If `bw` is not defined, 5-fold CV is performed.
  # if (is.null(bw)) {
  #   bw <- cv.local_kern_smooth(Lt = Lt, Ly = Ly, newt = NULL, 
  #                              kernel = kernel, loss = loss, K = 5, parallel = TRUE)
  # }
  
  w <- 1/length(Lt)
  mu_hat <- sapply(newt, function(t) {
    tmp <- (Lt - t) / bw
    
    if (kernel == "epanechnikov") {
      kern <- (3/4 * w) * (1 - tmp^2)   # Epanechnikov kernel
    } else if (kernel == "gauss") {
      kern <- w * 1/sqrt(2*pi) * exp(-1/2 * tmp^2)   # gaussian kernel
    }
    
    idx <- which(kern > 0)   # non-negative values
    W <- diag(kern[idx]) / bw
    X <- matrix(1, length(idx), 2)
    X[, 2] <- Lt[idx] - t
    Y <- Ly[idx]
    
    if (loss == "L2") {   # squared loss
      # Weighted least squares
      beta <- solve(t(X) %*% W %*% X) %*% t(X) %*% W %*% Y
      
      return(beta[1, ])
    } else if (loss == "Huber") {   # huber loss
      fit <- rlm(x = sqrt(W) %*% X,
                 y = sqrt(W) %*% Y,
                 maxit = 100,
                 scale.est = "Huber",
                 ...)
      # df <- data.frame(y = sqrt(W) %*% Ly[idx],
      #                  x = sqrt(W) %*% X)
      # fit <- rlm(y ~ .-1,
      #            data = df,
      #            method = "M",
      #            maxit = 100,
      #            scale.est = "Huber",
      #            k2 = 2)
      beta <- fit$coefficients
      
      return(beta[1])
      # return(fit$s)
    }
  })
  
  return( as.numeric(mu_hat) )
}

### K-fold cross validation to find optimal bandwidth for local polynomial kernel smoother
# Lt : a list of vectors containing time points for each curve
# Ly : a list of vectors containing observations for each curve
# K : the number of folds
# bw_cand : user defined bandwidth candidates for CV
# parallel : If parallel is TRUE, it implements `foreach()` in `doParallel` for CV.
# Other parameters are same with `local_kern_smooth()`.
cv.local_kern_smooth <- function(Lt, Ly, newt = NULL, kernel = "epanechnikov", loss = "L2", K = 5, 
                                 bw_cand = NULL, parallel = FALSE, k2 = 1.345, ...) {
  if (!(is.list(Lt) & is.list(Ly))) {
    stop("Lt and Ly should be only a list type.")
  }
  
  if (is.null(bw_cand)) {
    a <- min(unlist(Lt))
    b <- max(unlist(Lt))
    bw_cand <- 10^seq(-2, 0, length.out = 20) * (b - a)/3
  }
  
  # get index for each folds
  folds <- list()
  n <- length(Lt)   # the number of curves
  fold_num <- n %/% K   # the number of curves for each folds
  fold_sort <- sample(1:n, n)
  for (k in 1:K) {
    ind <- (fold_num*(k-1)+1):(fold_num*k)
    if (k == K) {
      ind <- (fold_num*(k-1)+1):n
    }
    folds[[k]] <- fold_sort[ind]
  }
  
  # K-fold cross validation
  if (parallel == TRUE) {
    require(doParallel)
    # Parallel computing setting
    ncores <- detectCores() - 3
    cl <- makeCluster(ncores)
    registerDoParallel(cl)
    
    cv_error <- foreach(i = 1:length(bw_cand), .combine = "c", 
                        .export = c("local_kern_smooth"), .packages = c("MASS")) %dopar% {
      err <- 0
      for (k in 1:K) {
        Lt_train <- Lt[ -folds[[k]] ]
        Ly_train <- Ly[ -folds[[k]] ]
        Lt_test <- Lt[ folds[[k]] ]
        Ly_test <- Ly[ folds[[k]] ]
        
        y_hat <- local_kern_smooth(Lt = Lt_train, Ly = Ly_train, newt = Lt_test, 
                                   bw = bw_cand[i], kernel = kernel, loss = loss, ...)
        y <- unlist(Ly_test)
        if (loss == "L2") {
          err <- err + sum((y - y_hat)^2)   # squared errors
        } else if (loss == "Huber") {
          a <- abs(y - y_hat)
          err_huber <- ifelse(a > k2, k2*(a - k2/2), a^2/2)
          err <- err + sum(err_huber)
        }
      }
      
      return(err)
    }
    stopCluster(cl)
  } else {
    cv_error <- rep(0, length(bw_cand))
    for (k in 1:K) {
      Lt_train <- Lt[ -folds[[k]] ]
      Ly_train <- Ly[ -folds[[k]] ]
      Lt_test <- Lt[ folds[[k]] ]
      Ly_test <- Ly[ folds[[k]] ]
      
      for (i in 1:length(bw_cand)) {
        y_hat <- local_kern_smooth(Lt = Lt_train, Ly = Ly_train, newt = Lt_test, 
                                   bw = bw_cand[i], kernel = kernel, loss = loss, ...)
        y <- unlist(Ly_test)
        # cv_error[i] <- cv_error[i] + sum((y - y_hat)^2)   # squared errors
        if (loss == "L2") {
          cv_error[i] <- cv_error[i] + sum((y - y_hat)^2)   # squared errors
        } else if (loss == "Huber") {
          a <- abs(y - y_hat)
          err_huber <- ifelse(a > k2, k2*(a - k2/2), a^2/2)
          cv_error[i] <- cv_error[i] + sum(err_huber)
        }
      }
    }
  }
  
  bw <- list(selected_bw = bw_cand[ which.min(cv_error) ],
             cv.error = data.frame(bw = bw_cand,
                                   error = cv_error))
  
  return(bw)
}



