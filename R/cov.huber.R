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
                      D=NULL)
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
        mu <- meanfunc(Lt, Ly, method = 'Huber')   # Huber option
    }
    
    if(is.function(mu)) mu.hat <- lapply(Lt,mu)
    else mu.hat <- predict(mu, unlist(Lt))       # Huber option
    
    if(is.null(sig2e)) sig2e <- sigma2(Lt,Ly)
    
    if(is.null(sig2x))
    {
        sig2x <- varfunc(Lt,Ly,mu=mu,sig2=sig2e)
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
                 method='SP')
    class(rslt) <- 'covfunc'
    
    if(!is.null(newt))
        rslt$fitted <- predict(rslt,newt)
    
    return(rslt)
}

