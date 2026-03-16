library(HDMT)
`%!in%` <- Negate(`%in%`)
MDACT_testing = function(p.M,p.Y,estws,significance_upper,control.method){
  sim.num <- length(p.M)
  
  pi.01.est = max(1e-3, estws$alpha01)
  pi.10.est = max(1e-3, estws$alpha10)
  pi.00.est = max(1e-3, estws$alpha00)
  c <- pi.01.est + pi.10.est + pi.00.est
  
  pi.01.est.sd = pi.01.est/c
  pi.10.est.sd = pi.10.est/c
  pi.00.est.sd = pi.00.est/c
  
  Ts <- pi.01.est*p.M + pi.10.est*p.Y + pi.00.est*pmax(p.M,p.Y)^2
  
  ss <- MDACT_thr_adjust(p.M, p.Y,pi.10.est,pi.01.est,pi.00.est,significance_upper,control.method)
  
  return(ss)
}

MDACT_thr_adjust <- function(p.M, p.Y,pi.10.est,pi.01.est,pi.00.est,significance_upper,control.method){
  sim.num <- length(p.M)
  c <- pi.01.est + pi.10.est + pi.00.est
  pi.11.est <- max(1 - c, 0)
  
  Ts <- pi.01.est*p.M + pi.10.est*p.Y + pi.00.est*pmax(p.M,p.Y)^2
  
  F00 <- function(t){
    F_0int <- function(p.y){
      c0 <- pi.00.est * p.y^2 + pi.01.est * p.y
      c1 <- t- pi.10.est * p.y
      c2 <- (t - pi.00.est * p.y^2 - pi.10.est * p.y)/pi.01.est
      m2 <- (-pi.01.est+sqrt(pmax(0,pi.01.est^2 + 4*pi.00.est*(t-pi.10.est*p.y))))/(2*pi.00.est)
      s1 <- pmin(pmax(c2,0), 1)
      s2 <- pmin(pmax(m2,0), 1)
      s <- rep(0, length(p.y))
      s[c0>=c1] <- s1[c0>=c1]
      s[c0<c1] <- s2[c0<c1]
      s
    }
    
    x <- tryCatch(integrate(F_0int, 0, 1, rel.tol = .Machine$double.eps), error = function(e) e)
    if(!inherits(x, "error")){
      res <- integrate(F_0int, 0, 1, rel.tol = .Machine$double.eps)
    }else{
      res <- integrate(F_0int, 0, 1)
    }
    res$value
  }
  
  F01 <- function(t){
    c0 <- pi.00.est * p.Y^2 + pi.01.est * p.Y
    c1 <- t- pi.10.est * p.Y
    c2 <- (t - pi.00.est * p.Y^2 - pi.10.est * p.Y)/pi.01.est
    m2 <- (-pi.01.est+sqrt(pmax(0,pi.01.est^2 + 4*pi.00.est*(t-pi.10.est*p.Y))))/(2*pi.00.est)
    s1 <- pmin(pmax(c2,0), 1)
    s2 <- pmin(pmax(m2,0), 1)
    
    (sum(s1[c0>=c1]) + sum(s2[c0<c1])) / length(p.Y)
  }
  
  F10 <- function(t){
    c0 <- pi.00.est * p.M^2 + pi.10.est * p.M
    c1 <- t- pi.01.est * p.M
    c2 <- (t - pi.00.est * p.M^2 - pi.01.est * p.M)/pi.10.est
    m2 <- (-pi.10.est+sqrt(pmax(0,pi.10.est^2 + 4*pi.00.est*(t-pi.01.est*p.M))))/(2*pi.00.est)
    s1 <- pmin(pmax(c2,0), 1)
    s2 <- pmin(pmax(m2,0), 1)
    
    (sum(s1[c0>=c1]) + sum(s2[c0<c1])) / length(p.M)
  }
  
  if(control.method == "FDR"){
    thr <- function(t) {
      c * (1/(1-pi.11.est)*(pi.10.est/(pi.10.est+pi.11.est)*F10(t) + 
                              pi.01.est/(pi.01.est+pi.11.est)*F01(t) + 
                              (pi.00.est - pi.10.est*(pi.00.est + pi.01.est)/(pi.10.est + pi.11.est) - 
                                 pi.01.est*(pi.00.est + pi.10.est)/(pi.01.est + pi.11.est))*F00(t))) / max(mean(Ts < t), 1 / sim.num) - significance_upper
    }
    
    #t_s <- uniroot(thr, c(min(c(p.M, p.Y)), 1))$root
    sorted_index <- order(Ts) # Get the indices that sort 'Ts'
    lower_bound <- 1
    upper_bound <- length(Ts)
    
    while (upper_bound - lower_bound > 1) {
      mid_index <- floor((lower_bound + upper_bound) / 2)
      mid_value <- Ts[sorted_index[mid_index]]
      
      x <- tryCatch(thr(mid_value), error = function(e) e)
      k <- 10
      while(inherits(x, "error")){
        x <- tryCatch(thr(round(mid_value,k)), error = function(e) e)
        k <- k-1
      }
      if (x < 0) {
        lower_bound <- mid_index
      } else {
        upper_bound <- mid_index
      }
    }
    t_s <- Ts[sorted_index[lower_bound]]
    
    if(thr(t_s) < 0){
      return(which(Ts <= t_s))
    }else{
      return(which(Ts < t_s))
    }
  }
  if(control.method == "FWER"){
    thr <- function(t) {
      c * (1/(1-pi.11.est)*(pi.10.est/(pi.10.est+pi.11.est)*F10(t) + 
                              pi.01.est/(pi.01.est+pi.11.est)*F01(t) + 
                              (pi.00.est - pi.10.est*(pi.00.est + pi.01.est)/(pi.10.est + pi.11.est) - 
                                 pi.01.est*(pi.00.est + pi.10.est)/(pi.01.est + pi.11.est))*F00(t)))* sim.num - significance_upper
    }
    
    
    sorted_index <- order(Ts) # Get the indices that sort 'Ts'
    lower_bound <- 1
    upper_bound <- length(Ts)
    
    while (upper_bound - lower_bound > 1) {
      mid_index <- floor((lower_bound + upper_bound) / 2)
      mid_value <- Ts[sorted_index[mid_index]]
      
      x <- tryCatch(thr(mid_value), error = function(e) e)
      k <- 10
      while(inherits(x, "error")){
        x <- tryCatch(thr(round(mid_value,k)), error = function(e) e)
        k <- k-1
      }
      if (x < 0) {
        lower_bound <- mid_index
      } else {
        upper_bound <- mid_index
      }
    }
    t_s <- Ts[sorted_index[lower_bound]]
    
    if(thr(t_s) < 0){
      return(which(Ts <= t_s))
    }else{
      return(which(Ts < t_s))
    }
  }
  
  if(control.method == "size"){
    thr <- function(t) {
      (1/(1-pi.11.est)*(pi.10.est/(pi.10.est+pi.11.est)*F10(t) + 
                          pi.01.est/(pi.01.est+pi.11.est)*F01(t) + 
                          (pi.00.est - pi.10.est*(pi.00.est + pi.01.est)/(pi.10.est + pi.11.est) - 
                             pi.01.est*(pi.00.est + pi.10.est)/(pi.01.est + pi.11.est))*F00(t))) - significance_upper
    }
    x <- tryCatch(uniroot(thr, c(min(c(p.M, p.Y)), 1))$root, error = function(e) e)
    if(!inherits(x, "error")){
      t_s <- x
    }else{
      F_emp <- function(t) {
        F_emp <- (1/(1-pi.11.est)*(pi.10.est/(pi.10.est+pi.11.est)*F10(t) + 
                                     pi.01.est/(pi.01.est+pi.11.est)*F01(t) + 
                                     (pi.00.est - pi.10.est*(pi.00.est + pi.01.est)/(pi.10.est + pi.11.est) - 
                                        pi.01.est*(pi.00.est + pi.10.est)/(pi.01.est + pi.11.est))*F00(t)))
        return(F_emp)
      }
      
      sorted_index <- order(Ts) # Get the indices that sort 'Ts'
      lower_bound <- 1
      upper_bound <- length(Ts)
      
      while (upper_bound - lower_bound > 1) {
        mid_index <- floor((lower_bound + upper_bound) / 2)
        mid_value <- Ts[sorted_index[mid_index]]
        
        if (F_emp(round(mid_value,6)) < significance_upper) {
          lower_bound <- mid_index
        } else {
          upper_bound <- mid_index
        }
      }
      t_s <- Ts[sorted_index[lower_bound]]
    }
    
    if(F_emp(round(t_s,6)) < significance_upper){
      return(which(Ts <= t_s))
    }else{
      return(which(Ts < t_s))
    }
  }
}

index_fpn_fdr_power_func <- function(sat_index,true_set){
  if(length(true_set)==0){
    fpn <- ifelse(length(sat_index)>0,1,0)
    return(fpn)
  }else{
    fdr <-length(which(sat_index %!in% true_set))/max(length(sat_index),1)
    power <- length(which(sat_index %in% true_set))/length(true_set)
    return(c(fdr,power))
  }
}
