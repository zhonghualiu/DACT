source("MDACT_functions.R")
options(warn = -1)

ws <- c(0.01,0.01,0.979,0.001) ### assumed four cases proportions
sim.num <- 1e5 ## testing number

combinations <- as.matrix(expand.grid(rep(list(c(0, 1)),2)))[c(2,3,1,4),]
sim.num_per_set <- round(ws * sim.num)
sets <- vector("list", length = length(ws))
current_index <- 1
for (i in 1:length(ws)) {
  if (sim.num_per_set[i] > 0) {
    sets[[i]] <- current_index:(current_index + sim.num_per_set[i] - 1)
    current_index <- current_index + sim.num_per_set[i]
  } else {
    sets[[i]] <- integer(0)  # Empty set
  }
}

true_set <- unlist(sets[4])

#### generate the p values from z scores from testing

Z.M = c()
Z.Y = c()
p.M = c()
p.Y = c()

mu_M <- sample(c(1,-1),length(c(unlist(sets[1]),unlist(sets[4]))),replace = T)*rnorm(length(c(unlist(sets[1]),unlist(sets[4]))),4,1)
mu_Y <- sample(c(1,-1),length(c(unlist(sets[2]),unlist(sets[4]))),replace = T)*rnorm(length(c(unlist(sets[2]),unlist(sets[4]))),4,1)
Z.M <- rnorm(sim.num,0,1); Z.M[c(unlist(sets[1]),unlist(sets[4]))] <- rnorm(length(c(unlist(sets[1]),unlist(sets[4]))),mu_M,1);
Z.Y <- rnorm(sim.num,0,1); Z.Y[c(unlist(sets[2]),unlist(sets[4]))] <- rnorm(length(c(unlist(sets[2]),unlist(sets[4]))),mu_Y,1);
p.M <- 2* (1 - pnorm(abs(Z.M))); p.Y <- 2* (1 - pnorm(abs(Z.Y)))
p.M[p.M==0] <- 1e-17 ; p.Y[p.Y==0] <- 1e-17

estws <- null_estimation(cbind(p.M,p.Y))

sig_level <- 0.05
control.method <- "FDR"

##MDACT_testing function is the MDACT method, return is the identified testing index

mdact_identify <- MDACT_testing(p.M,p.Y,estws,significance_upper = sig_level,control.method=control.method)

fdr <-length(which(mdact_identify %!in% true_set))/max(length(mdact_identify),1)
power <- length(which(mdact_identify %in% true_set))/length(true_set)
c(fdr,power) #FDR Power one simulation

