EfronCorrect = function(pval){
  z = stats::qnorm(1-pval)
  res <- locfdr(z,nulltype = 1)
  mean.emp = res$fp0["mlest","delta"]
  sd.emp = res$fp0["mlest","sigma"]
  pval.emp = stats::pnorm(z,mean = mean.emp,sd = sd.emp,lower.tail = F)
  return(pval.emp)
}
