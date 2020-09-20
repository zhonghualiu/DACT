# DACT
This R package is used for testing mediation effects in genome-wide epigenetic studies. 
# Description 
DACT is an R package for testing mediation effects in genome-wide epigenetic studies. DACT leverages the epigenome-wide data to estimate the relative proportions of the three null cases under the composite null hypothesis. The empirical null framework has been employed for large-scale inference. 
# Installation

```
library(devtools)
devtools::install_github("https://github.com/zhonghualiu/DACT")
```

# Usage 
```
## p_a and p_b are the p-value vectors for the exposure-mediator and mediator-outcome associations 
## in genome-wide epigenetic studies.
DACT(p_a,p_b,correction="JC")
```
# Example 
```
library(DACT)
p_a = runif(1e4) ## p-value vector for the exposure-mediator associations
p_b = runif(1e4) ## p-value vector for the mediator-outcome associations
res = DACT(p_a=p_a,p_b=p_b,correction="JC")
```

# License 
This software is licensed under MIT. 

