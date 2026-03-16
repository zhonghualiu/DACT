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
p_a = runif(1e4) ## a p-value vector for the exposure-mediator associations
p_b = runif(1e4) ## a p-value vector for the mediator-outcome associations
res = DACT(p_a=p_a,p_b=p_b,correction="JC") ## a p-value vector for mediation effect testing
```

# License 
This software is licensed under MIT. 

## References

1. Liu, Z., Shen, J., Barfield, R., Schwartz, J., Baccarelli, A., & Lin, X. (2022).  
   Large-Scale Hypothesis Testing for Causal Mediation Effects with Applications in Genome-Wide Epigenetic Studies.  
   *Journal of the American Statistical Association*, 117(537), 67–81.  
   https://doi.org/10.1080/01621459.2021.1914634

2. Yang, H., Liu, Z., Wang, R., Lai, E., Schwartz, J., Baccarelli, A., Huang, Y., & Lin, X. (2025).  
   Causal Mediation Analysis for Integrating Exposure, Genomic and Phenotype Data.  
   *Annual Review of Statistics and Its Application* (Invited Review Paper).  
   https://doi.org/10.1146/annurev-statistics-040622-031653
