library(tidyverse)
library(Rcpp)
library(RcppArmadillo)
library(foreach)
library(doParallel)

path <- "/Users/kevin-imac/Desktop/Github - Repo/"
if(! file.exists(path)){
  path <- "/Users/kevinkvp/Desktop/Github Repo/"
}

sourceCpp(paste0(path, "HW3EM/src/main.cpp"))

### User-defined functions
meanSD <- function(x, dplace = 5){
  mm <- round(mean(x), digits = dplace)
  ss <- round(sd(x), digits = dplace)
  paste0(mm, " (SD = ", ss, ")")
}

### Simulated the data and run the model
set.seed(31082, kind = "L'Ecuyer-CMRG")
registerDoParallel(5)
resultEM <- foreach(t = 1:100) %dopar% {
  
  ### Simulate the data
  clus_ind <- rbinom(100, 1, 0.25)
  y <- rexp(100, rate = ifelse(tt == 1, 1, 2))
  
  ### Model
  em_result <- EM_rcpp(y, p0 = 0.25, lambda0 = 1, mu0 = 2, eps = 1e-10)
  
  list(data = y, result = em_result)
  
}
stopImplicitCluster()

### Mean
estParam <- sapply(1:100, 
                   function(x){c(resultEM[[x]]$result$p, resultEM[[x]]$result$lambda, 
                                 resultEM[[x]]$result$mu)})
apply(estParam, 1, meanSD)

### Bias
apply(estParam - matrix(c(0.25, 1, 2), ncol = 100, nrow = 3), 1, meanSD)

