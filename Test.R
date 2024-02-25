library(tidyverse)
library(Rcpp)
library(RcppArmadillo)
library(foreach)
library(doParallel)
library(Deriv)

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

### Simulated the data
set.seed(31082, kind = "L'Ecuyer-CMRG")
registerDoParallel(5)
simDat <- foreach(t = 1:100) %dopar% {
  
  ### Simulate the data
  clus_ind <- rbinom(100, 1, 0.25)
  y <- rexp(100, rate = ifelse(clus_ind == 1, 1, 2))
  
  y
  
}
stopImplicitCluster()

### Run the model
set.seed(31082, kind = "L'Ecuyer-CMRG")
registerDoParallel(5)
resultEM <- foreach(t = 1:100, .combine = "rbind") %dopar% {
  
  em_result <- EM_rcpp(y = simDat[[t]], p0 = 0.25, lambda0 = 1, mu0 = 2, eps = 1e-10)
  c(em_result$p, em_result$lambda, em_result$mu)
  
}
stopImplicitCluster()

### Louis SE
set.seed(31082, kind = "L'Ecuyer-CMRG")
registerDoParallel(5)
VARLouis <- foreach(t = 1:100, .combine = "rbind") %dopar% {
  
  varMat <- iY(y = simDat[[t]], p = resultEM[t, 1], lb = resultEM[t, 2], mu = resultEM[t, 3]) - 
    iX(y = simDat[[t]], p = resultEM[t, 1], lb = resultEM[t, 2], mu = resultEM[t, 3])
  1/(diag(varMat))
  
}
stopImplicitCluster()

#### Bootstrap
set.seed(31082, kind = "L'Ecuyer-CMRG")
registerDoParallel(5)
VARBoots <- foreach(t = 1:100, .combine = "rbind") %dopar% {
  
  bootMat <- BootResult(y = simDat[[t]], p0 = 0.25, lambda0 = 1, mu0 = 2, 
                        eps = 1e-10, EMinit = resultEM[t, ], M = 1000)
  apply(bootMat, 2, var)
  
}
stopImplicitCluster()

### SEM
SEM(y = simDat[[3]], p0 = 0.24, lambda0 = 0.99, mu0 = 2.01, eps = 1e-10, EMfinal = resultEM[3, ])

### Derivative: iX
lfunction <- function(y, p, lb, mu){
  sum(log((p * lb * exp(-lb * y)) + ((1 - p) * mu * exp(-mu * y))))
}

Deriv(lfunction, c("p", "lb", "mu")) %>%
  Deriv(c("p", "lb", "mu"))

Qfunction <- function(y, db, p, lb, mu){
  sum((db * ((-lb * y) + log(p * lb))) + ((1 - db) * ((-mu * y) + log((1 - p) * mu))))
}

Deriv(Qfunction, c("p", "lb", "mu")) %>%
  Deriv(c("p", "lb", "mu"))
