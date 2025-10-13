
createTestData <- function(n, w, first_zero=FALSE, last_zero=FALSE, subset_present=FALSE, logical_strata=FALSE, na_strata=FALSE) {
  one <- rep(1, n)
  t <- c(1:(n/2), 1:(n/2))
  epsilon <- rep(1, n)
  epsilon[2] <- 2
  epsilon[3] <- 2
  if (first_zero==TRUE) {
    epsilon[1] <- 0
    epsilon[n/2+1] <- 0
  }
  if (last_zero==TRUE) {
    epsilon[n/2] <- 0
    epsilon[n] <- 0
  }
  w <- rep(w, n)
  if (logical_strata==TRUE) {
    strata <- (t %% 2 == 0)
  } else {
    strata <- as.factor((t %% 2 == 0))
  }
  if (na_strata==TRUE) {
    strata[1] <- NA
  }
  subset <- rep(1, n)
  if (subset_present==TRUE) {
    subset[1] <- 0
  }
  d <- as.numeric(epsilon>0)
  return(data.frame(id = 1:n, t = t, epsilon = epsilon, d = d, w = w, strata = strata, subset=subset))
}

fcif <- function(cif,rr,type=c("cif","cloglog","logistic")) {
  mcif <- max(cif[,2])
  if (type[1]=="cif") mcif <- mcif*rr
  if (type[1]=="cloglog") mcif <- 1- exp(-mcif*rr)
  if (type[1]=="logistic") mcif <- mcif* rr/(1 + mcif * rr)
  return(mcif)
}

simRR <- function(n,cif1,cif2,lrr1,lrr2,cens=NULL,type1=c("cif","cloglog","logistic"),type2=c("cif","cloglog","logistic")) {
  A <- rbinom(n,1,0.5)
  L <- rbinom(n,1,0.5)
  rr1 <- exp(cbind(A,L) %*% lrr1)
  rr2 <- exp(cbind(A,L) %*% lrr2)
  f1 <- fcif(cif1,max(rr1),type=type1[1])
  f2 <- fcif(cif2,max(rr2),type=type2[1])
  mmm<- f1+f2
  mcif1 <- fcif(cif1,rr1,type=type1[1])
  mcif2 <- fcif(cif2,rr2,type=type2[1])
  if (mmm>1) warning(" models not satisfying sum <=1\n")
  T1 <- simsubdist(cif1,rr1,type=type1[1])
  T2 <- simsubdist(cif2,rr2,type=type2[1])
  dies <- rbinom(n,1,mcif1+mcif2)
  sel1 <- rbinom(n,1,mcif2/(mcif1+mcif2))+1
  epsilon  <- dies*(sel1)
  T1$epsilon <- epsilon
  T1$A <- A
  T1$L <- L
  T1$time <- T1$timecause
  T1$time2 <- T2$timecause
  T1$status <- epsilon
  T1 <- dtransform(T1,time=100,epsilon==0)
  T1 <- dtransform(T1,status=0,epsilon==0)
  T1 <- dtransform(T1,time=time2,epsilon==2)
  T1 <- dtransform(T1,status=2,epsilon==2)
  data <- T1

  if (!is.null(cens))  {
    cc <- rexp(n)/cens
    data$status <- ifelse(data$time<cc,data$status,0)
    data$time <- pmin(data$time,cc)
  }
  return(data)
}

