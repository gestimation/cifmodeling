createTestData1 <- function(n, w, first_zero=FALSE, last_zero=FALSE, subset_present=FALSE, logical_strata=FALSE, na_strata=FALSE) {
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

createTestData2 <- function() {
  data.frame(
    t       = c(5,  7,  7,  9, 10, 12, NA,  4),
    status  = c(1,  0,  2,  0,  1,  2,  1,  0),
    x       = c(0,  1,  1,  0,  0,  1,  1,  0),
    strata  = factor(c("A","A","B","B","A","B","A","B")),
    w       = c(1,  2,  1,  1,  0.5, 1,  1,  1),
    stringsAsFactors = FALSE
  )
}
