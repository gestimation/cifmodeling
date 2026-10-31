# Direct working-Q construction for a single multinomial hazard. This does not
# use the production block builder or trace algorithm.
one_jump_inference <- function(event, weight, corrected = FALSE, equal = FALSE) {
  keep <- weight > 0; event <- event[keep]; weight <- weight[keep]
  n <- length(event); causes <- if (any(event == 2)) 1:2 else 1L
  h <- vapply(causes, function(k) sum(weight[event == k]) / sum(weight), numeric(1))
  V <- diag(h, length(h)) - tcrossprod(h)
  bread <- sum(weight) * V; inverse <- solve(bread); contrast <- V[, 1L]
  U <- vapply(causes, function(k) weight * ((event == k) - h[k]), numeric(n))
  H <- if (corrected) 1 / sqrt(1 - pmin(weight / sum(weight), .75)) else rep(1, n)
  T <- U * H; M <- crossprod(T)
  omega <- lapply(weight, function(w) w * V)
  rho <- vapply(omega, function(o) drop(contrast %*% (solve(bread - o) - inverse) %*% contrast), numeric(1))
  alpha <- if (equal) rep(1/n, n) else rho / sum(rho)
  Q <- matrix(0, n, n)
  for (j in seq_len(n)) {
    L <- do.call(rbind, lapply(seq_len(n), function(i)
      drop((contrast %*% inverse * H[i]) %*% ((i == j) * diag(length(h)) - omega[[i]] %*% inverse))))
    Q <- Q + alpha[j] * L %*% M %*% t(L)
  }
  psi <- drop(T %*% inverse %*% contrast)
  list(variance = sum(psi^2), df = sum(diag(Q))^2 / sum(Q^2))
}

test_that("Fay and reference-distribution API names are canonical", {
  expect_true("conf.dist" %in% names(formals(cifcurve)))
  expect_false("conf.method" %in% names(formals(cifcurve)))
  d <- data.frame(t = 1:8, e = c(1, 0, 2, 1, 0, 1, 0, 0), w = c(1, 2, 1, .5, 1, 2, 1, 1))
  for (dist in c("normal", "t")) {
    ref <- cifcurve(Event(t,e)~1, d, weights = w, outcome.type = "competing-risk",
                    error = "fay", conf.dist = dist)
    expect_identical(ref$error, "fay")
    expect_identical(ref$conf.dist, dist)
    expect_identical(ref$fay.bound, .75)
    expect_false(any(c("conf.method", "fg.bound") %in% names(ref)))
    for (alias in c("fg", "fay-graubard", " FAY ")) {
      fit <- cifcurve(Event(t,e)~1, d, weights = w, outcome.type = "competing-risk",
                      error = alias, conf.dist = dist)
      expect_identical(fit$error, "fay")
      expect_equal(fit$std.err, ref$std.err, tolerance = 0)
      expect_equal(fit$df, ref$df, tolerance = 0)
      expect_equal(fit$lower, ref$lower, tolerance = 0)
      expect_equal(fit$upper, ref$upper, tolerance = 0)
    }
  }
  automatic <- cifcurve(Event(t,e)~1, d, outcome.type = "competing-risk", conf.dist = "t")
  expect_identical(automatic$error, "if")
  expect_identical(automatic$conf.dist, "t")
  expect_identical(automatic$df.method, "deletion")
  expect_error(cifcurve(Event(t,e)~1, d, conf.method = "t"), "unused argument")
  for (alias in c("fay", "fg", "fay-graubard")) {
    native <- calculateAJ_Rcpp(d$t, as.integer(d$e), error = alias)
    expect_identical(native$error, "fay")
  }
})

test_that("the exact KM case derivative equals Greenwood without censoring", {
  for (k in c(2, 10, 18)) {
    d <- data.frame(t = c(seq_len(k), rep(21, 20-k)), e = c(rep(1,k),rep(0,20-k)))
    a <- cifcurve(Event(t,e)~1, d, outcome.type="survival", error="if", report.influence.function=TRUE)
    b <- cifcurve(Event(t,e)~1, d, outcome.type="survival", error="greenwood")
    expect_equal(a$std.err, b$std.err, tolerance=1e-12)
    expect_equal(a$std.err[k]^2, (k/20)*(1-k/20)/20, tolerance=1e-12)
    expect_equal(colSums(a$influence.function[[1]]^2), a$std.err^2, tolerance=1e-12)
  }
})

test_that("IF-t and FG inference match a directly assembled working Q", {
  for (competing in c(FALSE,TRUE)) {
    e <- if(competing)c(1,1,2,2,0,0,0,0) else c(1,1,1,0,0,0,0,0)
    d <- data.frame(t=ifelse(e>0,1,2),e=e,w=c(1,2,.5,3,1,2,1,0))
    type <- if(competing)"competing-risk" else "survival"
    for (corrected in c(FALSE,TRUE)) {
      error <- if(corrected)"fay" else "if"
      ref <- one_jump_inference(e,d$w,corrected)
      fit <- cifcurve(Event(t,e)~1,d,weights=w,outcome.type=type,error=error,
                      conf.type="plain",conf.dist="t")
      expect_equal(fit$std.err[1]^2,ref$variance,tolerance=1e-12)
      expect_equal(fit$df[1],ref$df,tolerance=1e-10)
      expect_equal(fit$critical.value[1],qt(.975,ref$df),tolerance=1e-12)
      expect_identical(fit$df.method,"deletion")
      normal <- cifcurve(Event(t,e)~1,d,weights=w,outcome.type=type,error=error,conf.type="plain")
      expect_equal(normal$std.err,fit$std.err)
      expect_lte(fit$lower[1],normal$lower[1])
      if(!corrected) {
        ref <- one_jump_inference(e,d$w,equal=TRUE)
        fit <- cifcurve(Event(t,e)~1,d,weights=w,outcome.type=type,error="if",
                        conf.dist="t",df.method="equal")
        expect_equal(fit$df[1],ref$df,tolerance=1e-10)
      }
    }
  }
})

test_that("multi-time AJ influence equals numerical case perturbations", {
  d <- data.frame(t=c(.3,.5,.5,.8,1,1.1,1.4,1.6),e=c(0,1,2,1,0,2,0,1),w=c(1,2,1,.5,3,1,2,0))
  estimate <- function(w,points) {
    S<-1;F<-0;z<-numeric(length(points));knots<-sort(unique(d$t[d$e>0 & w>0]));j<-1L
    for(a in seq_along(points)) {
      while(j<=length(knots) && knots[j]<=points[a]) {
        u<-knots[j];Y<-sum(w[d$t>=u]);h<-sum(w[d$t==u & d$e>0])/Y
        F<-F+S*sum(w[d$t==u & d$e==1])/Y;S<-S*(1-h);j<-j+1L
      };z[a]<-F
    };z
  }
  fit<-cifcurve(Event(t,e)~1,d,weights=w,error="if",outcome.type="competing-risk",
                report.influence.function=TRUE)
  numerical<-vapply(seq_len(nrow(d)),function(i) {
    plus<-minus<-d$w;plus[i]<-plus[i]*(1+1e-5);minus[i]<-minus[i]*(1-1e-5)
    (estimate(plus,fit$time)-estimate(minus,fit$time))/2e-5
  },numeric(length(fit$time)))
  expect_equal(fit$influence.function[[1]],t(numerical),tolerance=1e-8)
  expect_equal(fit$std.err^2,colSums(t(numerical)^2),tolerance=1e-8)
  fg<-cifcurve(Event(t,e)~1,d,weights=w,error="fay",outcome.type="competing-risk",
               report.influence.function=TRUE)
  expect_equal(fg$influence.function,fit$influence.function)
})

test_that("new intervals respect strata, zero weights, support and boundaries", {
  d<-data.frame(t=c(1,2,3,4,5,1,2,3,4,5),e=c(1,0,1,0,2,2,2,0,0,1),
                w=c(1,2,1,1,0,1,1,1,1,0),g=rep(c("a","b"),each=5))
  for(err in c("if","fay")) {
    fit<-cifcurve(Event(t,e)~g,d,weights=w,error=err,conf.dist="t",outcome.type="competing-risk")
    separate<-lapply(split(d,d$g),function(z)cifcurve(Event(t,e)~1,z,weights=w,
      error=err,conf.dist="t",outcome.type="competing-risk"))
    expect_equal(fit$std.err,unlist(lapply(separate,`[[`,"std.err"),use.names=FALSE))
    expect_equal(fit$df,unlist(lapply(separate,`[[`,"df"),use.names=FALSE))
    expect_equal(fit$lower,unlist(lapply(separate,`[[`,"lower"),use.names=FALSE))
    expect_true(all(is.na(fit$lower[fit$time==5])))
    expect_true(all(fit$conf.status[fit$time==5]=="outside_observed_followup"))
    b<-separate$b
    expect_equal(b$surv,rep(1,length(b$time)))
    expect_equal(b$lower[b$time<=4],rep(1,4))
  }
  complete<-data.frame(t=rep(1,4),e=rep(1,4))
  fit<-cifcurve(Event(t,e)~1,complete,outcome.type="survival",conf.dist="t")
  expect_equal(fit$lower,0);expect_equal(fit$upper,0)
  expect_identical(fit$conf.status,"endpoint_degenerate")
  no.target<-calculateAJ_Rcpp(c(1,2,3),c(2L,2L,0L),error="fay")
  expect_equal(no.target$surv,rep(1,3))
  expect_equal(no.target$std.err,rep(0,3))
  expect_equal(dim(no.target$influence.function[[1]]),c(3L,3L))
})

test_that("multiple-jump complete exits have structural boundary intervals", {
  options <- list(c("if", "normal", "deletion"), c("if", "t", "deletion"),
                  c("if", "t", "equal"), c("fay", "normal", "deletion"),
                  c("fay", "t", "deletion"))
  for (n in c(5L, 10L)) for (weighted in c(FALSE, TRUE)) {
    d <- data.frame(t = seq_len(n), e = 1L,
                    w = if (weighted) seq_len(n) / 3 else rep(1, n))
    for (opt in options) for (tr in c("plain", "log", "log-log", "logit", "arcsin")) {
      args <- list(formula = Event(t,e)~1, data = d, outcome.type = "survival",
                   error = opt[1], conf.dist = opt[2], df.method = opt[3],
                   conf.type = tr, report.influence.function = TRUE)
      if (weighted) args$weights <- d$w
      f <- do.call(cifcurve, args)
      expect_identical(tail(f$surv, 1), 0)
      expect_identical(tail(f$std.err, 1), 0)
      expect_identical(tail(f$lower, 1), 0)
      expect_identical(tail(f$upper, 1), 0)
      expect_identical(tail(f$conf.status, 1), "endpoint_degenerate")
      expect_equal(f$influence.function[[1]][, n], rep(0, n), tolerance = 0)
      expect_true(all(f$surv[-n] > 0 & f$surv[-n] < 1))
    }
  }

  # Earlier censoring and a tied final exit (using the engine's time groups)
  # also give a unit CIF, even with a later zero-weight censoring row.
  censored <- data.frame(t = c(1, 2, 3, 4, 4 + 1e-13, 5),
                         e = c(0L, 1L, 0L, 1L, 1L, 0L), w = c(1, 2, .5, 3, 4, 0))
  for (opt in options) {
    f <- cifcurve(Event(t,e)~1, censored, weights = w, outcome.type = "survival",
                  error = opt[1], conf.dist = opt[2], df.method = opt[3])
    end <- which(f$time == 4)
    expect_identical(f$surv[end], 0)
    expect_identical(f$std.err[end], 0)
    expect_identical(f$lower[end], 0)
    expect_identical(f$upper[end], 0)
    expect_identical(f$conf.status[end], "endpoint_degenerate")
    expect_identical(tail(f$conf.status, 1), "outside_observed_followup")
  }

  # Zero-weight competing events do not prevent a unit CIF. Zero-weight rows
  # after complete exit must still be outside positive-weight follow-up.
  d <- data.frame(t = c(1:5, 2, 6), e = c(rep(1L, 5), 2L, 0L),
                  w = c(1:5, 0, 0), g = "complete")
  mixed <- data.frame(t = 1:7, e = c(1L, 2L, rep(1L, 5)), w = 1, g = "mixed")
  d <- rbind(d, mixed)
  for (opt in options) {
    f <- cifcurve(Event(t,e)~g, d, weights = w, outcome.type = "competing-risk",
                  error = opt[1], conf.dist = opt[2], df.method = opt[3])
    a <- seq_len(f$strata[1]); b <- f$strata[1] + seq_len(f$strata[2])
    end <- a[f$time[a] == 5]
    expect_identical(f$surv[end], 0)
    expect_identical(f$std.err[end], 0)
    expect_identical(f$lower[end], 0)
    expect_identical(f$upper[end], 0)
    expect_identical(f$conf.status[end], "endpoint_degenerate")
    expect_identical(tail(f$conf.status[a], 1), "outside_observed_followup")
    expect_true(is.na(tail(f$lower[a], 1)))
    # Complete exit with genuinely competing events leaves an interior CIF.
    expect_equal(tail(f$surv[b], 1), 1/7, tolerance = 1e-14)
    expect_gt(tail(f$std.err[b], 1), 0)
    expect_identical(tail(f$conf.status[b], 1), "ok")
  }
})

test_that("near-boundary probabilities without complete exit stay interior", {
  for (w in list(c(1, 5e-12), c(5e-12, 1))) {
    d <- data.frame(t = 1:2, e = c(1L, 0L), w = w)
    for (err in c("if", "fay")) for (method in c("normal", "t")) {
      f <- cifcurve(Event(t,e)~1, d, weights = w, outcome.type = "survival",
                    error = err, conf.dist = method)
      expect_true(all(f$surv > 0 & f$surv < 1))
      expect_true(all(f$std.err > 0))
      expect_false(any(f$conf.status %in% c("endpoint_degenerate", "invalid_estimate")))
      if (method == "normal") expect_true(all(f$upper > f$lower))
    }
  }
})

test_that("critical values share all existing interval transformations", {
  d<-data.frame(t=1:12,e=c(1,0,1,2,0,1,0,2,0,0,0,0))
  for(tr in c("plain","log","log-log","logit","arcsin","a","p")) {
    f<-cifcurve(Event(t,e)~1,d,outcome.type="competing-risk",conf.dist="t",conf.type=tr)
    ref<-calculateCI(list(surv=f$surv,std.err=f$std.err),.95,f$conf.type,
                      critical_value=f$critical.value)
    interior<-f$surv>0 & f$surv<1
    expect_equal(f$lower[interior],ref$lower[interior])
    expect_equal(f$upper[interior],ref$upper[interior])
  }
  f<-cifcurve(Event(t,e)~1,d,outcome.type="competing-risk",conf.dist="t",conf.type="none")
  expect_null(f$lower);expect_null(f$upper)
  expect_error(cifcurve(Event(t,e)~1,d,conf.dist="t",error="delta"),"require")
  expect_error(cifcurve(Event(t,e)~1,d,conf.dist="t",error="fay",df.method="equal"),"IF-t")
  expect_error(cifcurve(Event(t,e)~1,d,conf.dist="t",engine="calculateKM"),"require")
  expect_error(cifcurve(Event(t,e)~1,d,conf.dist="t",conf.int=1),"conf.int")
})
