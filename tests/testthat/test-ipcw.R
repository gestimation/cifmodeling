test_that("shared censoring curves retain legacy estimates and event-first reverse KM", {
  t <- c(1,1,2,2,3,4,5,6)
  status <- c(1,0,2,0,1,0,2,0)
  old <- cifmodeling:::estimate_censoring_km(t,status)
  new <- cifmodeling:::estimate_censoring_km(t,status,tie.order="event_first")
  expect_equal(old$table$survival.right,
    cifmodeling:::calculateKM(t,as.integer(status==0),error="none",count.type="numeric")$surv)
  expect_equal(new$table$survival.right[1:2],c(6/7,24/35))
  expect_equal(cifmodeling:::predict_censoring_km(new,c(1,2),side="left"),c(1,6/7))
  expect_equal(cifmodeling:::predict_censoring_km(new,c(1,2),side="right"),c(6/7,24/35))
  expect_s3_class(new,"cif_censoring_km")
  expect_error(cifmodeling:::predict_censoring_km(new,NA_real_),"finite")
})

test_that("unknown outcomes have zero weight while zero survival for known outcomes fails", {
  g <- cifmodeling:::estimate_censoring_km(c(1,2),c(0,0),tie.order="event_first")
  w <- cifmodeling:::calculate_ipcw_weights(g,c(1,3),c(TRUE,FALSE))
  expect_equal(w$weights,c(1,0))
  expect_equal(w$survival.raw,c(1,0))
  expect_equal(w$diagnostics$n.zero,1L)
  expect_false(w$diagnostics$normalized)
  expect_error(cifmodeling:::calculate_ipcw_weights(g,c(1,3),c(TRUE,TRUE)),"positivity")
})

test_that("KM score correction agrees with explicit jump derivatives", {
  t <- c(1,1,1,2,2,3,4,5,6)
  status <- c(1,0,0,2,0,1,0,2,0)
  # One administrative completion shares a censoring-status code but has no jump.
  event <- status==0;event[3] <- FALSE
  v <- pmin(t,4)
  score <- cbind(seq_along(t)/10,cos(seq_along(t)))
  for (order in c("event_first","censoring_first")) {
    km <- cifmodeling:::estimate_censoring_km(t,status,censoring.event=event,tie.order=order)
    expected <- matrix(0,length(t),ncol(score))
    for (r in which(km$hazard.table$time < max(v))) {
      jump <- km$hazard.table[r,]
      risk <- if (order=="event_first") t>jump$time | (t==jump$time & status==0) else t>=jump$time
      dN <- t==jump$time & event
      expected <- expected + outer((dN-jump$hazard*risk)/
        (jump$n.risk.censoring*(1-jump$hazard)),colSums(score[v>jump$time,,drop=FALSE]))
    }
    actual <- cifmodeling:::censoring_km_score_correction(km,t,status,
      rep("pooled",length(t)),v,score)
    expect_equal(actual,expected,ignore_attr=TRUE,tolerance=1e-12)
    expect_equal(unname(colSums(actual)),rep(0,ncol(score)),tolerance=1e-12)
  }
  km <- cifmodeling:::estimate_censoring_km(t,status,tie.order="event_first")
  one <- cifmodeling:::censoring_km_score_correction(km,t,status,
    rep("pooled",length(t)),rep(1,length(t)),score[,1,drop=FALSE])
  expect_equal(one,matrix(0,length(t),1),ignore_attr=TRUE)
  expect_error(cifmodeling:::censoring_km_score_correction(km,t,status,
    rep("pooled",length(t)),v,score[-1,,drop=FALSE]),"aligned")
  weighted <- km;weighted$weights[1] <- 2
  expect_error(cifmodeling:::censoring_km_score_correction(weighted,t,status,
    rep("pooled",length(t)),v,score),"unit sampling")
  truncated <- km;truncated$prob.truncation <- .1
  expect_error(cifmodeling:::censoring_km_score_correction(truncated,t,status,
    rep("pooled",length(t)),v,score),"truncation")
})
