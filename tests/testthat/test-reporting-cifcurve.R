test_that("survfit summaries preserve probability SEs, limits and requested-time support", {
  d <- data.frame(t = rep(1:8, 2), e = c(1,0,2,1,0,1,0,0,0,1,2,0,1,0,1,0),
                  w = c(1,2,1,.5,1,2,1,0,1,1,2,1,.5,1,2,1), g = rep(c("a","b"), each = 8))
  times <- c(0, 2, 5, 7, 7.5, 8)
  for (stratified in c(FALSE, TRUE)) for (err in c("if", "fay", "delta")) {
    formula <- if (stratified) Event(t,e)~g else Event(t,e)~1
    args <- list(formula = formula, data = d, outcome.type = "competing-risk", error = err)
    if (err != "delta") args$weights <- d$w
    for (dist in if (err == "delta") "normal" else c("normal", "t")) {
      args$conf.dist <- dist
      args$report.survfit.std.err <- FALSE
      f <- do.call(cifcurve, args)
      expect_s3_class(f, "survfit")
      expect_identical(f$logse, FALSE)
      if (!stratified) expect_null(f$strata)
      s <- summary(survival::survfit0(f), times = times, extend = TRUE)
      for (logse in c(FALSE, TRUE)) {
        args$report.survfit.std.err <- logse
        fs <- do.call(cifcurve, args)
        ss <- summary(survival::survfit0(fs), times = times, extend = TRUE)
        expect_equal(ss$std.err, s$std.err, tolerance = 0)
        expect_equal(ss$lower, s$lower, tolerance = 0)
        expect_equal(ss$upper, s$upper, tolerance = 0)
      }
      counts <- if (length(f$strata)) f$strata else length(f$time)
      indices <- split(seq_along(f$time), rep(seq_along(counts), counts))
      for (g in seq_along(indices)) {
        ix <- indices[[g]]; j <- findInterval(times, f$time[ix]); has <- j > 0
        rows <- (g-1L) * length(times) + seq_along(times)
        expect_equal(s$std.err[rows[has]], f$std.err[ix[j[has]]], tolerance = 0)
        expect_equal(s$df[rows[has]], f$df[ix[j[has]]], tolerance = 0)
        support <- f$reporting$curves[[g]]$support
        supported <- has & times <= support
        expect_equal(s$lower[rows[supported]], f$lower[ix[j[supported]]], tolerance = 0)
        expect_true(all(is.na(s$lower[rows[times > support]])))
        expect_true(all(s$conf.status[rows[times > support]] == "outside_observed_followup"))
      }
    }
  }
})

test_that("curve subsetting retains the selected stratum's inference", {
  d <- data.frame(t = rep(1:6,2), e = c(1,0,1,0,1,0,0,2,1,0,2,1), g = rep(c("a","b"),each=6))
  f <- cifcurve(Event(t,e)~g, d, outcome.type="competing-risk", error="fay", conf.dist="t")
  for (i in list(2L, "g=b", c(2L,1L))) {
    child <- f[i]
    ss <- summary(child, times=c(2,4,6), extend=TRUE)
    labels <- if (is.character(i)) i else names(f$strata)[i]
    parent <- summary(f, times=c(2,4,6), extend=TRUE)
    rows <- unlist(lapply(labels, function(label) which(as.character(parent$strata)==label)))
    expect_equal(ss$std.err, parent$std.err[rows], tolerance=0)
    expect_equal(ss$df, parent$df[rows], tolerance=0)
    expect_equal(ss$lower, parent$lower[rows], tolerance=0)
  }
  native <- summary(f, times = c(2,4,6), extend = TRUE)
  scaled <- summary(f, times = c(2,4,6), scale = 2, extend = TRUE)
  expect_equal(scaled$time, native$time / 2, tolerance = 0)
  expect_equal(scaled$lower, native$lower, tolerance = 0)
  expect_equal(scaled$std.err, native$std.err, tolerance = 0)
  expect_equal(scaled$df, native$df, tolerance = 0)
  positional <- summary(f, c(2,4,6), FALSE, 2, extend = TRUE)
  expect_equal(positional$df, scaled$df, tolerance = 0)
  expect_equal(positional$lower, scaled$lower, tolerance = 0)
})

test_that("decimal times retain native inference when summary units change", {
  d <- data.frame(t=c(.1,1.1,2.1,3.1), e=c(1,1,0,0))
  times <- c(0, .1-.Machine$double.eps, .1, 1.1, 3.1,
             3.1+4*.Machine$double.eps)
  for (stratified in c(FALSE, TRUE)) for (dist in c("normal", "t")) {
    data <- if (stratified) rbind(transform(d,g="a"),transform(d,g="b")) else d
    formula <- if (stratified) Event(t,e)~g else Event(t,e)~1
    f <- cifcurve(formula, data, outcome.type="survival", error="if", conf.dist=dist)
    original <- serialize(f,NULL)
    native <- summary(f,times=times,extend=TRUE)
    scaled <- summary(f,times=times,extend=TRUE,scale=365.25)
    positional <- summary(f,times,FALSE,365.25,extend=TRUE)
    framed <- summary(f,times=times,extend=TRUE,scale=365.25,data.frame=TRUE)
    fields <- c("surv","std.err","lower","upper","df","df.reason",
                "critical.value","conf.status","n.risk","n.event","n.censor")
    for (nm in fields) {
      expect_identical(scaled[[nm]],native[[nm]])
      expect_identical(positional[[nm]],native[[nm]])
      expect_identical(framed[[nm]],native[[nm]])
    }
    expect_identical(scaled$time,native$time/365.25)
    # Without requested times, survival chooses the native event rows.
    events <- summary(f)
    scaled.events <- summary(f,scale=365.25)
    for (nm in fields) expect_identical(scaled.events[[nm]],events[[nm]])
    # Preserve survival's scaled aggregate table and its rmean conventions.
    parent <- utils::getS3method("summary","survfit")
    expected <- parent(f,times=times,extend=TRUE,scale=365.25)
    expect_identical(scaled$table,expected$table)
    expect_identical(scaled$rmean.endtime,expected$rmean.endtime)
    expect_identical(serialize(f,NULL),original)
  }
})

test_that("tbl_survfit uses saved CIF intervals and exposes no internal diagnostics", {
  skip_if_not_installed("gtsummary", "2.6.1")
  d <- data.frame(t=c(1,2,4,6,8), e=c(1L,2L,1L,0L,0L), w=c(1,1,1,1,0))
  for (err in c("if","fay")) for (dist in c("normal","t")) {
    f <- cifcurve(Event(t,e)~1,d,weights=w,outcome.type="competing-risk",error=err,conf.dist=dist)
    counts <- summary(f, times=c(1.5,2,5,7), extend=TRUE)
    expect_equal(as.numeric(counts$n.risk), c(3,3,1,0), tolerance=0)
    # Includes the competing event at time 2; counts refer to all causes.
    expect_equal(as.numeric(counts$n.event), c(1,1,1,0), tolerance=0)
    expect_silent(tbl <- gtsummary::tbl_survfit(f,times=c(2,5,7),type="risk"))
    ard <- cardx::ard_survival_survfit(f,times=c(2,5,7),type="risk")
    extract <- function(nm) as.numeric(unlist(ard$stat[ard$stat_name==nm]))
    expect_equal(extract("estimate")[1:2],1-f$surv[findInterval(c(2,5),f$time)],tolerance=0)
    expect_equal(extract("conf.low")[1:2],1-f$upper[findInterval(c(2,5),f$time)],tolerance=0)
    expect_equal(extract("conf.high")[1:2],1-f$lower[findInterval(c(2,5),f$time)],tolerance=0)
    expect_true(is.na(tail(extract("conf.low"),1)))
    expect_false(any(c("df","critical.value","reporting","fay.bound") %in% names(tbl$table_body)))
  }
  f <- cifcurve(Event(t,e)~1,d,weights=w,outcome.type="competing-risk",error="if",conf.type="none")
  expect_null(f$lower)
  s <- summary(f,times=c(2,5),extend=TRUE)
  expect_true(all(is.na(s$lower)))
  expect_true(all(s$conf.status=="not_requested"))
  expect_silent(gtsummary::tbl_survfit(f,times=c(2,5),type="risk",statistic="{estimate}"))
})
