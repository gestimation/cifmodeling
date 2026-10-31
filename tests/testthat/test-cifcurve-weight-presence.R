test_that("NULL weights forwarded by wrappers retain unweighted inference", {
  d <- data.frame(t=rep(1:6,2), e=rep(c(1,0,2,1,0,0),2), g=rep(0:1,each=6))
  forward <- function(data, outcome, method=NULL, weights=NULL)
    cifcurve(Event(t,e)~g,data,weights=weights,outcome.type=outcome,error=method)
  empty <- NULL
  fields <- c("time","surv","std.err","lower","upper","error","n.risk")
  for (outcome in c("survival","competing-risk")) {
    data <- d
    if (outcome == "survival") data$e <- as.integer(data$e > 0)
    methods <- if (outcome == "survival") c("greenwood","tsiatis") else c("aalen","delta")
    for (method in c(list(NULL),as.list(methods))) {
      direct <- cifcurve(Event(t,e)~g,data,outcome.type=outcome,error=method)
      expect_warning(wrapped <- forward(data,outcome,method),NA)
      expect_warning(named.null <- cifcurve(Event(t,e)~g,data,weights=empty,
        outcome.type=outcome,error=method),NA)
      expect_warning(computed.null <- cifcurve(Event(t,e)~g,data,weights={NULL},
        outcome.type=outcome,error=method),NA)
      for (nm in fields) {
        expect_identical(wrapped[[nm]],direct[[nm]])
        expect_identical(named.null[[nm]],direct[[nm]])
        expect_identical(computed.null[[nm]],direct[[nm]])
      }
    }
  }
})

test_that("actual supplied weights retain weighted inference and fallback", {
  d <- data.frame(t=1:8,e=c(1,0,2,1,0,1,0,0),w=1:8)
  direct <- cifcurve(Event(t,e)~1,d,weights="w",outcome.type="competing-risk")
  forward <- function(weights=NULL)
    cifcurve(Event(t,e)~1,d,weights=weights,outcome.type="competing-risk")
  for (f in list(forward("w"),forward(d$w),
                 cifcurve(Event(t,e)~1,d,weights=w,outcome.type="competing-risk"))) {
    expect_identical(f$error,"if")
    expect_identical(f$std.err,direct$std.err)
    expect_identical(f$surv,direct$surv)
  }
  # Explicitly supplied unit weights remain a weighted request.
  for (method in c("aalen","delta")) {
    expect_warning(f <- cifcurve(Event(t,e)~1,d,weights=rep(1,nrow(d)),
      outcome.type="competing-risk",error=method),"falling back to 'if'")
    expect_identical(f$error,"if")
  }
})

test_that("unweighted CIF plots retain the requested variance and default delta", {
  data(diabetes.complications)
  render <- function(plot) {
    grDevices::pdf(file=tempfile(fileext=".pdf"))
    on.exit(grDevices::dev.off())
    ggsurvfit::ggsurvfit_build(plot)
  }
  for (method in list(NULL,"aalen","delta")) {
    warnings <- character()
    p <- withCallingHandlers(
      cifplot(Event(t,epsilon)~fruitq1,data=diabetes.complications,
        outcome.type="competing-risk",label.x="Time (months)",label.y="CIF of cancer death",
        limits.x=c(0,8),breaks.x=seq(0,8,6),label.strata=c("Placebo","Estrogen"),
        level.strata=c(0,1),order.strata=c(0,1),add.risktable=FALSE,add.estimate.table=TRUE,
        error=method,conf.type="arcsine-square root",conf.int=.95,print.panel=FALSE),
      warning=function(w) {
        warnings <<- c(warnings,conditionMessage(w))
        invokeRestart("muffleWarning")
      })
    expect_length(warnings,1L)
    expect_match(warnings,"Max of `survfit_object\\$time`.*outside `limits.x`")
    fit <- p$survfit.info$survfit
    reference <- cifcurve(Event(t,epsilon)~fruitq1,data=diabetes.complications,
      outcome.type="competing-risk",error=method,conf.type="arcsine-square root",conf.int=.95)
    expect_identical(fit$error,if (is.null(method)) "delta" else method)
    for (nm in c("time","surv","std.err","lower","upper","n.risk"))
      expect_equal(fit[[nm]],reference[[nm]],tolerance=0)
    expect_warning(render(p$plot),NA)
  }
})

test_that("an unused weights column does not make a CIF plot weighted", {
  d <- data.frame(t=rep(1:6,2),e=rep(c(1,0,2,1,0,0),2),
                  g=rep(0:1,each=6),weights=seq_len(12))
  expect_warning(p <- cifplot(Event(t,e)~g,d,outcome.type="competing-risk",
    error="aalen",add.risktable=FALSE,print.panel=FALSE),NA)
  reference <- cifcurve(Event(t,e)~g,d,outcome.type="competing-risk",error="aalen")
  expect_identical(p$survfit.info$survfit$error,"aalen")
  expect_identical(p$survfit.info$survfit$std.err,reference$std.err)
  expect_warning(weighted <- cifplot(Event(t,e)~g,d,weights="weights",
    outcome.type="competing-risk",error="aalen",add.risktable=FALSE,print.panel=FALSE),
    "falling back to 'if'")
  expect_identical(weighted$survfit.info$survfit$error,"if")
})
