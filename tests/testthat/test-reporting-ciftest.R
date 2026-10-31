test_that("curve-test reporting uses the selected saved test, not a replacement", {
  skip_if_not_installed("gtsummary", "2.6.1")
  d <- data.frame(t=rep(1:8,2),e=c(1,0,2,1,0,1,2,0,2,1,0,1,2,0,1,0),
                  g=rep(c("a","b"),each=8))
  fits <- list(
    ciftest(Event(t,e)~g,d,test="gray",tau=6),
    ciftest(Event(t,e)~g,d,test="score",rho=1,gamma=0,tau=6),
    ciftest(Event(t,e)~g,d,test="score",rho=0,gamma=1,tau=6))
  result <- do.call(rbind,lapply(fits,function(f) {
    td<-generics::tidy(f)
    expect_equal(td$p.value,f$p.value,tolerance=0)
    expect_equal(td$statistic,unname(f$statistic),tolerance=0)
    expect_identical(td$method,f$method)
    data.frame(method=paste(td$method,td$weight),p.value=td$p.value)
  }))
  tbl <- gtsummary::as_gtsummary(result)
  expect_identical(tbl$table_body,result)
  expect_identical(vapply(fits,`[[`,"","test"),c("gray","score","score"))
})
