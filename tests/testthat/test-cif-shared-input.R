test_that("shared event preparation retains transformed-missing row provenance", {
  d <- data.frame(time=1:6,status=c(0,1,2,1,2,0),x=c(1,-1,2,3,4,5),
                  group=c("a","a",NA,"b","b","b"))
  prep <- suppressWarnings(cifmodeling:::cif_prepare_event(
    Event(time,status)~log(x), d, other.variables.analyzed="group",na.action=na.exclude))
  expect_identical(prep$row.index,c(1L,4L,5L,6L))
  expect_identical(which(!prep$included),c(2L,3L))
  expect_equal(prep$exclusion.reason[c(2,3)],c("missing","missing"))
  expect_s3_class(stats::model.response(prep$frame),"Event")
  expect_s3_class(prep$na.action,"exclude")
  expect_equal(prep$t,c(1,4,5,6))
})

test_that("shared response decoding supports arbitrary causes without global option changes", {
  d <- data.frame(time=1:4,status=c(9,0,7,8))
  prior <- getOption("cifmodeling.allowed")
  prep <- cifmodeling:::cif_prepare_event(Event(time,status)~1,d)
  expect_equal(prep$epsilon,d$status)
  expect_identical(getOption("cifmodeling.allowed"),prior)
  expect_error(cifmodeling:::cif_prepare_event(
    Event(time,status,allowed=c(0,1,2))~1,d),"codes")
})

test_that("legacy CIF preparation retains its RHS grouping and row meanings", {
  d <- data.frame(time=1:6,status=c(0,1,2,1,2,0),group=c("a","a",NA,"b","b","b"))
  prep <- cifmodeling:::cif_prepare_input(Event(time,status)~group,d,
                                        outcome.type="competing-risk")
  expect_equal(prep$t,c(1,2,4,5,6))
  expect_identical(prep$row.index,c(1L,2L,4L,5L,6L))
  expect_equal(as.character(prep$strata),c("a","a","b","b","b"))
  expect_equal(prep$strata_name,"group")
  expect_equal(prep$d0,as.integer(prep$epsilon==0))
})
