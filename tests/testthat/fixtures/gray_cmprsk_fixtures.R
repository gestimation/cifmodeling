# Frozen validation values generated with cmprsk 2.2.12.
#
# This file is data, not a runtime dependency on cmprsk. Regenerate or audit it
# with .github/scripts/validate-gray-cmprsk.R.
gray_cmprsk_fixture_metadata <- list(
  cmprsk.version = "2.2.12",
  censoring.code = 0L,
  event.of.interest.code = 1L,
  competing.event.code = 2L,
  group.order = "first appearance, then retained as factor levels",
  tie.convention = "cmprsk::cuminc()/crstm simultaneous-event convention",
  gamma = 0
)

gray_cmprsk_fixtures <- list(
  list(
    id = "two-group-ties-rho0",
    cmprsk.version = "2.2.12",
    time = c(1, 2, 2, 3, 4, 5, 5, 6, 7, 8, 8, 9),
    status = c(1L, 0L, 2L, 1L, 2L, 0L, 1L, 1L, 0L, 2L, 1L, 0L),
    group = c("A", "A", "B", "B", "A", "B", "A", "B", "A", "B", "A", "B"),
    rho = 0,
    statistic = 0.59226597289993776,
    p.value = 0.44154421636513463,
    df = 1L
  ),
  list(
    id = "two-group-ties-rho1",
    cmprsk.version = "2.2.12",
    time = c(1, 2, 2, 3, 4, 5, 5, 6, 7, 8, 8, 9),
    status = c(1L, 0L, 2L, 1L, 2L, 0L, 1L, 1L, 0L, 2L, 1L, 0L),
    group = c("A", "A", "B", "B", "A", "B", "A", "B", "A", "B", "A", "B"),
    rho = 1,
    statistic = 0.52239450203103921,
    p.value = 0.46982203912643472,
    df = 1L
  ),
  list(
    id = "three-group-ties-rho0",
    cmprsk.version = "2.2.12",
    time = c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 9, 9, 10),
    status = c(1L, 2L, 0L, 1L, 2L, 1L, 0L, 2L, 1L, 0L, 2L, 1L, 0L, 1L, 2L, 1L, 0L, 2L),
    group = c("A", "B", "C", "A", "B", "C", "A", "B", "C",
              "A", "B", "C", "A", "B", "C", "A", "B", "C"),
    rho = 0,
    statistic = 3.144885237704965,
    p.value = 0.20753762740642823,
    df = 2L
  ),
  list(
    id = "three-group-ties-rho-half",
    cmprsk.version = "2.2.12",
    time = c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 9, 9, 10),
    status = c(1L, 2L, 0L, 1L, 2L, 1L, 0L, 2L, 1L, 0L, 2L, 1L, 0L, 1L, 2L, 1L, 0L, 2L),
    group = c("A", "B", "C", "A", "B", "C", "A", "B", "C",
              "A", "B", "C", "A", "B", "C", "A", "B", "C"),
    rho = 0.5,
    statistic = 3.221288553958535,
    p.value = 0.19975887256400449,
    df = 2L
  )
)
