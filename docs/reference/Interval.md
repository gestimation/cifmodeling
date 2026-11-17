# Encode an interval-censored event history

Encode an interval-censored event history

## Usage

``` r
Interval(l, r, epsilon)
```

## Arguments

- l, r:

  Numeric vectors of time bounds. An event occurred in `(l, r]`. For
  events detected at visits, use the last event-free time and the first
  detection time. Set `l = r` only for a known exact event time. For
  right censoring, use the last event-free time and set `r` to `Inf`.

- epsilon:

  Event code. The defaults in
  [`cifcurve()`](https://gestimation.github.io/cifmodeling/reference/cifcurve.md)
  are 0 for right censoring, 1 for the event of interest, and 2 for the
  competing event. For an observed event, the code identifies which
  event occurred first.

## Value

An `Interval` matrix for use as the response in
[`cifcurve()`](https://gestimation.github.io/cifmodeling/reference/cifcurve.md).
