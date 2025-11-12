#' Common visual arguments for CIF functions
#'
#' These arguments are shared by \code{cifplot()}, \code{cifpanel()}, and
#' \code{cifcurve()}.
#' @param type.y Optional vector/list per panel: \code{"surv"} or \code{"risk"} (display 1 - survival i.e. CIF).
#' @param label.x Character x-axis label (default \code{"Time"}).
#' @param label.y Character y-axis label (default is chosen automatically from \code{outcome.type}
#'   and \code{type.y}, e.g. \dQuote{Survival} or \dQuote{Cumulative incidence}).
#'
#' @param limits.x Numeric length-2 vector specifying x-axis limits. If \code{NULL}, it is
#'   set from the fitted object (typically \code{c(0, max(time))}).
#' @param limits.y Numeric length-2 vector specifying y-axis limits. If \code{NULL}, it is
#'   set to \code{c(0, 1)} for probability-type outcomes.
#' @param breaks.x Numeric vector of x-axis breaks (default \code{NULL}).
#' @param breaks.y Numeric vector of y-axis breaks (default \code{NULL}).
#' @param use.coord.cartesian Logical; if \code{TRUE}, uses \code{coord_cartesian()} for zooming
#'   instead of changing the scale limits (default \code{FALSE}).
#'
#' @param add.conf Logical; if \code{TRUE}, adds a confidence-interval ribbon
#'   (via \code{ggsurvfit::add_confidence_interval()}). Default \code{TRUE}.
#'
#' @param add.censor.mark Logical; if \code{TRUE}, draws censoring marks on each curve
#'   (via \code{ggsurvfit::add_censor_mark()}). Default \code{TRUE}.
#' @param shape.censor.mark Integer point shape used for censoring marks (default \code{3}).
#' @param size.censor.mark Numeric point size used for censoring marks (default \code{2}).
#'
#' @param add.competing.risk.mark Logical; if \code{TRUE}, draws time marks for the competing event
#'   (event 2). If no times are supplied via \code{competing.risk.time}, the function tries to
#'   extract them automatically from the data. Default \code{FALSE}.
#' @param competing.risk.time A **named list** of numeric vectors. Each name must correspond to a
#'   strata label, and its numeric vector gives the times at which the competing event occurred
#'   in that stratum. Typically left as \code{list()} and filled internally.
#' @param shape.competing.risk.mark Integer point shape for competing-risk marks (default \code{16}).
#' @param size.competing.risk.mark Numeric point size for competing-risk marks (default \code{2}).
#'
#' @param add.intercurrent.event.mark Logical; if \code{TRUE}, overlays user-specified intercurrent-event
#'   times per stratum. Default \code{FALSE}.
#' @param intercurrent.event.time A **named list** of numeric vectors for intercurrent events
#'   (names must match strata labels).
#' @param shape.intercurrent.event.mark Integer point shape for intercurrent-event marks
#'   (default \code{1}).
#' @param size.intercurrent.event.mark Numeric point size for intercurrent-event marks
#'   (default \code{2}).
#'
#' @param add.quantile Logical; if \code{TRUE}, adds a quantile reference line (via
#'   \code{ggsurvfit::add_quantile()}). Default \code{FALSE}.
#' @param level.quantile Numeric quantile level to be shown (default \code{0.5} for the median).
#'
#' @param rows.columns.panel Optional integer vector \code{c(nrow, ncol)} controlling
#'   the layout of the panel returned by the panel modes. If \code{NULL}, an automatic
#'   layout is determined from the number of subplots.
#'
#' @param style Character choosing the base plot style: \code{"classic"}, \code{"bold"},
#' \code{"framed"}, \code{"grid"}, \code{"gray"} or \code{"ggsurvfit"} (default \code{"classic"}).
#'   Abbreviations such as \code{"C"}, \code{"B"}, \code{"F"}, or \code{"G"} are also accepted.
#' @param palette Optional character vector specifying the color palette to use across strata.
#' @param linewidth Optional numeric specifying the line width of curve (default \code{0.8}).
#' @param linetype Optional logical using different line types of curve (default \code{FALSE}).
#' @param font.family Character specifying the font family: \code{"sans"},  \code{"serif"}, or
#' \code{"mono"} (default \code{"sans"}).
#' @param font.size Integer specifying the base font size (default \code{12}).
#' @param print.panel Logical. When \code{TRUE}, panel displays created internally are
#'   printed automatically in interactive sessions; otherwise they are returned
#'   invisibly for further modification (default \code{FALSE}).
#' @param filename.ggsave Character; if non-\code{NULL}, save the plot to this file.
#' @param width.ggsave Numeric width passed to \code{ggplot2::ggsave()} (default \code{6}).
#' @param height.ggsave Numeric height passed to \code{ggplot2::ggsave()} (default \code{6}).
#' @param dpi.ggsave Numeric DPI passed to \code{ggplot2::ggsave()} (default \code{300}).
#'
#' @name cif-visual-arguments
#' @keywords internal
NULL
