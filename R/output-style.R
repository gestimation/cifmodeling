# Shared themes for cifplot and regression prediction curves.
plot_style_classic <- function(font.family = "sans", font.size = 12, legend.position = "top") {
  ggplot2::theme_classic(base_size = font.size, base_family = font.family) +
    ggplot2::theme(
      legend.position   = legend.position,
      axis.title        = ggplot2::element_text(size = font.size + 3, family = font.family),
      axis.text         = ggplot2::element_text(size = font.size + 1, family = font.family),
      legend.text       = ggplot2::element_text(size = font.size + 3, family = font.family),
      legend.background = ggplot2::element_rect(fill = "transparent", color = NA),
      panel.background  = ggplot2::element_rect(fill = "transparent", color = NA),
      panel.border      = ggplot2::element_blank()
    )
}

plot_style_bold <- function(font.family = "sans", font.size = 12, legend.position = "top") {
  ggplot2::theme_classic(base_size = font.size, base_family = font.family) +
    ggplot2::theme(
      legend.position   = legend.position,
      axis.title        = ggplot2::element_text(size = font.size + 2, face = "bold", family = font.family),
      axis.text         = ggplot2::element_text(size = font.size,     family = font.family),
      legend.text       = ggplot2::element_text(size = font.size + 2, family = font.family),
      legend.background = ggplot2::element_rect(fill = "transparent", color = NA),
      panel.background  = ggplot2::element_rect(fill = "transparent", color = NA),
      panel.border      = ggplot2::element_blank()
    )
}

plot_style_framed <- function(font.family = "sans", font.size = 12, legend.position = "top") {
  ggplot2::theme_bw(base_size = font.size, base_family = font.family) +
    ggplot2::theme(
      legend.position   = legend.position,
      axis.title        = ggplot2::element_text(size = font.size + 2, face = "bold", family = font.family),
      axis.text         = ggplot2::element_text(size = font.size,     family = font.family),
      legend.text       = ggplot2::element_text(size = font.size + 2, family = font.family),
      legend.background = ggplot2::element_rect(fill = "transparent", color = NA),
      panel.background  = ggplot2::element_rect(fill = "transparent", color = NA),
      panel.grid        = ggplot2::element_blank(),
      panel.border      = ggplot2::element_rect(color = "black", linewidth = 2)
    )
}

plot_style_grid <- function(font.family = "sans", font.size = 12, legend.position = "top") {
  ggplot2::theme_bw(base_size = font.size, base_family = font.family) +
    ggplot2::theme(
      legend.position   = legend.position,
      axis.title        = ggplot2::element_text(size = font.size + 2, face = "bold", family = font.family),
      axis.text         = ggplot2::element_text(size = font.size,     family = font.family),
      legend.text       = ggplot2::element_text(size = font.size + 2, family = font.family),
      legend.background = ggplot2::element_rect(fill = "transparent", color = NA),
      panel.background  = ggplot2::element_rect(fill = "transparent", color = NA),
      panel.border      = ggplot2::element_rect(color = "black", linewidth = 2)
    )
}

plot_style_gray <- function(font.family = "sans", font.size = 12, legend.position = "top") {
  ggplot2::theme_bw(base_size = font.size, base_family = font.family) +
    ggplot2::theme(
      legend.position   = legend.position,
      axis.title        = ggplot2::element_text(size = font.size + 2, face = "bold", family = font.family),
      axis.text         = ggplot2::element_text(size = font.size,     family = font.family),
      legend.text       = ggplot2::element_text(size = font.size + 2, family = font.family),
      legend.background = ggplot2::element_rect(fill = "transparent", color = NA),
      panel.background  = ggplot2::element_rect(fill = 'gray97'),
      panel.grid        = ggplot2::element_blank(),
      panel.border      = ggplot2::element_rect(color = "black", linewidth = 0.2)
    )
}

output_plot_style <- function(style = "classic", font.family = "sans",
                              font.size = 12, legend.position = "top") {
  style <- match.arg(style, c("classic", "bold", "framed", "grid", "gray", "minimal"))
  if (style == "minimal") return(ggplot2::theme_minimal(base_size = font.size,
    base_family = font.family) + ggplot2::theme(legend.position = legend.position))
  renderer <- switch(style, classic = plot_style_classic, bold = plot_style_bold,
    framed = plot_style_framed, grid = plot_style_grid, gray = plot_style_gray)
  renderer(font.family, font.size, legend.position)
}
