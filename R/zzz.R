.onLoad <- function(libname, pkgname) {
  ggplot2::theme_set(ggplot2::theme_light() %+replace% 
                     theme(panel.background = element_blank(), 
                           panel.grid.minor = element_blank(),
                           strip.text = element_text(size = 6, colour = "white", face = 'bold')))
  options(dplyr.width = 200)
  options(width=200)
}