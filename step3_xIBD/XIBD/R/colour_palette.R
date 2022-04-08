#' Internal Function
#'
#' XIBD Colour Palette
#'
#' Generate a spectrum colour palette with a specified number of colours.
#' @param number.colours numeric. The number of colours to return.
#' @return A character vector of length=\code{number.colours} containing a colour specturm.
#' @importFrom grDevices colorRampPalette
getColourPaletteMajor <- function(number.colours) {
  if (number.colours == 1 )
    return("#69B4FF")
  if (number.colours == 2)
    return(c("#69B4FF","#FF7575"))
  if (number.colours == 3)
    return(c("#69B4FF","#FFE455","#FF7575"))
  if (number.colours == 4)
    return(c("#69B4FF","#99DD55","#FFE455","#FF7575"))
  if (number.colours == 5)
    return(c("#F3ABF3","#69B4FF","#99DD55","#FFE455","#FF7575"))
  if (number.colours == 6)
    return(c("#F3ABF3","#69B4FF","#55DDDD","#99DD55","#FFE455","#FF7575"))
  if (number.colours >= 7) {
    my.palette <- colorRampPalette(c("#F3ABF3","#69B4FF","#55DDDD","#99DD55","#FFE455","#FFB255","#FF7575"))
    return(my.palette(number.colours))
  }
}

#' Generate a specified number of shades of a given colour.
#' @param major.colour character. The colour name or code.
#' @param number.colours numeric. The number of colours to return.
#' @return A character vector of length=\code{number.colours} containing colour shades, excluding white.
#' @importFrom grDevices colorRampPalette
getColourPaletteMinor <- function(major.colour, number.colours) {
  my.palette.1 <- colorRampPalette(c(major.colour, "white"))
  my.palette.2 <- my.palette.1(number.colours+1)
  return(my.palette.2[1:number.colours])
}
