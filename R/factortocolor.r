#' @title change a factor levels with color
#' @description change a factor levels with color
#' @param factor factor
#' @details use an intern define color palette
#' @return character of color value
#' @examples
#' # not run
#' # factortocolor(factor)
#' @author Florent Dumont <florent.dumont@univresite-paris-saclay.fr>
#' @importFrom magrittr %>%
#' @importFrom foreach %do%
#' @noRd
factortocolor <- function( factor = NULL )
{
  i=1
  c("red", "blue", "green","orange","violet","gray","yellow","cyan","magenta","green4","black","lightslategrey",
    "lightseagreen", "orangered4", "purple2", "rosybrown1","slateblue4" , "yellowgreen","darkred","steelblue",
    "palegreen1","lightcoral","lightgreen","thistle1","navyblue","olivedrab","wheat4","honeydew3",
    "gray28","chartreuse2","orchid3","grey32","plum1","darkorchid1","lightsalmon4","darkolivegreen3","chocolate3",
    "orangered3", "grey11","cyan1", "gray88","seagreen2","darkorange1","darkcyan", "gray95","chocolate1",
    "grey54","maroon4","sandybrown", "oldlace","mediumslateblue","wheat1","lightblue1", "grey30","snow4",
    "beige","red2","royalblue","mediumpurple2","palegreen4","gray67","chartreuse3","grey47","peachpuff4",
    "magenta2", "greenyellow","salmon","brown1","mediumorchid1","grey33","deepskyblue2","salmon2","cornsilk3",
    "slategray4","grey77","grey92","grey37","mistyrose1") -> Palette0
  factor %>% levels -> levels
  factor %>% as.character -> factorcol
  foreach(i=1:length(levels)) %do% { replace(factorcol,factor == levels[i],Palette0[i]) -> factorcol }
  factorcol %>% return()
}