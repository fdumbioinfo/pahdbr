#' @title ena barplot
#' @param dat data.frame enrichment analysis results with 3 column
#' @param top numeric top feature to display
#' @param labsize numeric feature size 
#' @param title character 
#' @return plot
#' @examples
#' # not run
#' # enabarplot( dat )
#' @author Florent Dumont <florent.dumont@univresite-paris-saclay.fr>
#' @import ggplot2
#' @importFrom gridExtra grid.arrange
#' @importFrom stats reorder
#' @noRd
enabarplot <- function( dat , top = 80 , labsize = 7 , title = "geneset enrichment" )
{
  if( nrow(dat) < top ){ top <- nrow(dat) }
  dat %>% dplyr::slice( 1:top ) %>% dplyr::select( c(1,2) ) -> Dat0
  Dat0$Name %>% substring(1,50) -> Dat0$Name
  # pval plot
  Dat0 %>% ggplot( aes( x= reorder( .[[1]] ,.[[2]]), y = .[[2]] ) ) -> p
  p + ggtitle( paste(title, "(log10pval)" , sep = "" ) ) -> p
  p + geom_bar( stat = "identity" , position="dodge2" , fill = "#E69F00" ) -> p
  p + coord_flip() -> p
  p + theme_minimal() -> p
  p + theme( 
    plot.title = element_text(angle=0, size=15, face="bold", vjust=1, hjust=1),
    axis.text = element_text(angle=0, size=labsize, face="bold", hjust=1.10),
    axis.title = element_blank(),
    legend.key = element_blank(),
    plot.background = element_blank(),
    axis.ticks = element_blank(),
    legend.background = element_blank(),
    strip.background = element_rect(fill="red"),
    panel.background = element_blank() ) -> p1
  # ena plot
  dat %>% dplyr::slice( 1:top ) %>% dplyr::select( c(1,3) ) -> Dat0
  Dat0$Name %>% substring(1,50) -> Dat0$Name
  Dat0 %>% ggplot( aes( x= reorder(.[[1]], .[[2]] ), y = .[[2]] ) ) -> p
  p + ggtitle( paste(title, "(Enrichment)" , sep = "" ) ) -> p
  p + geom_bar( stat = "identity" , position="dodge2" , fill = "#E69F00" ) -> p
  p + coord_flip() -> p
  p + theme_minimal() -> p
  p + theme( 
    plot.title = element_text(angle=0, size=12, face="bold", vjust=1 , hjust = 1 ),
    axis.text = element_text(angle=0, size=labsize, face="bold", hjust=1.10),
    axis.title = element_blank(),
    legend.key = element_blank(),
    plot.background = element_blank(),
    axis.ticks = element_blank(),
    legend.background = element_blank(),
    strip.background = element_rect(fill="red"),
    panel.background = element_blank() ) -> p2
  grid.arrange( p1, p2, nrow = 1 ) -> pp
  pp
}