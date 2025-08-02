#' @title ipabarplotup
#' @description parse and create plot from ipa all enrichment results
#' @param dat data.frame enrichment analysis results with 3 column
#' @param top numeric top feature to display
#' @param labsize numeric feature size 
#' @return plot
#' @examples
#' # not run
#' # enabarplot( dat )
#' @author Florent Dumont <florent.dumont@univresite-paris-saclay.fr>
#' @import ggplot2
#' @importFrom gridExtra grid.arrange
#' @noRd
ipabarplotup <- function( dat , top = 80 , labsize = 5 )
{
  if( nrow(dat) < top ){ top <- nrow(dat) }
  dat %>% dplyr::slice( 1:top ) %>% dplyr::select( c(1,2) ) -> Dat0
  # pval plot
  Dat0 %>% ggplot( aes( x= reorder(.data$Name,.data$log10pvalFDR), y = .data$log10pvalFDR ) ) -> p
  p + ggtitle( "Upstream regulators (log10pval)" ) -> p
  p + geom_bar( stat = "identity" , position="dodge2" , fill = "#E69F00" ) -> p
  p + coord_flip() -> p
  p + theme_minimal() -> p
  p + theme( 
    plot.title = element_text(angle=0, size=8, face="bold", vjust=1 , hjust = 1),
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
  Dat0 %>% ggplot( aes( x= reorder(.data$Name,.data$OverlapSize), y = .data$OverlapSize ) ) -> p
  p + ggtitle( "Upstream regulators (overlap)" ) -> p
  p + geom_bar( stat = "identity" , position="dodge2" , fill = "#E69F00" ) -> p
  p + coord_flip() -> p
  p + theme_minimal() -> p
  p + theme( 
    plot.title = element_text(angle=0, size=8, face="bold", vjust=1 , hjust = 1 ),
    axis.text = element_text(angle=0, size=labsize, face="bold", hjust=1.10),
    axis.title = element_blank(),
    legend.key = element_blank(),
    plot.background = element_blank(),
    axis.ticks = element_blank(),
    legend.background = element_blank(),
    strip.background = element_rect(fill="red"),
    panel.background = element_blank() ) -> p2
  gridExtra::grid.arrange( p1, p2, nrow = 1 ) -> pp
  pp
}
