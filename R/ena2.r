#' @title MSigDB ORA analysis for two list.
#' @description MSigDB ORA analysis for two list
#' @param list1 character Symbol list 1
#' @param list2 character Symbol list 2
#' @param list2 character additional Symbol list 2
#' @param species character hs mm rn dr ss
#' @param dirname character output directory name
#' @return directory
#' @examples
#' # not run
#' # ena2( list1 , list2 )
#' @author Florent Dumont <florent.dumont@universite-paris-saclay.fr>
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate arrange select
#' @importFrom rlang .data
#' @importFrom stats fisher.test setNames p.adjust
#' @importFrom foreach foreach %do% %:%
#' @importFrom utils capture.output
#' @noRd
ena2 <- function(
    list1 = NULL, list2 = NULL, addlist = NULL, listnames = NULL, species = "hs",
    dirname = ".", path = NULL)
{
  i=j=1
  ifelse( is.null(dirname), "ena2" -> DirName0, paste("ena2_",dirname,sep="") -> DirName0 )
  path %>% file.path(DirName0) -> Path0
  if(!dir.exists(Path0)){ Path0 %>% dir.create }
  
}
