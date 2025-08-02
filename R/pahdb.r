#' @title MSigDB enrichment analysis
#' @description MSigDB enrichment analysis
#' @param symbollist character Symbol or NCBI gene ID
#' @param foldchange data.frame Symbol list with fold-change values to display
#' @param species character hs mm rn dr ss
#' @param keywords data.frame keywords list for geneset selection 
#' @param background numeric
#' @param overlapmin numeric for minimum overlap between geneset and list
#' @param top numeric top features to plot
#' @param labsize numeric size of function in barplot
#' @param dpibarplot character barplot resolution
#' @param dolistnetwork logical create symbollist interaction network 
#' @param dogenesetnetwork logical create top geneset interaction network
#' @param intmaxdh numeric maximum number of interaction to use for Davidson and Harel algorithm layout
#' @param confidence numeric MSigDB interaction confidence score
#' @param intmax numeric maximum number of interaction to use
#' @param nodesize numeric change Symbol size
#' @param dotopgo logical do gene ontology enrichment analysis using topGO
#' @param ipafile character Ingenuity Pathways Analysis export all results file
#' @param nodegraph numeric number of top GO term to use for GO term nodes graph
#' @param path character for relative path of output directory
#' @param dirname character name for output
#' @return file with enrichment analysis results
#' @examples
#' # not run
#' # ena124( Symbollist , filtergeneset = "reactome")
#' @author Florent Dumont <florent.dumont@universite-paris-saclay.fr>
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate arrange select desc
#' @importFrom rlang .data
#' @importFrom stats fisher.test setNames p.adjust
#' @importFrom foreach foreach %do% %:%
#' @importFrom utils capture.output
#' @export
pahdb <- function(
    symbollist = NULL, foldchange = NULL, keywords = NULL, species = "hs",
    background = 25000, overlapmin = 2, top = 50, labsize = 11, dpibarplot = "screen",
    dolistnetwork = TRUE, dogenesetnetwork = FALSE, confidence = c(0,150,400,500,600,700,800,900,seq(950,990,10)),
    intmaxdh = 5000, intmax = 100000, nodesize = 0.39,
    dotopgo = FALSE, ipafile = NULL, nodegraph = c(5,10), path = "." , dirname = NULL )
{
  i=1
  ifelse( is.null(dirname), "ena" -> DirName0, paste("ena_",dirname,sep="") -> DirName0 )
  path %>% file.path(DirName0) -> Path0
  if(!dir.exists(Path0)){ Path0 %>% dir.create }
  if(!is.null(keywords)){
    Keywords0 <- foreach(i=1:nrow(keywords)) %do% { 
      list(
        c(keywords$Keywords[i] %>% strsplit(",") %>% unlist),
        c(paste(keywords$KeywordID[i],"-",keywords$Name[i] %>% toupper,sep="") ) ) }
  }else{Keywords0 <- NULL}
  # MSigDB enrichment analysis
  msigdb(
    symbollist=symbollist,foldchange=foldchange,species=species,top=top,
    dogenesetnetwork=dogenesetnetwork, confidence=confidence,
    intmax=intmax, intmaxdh=intmaxdh, nodesize=nodesize, path=Path0)
  # topgo
  if(dotopgo){ topgo(symbollist=symbollist,species=species,nodegraph=nodegraph,top=top,annot=F,path=Path0,dirname=NULL) }
  # ipa
  if(!is.null(ipafile)){ ipa(ipafile=ipafile,top=top,path=Path0,dirname=NULL) }
  # list network
  if(!is.null(keywords)){ TRUE -> dolistnetwork }
  #
  if(dolistnetwork)
  {
    Path0 %>% file.path("Networks") %>% dir.create
    Path0 %>% file.path("Networks") -> Path1
    foreach(i=1:length(confidence)) %do%
      {
        paste("conf",confidence[i] %>% as.character,sep="") -> DirName1
        network(nodelist=symbollist,species="hs",foldchange=foldchange,confidence=confidence[i],
                intmax=intmax,intmaxdh=intmaxdh,nodesize=nodesize,dirname=DirName1,path=Path1)
      }
  }
  # final network MSigDB
  Path0 %>% list.files -> Listena0 
  Path0 %>% file.path("Networks") %>% list.files -> ListNetwork0
  # MSigDB
  if( (length(ListNetwork0) > 0) & (Listena0 %>% grepl("MSigDB_log.txt",.) %>% any) ){
    msigdbfinalnetwork(
      resultpath=Path0,symbollist=symbollist,foldchange=foldchange,
      top=top,keywords=Keywords0,intmaxdh=intmaxdh)}
  # MSigDB + IPA + topGO
  if( (length(ListNetwork0) > 0) & dotopgo & (Listena0 %>% grepl("MSigDB_log.txt",.) %>% any) ){
    finalnetworkall(
      resultpath=Path0,symbollist=symbollist,foldchange=foldchange,
      top=top,keywords=Keywords0,intmaxdh=intmaxdh)}
  #
  Path0 %>% list.files(full.names=T) %>% grep(".tsv|.jpeg",.,value=T) -> rm0
  rm0 %>% file.remove
}