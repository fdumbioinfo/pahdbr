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
#' @param mings numeric minimal size of a gene set
#' @param maxgs numeric maximal size of a gene set
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
enabench <- function(
    symbollist = NULL, foldchange = NULL, keywords = NULL, species = "hs",
    background = 25000, overlapmin = 2, top = 50, labsize = 11, dpibarplot = "screen",
    dolistnetwork = TRUE, dogenesetnetwork = FALSE, confidence = 0,
    intmaxdh = 6000, intmax = 10000000, nodesize = 0.39, mings = 5, maxgs = 700,
    dotopgo = TRUE, ipafile = NULL, nodegraph = c(5,10), path = "." , dirname = NULL )
{
  i=1
  ifelse(is.null(dirname),"ena" -> DirName0,paste("ena_",dirname,sep="") -> DirName0)
  path %>% file.path(DirName0) -> Path0
  if(!dir.exists(Path0)){ Path0 %>% dir.create }
  # if(!is.null(keywords)){
  #   Keywords0 <- foreach(i=1:nrow(keywords)) %do% { 
  #     list(
  #       c(keywords$Keywords[i] %>% strsplit(",") %>% unlist),
  #       c(paste(keywords$KeywordID[i],"-",keywords$Name[i] %>% toupper,sep="") ) ) }
  # }else{Keywords0 <- NULL}
  # MSigDB enrichment analysis
  # msigdb(
  #   symbollist=symbollist,foldchange=foldchange,species=species,top=top,
  #   dogenesetnetwork=dogenesetnetwork, confidence=confidence,
  #   intmax=intmax, intmaxdh=intmaxdh, nodesize=nodesize, path=Path0)
  # -----
  # 1 - MSigDB
  # -----
  symbollist %>% moal::ena(species=species,dirname="msigdb",mings=mings,maxgs=maxgs,path = Path0)
  # -----
  # 2 - topGO
  # -----
  if(dotopgo){symbollist %>% topgo(species=species,nodegraph=nodegraph,top=top,annot=F,path=Path0,dirname=NULL)}
  # -----
  # 3 - IPA
  # -----
  if(!is.null(ipafile)){ipafile %>% ipa(top=top,path=Path0,dirname=NULL)}
  #
  # if(dolistnetwork)
  # {
  #   Path0 %>% file.path("Networks") %>% dir.create
  #   Path0 %>% file.path("Networks") -> Path1
  #   foreach(i=1:length(confidence)) %do%
  #     {
  #       paste("conf",confidence[i] %>% as.character,sep="") -> DirName1
  #       network(nodelist=symbollist,species="hs",foldchange=foldchange,confidence=confidence[i],
  #               intmax=intmax,intmaxdh=intmaxdh,nodesize=nodesize,dirname=DirName1,path=Path1)
  #     }
  # }
  # final network MSigDB
  # Path0 %>% list.files -> Listena0 
  # Path0 %>% file.path("Networks") %>% list.files -> ListNetwork0
  # # MSigDB
  # if( (length(ListNetwork0) > 0) & (Listena0 %>% grepl("MSigDB_log.txt",.) %>% any) ){
  #   msigdbfinalnetwork(
  #     resultpath=Path0,symbollist=symbollist,foldchange=foldchange,
  #     top=top,keywords=Keywords0,intmaxdh=intmaxdh)}
  # # MSigDB + IPA + topGO
  # if( (length(ListNetwork0) > 0) & dotopgo & (Listena0 %>% grepl("MSigDB_log.txt",.) %>% any) ){
  #   finalnetworkall(
  #     resultpath=Path0,symbollist=symbollist,foldchange=foldchange,
  #     top=top,keywords=Keywords0,intmaxdh=intmaxdh)}
  # #
  # Path0 %>% list.files(full.names=T) %>% grep(".tsv|.jpeg",.,value=T) -> rm0
  # rm0 %>% file.remove
  # -----
  # 4 - Occurence
  # -----
  Path0 %>% file.path("ena_msigdb","ena") %>% list.files(full.names = T,recursive = T) %>% grep(".tsv$",.,value=T) -> m0
  m0 %>% head
  Path0 %>% file.path("ipa") %>% list.files(full.names = T) %>% grep(".tsv$",.,value=T) %>% grep("InputList",.,value=T,invert = T) -> i0
  i0
  Path0 %>% file.path("topGO") %>% list.files(full.names = T) %>% grep(".tsv$",.,value=T) -> t0
  t0
  c(m0,i0,t0) -> l0
  rep(c("MSigDB","IPA","topGO"),c(m0 %>% length,i0 %>% length,t0 %>% length)) -> Collection0
  l0 %>% basename %>% strsplit("_") %>% lapply("[",1) %>% unlist %>% paste(Collection0,.,sep="_") -> Collection1
  # lpeAll$Collection %>% table %>% names -> Collection0
  # Symbol occurrence for significant enrich canonical pathways 
  f2 <- foreach(j=1:length(Collection1),.combine = "rbind") %do%
    {
      # pval rank
      l0[j] %>% moal::input(.) -> t
      t %>% dplyr::slice(1:top) -> tt
      tt %>% colnames %>% grep("SymbolList",.) %>% dplyr::select(tt,.) %>% unlist -> ttt
      ttt %>% strsplit("\\|") %>% unlist -> t1
      # lpeAll %>% dplyr::filter(.data$Collection == Collection0[j]) -> t0
      # t0$SymbolList %>% strsplit("\\|") %>% unlist -> t1
      if( table(t1) %>% length == 1 ){ data.frame(Symbol=t1,Occ=1) -> t2 }
      if( table(t1) %>% length > 1 ){ t1 %>% table %>% sort(decreasing = T) %>% data.frame %>% setNames(c("Symbol","Occ")) -> t2 }
      t2 %>% data.frame("Collection"=rep(Collection1[j],nrow(t2))) -> pvalr
      # NES rank
      t %>% dplyr::arrange(-.data$NES) %>% dplyr::slice(1:top) -> tt
      t %>% dplyr::slice(1:top) -> tt
      tt %>% head
      tt %>% colnames %>% grep("SymbolList",.) %>% dplyr::select(tt,.) %>% unlist -> ttt
      ttt %>% strsplit("\\|") %>% unlist -> t1
      # lpeAll %>% dplyr::filter(.data$Collection == Collection0[j]) -> t0
      # t0$SymbolList %>% strsplit("\\|") %>% unlist -> t1
      if( table(t1) %>% length == 1 ){ data.frame(Symbol=t1,Occ=1) -> t2 }
      if( table(t1) %>% length > 1 ){ t1 %>% table %>% sort(decreasing = T) %>% data.frame %>% setNames(c("Symbol","Occ")) -> t2 }
      t2 %>% data.frame("Collection"=rep(Collection1[j],nrow(t2))) -> nesr
      #
      nesr %>% rbind(pvalr) %>% group_by(Symbol) %>% summarise(Occ=sum(.data$Occ)) %>%
        data.frame("Collection"=rep(Collection1[j],nrow(t2)))
    }
  f2 %>% dim
    # f2 <- foreach(j=1:length(Collection0),.combine = "rbind") %do%
  #   {
  #     lpeAll %>% dplyr::filter(.data$Collection == Collection0[j]) -> t0
  #     t0$SymbolList %>% strsplit("\\|") %>% unlist -> t1
  #     if( table(t1) %>% length == 1 ){ data.frame(Symbol=t1,Occ=1) -> t2 }
  #     if( table(t1) %>% length > 1 ){ t1 %>% table %>% sort(decreasing = T) %>% data.frame %>% setNames(c("Symbol","Occ")) -> t2 }
  #     t2 %>% data.frame("Collection"=rep(Collection0[j],nrow(t2)))
  #   }
  f2$Symbol %>% moal::annot(.) -> a0
  a0 %>% dplyr::select(-.data$Symbol) %>% data.frame(f2,.) -> f3
  paste("Symbol-Occurences_",nrow(f3),".tsv",sep="") -> FileName0
  Path0 %>% file.path(FileName0) -> FileName1
  f3 %>% dplyr::arrange(-(.data$Occ)) -> f4
  f4 %>% output(FileName1)
  # Occurence Matrix
  f4 %>% head
  symbollist %>% unique %>% data.frame(Symbol=.) -> m0
  m1 <- foreach(j=1:length(Collection0),.combine="cbind") %do%
    {
      f4 %>% dplyr::filter(.data$Collection == Collection1[j]) %>% dplyr::select(.data$Symbol,.data$Occ) %>% unique -> t0
      t0 %>% head
      m0 %>% dplyr::left_join(t0) %>% dplyr::select(.data$Occ)
    }
  m1 %>% head
  m2 <- foreach(j=1:ncol(m1),.combine="cbind") %do% { m1[j] %>% unlist %>% is.na %>% which %>% replace(m1[j] %>% unlist,.,0) }
  m2 %>% head
  m2 %>% dim
  m2 %>% data.frame(m0$Symbol,.) %>% setNames(c("Symbol",Collection1)) -> m3
  m3$Symbol %>% moal::annot(.) -> a0
  m3 %>% data.frame(a0[,c(-1)]) %>% dplyr::arrange(desc(.[[2]])) -> m4
  m4 %>% head
  paste("Symbol-Occurences","_matrix_",nrow(m4),".tsv",sep="") -> FileName0
  # resultpath %>% file.path("Final-network-all",FileName0) -> FileName1
  Path0 %>% file.path(FileName0) -> FileName1
  m4 %>% output(FileName1)
  # Occurence Matrix hierarchical clustering
  m4 %>% dplyr::select(c(2:(length(Collection0)+1)))  -> m5
  paste("Symbol-Occurences","_matrix_HC_",nrow(m5),".pdf",sep="") -> FileName0
  Path0 %>% file.path(FileName0) -> FileName1
  if(ncol(m5)>2)
  {
    if(m5 %>% apply(2,sd) %>% "=="(0) %>% which %>% any )
    {
      m5 %>% apply(2,sd) %>% "=="(0) %>% which -> selsd
      m5[,-selsd] -> m5
    }
    pdf(FileName1)
    par(mar=c(6.1,4.1,4.1,2.1))
    m5 %>% moal:::hc(m5 %>% colnames %>% as.factor ,legendtitle="Collection",legend=F,cexlabel=0.7)
    graphics.off()
  } 
}