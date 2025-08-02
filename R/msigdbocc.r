#' @title MSigDB occurence analysis
#' @description MSigDB occurence analysis
#' @param resultpath character msigdb result directory path
#' @param symbollist character Symbol or NCBI gene ID
#' @param foldchange data.frame Symbol list with fold-change values to display
#' @param top numeric top enrichment score geneset selection
#' @param confidence numeric for minus interaction confidence score
#' @param dogenesetnetwork logical if TRUE (by default) do top geneset interaction network
#' @param dopar logical do parallelizatino for geneset networks
#' @return directory with top enrichment analysis and occurence results 
#' @examples
#' # not run
#' # msigdbocc( symbollist )
#' @author Florent Dumont <florent.dumont@universite-paris-saclay.fr>
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate arrange select desc
#' @importFrom foreach foreach %do%
#' @importFrom grDevices graphics.off pdf
#' @importFrom graphics par
#' @noRd
msigdbocc <- function( 
    resultpath = NULL, symbollist = NULL, foldchange = NULL, top = 50, confidence = 0, dogenesetnetwork = FALSE,
    intmax = 10000, intmaxdh = 2500, nodesize = 0.39, dopar = TRUE )
{
  i=j=k=1
  c("Pathways","Functions","Regulators","Patterns","Localization") -> CollDirName0
  c("c0|c2","c5","c3","c4|c6|c7","c1|c8") %>% paste("(",.,")",".*.tsv$",sep="") -> OccGrep0
  c("c0|c2","c5","c3","c4|c6|c7","c1|c8") %>% paste("(",.,")",".*.jpeg$",sep="") -> OccGrepjpeg0
  foreach(i=1:length(CollDirName0)) %do%
    {
      resultpath %>% file.path(CollDirName0[i]) -> OccPath0
      if(!dir.exists(OccPath0)){dir.create(OccPath0)} 
      c(2,top) -> Threshold
      c("01",top %>% as.character) -> ThresholdName
      resultpath %>% list.files(full.names = T) %>% grep(".tsv$",.,value=T) -> l0
      l0 %>% basename %>% grep(".tsv$",.,value=T) %>% grep(OccGrep0[i],.,value=F) -> selbasename
      # resultpath %>% list.files(full.names = T) %>% basename %>% grep(".tsv$",.,value=T) %>% grep(OccGrep0[i],.,value=F) -> selbasename
      l0[selbasename] -> l1
      # resultpath %>% list.files(full.names = T) %>% grep(".tsv$",.,value=T) %>% grep(OccGrep0[i],.,value=T) -> l0
      l1 %>% basename %>% strsplit("_") %>% lapply("[",3) %>% unlist -> Collection0
      #
      f0 <- foreach(j=1:length(l1),.combine="rbind") %do%
        {
          l1[j] %>% moal::input(.) -> t0
          t0 %>% dplyr::select(.data$Name,.data$SymbolList,.data$OverlapSize,.data$GenesetSize,.data$ENAScore,.data$log10pvalFDR) -> t1
          # top pval
          t1 %>% dplyr::filter( .data$log10pvalFDR > Threshold[1] ) %>%
            dplyr::arrange(desc(.data$log10pvalFDR)) %>% dplyr::slice(1:Threshold[2] ) -> toppval
          # top ENA
          t1 %>% dplyr::arrange( desc(.data$ENAScore) ) %>% dplyr::slice(1:Threshold[2] )-> topENA
          toppval %>% rbind(topENA) %>% unique -> t2
          rep(Collection0[j],nrow(t2)) -> Collection1
          t2 %>% data.frame("Collection"=Collection1) -> t3
          l1[j] %>% basename %>% sub("(.*_).*.tsv","\\1",.) -> top50FileName0
          paste("top",top,"_",top50FileName0,nrow(t3),".tsv",sep="") -> top50FileName1
          file.path(OccPath0,top50FileName1) -> top50FileName2
          t3 %>% moal::output(top50FileName2)
          t3
        }
      f0 %>% dplyr::arrange(desc(.data$log10pvalFDR)) -> f1
      paste("top",top,"_",CollDirName0[i],"_",nrow(f1),".tsv",sep="") -> FileName0
      file.path(OccPath0,FileName0) -> FileName1
      f1 %>% moal::output(FileName1)
      # Symbol occurrence for significant enrich canonical pathways 
      f2 <- foreach(j=1:length(Collection0),.combine = "rbind") %do%
        {
          f1 %>% dplyr::filter(.data$Collection == Collection0[j]) -> t0
          t0$SymbolList %>% strsplit("\\|") %>% unlist -> t1
          if( length(t1 %>% table) == 1 ){ data.frame(Symbol=t1 %>% table %>% names, Occ= t1 %>% table %>% as.numeric) -> t2 }else{
            t1 %>% table %>% sort(decreasing = T) %>% data.frame %>% setNames(c("Symbol","Occ")) -> t2  }
          t2 %>% data.frame("Collection"=rep(Collection0[j],nrow(t2)))
        }
      f2$Symbol %>% annot(.) -> a0
      a0 %>% dplyr::select(-.data$Symbol) %>% data.frame(f2,.) -> f3
      paste("top",top,"_","Symbol-Occurences_",CollDirName0[i],"_",nrow(f3),".tsv",sep="") -> FileName0
      file.path(OccPath0,FileName0) -> FileName1
      f3 %>% dplyr::arrange(-(.data$Occ)) -> f4
      f4 %>% moal::output(FileName1)
      # Occurence Matrix
      symbollist %>% unique %>% data.frame(Symbol=.) -> m0
      m1 <- foreach(j=1:length(Collection0),.combine="cbind") %do%
        {
          f4 %>% dplyr::filter(.data$Collection == Collection0[j]) %>% dplyr::select(.data$Symbol,.data$Occ) -> t0
          m0 %>% dplyr::left_join(t0) %>% dplyr::select(.data$Occ)
        }
      m2 <- foreach(j=1:ncol(m1),.combine="cbind") %do% { m1[j] %>% unlist %>% is.na %>% which %>% replace(m1[j] %>% unlist,.,0) }
      m2 %>% data.frame(m0$Symbol,.) %>% setNames(c("Symbol",Collection0)) -> m3
      m3$Symbol %>% annot(.) -> a0
      m3 %>% data.frame(a0[,c(-1)]) %>% dplyr::arrange(desc(.[[2]])) -> m4
      paste("top",top,"_","Symbol-Occurences_",CollDirName0[i],"_matrix_",nrow(m4),".tsv",sep="") -> FileName0
      file.path(OccPath0,FileName0) -> FileName1
      m4 %>% moal::output(FileName1)
      # Occurence Matrix hierarchical clustering
      m4 %>% dplyr::select(c(2:(length(Collection0)+1)))  -> m5
      paste("top",top,"_","Symbol-Occurences_",CollDirName0[i],"_matrix_HC_",nrow(m5),".pdf",sep="") -> FileName0
      file.path(OccPath0,FileName0) -> FileName1
      if(ncol(m5)>2)
      {
        pdf(FileName1)
        par(mar=c(6.1,4.1,4.1,2.1))
        m5 %>% moal:::hc(m5 %>% colnames %>% as.factor ,legendtitle="Collection",legend=F,cexlabel=0.9)
        graphics.off()
      }
      #
      # networks
      #
      if(dogenesetnetwork)
      {
        OccPath0 %>% file.path("networks") -> OccPath1
        if(!dir.exists(OccPath1)){dir.create(OccPath1)} 
        f1$SymbolList %>% strsplit("\\|") %>% lapply(length) %>% ">"(2) %>% which %>% f1[.,] -> f2
        if(dopar){ parallel::makeCluster(parallel::detectCores()) -> cl;registerDoParallel(cl);doParallel::registerDoParallel(cl) }
        foreach(k=1:nrow(f2),.packages=c("magrittr","dplyr","moal","foreach","stringr","igraph","Rgraphviz","moalstringdbhs")) %dopar%
          {
            f2[k,2] %>% strsplit("\\|") %>% unlist -> SymbolList0
            f2[k,1] %>% substr(1,20) %>%
              # gsub("\\/"," ",.) %>% gsub("\\.","",.) %>% gsub("-","",.) %>%
              paste("_",f2[k,3],"_conf",confidence[1] %>% as.character,sep="") %>% gsub("\\/"," ",.) -> DirNameGeneset0
            f2[k,1] %>% paste("_",f2[k,3],"_conf",confidence[1] %>% as.character,sep="") %>% gsub("\\/"," ",.) -> Title0
            network(nodelist=SymbolList0,foldchange=foldchange,species="hs",
                    confidence=confidence[1],title=Title0,path=OccPath1,dirname=DirNameGeneset0)
           
          }
      if(dopar){ parallel::stopCluster(cl) }
      }
      # move file
      l1 %>% file.copy(OccPath0)
      # l0 %>% file.remove
      resultpath %>% list.files(full.names = T) %>% grep(".jpeg$",.,value=T) -> l0
      l0
      l0 %>% basename %>% grep(OccGrepjpeg0[i],.,value=F) -> selbasename
      # l0 %>% basename %>% grep(".jpeg$",.,value=T) %>% grep(OccGrep0[i],.,value=F) -> selbasename
      # resultpath %>% list.files(full.names = T) %>% basename %>% grep(".tsv$",.,value=T) %>% grep(OccGrep0[i],.,value=F) -> selbasename
      l0[selbasename] -> l1
      resultpath %>% list.files(full.names=T) %>% grep(OccGrep0[i] %>% gsub("tsv","jpeg",.),.,value=T) -> l0
      l1 %>% file.copy(OccPath0)
      # l0 %>% file.remove
    }
}
