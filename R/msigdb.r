#' @title MSigDB enrichment analysis
#' @description MSigDB enrichment analysis
#' @param symbollist character Symbol or NCBI gene ID
#' @param geneannot data.frame for omic function only
#' @param species character hs mm rn dr ss
#' @param background numeric
#' @param overlapmin numeric for minimum overlap between geneset and list
#' @param enascoremin numeric for minimum ratio ena
#' @param top numeric top features to plot
#' @param labsize numeric size of function in barplot
#' @param dpibarplot character barplot resolution
#' @param foldchange data.frame Symbol list with fold-change values to display
#' @param top numeric top enrichment score geneset selection
#' @param dogenesetnetwork logical if TRUE (by default) do top geneset interaction network
#' @param confidence numeric for minus interaction confidence score
#' @param dopar logical do parallelizatino for geneset networks
#' @param path character for relative path of output directory
#' @return file with enrichment analysis results
#' @examples
#' # not run
#' # msigdb( Symbollist , filtergeneset = "reactome")
#' @author Florent Dumont <florent.dumont@universite-paris-saclay.fr>
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate arrange select
#' @importFrom rlang .data
#' @importFrom stats fisher.test setNames p.adjust
#' @importFrom foreach foreach %do% %:%
#' @importFrom utils capture.output
#' @noRd
msigdb <- function(
    symbollist = NULL, geneannot = NULL, foldchange = NULL, species = "hs", 
    background = 25000 , overlapmin = 1 , enascoremin = 1, dogenesetnetwork = FALSE, confidence = 400,
    intmax = 10000 , intmaxdh = 2500, nodesize = 0.39,
    dopar = TRUE, top = 25, labsize = 11, dpibarplot = "screen", path = ".")
{
  i=j=1
  data.frame() -> a0
  FALSE -> OK1; FALSE -> OK2 ; FALSE -> OK3
  ifelse(species != "hs", Ortholog <- T, Ortholog <- FALSE)
  #
  if( is.null(geneannot) ) {
    symbollist %>% as.character %>% unique -> symbollist
    symbollist %>% is.na %>% "!"(.) %>% which %>% symbollist[.] -> symbollist
    symbollist %>% grep("^$",.,invert=T) %>% symbollist[.] -> symbollist
    # symbollist %>% grep("NA",.,invert=T) %>% symbollist[.] -> symbollist
    symbollist %>% "=="("NA") %>% "!"(.) %>% which %>% symbollist[.] -> symbollist
    if(length(symbollist) > overlapmin ){ TRUE -> OK1 }
    if( OK1 ){ symbollist %>% annot( species=species, ortholog=Ortholog ) -> a0 }
  }else{ 
    if(nrow(geneannot) > overlapmin){ TRUE -> OK1 }
    if( OK1 ){ geneannot -> a0 } 
  }
  #
  if( nrow(a0) > overlapmin ){ TRUE -> OK2 }
  #
  if(OK2){
    a0$GeneID -> GeneidList0
    moalannotgene::genesetdb -> GenesetDb0
    # 
    Ena0 <- foreach( i=1:length(GenesetDb0) ) %:%
      foreach( j=1:length(GenesetDb0[[i]]), .combine="rbind" ) %do%
      {
        GenesetDb0[[i]][[j]][[1]] -> GenesetName
        GeneidList0 %in% GenesetDb0[[i]][[j]][[4]] %>% which %>% a0$GeneID[.] -> OverlapGeneIDList0
        GeneidList0 %in% GenesetDb0[[i]][[j]][[4]] %>% which %>% a0$Symbol[.] -> OverlapSymbolList0
        a0$Symbol %>% length -> ListSize0
        OverlapSymbolList0 %>% length -> OverlapSize0
        GenesetDb0[[i]][[j]][[4]] %>% length -> GenesetSize0
        ( OverlapSize0 / GenesetSize0 ) %>% round(2) -> OverlapRatio0
        ( ( OverlapSize0 / ListSize0 ) / (  GenesetSize0 / background ) ) %>% round(2) -> EnaScore0
        if(  OverlapSize0 >= overlapmin & EnaScore0 > enascoremin )
        {
          c(OverlapSize0,ListSize0,GenesetSize0,background) %>% matrix(ncol=2) -> ContTable0
          ContTable0 %>% fisher.test %>% unlist %>% "["(1) %>% as.numeric -> Pval0
          c(
            GenesetName,
            paste0( OverlapSymbolList0, collapse = "|" ),
            paste0( OverlapGeneIDList0, collapse = "|" ),
            OverlapSize0,GenesetSize0,OverlapRatio0,EnaScore0,Pval0)
        }
      }
    Ena0 %>% setNames( GenesetDb0 %>% names ) -> Ena0
    Ena0 %>% sapply(is.null) %>% "!"(.) %>% which %>% Ena0[.] -> Ena0
    if(length(Ena0) > 0){ TRUE -> OK3 }
  }
  #
  if(OK3)
  {
    foreach(i=1:length(Ena0)) %do%
      {
        Ena0[[i]] %>% data.frame(stringsAsFactors=F) %>%
          stats::setNames(c("Name","SymbolList","GeneIDList","OverlapSize","GenesetSize","OverlapRatio","ENAScore","pval")) -> Ena1
        Ena1$OverlapSize %>% as.numeric -> Ena1$OverlapSize
        Ena1$GenesetSize %>% as.numeric -> Ena1$GenesetSize
        Ena1$OverlapRatio %>% as.numeric -> Ena1$OverlapRatio
        Ena1$ENAScore %>% as.numeric -> Ena1$ENAScore
        Ena1$pval %>% as.numeric -> Ena1$pval
        Ena1 %>% mutate( pvalFDR = .data$pval %>% p.adjust(method = "fdr") ) %>%
          mutate( log10pvalFDR = .data$pvalFDR %>% log10 %>% '*'( . , -1 ) %>% round( . , 4) ) %>%
          arrange( .data$pvalFDR  ) -> Ena1
        # output
        paste(
          Ena0[i] %>% names %>% strsplit("\\|") %>% unlist %>% "[["(1),"_",
          Ena0[i] %>% names %>% strsplit("\\|") %>% unlist %>% "[["(2),"_",
          dim(Ena1)[1],".tsv", sep="" ) -> FileName0
        path %>% file.path(FileName0) %>% output(Ena1,.)
        # barplot
        Ena1 %>% dplyr::select(c(1,10,7)) -> Ena2
        if( nrow(Ena2) > 0 )
        {
          paste( Ena0[i] %>% names %>% strsplit("\\|") %>% unlist %>% "[["(1),
                 Ena0[i] %>% names %>% strsplit("\\|") %>% unlist %>% "[["(2) , sep = " "  ) -> title
          enabarplot(dat=Ena2,top=top,labsize=labsize,title=title) -> p
          # output
          ifelse( nrow(Ena2) < top, nrow(Ena2) -> Topfilename0, top -> Topfilename0 )
          paste(
            "top",top,"_",Ena0[i] %>% names %>% strsplit("\\|") %>% unlist %>% "[["(1),"_",
            Ena0[i] %>% names %>% strsplit("\\|") %>% unlist %>% "[["(2),"_",
            Topfilename0, ".jpeg" , sep = "" ) -> FileName0
          path %>% file.path(FileName0) -> FileName1
          ggsave(filename=FileName1,plot=p,width=12,height=15,dpi=dpibarplot )
        }
    }
    # log
    Enainfo0 <- list()
    path %>% basename %>% gsub("(^ena_)","",.) %>% paste("List name: ",.,sep="") -> Enainfo0[[1]]
    ListSize0 %>% paste("List size: ",.,sep="") -> Enainfo0[[2]]
    GenesetDb0 %>% unlist(recursive=F) %>% length %>% paste("Genesets used: ",.,sep="") -> Enainfo0[[3]]
    background %>% paste("background used: ",.,sep="") -> Enainfo0[[4]]
    a0$Symbol %>% paste0(collapse="|") %>% paste("Symbols : ",.,sep="") -> Enainfo0[[5]]
    path %>% file.path("MSigDB_log.txt") %>% file("w") -> f
    foreach(i=1:length(Enainfo0)) %do% { Enainfo0[[i]] %>% paste(.,"\n",sep="") %>% writeLines(f) }
    close(f)
    # Occ
    msigdbocc(
      resultpath=path, symbollist=symbollist, foldchange=foldchange, top=top, confidence=confidence,
      dogenesetnetwork=dogenesetnetwork, intmax=intmax, intmaxdh=intmaxdh, nodesize=nodesize, dopar=dopar)
  }
}
