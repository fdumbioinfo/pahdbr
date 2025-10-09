#' @title topGO analysis
#' @description do topgo analysis
#' @param symbollist character list of Symbol
#' @param species character
#' @param overlapmin numeric min overlap between Symbol and GO term
#' @param nodegraph numeric number of top GO term to use for GO term nodes graph
#' @param topena numeric number of GOTerm to display on barplot
#' @param annot logical if TRUE (by default) annotate list before analysis
#' @param path character file path
#' @param dirname character results directory name
#' @details topGO analysis using classic fisher method
#' @return directory with barplot and graph node of GO term
#' @examples
#' # not run
#' # topgo( symbollist )
#' @author Florent Dumont <florent.dumont@univresite-paris-saclay.fr>
#' @importFrom magrittr %>%
#' @importFrom dplyr select
#' @importFrom ggplot2 ggsave
#' @importFrom topGO annFUN.org printGraph GenTable genesInTerm runTest annFUN.GO2genes usedGO
#' @importFrom methods new
#' @importFrom stringr str_to_title str_to_upper
#' @import AnnotationDbi
#' @import GO.db
#' @import org.Hs.eg.db
#' @import moalannotgene
#' @export
topgo <- function(
    symbollist = NULL, species = "hs", overlapmin = 1,
    nodegraph = c(5,10), topena = 50, annot = FALSE,
    path = ".", dirname = NULL )
{
  i=j=k=1
  symbollist %>% as.character -> SymbolList0
  moal:::orthoinfo -> orthoinfo 
  GeneDb0 <- NULL
  if(annot){
    SymbolList0 %>% annot(species=species) %>% dplyr::select(.data$Symbol) %>% unlist %>% as.character -> SymbolList1 }else{
      SymbolList0 -> SymbolList1 }
  # background
  ifelse( is.null(species),orthoinfo[[1]] -> Species0,
          orthoinfo %>% sapply("[[",1) %>% grep(species,.) %>% "[["(orthoinfo,.) -> Species0 )
  paste("moalannotgene::",paste0("genedb", Species0[1])," -> GeneDb0",sep="") -> text
  eval(expr = parse(text = text ) )
  GeneDb0 %>% dplyr::filter(.data$Species == Species0[1] & .data$GeneType == "protein-coding") -> GeneDb1
  GeneDb1 %>% dplyr::select(.data$Symbol) %>% unlist %>% unique -> Background0
  # output
  paste("topGO",sep="") -> DirName
  if(!is.null(dirname)){ paste(DirName,"_",dirname,sep="") -> DirName }
  path %>% file.path( DirName ) -> Path
  Path %>% dir.create
  Background0 %in% SymbolList1 %>% which -> sel
  if(length(sel) > overlapmin)
  {
    # background
    rep(c("0"),length(Background0)) %>% setNames(Background0) -> Background1
    replace(Background1,sel,c("1")) %>% as.factor -> Background2
    # GO annotation
    c("BP","MF","CC") -> db0
    c("FUNCTIONS","MOLECULES","LOCALIZATION") -> dbDirName0
    #
    foreach(i=1:length(db0)) %do% 
    {
      # Path %>% file.path(dbDirName0[i]) %>% dir.create
      if(Species0[1] == "hs"){ topGO::annFUN.org(whichOnto=db0[i],feasibleGenes=NULL,mapping="org.Hs.eg.db",ID="symbol") -> allGO2genes }
      # run topgo
      GOdata <- new("topGOdata",ontology=db0[i],allGenes=Background2,annot=annFUN.GO2genes,GO2genes=allGO2genes,nodeSize=overlapmin)
      results <- runTest(GOdata, algorithm="classic", statistic="fisher")
      GenTable(GOdata,classicFisher=results,topNodes=length(usedGO(GOdata)),numChar=1000) -> goEnrichment0
      goEnrichment0$Term -> Name
      goEnrichment0$GO.ID -> GOID
      goEnrichment0$classicFisher %>% "=="("< 1e-30") %>% replace(goEnrichment0$classicFisher,.,1e-30) %>% as.numeric -> pval
      pval %>% p.adjust(method = "fdr") -> pvalFDR
      pvalFDR %>% log10(.) %>% "*"(-1) -> log10pvalFDR
      goEnrichment0$Annotated -> GenesetSize
      goEnrichment0$Significant -> OverlapSize
      OverlapSize %>% "/"(goEnrichment0$Expected) -> NES
      data.frame(Name,GOID,OverlapSize,GenesetSize,NES,pval,pvalFDR,log10pvalFDR) -> goEnrichment1
      goEnrichment1 %>% dim
      goEnrichment1 %>% dplyr::filter(OverlapSize > 1 & NES > 1) -> goEnrichment2
      goEnrichment2 %>% dim
      goEnrichment2 %>% colnames
      GOdata %>% genesInTerm -> gt0
      gt0 %>% names %in% goEnrichment2$GOID %>% which -> sel
      if(length(sel) > nodegraph[length(nodegraph)])
      {
        gt0[sel] -> goLists1
        goLists1 %>% names -> GOID
        goLists1 %>% lapply(paste0,collapse="|") %>% unlist -> SymbolList
        goEnrichment2 %>% dplyr::inner_join(data.frame(GOID,SymbolList,stringsAsFactors=FALSE)) -> goEnrichment3
        goEnrichment3 %>% colnames
        goEnrichment3 %>% dplyr::select(c(1,2,ncol(goEnrichment3),3:(ncol(goEnrichment3)-1))) -> goEnrichment3
        goEnrichment3 %>% colnames
        # goEnrichment3 %>% dplyr::select(.data$Name,.data$GOID,.data$SymbolList,.data$OverlapSize,.data$GenesetSize,.data$NES,.data$log10pvalFDR) -> goEnrichment3
        #
        foreach(k=1:nrow(goEnrichment3)) %do%
          {
            goEnrichment3$SymbolList[k] %>% strsplit("\\|") %>% unlist -> t
            t %in% symbollist %>% which -> sel
            if(length(sel)>0){ t[sel] %>% paste0(collapse = "|") -> goEnrichment3$SymbolList[k] }
          }
        # output
        paste(dbDirName0[i],"_",nrow(goEnrichment3),".tsv",sep = "") -> FileName0
        Path %>% file.path(FileName0) -> FileName1
        goEnrichment3 %>% output(FileName1)
        # top
        # goEnrichment3 %>% dplyr::slice(1:top) -> toppval
        # goEnrichment3 %>% dplyr::arrange(-.data$ENAScore) %>% dplyr::slice(1:top) -> topena
        # rbind(toppval,topena) %>% unique -> goEnrichment4
        # output
        # paste("top",top,"_",dbDirName0[i],"_",nrow(goEnrichment4),".tsv",sep = "") -> FileName0
        # Path %>% file.path(dbDirName0[i], FileName0) -> FileName1
        # goEnrichment4 %>% output(FileName1)
        # barplot
        # pval ranking
        goEnrichment3 %>% colnames
        goEnrichment3 %>% dplyr::slice(1:topena) -> fgseapval1plot
        fgseapval1plot %>% colnames
        fgseapval1plot$pval %>% log10 %>% "*"(-1) -> Log10Pval
        fgseapval1plot %>% data.frame(Log10Pval) -> fgseapval1plot2
        # reduce long geneset name
        fgseapval1plot2$Name %>% gsub(" [a-z]*\\.\\.\\.$","",.) %>% gsub("\\.\\.\\.$","",.) -> fgseapval1plot2$Name
        fgseapval1plot2$Name %>% as.character %>% nchar -> Nchar0
        Nchar0 %>% ">"(50) %>% which -> selNchar
        if(length(selNchar)>0)
        {
          fgseapval1plot2$Name[selNchar] %>% substr(1,45) -> Head0
          fgseapval1plot2$Name[selNchar] %>% substr(Nchar0-10,Nchar0)-> Tail0
          fgseapval1plot2$Name %>% as.character -> NameNchar0
          fgseapval1plot2$Name[selNchar] <- paste(Head0,Tail0,sep="...")
        }
        fgseapval1plot2 %>% dplyr::mutate(Name=forcats::fct_reorder(.data$Name,Log10Pval)) -> fgseapval1plot3
        fgseapval1plot3$OverlapSize %>% paste("/",fgseapval1plot3$GenesetSize) %>% 
          lapply(paste0,collapse="") %>% unlist -> Overlap
        fgseapval1plot3 %>% data.frame(Overlap) -> fgseapval1plot3
        fgseapval1plot3 -> datp
        # NES ranking
        goEnrichment3 %>% colnames
        goEnrichment3 %>% dplyr::arrange( -abs(.data$NES) ) %>% dplyr::slice(1:topena) -> fgseapval1plot
        # rl5 %>% dplyr::slice(1:topena) -> fgseapval1plot
        fgseapval1plot %>% colnames
        fgseapval1plot$pval %>% log10 %>% "*"(-1) -> Log10Pval
        fgseapval1plot %>% data.frame(Log10Pval) -> fgseapval1plot2
        # reduce long geneset name
        fgseapval1plot2$Name %>% as.character %>% nchar -> Nchar0
        Nchar0 %>% ">"(50) %>% which -> selNchar
        if(length(selNchar)>0)
        {
          fgseapval1plot2$Name[selNchar] %>% substr(1,45) -> Head0
          fgseapval1plot2$Name[selNchar] %>% substr(Nchar0-10,Nchar0)-> Tail0
          fgseapval1plot2$Name %>% as.character -> NameNchar0
          fgseapval1plot2$Name[selNchar] <- paste(Head0,Tail0,sep="...")
        }
        fgseapval1plot2 %>% dplyr::mutate(Name=forcats::fct_reorder(.data$Name,Log10Pval)) -> fgseapval1plot3
        fgseapval1plot3$OverlapSize %>% paste("/",fgseapval1plot3$GenesetSize) %>% 
          lapply(paste0,collapse="") %>% unlist -> Overlap
        fgseapval1plot3 %>% data.frame(Overlap) -> fgseapval1plot3
        fgseapval1plot3 -> datnes
        datnes %>% head
        
        datp %>% enabarplot(datnes=datnes,title=dbDirName0[i]) -> p
        
        
        
        # fgseapval1plot3 %>% enabarplot(title=dbDirName0[i]) -> p
        paste(dbDirName0[i],"_",nrow(fgseapval1plot3),".pdf",sep = "") -> FileName0
        Path %>% file.path(FileName0) -> FileName1
        ggsave(plot=p,filename=FileName1,width=40,height=30,scale=1,units="cm")
        
        # goEnrichment3$Name %>% gsub(" [a-z]*\\.\\.\\.$","",.) %>% gsub("\\.\\.\\.$","",.) -> GOTerm0
        # paste(goEnrichment3$GOID,GOTerm0,sep=", ") %>% substr(1,50) -> GOTerm1
        # GOTerm1 %>% factor(levels=rev(GOTerm1)) -> GOTerm2
        # goEnrichment3$log10pvalFDR -> Values0
        # paste(dbDirName0[i],"_Barplot_top",top,"_GOTerms", sep = "") -> FileName0
        # if(length(GOTerm2) < top){ top <- length(GOTerm2) }
        # data.frame(GOTerm2,Values0) %>% dplyr::slice(1:top) %>% moal:::barplot(title=paste(dbDirName0[i],"_Barplot top ",top," GO Terms",sep="")) -> p
        # paste(FileName0,".pdf", sep = "") -> FileName1
        # Path %>% file.path(dbDirName0[i],FileName1) -> FileName2
        # ggsave(FileName2,p)
        # graph node
        foreach( j=nodegraph ) %do%
          {
            paste(dbDirName0[i],"_NodeGraph",sep="") -> FileName0
            Path %>% file.path(FileName0) -> FileName1
            printGraph(GOdata, results, firstSigNodes=j, fn.prefix=FileName1, useInfo="all", pdfSW=TRUE)
          }
      }
    }
  }
}