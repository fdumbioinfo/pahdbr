#' @title ipa
#' @description parse and create plot from ipa all enrichment results
#' @param ipafile character path to the file
#' @param top numeric top feature to display
#' @param labsize numeric feature size 
#' @param background numeric
#' @param overlapmin numeric for minimum overlap between geneset and list
#' @param enascoremin numeric for minimum ratio ena
#' @param dpi character jpeg resolution between retina print screen 
#' @param path character for relative path of output directory
#' @param dirname character name for output
#' @return no values
#' @examples
#' # not run
#' # ipa( file , path "" )
#' @author Florent Dumont <florent.dumont@univresite-paris-saclay.fr>
#' @import ggplot2
#' @importFrom rlang .data
#' @export
ipa <- function( 
    ipafile = "", top = 50, labsize = 8 , dpi = "retina",
    overlapmin = 2 , enascoremin = 1, background = 25000,
    path = "." , dirname = NULL )
{
  i=1
  # output directory
  ifelse( is.null(dirname), "ipa" -> DirName, paste("ipa_",dirname,sep="") -> DirName )
  path %>% file.path( DirName ) -> Path
  Path %>% dir.create
  #
  # read file
  #
  ipafile %>% file("r") -> f
  f %>% readLines() -> rl0
  close(f)
  #
  rl0 %>% grep("^$",.) -> sel
  c("Analysis Details","Ingenuity Canonical Pathways" , "Upstream Regulators",
    "Diseases and Bio Functions","Tox Functions","Regulator Effects","Analysis Ready Molecules") -> IPAGrep0
  #
  # Analysis setting
  #
  rl0 %>% grep(IPAGrep0[1],.) -> IPAGrep1
  rl0[IPAGrep1:(sel[3]-1)] -> rl1
  "IPACoreAnalysisSettings.txt" %>% file.path(Path , . ) %>% file("w") -> f
  rl1 %>% writeLines(f)
  close(f)
  #
  # input List
  #
  rl0 %>% grep(IPAGrep0[7],.) -> IPAGrep1
  rl0[(IPAGrep1+1):(sel[length(sel)]-1)] -> rl1
  rl1[-1] %>% strsplit("\t") %>% unlist -> rl2
  length(rl1) -> listsize
  rl1[1] %>% strsplit("\t") %>% unlist %>% gsub(" ","",.) -> header
  header[1] <- "rowID"
  rl2 %>% matrix( ncol = header %>% length , byrow = T ) %>% data.frame %>% setNames(header) -> rl3
  rl3$rowID %>% strsplit("\\/") %>% lapply("[",1) %>% unlist -> rl3$rowID
  rl3$Symbol %>% strsplit("\\/") %>% lapply("[",1) %>% unlist -> rl3$Symbol
  paste("InputList_",listsize,".tsv",sep="") %>% file.path(Path,.) -> FileName0
  rl3$Symbol -> InputList0
  rl3 %>% output(FileName0)
  #
  # Pathways
  #
  # create results table
  Path %>% file.path("Pathways") %>% dir.create
  rl0 %>% grep(IPAGrep0[2],.) -> IPAGrep2
  rl0[IPAGrep2:(sel[4]-1)] -> rl1
  rl1[-1] %>% strsplit("\t") %>% unlist -> rl2
  rl1[1] %>% strsplit("\t") %>% unlist %>% gsub(" ","",.) -> header
  header[1] <- "Pathways"
  header[2] <- "log10pval"
  rl2 %>% matrix(ncol = header %>% length , byrow = T) %>% data.frame %>% setNames(header) -> rl3
  fO <- foreach(i=1:nrow(rl3),.combine="c") %do% 
    { 
      rl3$Molecules[i] %>% gsub(" \\(includes others\\)","",.) %>% strsplit(",") %>% unlist %>%
        strsplit("\\/") %>% lapply("[",1) %>% unlist %>% paste0(collapse = "|")
    }
  fO -> rl3$Molecules
  rl3$Ratio %>% as.numeric -> rl3$Ratio
  rl3$log10pval %>% as.numeric -> rl3$log10pval
  # compute enrichment score
  f0 <- foreach( i=1:nrow(rl3), .combine = "rbind" ) %do% 
    {
      rl3[i,5] %>% strsplit("\\|") %>% unlist %>% length %>% as.numeric -> OverlapSize
      OverlapSize %>% "/"( rl3[i,4] %>% as.numeric ) %>% round -> GenesetSize
      OverlapSize / GenesetSize -> OverlapRatio
      listsize * ( GenesetSize / background ) -> exp
      OverlapSize / exp -> ENAScore
      c(OverlapSize,listsize-OverlapSize,GenesetSize-OverlapSize,background-OverlapSize-GenesetSize) %>%
        matrix(ncol = 2) %>% fisher.test %>% unlist %>% "["(1) %>% as.numeric -> pval
      c(OverlapSize, GenesetSize, OverlapRatio, ENAScore, pval)
    }
  f0[,5] %>% p.adjust(method = "fdr") -> pvalFDR
  pvalFDR %>% log10 %>% '*'(.,-1) %>% round(.,4) -> log10pvalFDR
  c("Name","SymbolList","OverlapSize","GenesetSize","OverlapRatio","ENAScore","pval","pvalFDR","log10pvalFDR") -> Header0
  rl3 %>% data.frame(f0,pvalFDR,log10pvalFDR) %>% dplyr::select( c(1,5,6:10,11,12) ) %>% setNames(Header0) -> rl4
  rl4 %>% dplyr::filter(OverlapSize > overlapmin & ENAScore > enascoremin) %>% dplyr::arrange( desc(log10pvalFDR ) ) -> rl5
  # output
  paste("Pathways_",nrow(rl5),".tsv",sep="") -> FileName0
  Path %>% file.path("Pathways",FileName0) -> FileName1
  rl5 %>% output(FileName1)
  # output top
  rl5 %>% dplyr::slice(1:top) -> toppval
  rl5 %>% dplyr::arrange(-.data$ENAScore) %>% dplyr::slice(1:top) -> topena
  rbind(toppval,topena) %>% unique -> rl6
  paste("top",top,"_Pathways_",nrow(rl6),".tsv",sep="") -> FileName0
  Path %>% file.path("Pathways",FileName0) -> FileName1
  rl6 %>% output(FileName1)
  # barplot
  rl5 %>% dplyr::select(c(1,9,6)) -> rl6
  ipabarplotcp(dat=rl6, top=top, labsize=labsize) -> p
  # output
  paste("Pathways_barplot_top",top,".jpeg",sep="") -> FileName0
  Path %>% file.path("Pathways",FileName0) -> FileName1
  ggsave(filename=FileName1, plot=p, width=12, height=10, dpi=dpi)
  #
  # Regulators
  #
  # create results table
  Path %>% file.path("Regulators") %>% dir.create
  rl0 %>% grep(IPAGrep0[3],.) -> IPAGrep3
  rl0[IPAGrep3:(sel[5]-1)] -> rl1
  rl1[-c(1:2)] %>% strsplit("\t") %>% unlist -> rl2
  rl1[2] %>% strsplit("\t") %>% unlist %>% gsub(" ","",.) -> header
  rl2 %>% matrix(ncol = header %>% length , byrow = T) %>% data.frame %>% setNames(header) -> rl3
  rl3$UpstreamRegulator -> Name
  rl3$MoleculeType -> MoleculeType
  rl3$`p-valueofoverlap` %>% as.numeric -> pval
  pval %>% p.adjust(method = "fdr") -> pvalFDR
  pvalFDR %>% log10(.) %>% "*"(-1) -> log10pvalFDR
  rl3$Targetmoleculesindataset %>% lapply(FUN=gsub,pattern=",",replacement="\\|") %>% unlist -> SymbolList
  f0 <- foreach(i=1:nrow(rl3),.combine="c") %do% 
    { 
      rl3$Targetmoleculesindataset[i] %>% gsub(" \\(includes others\\)","",.) %>% strsplit(",") %>% unlist %>%
        strsplit("\\/") %>% lapply("[",1) %>% unlist %>% paste0(collapse = "|")
    }
  f0 -> SymbolList
  SymbolList %>% strsplit("\\|") %>% lapply(length) %>% unlist -> OverlapSize
  OverlapSize -> ENAScore
  data.frame(Name,MoleculeType,SymbolList,OverlapSize,ENAScore,pval,pvalFDR,log10pvalFDR) -> rl4
  rl4 %>% dplyr::filter(OverlapSize > overlapmin & ENAScore > enascoremin) %>% dplyr::arrange( desc(log10pvalFDR ) ) -> rl5
  # output
  paste("Regulators_",nrow(rl5),".tsv",sep="") -> FileName0
  Path %>% file.path("Regulators",FileName0) -> FileName1
  rl5 %>% output(FileName1)
  # output top
  rl5 %>% dplyr::slice(1:top) -> rl6
  paste("top",top,"_Regulators_",nrow(rl6),".tsv",sep="") -> FileName0
  Path %>% file.path("Regulators",FileName0) -> FileName1
  rl6 %>% output(FileName1)
  # barplot
  rl5 %>% dplyr::select(Name, log10pvalFDR, OverlapSize) -> rl6
  ipabarplotup(dat=rl6, top=top, labsize=labsize) -> p
  # output
  paste("Regulators_barplot_",top,".jpeg",sep="") -> FileName0
  Path %>% file.path("Regulators",FileName0) -> FileName1
  ggsave(filename=FileName1, plot=p, width=12, height=10, dpi="retina")
  #
  # Functions
  #
  # create results table
  Path %>% file.path("Functions") %>% dir.create
  rl0 %>% grep(IPAGrep0[4],.) -> IPAGrep4
  rl0[IPAGrep4:(sel[7]-1)] -> rl1
  rl1[-c(1:2)] %>% strsplit("\t") %>% unlist -> rl2
  rl1[2] %>% strsplit("\t") %>% unlist %>% gsub(" ","",.) -> header
  rl2 %>% matrix(ncol = header %>% length , byrow = T) %>% data.frame %>% setNames(header) -> rl3
  rl3$DiseasesorFunctionsAnnotation -> Name
  rl3$Categories -> IPACategories
  rl3$Functions -> IPAFunctions
  rl3$Molecules %>% lapply(FUN=gsub,pattern=",",replacement="\\|") %>% unlist -> SymbolList
  SymbolList <- foreach(i=1:nrow(rl3),.combine="c") %do% 
    { 
      rl3$Molecules[i] %>% gsub(" \\(includes others\\)","",.) %>% strsplit(",") %>% unlist %>%
        strsplit("\\/") %>% lapply("[",1) %>% unlist %>% paste0(collapse = "|")
    }
  SymbolList %>% strsplit("\\|") %>% lapply(length) %>% unlist %>% as.numeric -> OverlapSize
  OverlapSize -> ENAScore
  rl3$`p-Value` %>% as.numeric -> pval
  pval %>% p.adjust(method = "fdr") -> pvalFDR
  pvalFDR %>% log10(.) %>% "*"(-1) -> log10pvalFDR
  data.frame( Name,IPACategories,IPAFunctions,SymbolList,OverlapSize,ENAScore,pval,pvalFDR,log10pvalFDR ) -> rl4
  rl4 %>% dplyr::filter(OverlapSize > overlapmin & ENAScore > enascoremin) %>% dplyr::arrange( desc(log10pvalFDR ) ) -> rl5
  # output
  paste("Functions_",nrow(rl5),".tsv",sep="") -> FileName0
  Path %>% file.path("Functions",FileName0) -> FileName1 
  rl5 %>% output(FileName1)
  # output top
  rl5 %>% dplyr::slice(1:top) -> rl6
  paste("top",top,"_Functions_",nrow(rl6),".tsv",sep="") -> FileName0
  Path %>% file.path("Functions",FileName0) -> FileName1
  rl6 %>% output(FileName1)
  # barplot
  rl5 %>% dplyr::select(Name,log10pvalFDR,OverlapSize) -> rl6
  ipabarplotfun(dat=rl6, top=top, labsize=labsize) -> p
  # output
  paste("Functions_barplot_",top,".jpeg",sep="") -> FileName0
  Path %>% file.path("Functions",FileName0) -> FileName1
  ggsave(filename=FileName1, plot=p, width=12, height=10, dpi="retina")
  #
  # Diseases
  #
  # create results table
  Path %>% file.path("Diseases") %>% dir.create
  rl0 %>% grep(IPAGrep0[5],.) -> IPAGrep5
  rl0[IPAGrep5:(sel[8]-1)] -> rl1
  rl1[-c(1:2)] %>% strsplit("\t") %>% unlist -> rl2
  rl1[2] %>% strsplit("\t") %>% unlist %>% gsub(" ","",.) -> header
  rl2 %>% matrix(ncol = header %>% length , byrow = T) %>% data.frame %>% setNames(header) -> rl3
  rl3$DiseasesorFunctionsAnnotation -> Name
  rl3$Categories -> IPACategories
  rl3$Functions -> IPAFunctions
  SymbolList <- foreach(i=1:nrow(rl3),.combine="c") %do% 
    { 
      rl3$Molecules[i] %>% gsub(" \\(includes others\\)","",.) %>% strsplit(",") %>% unlist %>%
        strsplit("\\/") %>% lapply("[",1) %>% unlist %>% paste0(collapse = "|")
    }
  SymbolList %>% strsplit("\\|") %>% lapply(length) %>% unlist %>% as.numeric -> OverlapSize
  OverlapSize -> ENAScore
  rl3$`p-Value` %>% as.numeric -> pval
  pval %>% p.adjust(method = "fdr") -> pvalFDR
  pvalFDR %>% log10(.) %>% "*"(-1) -> log10pvalFDR
  data.frame( Name,IPACategories,IPAFunctions,SymbolList,OverlapSize,ENAScore,pval,pvalFDR,log10pvalFDR ) -> rl4
  rl4 %>% dplyr::filter(OverlapSize > overlapmin & ENAScore > enascoremin) %>% dplyr::arrange( desc(log10pvalFDR ) ) -> rl5
  # output
  paste("Diseases_",nrow(rl5),".tsv",sep="") -> FileName0
  Path %>% file.path("Diseases",FileName0) -> FileName1
  rl5 %>% output(FileName1)
  # output top
  rl5 %>% dplyr::slice(1:top) -> rl6
  paste("top",top,"_Diseases_",nrow(rl6),".tsv",sep="") -> FileName0
  Path %>% file.path("Diseases",FileName0) -> FileName1
  rl6 %>% output(FileName1)
  # barplot
  rl5 %>% dplyr::select(Name,log10pvalFDR,OverlapSize) -> rl6
  ipabarplottox(dat=rl6, top=top, labsize=labsize) -> p
  # output
  paste("Diseases_barplot_",top,".jpeg",sep="") -> FileName0
  Path %>% file.path("Diseases",FileName0) -> FileName1
  ggsave(filename=FileName1, plot=p, width=12, height=10, dpi="retina")
}

