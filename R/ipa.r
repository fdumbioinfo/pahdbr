#' @title ipa
#' @description parse and create plot from ipa all enrichment results
#' @param ipafile character path to the file
#' @param topena numeric top feature to display
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
    ipafile = "", topena = 50, labsize = 8 , dpi = "retina",
    overlapmin = 2 , enascoremin = 1, background = 25000,
    path = "." , dirname = NULL )
{
  i=1
  # create output directory
  ifelse(is.null(dirname), "ipa" -> DirName, paste("ipa_",dirname,sep="") -> DirName)
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
  # ----
  # 1 - Analysis setting
  # ----
  rl0 %>% grep(IPAGrep0[1],.) -> IPAGrep1
  rl0[IPAGrep1:(sel[3]-1)] -> rl1
  "IPACoreAnalysisSettings.txt" %>% file.path(Path , . ) %>% file("w") -> f
  rl1 %>% writeLines(f)
  close(f)
  # ----
  # 2 - input List
  # -----
  rl0 %>% grep(IPAGrep0[7],.) -> IPAGrep1
  rl0[(IPAGrep1+1):(sel[length(sel)]-1)] -> rl1
  rl1[-1] %>% strsplit("\t") %>% unlist -> rl2
  rl2 %>% unique %>% moal::annot(.) -> a0
  length(rl1) -> listsize
  # rl1[1] %>% strsplit("\t") %>% unlist %>% gsub(" ","",.) -> header
  # header[1] <- "rowID"
  # rl2 %>% matrix( ncol = header %>% length , byrow = T ) %>% data.frame %>% setNames(header) -> rl3
  # rl3 %>% head
  # rl3$rowID %>% strsplit("\\/") %>% lapply("[",1) %>% unlist -> rl3$rowID
  # rl3$Symbol %>% strsplit("\\/") %>% lapply("[",1) %>% unlist -> rl3$Symbol
  paste("InputList_",nrow(a0),".tsv",sep="") %>% file.path(Path,.) -> FileName0
  a0 %>% output(FileName0)
  # -----
  # 3 - Pathways
  # -----
  # create results table
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
  f0 <- foreach(i=1:nrow(rl3),.combine="rbind") %do% 
    {
      rl3[i,5] %>% strsplit("\\|") %>% unlist %>% length %>% as.numeric -> OverlapSize
      OverlapSize %>% "/"( rl3[i,4] %>% as.numeric ) %>% round -> GenesetSize
      # OverlapSize / GenesetSize -> OverlapRatio
      listsize * ( GenesetSize / background ) -> exp
      OverlapSize / exp -> ENAScore
      c(OverlapSize,listsize-OverlapSize,GenesetSize-OverlapSize,background-OverlapSize-GenesetSize) %>%
        matrix(ncol = 2) %>% fisher.test %>% unlist %>% "["(1) %>% as.numeric -> pval
      c(OverlapSize,GenesetSize,ENAScore,pval)
    }
  f0[,4] %>% p.adjust(method = "fdr") -> pvalFDR
  pvalFDR %>% log10 %>% '*'(.,-1) %>% round(.,4) -> log10pvalFDR
  c("Name","SymbolList","OverlapSize","GenesetSize","NES","pval","pvalFDR","log10pvalFDR") -> Header0
  rl3 %>% data.frame(f0,pvalFDR,log10pvalFDR) %>% dplyr::select( c(1,5,6:9,10,11) ) %>% setNames(Header0) -> rl4
  rl4 %>% dim
  # rl4 %>% dplyr::filter(OverlapSize > overlapmin & ENAScore > enascoremin) %>% dplyr::arrange( desc(log10pvalFDR ) ) -> rl5
  rl4 %>% dplyr::arrange(desc(log10pvalFDR)) -> rl5
  # output
  paste("PATHWAYS_",nrow(rl5),".tsv",sep="") -> FileName0
  Path %>% file.path(FileName0) -> FileName1
  rl5 %>% output(FileName1)
  rl5 -> rlp
  # output top
  # rl5 %>% dplyr::slice(1:top) -> toppval
  # rl5 %>% dplyr::arrange(-.data$ENAScore) %>% dplyr::slice(1:top) -> topena
  # rbind(toppval,topena) %>% unique -> rl6
  # paste("top",top,"_Pathways_",nrow(rl6),".tsv",sep="") -> FileName0
  # Path %>% file.path("Pathways",FileName0) -> FileName1
  # rl6 %>% output(FileName1)
  # barplot
  # pval ranking
  rl5 %>% colnames
  rl5 %>% dplyr::slice(1:topena) -> fgseapval1plot
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
  fgseapval1plot3 -> datp
  # NES ranking
  rl5 %>% colnames
  rl5 %>% dplyr::arrange( -abs(.data$NES) ) %>% dplyr::slice(1:topena) -> fgseapval1plot
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
  datp %>% enabarplot(datnes=datnes,title="Pathways") -> p
  
  # # barplot 
  # fgseapval1plot3 %>% ggplot( aes(x=Log10Pval,y=.data$Name,fill=.data$NES)) -> p
  # p + geom_bar(stat="identity",width = 0.95) -> p
  # p + scale_color_gradient2(low="blue",mid="white",high="red",aesthetics="fill") -> p
  # # p + theme_minimal() -> p
  # # p + theme_light() -> p
  # # p + theme_linedraw() -> p
  # p + theme_bw() -> p
  # p + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1,size = 5),
  #           axis.text.y = element_text(angle = 0, vjust = 0.4, hjust=1,size = 6)) -> p
  # p + labs(x="log10(p-value)",y="Genesets") -> p
  # p + ggtitle("Canonical Pathways") -> p
  # #
  # T -> addratioena
  # if(addratioena){p + geom_text(data=fgseapval1plot3,aes(label=.data$Overlap),hjust=1.1,size=2) -> p}
  # p + theme(plot.title=element_text(size=10,hjust=0.5),
  #           plot.subtitle=element_text(size=8,hjust=0.5),
  #           axis.title.x=element_text(size=6,face="bold"),
  #           axis.text.x = element_text(size=6,face="bold"),
  #           axis.title.y=element_text(size=6,face="bold"),
  #           axis.text.y=element_text(size=8,face="bold"),
  #           legend.position="bottom",
  #           legend.title=element_text(size = 5),
  #           legend.text=element_text(size = 5)) -> p
  # p + geom_vline(xintercept=-log10(0.05),color="darkgoldenrod1",size=0.4,linetype="dashed") -> p
  # p -> p1
  # T -> addenarankbarplot
  # if(addenarankbarplot)
  # {
  #   fgseapval1plot2 %>% dplyr::mutate(Name=forcats::fct_reorder(.data$Name,abs(.data$NES))) -> fgseapval1plot3
  #   fgseapval1plot3$OverlapSize %>% paste("/",fgseapval1plot3$GenesetSize) %>% 
  #     lapply(paste0,collapse="") %>% unlist -> Overlap
  #   fgseapval1plot3 %>% data.frame(Overlap) -> fgseapval1plot3
  #   fgseapval1plot3 %>% ggplot( aes(x=.data$Log10Pval,y=.data$Name,fill=.data$NES)) -> p
  #   p + geom_bar(stat="identity") -> p
  #   p + scale_color_gradient2(low="blue",mid="white",high="red",aesthetics="fill") -> p
  #   p + theme_bw() -> p
  #   if(addratioena){p + geom_text(data = fgseapval1plot3, aes(label = .data$Overlap),hjust=1.1,size = 2) -> p}
  #   p + theme(plot.title=element_text(size=10,hjust=0.5),
  #             plot.subtitle=element_text(size=8,hjust=0.5),
  #             axis.title.x=element_text(size=6,face="bold"),
  #             axis.text.x=element_text(size=6,face="bold"),
  #             axis.title.y=element_text(size=6,,face="bold"),
  #             axis.text.y=element_text(size=8,face="bold"),
  #             legend.position="bottom",
  #             legend.title=element_text(size = 5), 
  #             legend.text=element_text(size = 5)) -> p
  #   p + labs(x="log10(p-value)",y="Genesets") -> p
  #   p + ggtitle("Canonical Pathways") -> p
  #   p + geom_vline(xintercept=-log10(0.05),color="darkgoldenrod1",size=0.4,linetype="dashed") -> p
  #   p -> p2
  #   # gridExtra::grid.arrange( p1, p2, nrow = 1 ) -> p
  #   # gridExtra::grid.arrange(p1,p2,heights = c(5, 5),nrow=2,) -> p
  #   # gridExtra::grid.arrange(p1,p2,heights = c(10, 10),nrow=2) -> p
  #   # gridExtra::grid.arrange(p1,p2,heights = c(10, 10),nrow=2) -> p
  #   # gridExtra::grid.arrange(p1,p2,heights = c(10, 10),ncol=2) -> p
  #   # gridExtra::grid.arrange(p1,p2,heights = c(5,5),ncol=2) -> p
  #   gridExtra::grid.arrange(p1,p2,ncol=2) -> p
  #   # gridExtra::grid.arrange(p1,p2,heights=c(10, 10)) -> p
  #   # ggpubr::ggarrange(p1, p2, widths = c(5,5)) -> p
  #   # plot_grid(p1, p2, labels = c('A', 'B'), label_size = 12)
  # }
  paste("PATHWAYS_",topena,".pdf",sep="") -> FileName0
  Path %>% file.path(FileName0) -> FileName1
  ggsave(plot=p,filename=FileName1,width=40,height=30,scale=1,units="cm")
  
  
  # rl5 %>% dplyr::select(c(1,9,6)) -> rl6
  # ipabarplotcp(dat=rl6, top=top, labsize=labsize) -> p
  # output
  # paste("Pathways_barplot_top",top,".jpeg",sep="") -> FileName0
  # Path %>% file.path("Pathways",FileName0) -> FileName1
  # ggsave(filename=FileName1, plot=p, width=12, height=10, dpi=dpi)
  # -----
  # 4 - REGULATORS
  # -----
  # create results table
  # Path %>% file.path("Regulators") %>% dir.create
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
  # rl3$Targetmoleculesindataset %>% lapply(FUN=gsub,pattern=",",replacement="\\|") %>% unlist -> SymbolList
  rlp %>% dplyr::arrange(-.data$OverlapSize) -> rlp
  rlp$OverlapSize
  # f0 <- foreach(i=1:nrow(rl3),.combine="rbind") %do% 
  #   { 
  #     rl3$Targetmoleculesindataset[i] %>% gsub(" \\(includes others\\)","",.) %>% strsplit(",") %>% unlist %>%
  #       strsplit("\\/") %>% lapply("[",1) %>% unlist %>% paste0(collapse = "|") -> SymbolList
  #     rl3$Targetmoleculesindataset[i] %>% gsub(" \\(includes others\\)","",.) %>% strsplit(",") %>% unlist %>%
  #       strsplit("\\/") %>% lapply("[",1) %>% unlist %>% length -> OverlapSize
  #     OverlapSize %>% "<"(rlp$OverlapSize) %>% which %>% max %>% rlp$GenesetSize[.] -> GenesetSize
  #     OverlapSize %>% "<"(rlp$OverlapSize) %>% which %>% max %>% rlp$NES[.] -> NES
  #     c(SymbolList,OverlapSize,GenesetSize,NES)
  #   }
  f0 <- foreach(i=1:nrow(rl3),.combine="rbind") %do% 
    { 
      rl3$Targetmoleculesindataset[i] %>% gsub(" \\(includes others\\)","",.) %>% strsplit(",") %>% unlist %>%
        strsplit("\\/") %>% lapply("[",1) %>% unlist %>% paste0(collapse = "|") -> SymbolList
      rl3$Targetmoleculesindataset[i] %>% gsub(" \\(includes others\\)","",.) %>% strsplit(",") %>% unlist %>%
        strsplit("\\/") %>% lapply("[",1) %>% unlist %>% length -> OverlapSize
      OverlapSize -> NES
      NA -> GenesetSize
      c(SymbolList,OverlapSize,GenesetSize,NES)
    }
  f0 %>% data.frame %>% head
  f0 %>% data.frame %>% head
  # f0[,4] %>% p.adjust(method = "fdr") -> pvalFDR
  # pvalFDR %>% log10 %>% '*'(.,-1) %>% round(.,4) -> log10pvalFDR
  c("Name","MoleculeType","SymbolList","OverlapSize","GenesetSize","NES","pval","pvalFDR","log10pvalFDR") -> Header0
  data.frame(Name,MoleculeType,f0 %>% data.frame,pval,pvalFDR,log10pvalFDR) %>% setNames(Header0) -> rl4
  # rl3 %>% data.frame(f0,pvalFDR,log10pvalFDR) %>% dplyr::select( c(1,5,6:9,10,11) ) %>% setNames(Header0) -> rl4
  rl4 %>% head
  rl4 %>% dim
  rl4 %>% str
  mode(rl4$OverlapSize) <- "numeric"
  mode(rl4$GenesetSize) <- "numeric"
  mode(rl4$NES) <- "numeric"
  # rl4 %>% dplyr::filter(OverlapSize > overlapmin & ENAScore > enascoremin) %>% dplyr::arrange( desc(log10pvalFDR ) ) -> rl5
  rl4 %>% dplyr::arrange(desc(log10pvalFDR)) -> rl5
  # f0 -> SymbolList
  # f0 %>% head
  # SymbolList %>% strsplit("\\|") %>% lapply(length) %>% unlist -> OverlapSize
  # OverlapSize -> NES
  # data.frame(Name,MoleculeType,SymbolList,Overlapsize,NES,pval,pvalFDR,log10pvalFDR) -> rl4
  # rl4 %>% dplyr::filter(OverlapSize > overlapmin & ENAScore > enascoremin) %>% dplyr::arrange( desc(log10pvalFDR ) ) -> rl5
  # rl5 %>% colnames
  # output
  paste("REGULATORS_",nrow(rl5),".tsv",sep="") -> FileName0
  Path %>% file.path(FileName0) -> FileName1
  rl5 %>% output(FileName1)
  # output top
  # rl5 %>% dplyr::slice(1:top) -> rl6
  # paste("top",top,"_Regulators_",nrow(rl6),".tsv",sep="") -> FileName0
  # Path %>% file.path("Regulators",FileName0) -> FileName1
  # rl6 %>% output(FileName1)
  # barplot
  # rl5 %>% dplyr::select(Name, log10pvalFDR, OverlapSize) -> rl6
  # ipabarplotup(dat=rl6, top=top, labsize=labsize) -> p
  # barplot
  # pval ranking
  rl5 %>% dplyr::slice(1:topena) -> fgseapval1plot
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
  fgseapval1plot3 %>% head
  fgseapval1plot3 -> datp
  # NES ranking
  rl5 %>% colnames
  rl5 %>% dplyr::arrange( -abs(.data$NES) ) %>% dplyr::slice(1:topena) -> fgseapval1plot
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
  datp %>% enabarplot(datnes=datnes,title="Regulators") -> p
  # fgseapval1plot3 %>% enabarplot(title="Regulators") -> p
  paste("REGULATORS_",topena,".pdf",sep="") -> FileName0
  Path %>% file.path(FileName0) -> FileName1
  ggsave(filename=FileName1, plot=p, width=12, height=10, dpi="retina")
  # -----
  # 6 - Functions
  # -----
  # create results table
  # Path %>% file.path("Functions") %>% dir.create
  rl0 %>% grep(IPAGrep0[4],.) -> IPAGrep4
  rl0[IPAGrep4:(sel[7]-1)] -> rl1
  rl1[-c(1:2)] %>% strsplit("\t") %>% unlist -> rl2
  rl1[2] %>% strsplit("\t") %>% unlist %>% gsub(" ","",.) -> header
  rl2 %>% matrix(ncol = header %>% length , byrow = T) %>% data.frame %>% setNames(header) -> rl3
  rl3$DiseasesorFunctionsAnnotation -> Name
  rl3$Categories -> IPACategories
  rl3$Functions -> IPAFunctions
  # rl3$Molecules %>% lapply(FUN=gsub,pattern=",",replacement="\\|") %>% unlist -> SymbolList
  # SymbolList <- foreach(i=1:nrow(rl3),.combine="c") %do% 
  #   { 
  #     rl3$Molecules[i] %>% gsub(" \\(includes others\\)","",.) %>% strsplit(",") %>% unlist %>%
  #       strsplit("\\/") %>% lapply("[",1) %>% unlist %>% paste0(collapse = "|")
  #   }
  # SymbolList %>% strsplit("\\|") %>% lapply(length) %>% unlist %>% as.numeric -> OverlapSize
  # OverlapSize -> ENAScore
  rl3$`p-Value` %>% as.numeric -> pval
  pval %>% p.adjust(method = "fdr") -> pvalFDR
  pvalFDR %>% log10(.) %>% "*"(-1) -> log10pvalFDR
  f0 <- foreach(i=1:nrow(rl3),.combine="rbind") %do% 
    { 
      rl3$Molecules[i] %>% gsub(" \\(includes others\\)","",.) %>% strsplit(",") %>% unlist %>%
        strsplit("\\/") %>% lapply("[",1) %>% unlist %>% paste0(collapse = "|") -> SymbolList
      rl3$Molecules[i] %>% gsub(" \\(includes others\\)","",.) %>% strsplit(",") %>% unlist %>%
        strsplit("\\/") %>% lapply("[",1) %>% unlist %>% length -> OverlapSize
      OverlapSize -> NES
      NA -> GenesetSize
      c(SymbolList,OverlapSize,GenesetSize,NES)
    }
  f0 %>% dim
  c("Name","IPACategories","IPAFunctions","SymbolList","OverlapSize","GenesetSize","NES","pval","pvalFDR","log10pvalFDR") -> Header0
  # data.frame(Name,MoleculeType,f0 %>% data.frame,pval,pvalFDR,log10pvalFDR) %>% setNames(Header0) -> rl4
  data.frame( Name,IPACategories,IPAFunctions,f0,pval,pvalFDR,log10pvalFDR ) %>% setNames(Header0) -> rl4
  mode(rl4$OverlapSize) <- "numeric"
  mode(rl4$GenesetSize) <- "numeric"
  mode(rl4$NES) <- "numeric"
  
  # rl4 %>% dplyr::filter(OverlapSize > overlapmin & ENAScore > enascoremin) %>% dplyr::arrange( desc(log10pvalFDR ) ) -> rl5
  rl4 %>% dplyr::arrange(desc(log10pvalFDR )) -> rl5
  # output
  paste("FUNCTIONS_",nrow(rl5),".tsv",sep="") -> FileName0
  Path %>% file.path(FileName0) -> FileName1 
  rl5 %>% output(FileName1)
  # output top
  # rl5 %>% dplyr::slice(1:top) -> rl6
  # paste("top",top,"_Functions_",nrow(rl6),".tsv",sep="") -> FileName0
  # Path %>% file.path("Functions",FileName0) -> FileName1
  # rl6 %>% output(FileName1)
  # barplot
  # rl5 %>% dplyr::select(Name,log10pvalFDR,OverlapSize) -> rl6
  # ipabarplotfun(dat=rl6, top=top, labsize=labsize) -> p
  # pval ranking
  rl5 %>% dplyr::slice(1:topena) -> fgseapval1plot
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
  fgseapval1plot3 %>% head
  fgseapval1plot3 -> datp
  # NES ranking
  rl5 %>% colnames
  rl5 %>% dplyr::arrange( -abs(.data$NES) ) %>% dplyr::slice(1:topena) -> fgseapval1plot
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
  datp %>% enabarplot(datnes=datnes,title="Functions") -> p
  
  
  
  # fgseapval1plot3 %>% enabarplot(title="Functions") -> p
  # output
  paste("FUNCTIONS_",topena,".pdf",sep="") -> FileName0
  Path %>% file.path(FileName0) -> FileName1
  ggsave(filename=FileName1, plot=p, width=12, height=10, dpi="retina")
  # -----
  # 5 - Diseases
  # -----
  # create results table
  rl0 %>% grep(IPAGrep0[5],.) -> IPAGrep5
  rl0[IPAGrep5:(sel[8]-1)] -> rl1
  rl1[-c(1:2)] %>% strsplit("\t") %>% unlist -> rl2
  rl1[2] %>% strsplit("\t") %>% unlist %>% gsub(" ","",.) -> header
  rl2 %>% matrix(ncol = header %>% length , byrow = T) %>% data.frame %>% setNames(header) -> rl3
  rl3$DiseasesorFunctionsAnnotation -> Name
  rl3$Categories -> IPACategories
  rl3$Functions -> IPAFunctions
  
  rl3$`p-Value` %>% as.numeric -> pval
  pval %>% p.adjust(method = "fdr") -> pvalFDR
  pvalFDR %>% log10(.) %>% "*"(-1) -> log10pvalFDR
  f0 <- foreach(i=1:nrow(rl3),.combine="rbind") %do% 
    { 
      rl3$Molecules[i] %>% gsub(" \\(includes others\\)","",.) %>% strsplit(",") %>% unlist %>%
        strsplit("\\/") %>% lapply("[",1) %>% unlist %>% paste0(collapse = "|") -> SymbolList
      rl3$Molecules[i] %>% gsub(" \\(includes others\\)","",.) %>% strsplit(",") %>% unlist %>%
        strsplit("\\/") %>% lapply("[",1) %>% unlist %>% length -> OverlapSize
      OverlapSize -> NES
      NA -> GenesetSize
      c(SymbolList,OverlapSize,GenesetSize,NES)
    }
  f0 %>% dim
  c("Name","IPACategories","IPAFunctions","SymbolList","OverlapSize","GenesetSize","NES","pval","pvalFDR","log10pvalFDR") -> Header0
  # data.frame(Name,MoleculeType,f0 %>% data.frame,pval,pvalFDR,log10pvalFDR) %>% setNames(Header0) -> rl4
  data.frame( Name,IPACategories,IPAFunctions,f0,pval,pvalFDR,log10pvalFDR ) %>% setNames(Header0) -> rl4
  mode(rl4$OverlapSize) <- "numeric"
  mode(rl4$GenesetSize) <- "numeric"
  mode(rl4$NES) <- "numeric"
  
  # rl4 %>% dplyr::filter(OverlapSize > overlapmin & ENAScore > enascoremin) %>% dplyr::arrange( desc(log10pvalFDR ) ) -> rl5
  rl4 %>% dplyr::arrange(desc(log10pvalFDR )) -> rl5
  
  
  # SymbolList <- foreach(i=1:nrow(rl3),.combine="c") %do% 
  #   { 
  #     rl3$Molecules[i] %>% gsub(" \\(includes others\\)","",.) %>% strsplit(",") %>% unlist %>%
  #       strsplit("\\/") %>% lapply("[",1) %>% unlist %>% paste0(collapse = "|")
  #   }
  # SymbolList %>% strsplit("\\|") %>% lapply(length) %>% unlist %>% as.numeric -> OverlapSize
  # OverlapSize -> ENAScore
  # rl3$`p-Value` %>% as.numeric -> pval
  # pval %>% p.adjust(method = "fdr") -> pvalFDR
  # pvalFDR %>% log10(.) %>% "*"(-1) -> log10pvalFDR
  # data.frame( Name,IPACategories,IPAFunctions,SymbolList,OverlapSize,ENAScore,pval,pvalFDR,log10pvalFDR ) -> rl4
  # rl4 %>% dplyr::filter(OverlapSize > overlapmin & ENAScore > enascoremin) %>% dplyr::arrange( desc(log10pvalFDR ) ) -> rl5
  # output
  paste("DISEASES_",nrow(rl5),".tsv",sep="") -> FileName0
  Path %>% file.path(FileName0) -> FileName1
  rl5 %>% output(FileName1)
  # output top
  # rl5 %>% dplyr::slice(1:topena) -> rl6
  # paste("top",top,"_Diseases_",nrow(rl6),".tsv",sep="") -> FileName0
  # Path %>% file.path("Diseases",FileName0) -> FileName1
  # rl6 %>% output(FileName1)
  # barplot
  # rl5 %>% dplyr::select(Name,log10pvalFDR,OverlapSize) -> rl6
  # ipabarplottox(dat=rl6, top=top, labsize=labsize) -> p
  # pval ranking
  rl5 %>% dplyr::slice(1:topena) -> fgseapval1plot
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
  fgseapval1plot3 %>% head
  fgseapval1plot3 -> datp
  # NES ranking
  rl5 %>% colnames
  rl5 %>% dplyr::arrange( -abs(.data$NES) ) %>% dplyr::slice(1:topena) -> fgseapval1plot
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
  datnes %>% dim
  datp %>% enabarplot(datnes=datnes,title="Diseases") -> p
  
  
  # fgseapval1plot3 %>% enabarplot(title="Diseases") -> p
  # output
  paste("DISEASES_",topena,".pdf",sep="") -> FileName0
  Path %>% file.path(FileName0) -> FileName1
  ggsave(filename=FileName1, plot=p, width=12, height=10, dpi="retina")
}

