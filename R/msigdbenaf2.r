#' @title MSigDB Enrichment Analysis For Two list.
#' @description MSigDB Enrichment Analysis For Two list
#' @param list1 character Symbol list 1
#' @param list2 character Symbol list 2
#' @param addlist character additional Symbol list 2
#' @param listnames character vector with 3 list names respectively
#' @param species character hs mm rn dr ss
#' @param collection character MSigDB collection name to extract top pathways
#' @param top numeric number of top patways to extract
#' @param intmaxdh numeric maximum number of interaction to use for Davidson and Harel algorithm layout
#' @param mings numeric minimal size of a gene set
#' @param maxgs numeric maximal size of a gene set
#' @param overlapmin numeric minimal overlap to keep for gene set analysis
#' @param dirname character output directory name
#' @param path character path
#' @return directory
#' @examples
#' # not run
#' # menaf2( list1 , list2 )
#' @author Florent Dumont <florent.dumont@universite-paris-saclay.fr>
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate arrange select
#' @importFrom rlang .data
#' @importFrom stats fisher.test setNames p.adjust
#' @importFrom foreach foreach %do% %:%
#' @importFrom utils capture.output
#' @export
msigdbenaf2 <- function(
    list1 = NULL, list2 = NULL, addlist = NULL, listnames = NULL,species = "hs",
    top = 50, collection = NULL, intmaxdh = 5000,mings = 5, maxgs = 500, overlapmin = 2,
    dirname = NULL, path = ".")
{
  i=j=1
  if(is.null(collection)){"reactome" -> collection}
  if(length(collection) == 1 & collection[1] == "pathways"){ c("reactome","kegg","pid","biocarta","cgp","hallmark","kegg") -> collection}
  if((length(collection) == 1) & collection[1] == "c2"){ c("reactome","kegg","pid","biocarta","cgp","hallmark","kegg") -> collection}
  if(is.null(listnames)){ c("list1","list2","addlist") -> listnames }
  paste(listnames[1],"x",listnames[2],sep="") -> DirName0
  ifelse(is.null(dirname),paste("msigdbenaf2_",DirName0,sep="") -> DirName1,paste("msigdbenaf2_",dirname,sep="") -> DirName1 )
  path %>% file.path(DirName1) -> Path0
  if(!dir.exists(Path0)){ Path0 %>% dir.create }
  # addlist = NULL
  list(list1,list2) %>% moal::venn(listnames=listnames[c(1,2)],export=T,dirname="venn2",path = Path0)
  Path0 %>% list.files(full.names = T) %>% grep("venn2",.,value=T) %>% list.files(full.names=T) %>% grep("IntAB",.,value=T) -> i0
  i0 %>% moal::input(.) -> i1
  if(length(i1$rowID)>0){ i1$rowID -> i2 }else{ NULL -> i2}
  # Venn diagramm with addlist and add addlist to overlap list
  if(!is.null(addlist) & length(addlist)>0)
  {
    list(list1,list2,addlist) %>% 
      moal::venn(listnames=listnames %>% c("addlist"),export=T,dirname="venn3",path=Path0)
    if(!is.null(i2)){ c(i2,addlist) %>% unique -> i2 }else{addlist %>% unique -> i2}
  }
  # output overlap (with addlist if not NULL)
  FALSE -> OVERLAP 
  if(length(i2)>0)
  { 
    i2 %>% moal::annot(species=species) -> a0
    paste(DirName0,"_",length(i2),".tsv",sep="") -> FileNamei
    a0 %>% moal::output(Path0 %>% file.path(FileNamei))
    TRUE -> OVERLAP
  }
  #
  # ORA
  #
  # list1 ORA
  if(length(list1)>0)
  {
    if(length(list1)>200){ 1 -> layout }else{ 2 -> layout }
    paste(listnames[1],"_",length(list1),sep="") -> DirNamel1
    # data.frame(Symbol=list1,FC=rep(5,length(list1))) -> foldchange
    # list1 %>% pahdb(foldchange=foldchange,species=species,top=top,dirname=DirNamel1,path=Path0,dotopgo=F)
    list1 %>% moal::ena(species=species,topdeg=length(list1),layout=layout,intmaxdh=intmaxdh,
                        mings=mings,maxgs=maxgs,overlapmin=overlapmin,dirname=DirNamel1,path=Path0)
  }
  # list1 ORA
  if(length(list2)>0)
  {
    if(length(list2)>200){ 1 -> layout }else{ 2 -> layout }
    paste(listnames[2],"_",length(list2),sep="") -> DirNamel2
    # data.frame(Symbol=list2,FC=rep(5,length(list2))) -> foldchange
    # list2 %>% pahdb(foldchange=foldchange,species=species,top=top,dirname=DirNamel2,path=Path0,dotopgo=F)
    list2 %>% moal::ena(species=species,topdeg=length(list2),layout=layout,intmaxdh=intmaxdh,
                        mings=mings,maxgs=maxgs,overlapmin=overlapmin,dirname=DirNamel2,path=Path0)
  }
  #
  # if overlap or addlist not NULL
  #
  if(OVERLAP)
  {
    # ORA overlap
    if(length(i2)>200){ 1 -> layout }else{ 2 -> layout }
    paste(DirName0,"_",length(i2),sep="") -> DirNamei2
    i2 %>% moal::ena(species=species,topdeg=length(i2),layout=layout,intmaxdh=intmaxdh,
                     mings=mings,maxgs=maxgs,overlapmin=overlapmin,dirname=DirNamei2,path=Path0)
    #
    Path0 %>% list.files(full.names=T,recursive=T) -> lp0
    # list1
    p0 <- foreach(i=1:length(collection),.combine = "rbind") %do%
      {
        lp0 %>% grep(collection[i],.,ignore.case=T) %>% lp0[.] -> lp1
        lp1 %>% basename %>% grep(".*.tsv$",.) %>% lp1[.] -> lp2
        paste("_",listnames[1],"_",sep="") -> Grep0
        lp2 %>% grep(Grep0,.) %>% lp2[.] -> lp3
        lp3 %>% moal::input(.) %>% dplyr::slice(1:top)
      }
    # list2
    pp0 <- foreach(i=1:length(collection),.combine = "rbind") %do%
      {
        lp0 %>% grep(collection[i],.,ignore.case=T) %>% lp0[.] -> lp1
        lp1 %>% basename %>% grep(".*.tsv$",.) %>% lp1[.] -> lp2
        paste("_",listnames[2],"_",sep="") -> Grep0
        lp2 %>% grep(Grep0,.) %>% lp2[.] -> lp3
        lp3 %>% moal::input(.) %>% dplyr::slice(1:top)
      }
    # addlist
    ppp0 <- foreach(i=1:length(collection),.combine = "rbind") %do%
      {
        lp0 %>% grep(collection[i],.,ignore.case=T) %>% lp0[.] -> lp1
        lp1 %>% basename %>% grep(".*.tsv$",.) %>% lp1[.] -> lp2
        paste("_",DirName0,"_",sep="") -> Grep0
        lp2 %>% grep(Grep0,.) %>% lp2[.] -> lp3
        lp3 %>% moal::input(.) %>% dplyr::slice(1:top)
      }
    # p0 %>% grep(collection,.,ignore.case=T) %>% p0[.] -> p1
    # p1 %>% grep("top.*.tsv$",.) %>% p1[.] -> p2
    # p1 %>% basename %>% grep("^c[0-9]_.*.tsv$",.) %>% p1[.] -> p2
    # p1 %>% basename %>% grep(".*.tsv$",.) %>% p1[.] -> p2
    # p2 %>% grep(paste("_",listnames[1],"_",sep=""),.) %>% p2[.] %>% moal::input(.) -> p0
    # p0 %>% dplyr::slice(1:top) -> p0
    # p2 %>% grep(paste("_",listnames[2],"_",sep=""),.) %>% p2[.] %>% moal::input(.) -> pp0
    # pp0 %>% dplyr::slice(1:top) -> pp0
    # p2 %>% grep(paste("_",DirName0,"_",sep=""),.) %>% p2[.] %>% moal::input(.) -> ppp0
    # ppp0 %>% dplyr::slice(1:top) -> ppp0
    if(length(collection)>0){ paste0(collection,collapse="-") -> Collection0 }else{ collection -> Collection0}
    paste("venn3_pathways_",Collection0,sep="") -> DirNameP
    c(DirName0,listnames[1],listnames[2]) -> listnamesP
    list(ppp0$Name,p0$Name,pp0$Name) %>% moal::venn(listnames=listnamesP,dirname=DirNameP,path=Path0,export=T)
    #
    foreach(i=1:length(collection)) %do%
      {
        lp0 %>% grep(collection[i],.,ignore.case=T) %>% lp0[.] -> lp1
        lp1 %>% basename %>% grep(".*.tsv$",.) %>% lp1[.] -> lp2
        paste("_",DirName0,"_",sep="") -> Grep0
        lp2 %>% grep(Grep0,.) %>% lp2[.] -> lp3
        lp3 %>% moal::input(.) %>% dplyr::slice(1:top) -> lp4
        paste("top",top,"_",collection[i],"_",DirName0,"_",length(i2),"_",nrow(lp4),".tsv",sep="") -> FileNameppp1
        lp4 %>% moal::output(Path0 %>% file.path(FileNameppp1))
      }
    #
    # add annotation
    p0 %>% setNames(c("Name",paste(colnames(p0)[-1],listnames[1],sep=""))) -> p1
    pp0 %>% setNames(c("Name",paste(colnames(pp0)[-1],listnames[2],sep=""))) -> pp1
    ppp0 %>% setNames(c("Name",paste(colnames(pp0)[-1],DirName0,sep=""))) -> ppp1
    Path0 %>% list.files(full.names = T) %>% grep("venn3_pathways_",.,value=T) %>%
      list.files(full.names = T) %>% grep(".tsv$",.,value=T) -> v0
    tt <- foreach(i=1:length(v0)) %do% { v0[i] %>% input -> t ; t$rowID }
    tt %>% unlist %>% unique %>% data.frame(Name=.) %>% left_join(ppp1) %>% left_join(pp1) %>% left_join(p1) -> pAll
    # add annotation to canonical pathway lists
    foreach(i=1:length(v0)) %do%
      {
        v0[i] %>% input %>% setNames(c("Name")) %>% inner_join(pAll) -> t
        t %>% apply(2,is.na) %>% apply(2,all) %>% "!"(.) %>% which %>% dplyr::select(t,.) -> tt
        tt %>% moal::output(v0[i])
      }
    #
    # top pathways
    #
    # top %>% paste("top",.,"_pathways_",collection,"_overlap",sep="") -> DirNameTP
    # Path0 %>% file.path(DirNameTP) %>% dir.create
    # v0 %>% grep("SpeA|SpeIntAB|IntABC|SpeIntAC",.,value=T) %>% file.copy(Path0 %>% file.path(DirNameTP))
    # Path0 %>% file.path(DirNameTP) %>% list.files(full.names = T) -> l0
    # l0 %>% basename %>% strsplit("_") %>% lapply("c",paste(DirName0,"_",length(i2),"_pathways",sep="")) %>% 
    #   lapply("[",c(1,4,2,3)) %>% lapply("paste0",collapse="_") %>% file.path(l0 %>% dirname,.) %>% file.rename(l0,.)
    #
    # gmt
    #
    Path0 %>% file.path("gmt") %>% dir.create
    moalannotgene::genesetdb -> Genesetdb0
    Genesetdb0 %>% names
    # Genesetdb0 %>% names %>% grep(collection,.,ignore.case=T) %>% Genesetdb0[.] %>% unlist(recursive = F) -> Genesetdb1
    Genesetdb1 <- foreach(i=1:length(collection),.combine = "c") %do%
      {
        Genesetdb0 %>% names %>% grep(collection[i],.,ignore.case=T) %>% Genesetdb0[.] %>% unlist(recursive=F)
      }
    Genesetdb1 %>% names
    t %>% length
    # Genesetdb1 %>% length
    # Genesetdb1 %>% names
    # Genesetdb1 %>% head
    # Genesetdb1[[1]]
    Genesetdb1 %>% lapply("[",1) %>% unlist -> GeneSetName0
    Genesetdb1 %>% lapply("[",4) -> GeneSetList0
    # GeneSetList0 %>% head
    # GeneSetName0 %>% head
    # GeneSetName0 %>% length
    # Path0 %>% file.path("gmt") %>% dir.create 
    # create C1_NMDARxPAH-KB_152 gmt file
    # "NMDARxPAH-KB/C1_NMDARxPAH-KB_152.tsv" %>% moal::input(.) -> l0
    #
    # paste(DirName0,"_",length(i2),".gmt",sep="") -> gmtFileName
    # gmt overlap
    paste("top",top,"_",DirName0,"_",length(i2),".gmt",sep="") -> gmtFileName
    paste(DirName0,"_",length(i2),sep="") -> gmtName
    # Path0 %>% file.path("gmt",gmtFileName) %>% file("w") -> f
    Path0 %>% file.path("gmt",gmtFileName) %>% file("w") -> f
    i2 %>% strsplit("\\|") %>% unlist %>% paste0(collapse="\t" ) -> t
    paste(gmtName,"\t",gmtName,"\t",t,"\n",sep="") %>% writeLines(f,sep="") 
    close(f)
    # gmt by collection
    foreach(i=1:length(collection)) %do%
      {
        lp0 %>% grep(collection[i],.,ignore.case=T) %>% lp0[.] -> lp1
        lp1 %>% basename %>% grep(".*.tsv$",.) %>% lp1[.] -> lp2
        paste("_",DirName0,"_",sep="") -> Grep0
        lp2 %>% grep(Grep0,.) %>% lp2[.] -> lp3
        lp3 %>% moal::input(.) %>% dplyr::slice(1:top) -> lp4
        paste("top",top,"_",collection[i],"_",DirName0,"_",length(i2),"_",nrow(ppp0),".gmt",sep="") -> gmtFileName
        paste(DirName0,"_",length(i2),sep="") -> gmtName
        # Path0 %>% file.path("gmt",gmtFileName) %>% file("w") -> f
        Path0 %>% file.path("gmt",gmtFileName) %>% file("w") -> f
        # i2 %>% strsplit("\\|") %>% unlist %>% paste0(collapse="\t" ) -> t
        # paste(gmtName,"\t",gmtName,"\t",t,"\n",sep="") %>% writeLines(f,sep="")
        ppp0
        foreach(j=1:nrow(lp4)) %do%
          {
            GeneSetName0 %>% grep(lp4[j,1],.) %>% GeneSetList0[.] %>% unlist %>% moal::annot(.) %>% dplyr::select(Symbol) %>%
              unlist %>% strsplit("\\|") %>% unlist %>% paste0(collapse="\t" ) -> t
            # ll0[j,2] %>% strsplit("\\|") %>% unlist %>% paste0(collapse="\t" ) -> t
            paste(lp4[j,1],"\t",lp4[j,1],"\t",t,"\n",sep="") %>% writeLines(f,sep="") 
          }
        close(f)
        
      }
    
    
    
    
    
    
    
    
    
    # create top pathways gmt files
    # Path0 %>% file.path(DirNameTP) %>% list.files(full.names = T) -> l0
    # foreach(i=1:length(l0)) %do%
    #   {
    #     l0[i] %>% basename %>% sub("(.*).tsv","\\1",.) %>% paste(".gmt",sep="") -> FileName
    #     Path0 %>% file.path("gmt",FileName) %>% file("w") -> f
    #     l0[i] %>% moal::input(.) -> ll0
    #     foreach(j=1:nrow(ll0)) %do%
    #       {
    #         GeneSetName0 %>% grep(ll0[j,1],.) %>% GeneSetList0[.] %>% unlist %>% moal::annot(.) %>% dplyr::select(Symbol) %>%
    #           unlist %>% strsplit("\\|") %>% unlist %>% paste0(collapse="\t" ) -> t
    #         # ll0[j,2] %>% strsplit("\\|") %>% unlist %>% paste0(collapse="\t" ) -> t
    #         paste(ll0[j,1],"\t",ll0[j,1],"\t",t,"\n",sep="") %>% writeLines(f,sep="") 
    #       }
    #     close(f)
    #   }
  }
}
