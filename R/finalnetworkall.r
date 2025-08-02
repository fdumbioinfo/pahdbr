#' @title MSigDB summary network
#' @description MSigDB summary network
#' @param resultpath character msigdb result directory path
#' @param symbollist character Symbol or NCBI gene ID
#' @param foldchange data.frame Symbol list with fold-change values to display
#' @param top numeric top enrichment score geneset selection
#' @param intmaxdh numeric maximum number of interaction to use for Davidson and Harel algorithm layout
#' @return directory
#' @author Florent Dumont <florent.dumont@universite-paris-saclay.fr>
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate arrange select desc
#' @importFrom foreach foreach %do%
#' @importFrom grDevices graphics.off pdf
#' @importFrom graphics par
#' @noRd
finalnetworkall <- function(resultpath = ".", symbollist = NULL, top = 50, keywords = NULL, foldchange = NULL, intmaxdh = 2500)
{
  i=j=k=l=m=p=1
  moal:::palette0 -> Palette0
  resultpath %>% file.path("Final-network-all") %>% dir.create
  if(is.null(foldchange)){ data.frame( Symbol=symbollist, FC=rep(5,length(symbollist)) ) -> FoldChange0 }else{ foldchange -> FoldChange0  }
  #
  # top geneset and keywords filtering
  #
  if(!is.null(keywords))
  {
    keywords -> Keywords0
    #
    parallel::makeCluster(parallel::detectCores()) -> cl;registerDoParallel(cl);doParallel::registerDoParallel(cl)
    foreach(i=1:length(Keywords0),.packages=c("moal","magrittr","dplyr","gplots","foreach","stringr","igraph","Rgraphviz","moalstringdbhs")) %dopar%
      {
        Keywords0[[i]][1] %>% unlist %>% paste0(collapse = "|") -> Keywords1
        Keywords0[[i]][2] %>% unlist -> GrepName
        lpmsigdb <- NULL
        lemsigdb <- NULL
        lpemsigdb <- NULL
        # MSigDB
        resultpath %>% list.files(full.names = T, recursive = T) -> l0
        l0 %>% basename %>% grep("^c[0-9]_.*.tsv$",.) %>% l0[.] -> l1
        lpmsigdb <- foreach(j=1:length(l1),.combine = "rbind") %do%
          {
            l1[j] %>% input(.) -> ll0
            l1[j] %>% basename %>% strsplit("_") %>% unlist %>% "["(3) -> db
            ll0$Name %>% gsub("_"," ",.) %>% grep(Keywords1,.,value = F, ignore.case=T) -> sel
            if(length(sel) > 0){ 
              ll0[sel,] %>% dplyr::select(.data$Name,.data$SymbolList,.data$OverlapSize,.data$ENAScore,.data$log10pvalFDR) %>% data.frame(Collection=db, rankpval=sel) } 
          }
        lpmsigdb %>% head
        lpmsigdb %>% dim
        lemsigdb <- foreach(j=1:length(l1),.combine = "rbind") %do%
          {
            l1[j] %>% input(.) -> ll0
            ll0 %>% arrange(desc(.data$ENAScore)) -> ll0
            l1[j] %>% basename %>% strsplit("_") %>% unlist %>% "["(3) -> db
            ll0$Name %>% gsub("_"," ",.) %>% grep(Keywords1,.,value = F, ignore.case=T) -> sel
            if(length(sel) > 0){ 
              ll0[sel,] %>% dplyr::select(.data$Name,.data$SymbolList,.data$OverlapSize,.data$ENAScore,.data$log10pvalFDR) %>% data.frame(Collection=db, rankena=sel) } 
          }
        if( is.null(lpmsigdb) )
        { data.frame(Name=character(),SymbolList=character(),ENAScore=numeric(),log10pvalFDR=numeric(),
                     Collection=character(),rankpval=character(),rankena=character() ) -> lpmsigdb }
        if( is.null(lemsigdb) )
        { data.frame(Name=character(),SymbolList=character(),ENAScore=numeric(),log10pvalFDR=numeric(),
                     Collection=character(),rankpval=character(),rankena=character() ) -> lemsigdb }
        if( !nrow(lpmsigdb) > 0 & !nrow(lemsigdb) > 0 ){ lpmsigdb %>% inner_join(lemsigdb) -> lpemsigdb }
        if( nrow(lpmsigdb) > 0 & nrow(lemsigdb) > 0 ){ lpmsigdb %>% inner_join(lemsigdb) -> lpemsigdb }
        # lpemsigdb %>% unique -> lpeAll
        # IPA
        lpipa <- NULL
        leipa <- NULL
        lpeipa <- NULL
        resultpath %>% file.path("ipa") %>% list.files(full.names=T, recursive=T) -> l0
        l0 %>% basename %>% grep(".*.tsv$",.) %>% l0[.] -> l1
        l1 %>% basename %>% grep("Input",.,invert=T) %>% l1[.] -> l2
        l1 %>% basename %>% grep("top",.,invert=T) %>% l1[.] -> l2
        lpipa <- foreach(j=1:length(l2),.combine = "rbind") %do%
          {
            l2[j] %>% input(.) -> ll0
            l2[j] %>% basename %>% strsplit("_") %>% unlist %>% "["(1) %>% paste("IPA",.,sep="") -> db
            ll0$Name %>% grep(Keywords1,.,value = F, ignore.case=T) -> sel
            if(length(sel) > 0){ 
              ll0[sel,] %>% dplyr::select(.data$Name,.data$SymbolList,.data$OverlapSize,.data$ENAScore,.data$log10pvalFDR) %>%
                data.frame(Collection=db, rankpval=sel) } 
          }
        lpipa %>% head
        lpipa %>% dim
        leipa <- foreach(j=1:length(l2),.combine = "rbind") %do%
          {
            l2[j] %>% input(.) -> ll0
            l2[j] %>% basename %>% strsplit("_") %>% unlist %>% "["(1) %>% paste("IPA",.,sep="") -> db
            ll0$Name %>% grep(Keywords1,.,value = F, ignore.case=T) -> sel
            if(length(sel) > 0){ 
              ll0[sel,] %>% dplyr::select(.data$Name,.data$SymbolList,.data$OverlapSize,.data$ENAScore,.data$log10pvalFDR) %>%
                data.frame(Collection=db, rankena=sel) } 
          }
        leipa %>% head
        leipa %>% dim
        if( is.null(lpipa) )
        { data.frame(Name=character(),SymbolList=character(),ENAScore=numeric(),log10pvalFDR=numeric(),
                     Collection=character(),rankpval=character(),rankena=character() ) -> lpipa }
        if( is.null(leipa) )
        { data.frame(Name=character(),SymbolList=character(),ENAScore=numeric(),log10pvalFDR=numeric(),
                     Collection=character(),rankpval=character(),rankena=character() ) -> leipa }
        if( !nrow(lpipa) > 0 & !nrow(leipa) > 0 ){ lpipa %>% inner_join(leipa) -> lpeipa }
        if( nrow(lpipa) > 0 & nrow(leipa) > 0 ){ lpipa %>% inner_join(leipa) -> lpeipa }
        # topGO
        lptopgo <- NULL
        letopgo <- NULL
        lpetopgo <- NULL
        resultpath %>% file.path("topGO") %>% list.files(full.names=T,recursive=T) -> l0
        l0 %>% basename %>% grep(".tsv$",.) %>% l0[.] -> l1
        l1 %>% basename %>% grep("top",.,invert = T) %>% l1[.] -> l1
        lptopgo <- foreach(j=1:length(l1),.combine="rbind") %do%
          {
            l1[j] %>% input(.) -> ll0
            l1[j] %>% basename %>% strsplit("_") %>% unlist %>% "["(1) %>% paste("topGO",.,sep="") -> db
            ll0$Name %>% grep(Keywords1,.,value=F,ignore.case=T) -> sel
            if(length(sel) > 0){ 
              ll0[sel,] %>% dplyr::select(.data$Name,.data$SymbolList,.data$OverlapSize,.data$ENAScore,.data$log10pvalFDR) %>%
                data.frame(Collection=db, rankpval=sel) } 
          }
        lptopgo %>% head
        lptopgo %>% dim
        letopgo <- foreach(j=1:length(l1),.combine = "rbind") %do%
          {
            l1[j] %>% input(.) -> ll0
            l1[j] %>% basename %>% strsplit("_") %>% unlist %>% "["(1) %>% paste("topGO",.,sep="") -> db
            ll0$Name %>% grep(Keywords1,.,value = F, ignore.case=T) -> sel
            if(length(sel) > 0){ 
              ll0[sel,] %>% dplyr::select(.data$Name,.data$SymbolList,.data$OverlapSize,.data$ENAScore,.data$log10pvalFDR) %>%
                data.frame(Collection=db, rankena=sel) } 
          }
        letopgo %>% head
        letopgo %>% dim
        letopgo %>% colnames
        if( is.null(lptopgo) )
        { data.frame(Name=character(),SymbolList=character(),ENAScore=numeric(),log10pvalFDR=numeric(),
                     Collection=character(),rankpval=character(),rankena=character() ) -> lptopgo }
        if( is.null(letopgo) )
        { data.frame(Name=character(),SymbolList=character(),ENAScore=numeric(),log10pvalFDR=numeric(),
                     Collection=character(),rankpval=character(),rankena=character() ) -> letopgo }
        if( !nrow(lptopgo) > 0 & !nrow(letopgo) > 0 ){ lptopgo %>% inner_join(letopgo) -> lpetopgo }
        if( nrow(lptopgo) > 0 & nrow(letopgo) > 0 ){ lptopgo %>% inner_join(letopgo) -> lpetopgo }
        # merging all database results
        lpemsigdb %>% colnames
        lpeipa %>% colnames
        lpetopgo %>% colnames
        rbind(lpemsigdb,lpeipa,lpetopgo) %>% unique -> lpeAll
        lpeAll %>% head
        lpeAll %>% dim
        paste(GrepName,"_results_",nrow(lpeAll),".tsv",sep="") -> FileName0
        resultpath %>% file.path("Final-network-all",FileName0) -> FileName1
        lpeAll %>% output(FileName1)
      }
    parallel::stopCluster(cl)
    #
    # Occurences and Final network
    #
    resultpath %>% file.path("Final-network-all") %>% list.files(full.names = T) %>% grep("results",.,value = T) -> l0
    l0 %>% basename %>% sub(".*_(.*)\\.tsv","\\1",.) %>% as.numeric %>% ">"(0) %>% which %>% l0[.] -> l1
    # Path0 %>% list.files(full.names=T,recursive=F) %>% grep("nw_",.,value = T) %>%
    # list.files(full.names=T,recursive=F) %>% grep("nw_",.,value = T) -> conf0
    #
    parallel::makeCluster(parallel::detectCores()) -> cl;registerDoParallel(cl);doParallel::registerDoParallel(cl)
    foreach(i=1:length(l1),.packages=c("moal","magrittr","dplyr","gplots","foreach","stringr","igraph","Rgraphviz","moalstringdbhs")) %dopar%
      {
        # occurrence
        l1[i] %>% basename %>% strsplit("_") %>% unlist %>% "["(1) -> GrepDirName0
        resultpath %>% file.path("Final-network-all",GrepDirName0) %>% dir.create
        l1[i] %>% file.copy(resultpath %>% file.path("Final-network-all",GrepDirName0))
        l1[i] %>% input(.) -> ll0
        ll0 %>% head
        # l1[i] %>% file.remove
        ll0$Collection %>% table %>% names -> Collection0
        # Symbol occurrence for significant enrich canonical pathways 
        f2 <- foreach(j=1:length(Collection0),.combine = "rbind") %do%
          {
            ll0 %>% dplyr::filter(.data$Collection == Collection0[j]) -> t0
            t0$SymbolList %>% strsplit("\\|") %>% unlist -> t1
            if( table(t1) %>% length == 1 ){ data.frame(Symbol=t1,Occ=1) -> t2 }
            if( table(t1) %>% length > 1 ){ t1 %>% table %>% sort(decreasing = T) %>% data.frame %>% setNames(c("Symbol","Occ")) -> t2 }
            t2 %>% data.frame("Collection"=rep(Collection0[j],nrow(t2)))
          }
        f2$Symbol %>% annot(.) -> a0
        a0 %>% dplyr::select(-.data$Symbol) %>% data.frame(f2,.) -> f3
        paste("Symbol-Occurences_",GrepDirName0,"_",nrow(f3),".tsv",sep="") -> FileName0
        resultpath %>% file.path("Final-network-all",GrepDirName0,FileName0) -> FileName1
        f3 %>% dplyr::arrange(-(.data$Occ)) -> f4
        f4 %>% output(FileName1)
        # Occurence Matrix
        symbollist %>% unique %>% data.frame(Symbol=.) -> m0
        m1 <- foreach(j=1:length(Collection0),.combine="cbind") %do%
          {
            f4 %>% dplyr::filter(.data$Collection == Collection0[j]) %>% dplyr::select(.data$Symbol,.data$Occ) %>% unique -> t0
            m0 %>% dplyr::left_join(t0) %>% dplyr::select(.data$Occ)
          }
        m2 <- foreach(j=1:ncol(m1),.combine="cbind") %do% { m1[j] %>% unlist %>% is.na %>% which %>% replace(m1[j] %>% unlist,.,0) }
        m2 %>% data.frame(m0$Symbol,.) %>% setNames(c("Symbol",Collection0)) -> m3
        m3$Symbol %>% annot(.) -> a0
        m3 %>% data.frame(a0[,c(-1)]) %>% dplyr::arrange(desc(.[[2]])) -> m4
        m4 %>% head
        paste("Symbol-Occurences_",GrepDirName0,"_matrix_",nrow(m4),".tsv",sep="") -> FileName0
        resultpath %>% file.path("Final-network-all",GrepDirName0,FileName0) -> FileName1
        m4 %>% output(FileName1)
        # Occurence Matrix hierarchical clustering
        m4 %>% dplyr::select(c(2:(length(Collection0)+1)))  -> m5
        paste("Symbol-Occurences_",GrepDirName0,"_matrix_HC_",nrow(m5),".pdf",sep="") -> FileName0
        resultpath %>% file.path("Final-network-all",GrepDirName0,FileName0) -> FileName1
        if(ncol(m5)>2)
        {
          if(m5 %>% apply(2,sd) %>% "=="(0) %>% which %>% any )
          {
            m5 %>% apply(2,sd) %>% "=="(0) %>% which -> selsd
            m5[,-selsd] -> m5
          }
          pdf(FileName1)
          par(mar=c(6.1,4.1,4.1,2.1))
          m5 %>% moal:::hc(m5 %>% colnames %>% as.factor ,legendtitle="Collection",legend=F,cexlabel=0.9)
          graphics.off()
        }
      }
    parallel::stopCluster(cl)
    #
    # network
    #
    resultpath %>% file.path("Networks") %>% list.files(full.names=T,recursive=F) %>% grep("nw",.,value = T) -> conf0
    #
    parallel::makeCluster(parallel::detectCores()) -> cl;registerDoParallel(cl);doParallel::registerDoParallel(cl)
    foreach(i=1:length(l1),.packages=c("moal","magrittr","dplyr","gplots","foreach","stringr","igraph","Rgraphviz","moalstringdbhs")) %dopar%
      {
        l1[i] %>% basename %>% strsplit("_") %>% unlist %>% "["(1) -> GrepDirName0
        l1[i] %>% input(.) -> ll0
        ll0 %>% dplyr::filter(.data$rankpval < top | .data$rankena < top) -> ll1
        ll1 %>% head
        ll1 %>% dim
        l1[i] %>% basename %>% sub("(.*)_.*.tsv$","\\1",.) %>% paste("top",top,"_",.,"_",nrow(ll1),".tsv",sep="") -> FileName0
        resultpath %>% file.path("Final-network-all",GrepDirName0,FileName0) -> FileName1
        ll1 %>% output(FileName1)
        foreach(j=1:length(conf0)) %do%
          {
            # conf0[j] %>% basename %>% strsplit("_") %>% unlist %>% "["(7) %>% gsub("c","conf",.) -> ConfidenceList0
            conf0[j] %>% basename %>% strsplit("_") %>% unlist %>% "["(2) -> ConfidenceList0
            conf0[j] %>% list.files(full.names=T) %>% grep("Nodes",.,value=T) %>% input(.) -> Nodes0
            Nodes0 %>% head
            Nodes0 %>% dim
            Nodes0 %>% dplyr::select(.data$ENSPID,.data$Symbol) -> Nodes1
            conf0[j]  %>% list.files(full.names=T) %>% grep("Edges",.,value=T) %>% input(.) -> Edges0
            Edges0 %>% head
            Edges0 %>% dim
            Edges0 %>% dplyr::select(.data$NodeA,.data$NodeB) -> Edges1
            pedges0 <- foreach(k=1:nrow(ll1),.combine="rbind") %do% { data.frame(ll1$Name[k],ll1$SymbolList[k] %>% strsplit("\\|") %>% unlist) }
            pedges0 %>% setNames(c("NodeA","NodeB")) -> pedges0
            pedges0 %>% head
            # add ENSPID
            Edges0 %>% dplyr::select(.data$NodeA,.data$SymbolA) %>% rbind(data.frame("NodeA"=Edges0$NodeB,"SymbolA"=Edges0$SymbolB)) %>% 
              setNames(c("ENSPID","Symbol")) %>% unique -> t
            pedges0 %>% setNames(c("NodeA","Symbol")) %>% inner_join(t) -> pedges1
            pedges1 %>% head
            pedges1 %>% dplyr::select(.data$NodeA,.data$ENSPID) %>% setNames(c("NodeA","NodeB")) -> pedges2
            pedges2 %>% head
            rbind(Edges1,pedges2) %>% unique -> edges
            # geneset edges color 
            rep("gray",nrow(edges)) -> colorEdges
            edges$NodeA %>% grep("ENSP",.,value=F,invert=T) %>% replace(colorEdges,.,"orange") -> colorEdges
            data.frame(edges,colorEdges) -> edges
            #
            # Nodes
            data.frame(ENSPID=pedges1$NodeA,Symbol=pedges1$NodeA) %>% unique %>%
              setNames(c("ENSPID","Symbol")) %>% rbind(Nodes1) %>% unique -> Nodes2
            # Geneset ColorNode
            rep("orange",length(pedges1$NodeA %>% unique)) -> FunctionColor0
            rep("gray",nrow(Nodes2)-length(pedges1$NodeA %>% unique)) -> SymbolColor0
            c(FunctionColor0,SymbolColor0) -> colorNode
            Nodes2 %>% data.frame(colorNode=colorNode) -> Nodes3
            # ColorName
            rep("orangered4",length(pedges1$NodeA %>% unique)) -> FunctionColorName0
            rep("blue",nrow(Nodes2)-length(pedges1$NodeA %>% unique)) -> SymbolColorName0
            c(FunctionColorName0,SymbolColorName0) -> colorName
            Nodes3 %>% data.frame(colorName=colorName) -> Nodes4
            # vertex.label.cex
            rep(0.7,length(pedges1$NodeA %>% unique)) -> FunctionLabelCex0
            rep(0.4,nrow(Nodes2)-length(pedges1$NodeA %>% unique)) -> SymbolLabelCex0
            c(FunctionLabelCex0,SymbolLabelCex0) -> LabelCex
            Nodes4 %>% data.frame(LabelCex=LabelCex) -> Nodes5
            #
            igraph::graph_from_data_frame(d=edges, vertices=Nodes4, directed=F) -> NetWork0
            #
            foreach(p=2:ncol(FoldChange0)) %do%
              {
                moal:::networkcolor(foldchange=FoldChange0[,c(1,p)], colornode=Nodes4$colorNode, networkobj=NetWork0 ) -> NetworkColor0
                NetworkColor0[[1]] -> Nodes4$colorNode
                Nodes4 -> nodes
                igraph::graph_from_data_frame(d=edges, vertices=nodes, directed=F) -> NetWork0
                #
                # display
                igraph::degree(NetWork0, mode="all") -> Deg0
                list( c("layout_as_tree","tree"),c("layout_in_circle","circle"),c("layout_on_grid","grid"),
                      c("layout_on_sphere","sphere"),c("layout_with_fr","fr"),c("layout_with_dh","dh") ) -> NetworkLayout0
                if( nrow(edges) < intmaxdh ){ NetworkLayout0 -> NetworkLayout1 }else{ NetworkLayout0[-6] -> NetworkLayout1 }
                # igraph::layout_with_dh(NetWork0) -> LayOut
                paste("nw_",ConfidenceList0,sep="") -> DirNameNw0
                resultpath %>% file.path("Final-network-all",GrepDirName0,DirNameNw0) %>% dir.create
                # output edges
                edges$NodeA %>% annot(idtype="ENSP",dboutput="EBI") %>% dplyr::select(.data$ENSPID,.data$Symbol) %>% setNames(c("NodeA","SymbolA")) -> t
                edges$NodeB %>% annot(idtype="ENSP",dboutput="EBI") %>% dplyr::select(.data$ENSPID,.data$Symbol) %>% setNames(c("NodeB","SymbolB"))  -> tt
                cbind(t,tt) -> edgeso
                paste("Edges_",nrow(edgeso),".tsv",sep="") -> FileName0
                resultpath %>% file.path("Final-network-all",GrepDirName0,DirNameNw0,FileName0) -> FileName1
                edgeso %>% output(FileName1)
                # output nodes
                Nodes5$ENSPID %>% data.frame(ENSPID=.) %>% dplyr::left_join(Nodes0) -> nodeso
                paste("Nodes_",nrow(nodeso),".tsv",sep="") -> FileName0
                resultpath %>% file.path("Final-network-all",GrepDirName0,DirNameNw0,FileName0) -> FileName1
                nodeso %>% output(FileName1)
                #
                foreach(l=1:length(NetworkLayout1) ) %do% {
                  set.seed(1234567)
                  eval( parse( text = paste(NetworkLayout1[[l]][1],"(NetWork0)",sep ="") ) ) -> LayOut
                  NetworkLayout1[[l]][1] %>% strsplit("_") %>% unlist %>% "["(3) -> LayOutName
                  paste(ConfidenceList0,"_",LayOutName,"_",colnames(FoldChange0)[p],"_",nrow(edges),"_",nrow(nodes),".pdf", sep="") -> FileName0
                  resultpath %>% file.path("Final-network-all",GrepDirName0,DirNameNw0,FileName0) -> FileName1
                  FileName0 %>% paste(GrepDirName0,"_",.,sep="") %>% sub(".pdf","",.) -> Title0
                  pdf(FileName1, width=20, height=20 )
                  plot.igraph(
                    NetWork0, main=Title0, cex.main=2,
                    edge.width=0.2, edge.weight=0.1, edge.color=E(NetWork0)$colorEdges,
                    layout=LayOut, edge.curved=0.1, arrow.mode=3, edge.arrow.size=0.1, edge.arrow.width=0.5,
                    vertex.label=V(NetWork0)$Symbol %>% as.character, vertex.size=(log2(Deg0)+4)*0.31, vertex.label.cex=V(NetWork0)$LabelCex,
                    vertex.frame.color=V(NetWork0)$colorNode, vertex.label.font=2, vertex.label.color=V(NetWork0)$colorName, vertex.color=V(NetWork0)$colorNode )
                  graphics.off()
                }
              }
          }
      }
    parallel::stopCluster(cl)
    #
    # summary-network split
    #
    resultpath %>% file.path("Final-network-MSigDB") %>% list.files(recursive = T,full.names =T) -> l0
    l0
    l0 %>% basename %>% grep("^top.*_[0-9]{1,2}-(.*)_results_.*.tsv",.) %>% l0[.] -> l1
    l1 %>% sub(".*_(.*).tsv","\\1",.) %>% as.numeric %>% ">"(0) %>% which %>% l1[.] -> l2
    l2 %>% basename %>% sub("^top.*_[0-9]{1,2}-(.*)_results_.*.tsv","\\1",.) -> KeywordName0
    l2 %>% basename %>% sub("^top.*_([0-9]{1,2}-.*)_results_.*.tsv","\\1",.) -> KeywordDirNameName0
    resultpath %>% file.path("Final-network-all") %>% list.files(recursive = T,full.names =T) -> l0
    l0
    l0 %>% basename %>% grep("^top.*_[0-9]{1,2}-(.*)_results_.*.tsv",.) %>% l0[.] -> l1
    l1 %>% sub(".*_(.*).tsv","\\1",.) %>% as.numeric %>% ">"(0) %>% which %>% l1[.] -> l2
    l2 %>% basename %>% sub("^top.*_[0-9]{1,2}-(.*)_results_.*.tsv","\\1",.) -> KeywordName0
    l2 %>% basename %>% sub("^top.*_([0-9]{1,2}-.*)_results_.*.tsv","\\1",.) -> KeywordDirNameName0
    # gene set color network
    moal:::palette0 -> Palette0
    Palette0[-6] -> Palette0
    foreach(i=1:length(KeywordName0)) %do%
      {
        # rep(c("orange","tan","sienna1"),(length(KeywordName0)/2+1)) -> Palette0
        resultpath %>% file.path("Final-network-all",KeywordDirNameName0[i],"Summary-network") %>% dir.create
        l2[i]
        l2[i] %>% input(.) -> ll0 
        ll0$SymbolList %>% strsplit("\\|") %>% unlist %>% unique -> NodeB 
        data.frame(NodeA=KeywordName0[i],NodeB) -> pedges0
        pedges0$NodeA %>% table
        #
        resultpath %>% list.files(full.names=T,recursive=F) %>% grep("Networks",.,value = T) %>%
          list.files(full.names=T,recursive=F) %>% grep("nw_",.,value = T) -> conf0
        #
        foreach(j=1:length(conf0)) %do%
          {
            conf0[j] %>% basename %>% strsplit("_") %>% unlist %>% "["(2) -> ConfidenceList0
            conf0[j] %>% list.files(full.names=T) %>% grep("Nodes",.,value=T) %>% input(.) -> Nodes0
            Nodes0 %>% head
            Nodes0 %>% dim
            Nodes0 %>% dplyr::select(.data$ENSPID,.data$Symbol) -> Nodes1
            conf0[j]  %>% list.files(full.names=T) %>% grep("Edges",.,value=T) %>% input(.) -> Edges0
            Edges0 %>% head
            Edges0 %>% dim
            Edges0 %>% dplyr::select(.data$NodeA,.data$NodeB) -> Edges1
            # add ENSPID
            Edges0 %>% dplyr::select(.data$NodeA,.data$SymbolA) %>% rbind(data.frame("NodeA"=Edges0$NodeB,"SymbolA"=Edges0$SymbolB)) %>% 
              setNames(c("ENSPID","Symbol")) %>% unique -> t
            pedges0 %>% setNames(c("NodeA","Symbol")) %>% inner_join(t) -> pedges1
            pedges1 %>% head
            pedges1 %>% dplyr::select(.data$NodeA,.data$ENSPID) %>% setNames(c("NodeA","NodeB")) -> pedges2
            pedges2 %>% head
            rbind(Edges1,pedges2) %>% unique -> edges
            # geneset edges color 
            rep("gray",nrow(edges)) -> colorEdges
            colorEdges %>% table
            edges$NodeA %>% grep("ENSP",.,value=T,invert=T) %>% table
            # edges$NodeA %>% grep("ENSP",.,value=T,invert=T) %>% table %>% names -> KeywordName0
            edges$NodeA %>% grep("ENSP",.,value=T,invert=T) %>% table %>% names -> KeywordName1
            foreach(k=1:length(KeywordName0[i])) %do% {
              edges$NodeA %>% "=="(KeywordName1[k]) %>% which -> sel
              colorEdges[sel] <- Palette0[k] }
            data.frame(edges,colorEdges) -> edges
            edges$colorEdges %>% table
            # Nodes
            data.frame(ENSPID=pedges1$NodeA,Symbol=pedges1$NodeA) %>% unique %>%
              setNames(c("ENSPID","Symbol")) %>% rbind(Nodes1) %>% unique -> Nodes2
            # Geneset ColorNode
            rep("gray",nrow(Nodes2)) -> colorNode0
            foreach(k=1:length(KeywordName0[i])) %do% {
              Nodes2$Symbol %>% "=="(KeywordName1[k]) %>% which -> sel
              colorNode0[sel] <- Palette0[k] }
            colorNode0 %>% table
            Nodes2 %>% data.frame(colorNode=colorNode0) -> Nodes3
            # ColorName
            rep("blue",nrow(Nodes2)) -> colorName0
            foreach(k=1:length(KeywordName0[i])) %do% {
              Nodes2$Symbol %>% "=="(KeywordName1[k]) %>% which -> sel
              colorName0[sel] <- "black" }
            colorName0 %>% table
            # rep("orangered4",length(pedges1$NodeA %>% unique)) -> FunctionColorName0
            # rep("blue",nrow(Nodes2)-length(pedges1$NodeA %>% unique)) -> SymbolColorName0
            # c(FunctionColorName0,SymbolColorName0) -> colorName
            # Nodes3 %>% data.frame(colorName=colorName) -> Nodes4
            Nodes3 %>% data.frame(colorName=colorName0) -> Nodes4
            # vertex.label.cex
            rep(0.7,nrow(Nodes2)) -> LabelCex0
            foreach(k=1:length(KeywordName0[i])) %do% { 
              Nodes2$Symbol %>% "=="(KeywordName1[k]) %>% which -> sel
              LabelCex0[sel] <- 0.9 }
            LabelCex0 %>% table
            # rep(0.7,length(pedges1$NodeA %>% unique)) -> FunctionLabelCex0
            # rep(0.4,nrow(Nodes2)-length(pedges1$NodeA %>% unique)) -> SymbolLabelCex0
            # c(FunctionLabelCex0,SymbolLabelCex0) -> LabelCex
            # Nodes4 %>% data.frame(LabelCex=LabelCex) -> Nodes5
            Nodes4 %>% data.frame(LabelCex=LabelCex0) -> Nodes5
            # add Connectivity geneset
            edges$NodeA %>% grep("ENSP",.,value=T,invert=T) %>% table -> GenesetConnectivity0
            foreach(m=1:length(GenesetConnectivity0)) %do%
              {
                Nodes0$ENSPID %>% "=="(GenesetConnectivity0 %>% names %>% "["(m) ) %>% which -> sel
                if(length(sel) > 0){GenesetConnectivity0[m] %>% as.numeric -> Nodes5$Connectivity[sel]}
              }
            #
            igraph::graph_from_data_frame(d=edges, vertices=Nodes5, directed=F) -> NetWork0
            # Node Color fold-change
            foreach(p=2:ncol(FoldChange0)) %do%
              {
                moal:::networkcolor(foldchange=FoldChange0[,c(1,p)], colornode=Nodes5$colorNode, networkobj=NetWork0 ) -> NetworkColor0
                NetworkColor0[[1]] -> Nodes5$colorNode
                Nodes5 -> nodes
                igraph::graph_from_data_frame(d=edges, vertices=nodes, directed=F) -> NetWork0
                # display
                igraph::degree(NetWork0, mode="all") -> Deg0
                list( c("layout_as_tree","tree"),c("layout_in_circle","circle"),c("layout_on_grid","grid"),
                      c("layout_on_sphere","sphere"),c("layout_with_fr","fr"),c("layout_with_dh","dh") ) -> NetworkLayout0
                if( nrow(edges) < intmaxdh ){ NetworkLayout0 -> NetworkLayout1 }else{ NetworkLayout0[-6] -> NetworkLayout1 }
                # igraph::layout_with_dh(NetWork0) -> LayOut
                resultpath %>% file.path("Final-network-all",KeywordDirNameName0[i],"Summary-network",ConfidenceList0) %>% dir.create
                # output edges
                edges$NodeA %>% annot(idtype="ENSP",dboutput="EBI") %>% dplyr::select(.data$ENSPID,.data$Symbol) %>% setNames(c("NodeA","SymbolA")) -> t
                edges$NodeB %>% annot(idtype="ENSP",dboutput="EBI") %>% dplyr::select(.data$ENSPID,.data$Symbol) %>% setNames(c("NodeB","SymbolB"))  -> tt
                cbind(t,tt) -> edgeso
                edgeso %>% dim
                edgeso %>% data.frame(colorEdges=edges$colorEdges) -> edgesoo
                paste("Edges_",nrow(edgesoo),".tsv",sep="") -> FileName0
                resultpath %>% file.path("Final-network-all",KeywordDirNameName0[i],"Summary-network",ConfidenceList0,FileName0) -> FileName1
                edgesoo %>% output(FileName1)
                # output nodes
                paste("Nodes_",nrow(nodes),".tsv",sep="") -> FileName0
                resultpath %>% file.path("Final-network-all",KeywordDirNameName0[i],"Summary-network",ConfidenceList0,FileName0) -> FileName1
                nodes %>% output(FileName1)
                #
                foreach(l=1:length(NetworkLayout1) ) %do% {
                  set.seed(1234567)
                  eval( parse( text = paste(NetworkLayout1[[l]][1],"(NetWork0)",sep ="") ) ) -> LayOut
                  NetworkLayout1[[l]][1] %>% strsplit("_") %>% unlist %>% "["(3) -> LayOutName
                  paste(ConfidenceList0,"_",LayOutName,"_",colnames(FoldChange0)[p],"_",nrow(edges),"_",nrow(nodes),".pdf", sep="") -> FileName0
                  resultpath %>% file.path("Final-network-all",KeywordDirNameName0[i],"Summary-network",ConfidenceList0,FileName0) -> FileName1
                  FileName0 %>% paste("summary_",.,sep="") %>% sub(".pdf","",.) -> Title0
                  pdf(FileName1, width=20, height=20 )
                  plot.igraph(
                    NetWork0, main=Title0, cex.main=2,
                    edge.width=0.2, edge.weight=0.1, edge.color=E(NetWork0)$colorEdges,
                    layout=LayOut, edge.curved=0.1, arrow.mode=3, edge.arrow.size=0.1, edge.arrow.width=0.5,
                    vertex.label=V(NetWork0)$Symbol %>% as.character, vertex.size=(log2(Deg0)+4)*0.31, vertex.label.cex=V(NetWork0)$LabelCex,
                    vertex.frame.color=V(NetWork0)$colorNode, vertex.label.font=2, vertex.label.color=V(NetWork0)$colorName, vertex.color=V(NetWork0)$colorNode )
                  graphics.off()
                }
              } 
          }
      }
    resultpath %>% file.path("Final-network-all") %>% list.files(full.names = T) %>% grep("results",.,value = T) -> l0
    l0 %>% basename %>% sub(".*_(.*)\\.tsv","\\1",.) %>% as.numeric %>% ">"(0) %>% which %>% l0[.] -> l1
    l1 %>% file.remove
  }
  #
  # without keywords
  #
  # MSigDB
  resultpath %>% file.path("Final-network-MSigDB") %>% list.files(full.names=T,recursive=T) -> l0
  paste("top",top,"_MSigDB_[0-9]*.tsv",sep="") -> Grep0
  l0 %>% basename %>% grep(Grep0,.) %>% l0[.] -> l1
  l1 %>% head
  l1 %>% input(.) -> lpemsigdb
  lpemsigdb %>% dplyr::select(-.data$GenesetSize) -> lpemsigdb
  lpemsigdb %>% head
  # IPA
  resultpath %>% file.path("ipa") %>% list.files(full.names=T,recursive=T) -> l0
  paste("top",top,"_.*.tsv",sep="") -> Grep0
  l0 %>% basename %>% grep(Grep0,.) %>% l0[.] -> lipa
  lipa %>% head
  lpeipa <- foreach(j=1:length(lipa), .combine = "rbind") %do%
    {
      lipa[j] %>% input(.) -> ll0
      lipa[j] %>% basename %>% strsplit("_") %>% unlist %>% "["(2) %>% paste("IPA",.,sep="") -> db
      # ll0 %>% dplyr::slice(1:top) -> ll1
      ll0 %>% dplyr::select(.data$Name,.data$SymbolList,.data$OverlapSize,.data$ENAScore,.data$log10pvalFDR) %>% data.frame(Collection=db) -> t
      t %>% dplyr::arrange(-.data$log10pvalFDR) %>% data.frame(rankpval=1:nrow(t)) -> tt
      tt %>% dplyr::arrange(-.data$ENAScore) %>% data.frame(rankena=1:nrow(t)) -> ttt
      ttt %>% dplyr::arrange(-.data$log10pvalFDR) %>% dplyr::slice(1:top) -> tp
      ttt %>% dplyr::arrange(-.data$ENAScore) %>% dplyr::slice(1:top) -> te
      rbind(tp,te) %>% unique
    }
  # lpeipa <- foreach(i=1:length(l1), .combine = "rbind") %do%
  #   {
  #     l1[i] %>% input(.) -> ll0
  #     ll0 %>% head
  #     l1[i] %>% basename %>% strsplit("_") %>% unlist %>% "["(2) %>% paste("IPA",.,sep="") -> db
  #     ll0 %>% dplyr::select(.data$Name,.data$SymbolList,.data$OverlapSize,.data$ENAScore,.data$log10pvalFDR) %>% data.frame(Collection=db)
  #   }
  lpeipa %>% head
  lpeipa %>% colnames
  # topGO
  resultpath %>% file.path("topGO") %>% list.files(full.names=T,recursive=T) -> l0
  paste("top",top,"_.*.tsv",sep="") -> Grep0
  l0 %>% basename %>% grep(Grep0,.) %>% l0[.] -> ltopgo
  l1 %>% head
  lpetopgo <- foreach(j=1:length(ltopgo), .combine = "rbind") %do%
    {
      ltopgo[j] %>% input(.) -> ll0
      ltopgo[j] %>% basename %>% strsplit("_") %>% unlist %>% "["(2) %>% paste("topGO",.,sep="") -> db
      # ll0 %>% dplyr::slice(1:top) -> ll1
      ll0 %>% dplyr::select(.data$Name,.data$SymbolList,.data$OverlapSize,.data$ENAScore,.data$log10pvalFDR) %>% data.frame(Collection=db) -> t
      t %>% dplyr::arrange(-.data$log10pvalFDR) %>% data.frame(rankpval=1:nrow(t)) -> tt
      tt %>% dplyr::arrange(-.data$ENAScore) %>% data.frame(rankena=1:nrow(t)) -> ttt
      ttt %>% dplyr::arrange(-.data$log10pvalFDR) %>% dplyr::slice(1:top) -> tp
      ttt %>% dplyr::arrange(-.data$ENAScore) %>% dplyr::slice(1:top) -> te
      rbind(tp,te) %>% unique
    }
  lpetopgo %>% head
  lpetopgo %>% colnames
  lpemsigdb %>% colnames
  lpeipa %>% colnames
  lpetopgo %>% colnames
  rbind(lpemsigdb,lpeipa,lpetopgo) %>% unique -> lpeAll
  paste("top",top,"_all_genesets_",nrow(lpeAll),".tsv",sep="") -> FileName0
  resultpath %>% file.path("Final-network-all",FileName0) -> FileName1
  lpeAll %>% output(FileName1)
  #
  # top categories
  resultpath %>% list.files(full.names=T,recursive=T) -> l0
  l0 %>% basename %>% grep("^c[0-9]_.*.tsv$",.) %>% l0[.] -> l1
  c(l1,lipa,ltopgo) -> l1
  l1 %>% length
  c("Localization","Functions","Pathways","Patterns","Regulators","Diseases","Molecules") -> CollCat0
  foreach(q=1:length(CollCat0)) %do%
    {
      l1 %>% grep(CollCat0[q],.) %>% l1[.] -> l2
      lpemsigdbcollcat <- foreach(j=1:length(l2),.combine = "rbind") %do%
        {
          l2[j] %>% input(.) -> ll0
          if(l2[j] %>% grepl("topGO",.)){ l2[j] %>% basename %>% strsplit("_") %>% unlist %>% "["(2) %>% paste("topGO",.,sep="")-> db }
          if(l2[j] %>% grepl("ipa",.)){ l2[j] %>% basename %>% strsplit("_") %>% unlist %>% "["(2) %>% paste("IPA",.,sep="")-> db }
          if(l2[j] %>% grepl("ipa|topGO",.) %>% "!"(.)){ l2[j] %>% basename %>% strsplit("_") %>% unlist %>% "["(3) %>% paste("MSigDB",.,sep="")-> db }
          # l2[j] %>% basename %>% strsplit("_") %>% unlist %>% "["(3)  -> db
          # ll0 %>% dplyr::slice(1:top) -> ll1
          ll0 %>% dplyr::select(.data$Name,.data$SymbolList,.data$OverlapSize,.data$ENAScore,.data$log10pvalFDR) %>% data.frame(Collection=db) -> t
          t %>% dplyr::arrange(-.data$log10pvalFDR) %>% data.frame(rankpval=1:nrow(t)) -> tt
          tt %>% dplyr::arrange(-.data$ENAScore) %>% data.frame(rankena=1:nrow(t)) -> ttt
          ttt %>% dplyr::arrange(-.data$log10pvalFDR) %>% dplyr::slice(1:top) -> tp
          ttt %>% dplyr::arrange(-.data$ENAScore) %>% dplyr::slice(1:top) -> te
          rbind(tp,te) %>% unique
        }
      lpemsigdbcollcat %>% unique %>% dplyr::arrange(-.data$log10pvalFDR) -> lpemsigdbcollcat1
      paste("top",top,"_all_",CollCat0[q],"_",nrow(lpemsigdbcollcat1),".tsv",sep="") -> FileName0
      resultpath %>% file.path("Final-network-all",FileName0) -> FileName1
      lpemsigdbcollcat1 %>% output(FileName1)
    }
  #
  # Occurences
  #
  # l1[i] %>% file.remove
  lpeAll$Collection %>% table %>% names -> Collection0
  # Symbol occurrence for significant enrich canonical pathways 
  f2 <- foreach(j=1:length(Collection0),.combine = "rbind") %do%
    {
      lpeAll %>% dplyr::filter(.data$Collection == Collection0[j]) -> t0
      t0$SymbolList %>% strsplit("\\|") %>% unlist -> t1
      if( table(t1) %>% length == 1 ){ data.frame(Symbol=t1,Occ=1) -> t2 }
      if( table(t1) %>% length > 1 ){ t1 %>% table %>% sort(decreasing = T) %>% data.frame %>% setNames(c("Symbol","Occ")) -> t2 }
      t2 %>% data.frame("Collection"=rep(Collection0[j],nrow(t2)))
    }
  f2$Symbol %>% annot(.) -> a0
  a0 %>% dplyr::select(-.data$Symbol) %>% data.frame(f2,.) -> f3
  paste("Symbol-Occurences_",nrow(f3),".tsv",sep="") -> FileName0
  resultpath %>% file.path("Final-network-all",FileName0) -> FileName1
  f3 %>% dplyr::arrange(-(.data$Occ)) -> f4
  f4 %>% output(FileName1)
  # Occurence Matrix
  symbollist %>% unique %>% data.frame(Symbol=.) -> m0
  m1 <- foreach(j=1:length(Collection0),.combine="cbind") %do%
    {
      f4 %>% dplyr::filter(.data$Collection == Collection0[j]) %>% dplyr::select(.data$Symbol,.data$Occ) %>% unique -> t0
      m0 %>% dplyr::left_join(t0) %>% dplyr::select(.data$Occ)
    }
  m1 %>% head
  m2 <- foreach(j=1:ncol(m1),.combine="cbind") %do% { m1[j] %>% unlist %>% is.na %>% which %>% replace(m1[j] %>% unlist,.,0) }
  m2 %>% head
  m2 %>% dim
  m2 %>% data.frame(m0$Symbol,.) %>% setNames(c("Symbol",Collection0)) -> m3
  m3$Symbol %>% annot(.) -> a0
  m3 %>% data.frame(a0[,c(-1)]) %>% dplyr::arrange(desc(.[[2]])) -> m4
  m4 %>% head
  paste("Symbol-Occurences","_matrix_",nrow(m4),".tsv",sep="") -> FileName0
  resultpath %>% file.path("Final-network-all",FileName0) -> FileName1
  m4 %>% output(FileName1)
  # Occurence Matrix hierarchical clustering
  m4 %>% dplyr::select(c(2:(length(Collection0)+1)))  -> m5
  paste("Symbol-Occurences","_matrix_HC_",nrow(m5),".pdf",sep="") -> FileName0
  resultpath %>% file.path("Final-network-all",FileName0) -> FileName1
  if(ncol(m5)>2)
  {
    if(m5 %>% apply(2,sd) %>% "=="(0) %>% which %>% any )
    {
      m5 %>% apply(2,sd) %>% "=="(0) %>% which -> selsd
      m5[,-selsd] -> m5
    }
    pdf(FileName1)
    par(mar=c(6.1,4.1,4.1,2.1))
    m5 %>% moal:::hc(m5 %>% colnames %>% as.factor ,legendtitle="Collection",legend=F,cexlabel=0.9)
    graphics.off()
  } 
  #
  # summary-network all keyword
  #
  if(!is.null(keywords))
  {
    
    resultpath %>% file.path("Final-network-all") %>% list.files(recursive = T,full.names =T) -> l0
    l0
    l0 %>% basename %>% grep("^top.*_[0-9]{1,2}-(.*)_results_.*.tsv",.) %>% l0[.] -> l1
    l1 %>% sub(".*_(.*).tsv","\\1",.) %>% as.numeric %>% ">"(0) %>% which %>% l1[.] -> l2
    l2 %>% basename %>% sub("^top.*_[0-9]{1,2}-(.*)_results_.*.tsv","\\1",.) -> KeywordName0
    # gene set color network
    Palette0[-6] -> Palette0
    # rep(c("orange","tan","sienna1"),(length(KeywordName0)/2+1)) -> Palette0
    resultpath %>% file.path("Final-network-all","Summary-network") %>% dir.create
    pedges0 <- foreach(i=1:length(l2),.combine = "rbind") %do%
      {
        l2[i] %>% input(.) -> ll0 
        ll0$SymbolList %>% strsplit("\\|") %>% unlist %>% unique -> NodeB 
        data.frame(NodeA=KeywordName0[i],NodeB)
      }
    pedges0$NodeA %>% table
    #
    resultpath %>% list.files(full.names=T,recursive=F) %>% grep("Networks",.,value = T) %>%
      list.files(full.names=T,recursive=F) %>% grep("nw_",.,value = T) -> conf0
    #
    foreach(j=1:length(conf0)) %do%
      {
        # conf0[j] %>% basename %>% strsplit("_") %>% unlist %>% "["(7) %>% gsub("c","conf",.) -> ConfidenceList0
        conf0[j] %>% basename %>% strsplit("_") %>% unlist %>% "["(2) -> ConfidenceList0
        conf0[j] %>% list.files(full.names=T) %>% grep("Nodes",.,value=T) %>% input(.) -> Nodes0
        Nodes0 %>% head
        Nodes0 %>% dim
        Nodes0 %>% dplyr::select(.data$ENSPID,.data$Symbol) -> Nodes1
        conf0[j]  %>% list.files(full.names=T) %>% grep("Edges",.,value=T) %>% input(.) -> Edges0
        Edges0 %>% head
        Edges0 %>% dim
        Edges0 %>% dplyr::select(.data$NodeA,.data$NodeB) -> Edges1
        # pedges0 <- foreach(k=1:nrow(ll1), .combine="rbind") %do% {
        #   data.frame(ll1$Name[k],ll1$SymbolList[k] %>% strsplit("\\|") %>% unlist) }
        # pedges0 %>% setNames(c("NodeA","NodeB")) -> pedges0
        # pedges0 %>% head
        # add ENSPID
        Edges0 %>% dplyr::select(.data$NodeA,.data$SymbolA) %>% rbind(data.frame("NodeA"=Edges0$NodeB,"SymbolA"=Edges0$SymbolB)) %>% 
          setNames(c("ENSPID","Symbol")) %>% unique -> t
        pedges0 %>% setNames(c("NodeA","Symbol")) %>% inner_join(t) -> pedges1
        pedges1 %>% head
        pedges1 %>% dplyr::select(.data$NodeA,.data$ENSPID) %>% setNames(c("NodeA","NodeB")) -> pedges2
        pedges2 %>% head
        rbind(Edges1,pedges2) %>% unique -> edges
        # geneset edges color 
        rep("gray",nrow(edges)) -> colorEdges
        colorEdges %>% table
        edges$NodeA %>% grep("ENSP",.,value=T,invert=T) %>% table
        edges$NodeA %>% grep("ENSP",.,value=T,invert=T) %>% table %>% names -> KeywordName0
        foreach(k=1:length(KeywordName0)) %do% {
          edges$NodeA %>% "=="(KeywordName0[k]) %>% which -> sel
          colorEdges[sel] <- Palette0[k] }
        # foreach(k=1:length(KeywordName0)) %do% { 
        #   edges$NodeA %>% "=="(KeywordName0[k]) %>% which -> sel 
        #   colorEdges[sel] <- "orange" }
        data.frame(edges,colorEdges) -> edges
        edges$colorEdges %>% table
        #
        # Nodes
        data.frame(ENSPID=pedges1$NodeA,Symbol=pedges1$NodeA) %>% unique %>%
          setNames(c("ENSPID","Symbol")) %>% rbind(Nodes1) %>% unique -> Nodes2
        # Geneset ColorNode
        rep("gray",nrow(Nodes2)) -> colorNode0
        foreach(k=1:length(KeywordName0)) %do% {
          Nodes2$Symbol %>% "=="(KeywordName0[k]) %>% which -> sel
          colorNode0[sel] <- Palette0[k] }
        # foreach(k=1:length(KeywordName0)) %do% { 
        #   Nodes2$Symbol %>% "=="(KeywordName0[k]) %>% which -> sel
        #   colorNode0[sel] <- "orange" }
        colorNode0 %>% table
        # rep("orange",length(pedges1$NodeA %>% unique)) -> FunctionColor0
        # rep("gray",nrow(Nodes2)-length(pedges1$NodeA %>% unique)) -> SymbolColor0
        # c(FunctionColor0,SymbolColor0) -> colorNode
        # Nodes2 %>% data.frame(colorNode=colorNode) -> Nodes3
        Nodes2 %>% data.frame(colorNode=colorNode0) -> Nodes3
        # ColorName
        rep("blue",nrow(Nodes2)) -> colorName0
        # foreach(k=1:length(KeywordName0)) %do% {
        #   Nodes2$Symbol %>% "=="(KeywordName0[k]) %>% which -> sel
        #   colorName0[sel] <- Palette0[k] }
        foreach(k=1:length(KeywordName0)) %do% {
          Nodes2$Symbol %>% "=="(KeywordName0[k]) %>% which -> sel
          colorName0[sel] <- "black" }
        colorName0 %>% table
        # rep("orangered4",length(pedges1$NodeA %>% unique)) -> FunctionColorName0
        # rep("blue",nrow(Nodes2)-length(pedges1$NodeA %>% unique)) -> SymbolColorName0
        # c(FunctionColorName0,SymbolColorName0) -> colorName
        # Nodes3 %>% data.frame(colorName=colorName) -> Nodes4
        Nodes3 %>% data.frame(colorName=colorName0) -> Nodes4
        # vertex.label.cex
        rep(0.7,nrow(Nodes2)) -> LabelCex0
        foreach(k=1:length(KeywordName0)) %do% { 
          Nodes2$Symbol %>% "=="(KeywordName0[k]) %>% which -> sel
          LabelCex0[sel] <- 0.9 }
        LabelCex0 %>% table
        # rep(0.7,length(pedges1$NodeA %>% unique)) -> FunctionLabelCex0
        # rep(0.4,nrow(Nodes2)-length(pedges1$NodeA %>% unique)) -> SymbolLabelCex0
        # c(FunctionLabelCex0,SymbolLabelCex0) -> LabelCex
        # Nodes4 %>% data.frame(LabelCex=LabelCex) -> Nodes5
        Nodes4 %>% data.frame(LabelCex=LabelCex0) -> Nodes5
        # add Connectivity geneset
        edges$NodeA %>% grep("ENSP",.,value=T,invert=T) %>% table -> GenesetConnectivity0
        foreach(m=1:length(GenesetConnectivity0)) %do%
          {
            Nodes0$ENSPID %>% "=="(GenesetConnectivity0 %>% names %>% "["(m) ) %>% which -> sel
            if(length(sel) > 0){GenesetConnectivity0[m] %>% as.numeric -> Nodes5$Connectivity[sel]}
          }
        #
        igraph::graph_from_data_frame(d=edges, vertices=Nodes5, directed=F) -> NetWork0
        # Node Color fold-change
        foreach(p=2:ncol(FoldChange0)) %do%
          {
            moal:::networkcolor(foldchange=FoldChange0[,c(1,p)], colornode=Nodes5$colorNode, networkobj=NetWork0 ) -> NetworkColor0
            NetworkColor0[[1]] -> Nodes5$colorNode
            Nodes5 -> nodes
            igraph::graph_from_data_frame(d=edges, vertices=nodes, directed=F) -> NetWork0
            # display
            igraph::degree(NetWork0, mode="all") -> Deg0
            list( c("layout_as_tree","tree"),c("layout_in_circle","circle"),c("layout_on_grid","grid"),
                  c("layout_on_sphere","sphere"),c("layout_with_fr","fr"),c("layout_with_dh","dh") ) -> NetworkLayout0
            if( nrow(edges) < intmaxdh ){ NetworkLayout0 -> NetworkLayout1 }else{ NetworkLayout0[-6] -> NetworkLayout1 }
            # igraph::layout_with_dh(NetWork0) -> LayOut
            resultpath %>% file.path("Final-network-all","Summary-network",ConfidenceList0) %>% dir.create
            # output edges
            edges$NodeA %>% annot(idtype="ENSP",dboutput="EBI") %>% dplyr::select(.data$ENSPID,.data$Symbol) %>% setNames(c("NodeA","SymbolA")) -> t
            edges$NodeB %>% annot(idtype="ENSP",dboutput="EBI") %>% dplyr::select(.data$ENSPID,.data$Symbol) %>% setNames(c("NodeB","SymbolB"))  -> tt
            cbind(t,tt) -> edgeso
            edgeso %>% dim
            edgeso %>% data.frame(colorEdges=edges$colorEdges) -> edgesoo
            paste("Edges_",nrow(edgesoo),".tsv",sep="") -> FileName0
            resultpath %>% file.path("Final-network-all","Summary-network",ConfidenceList0,FileName0) -> FileName1
            edgesoo %>% output(FileName1)
            # output nodes
            paste("Nodes_",nrow(nodes),".tsv",sep="") -> FileName0
            resultpath %>% file.path("Final-network-all","Summary-network",ConfidenceList0,FileName0) -> FileName1
            nodes %>% output(FileName1)
            #
            foreach(l=1:length(NetworkLayout1) ) %do% {
              set.seed(1234567)
              eval( parse( text = paste(NetworkLayout1[[l]][1],"(NetWork0)",sep ="") ) ) -> LayOut
              NetworkLayout1[[l]][1] %>% strsplit("_") %>% unlist %>% "["(3) -> LayOutName
              paste(ConfidenceList0,"_",LayOutName,"_",colnames(FoldChange0)[p],"_",nrow(edges),"_",nrow(nodes),".pdf", sep="") -> FileName0
              resultpath %>% file.path("Final-network-all","Summary-network",ConfidenceList0,FileName0) -> FileName1
              FileName0 %>% paste("summary_",.,sep="") %>% sub(".pdf","",.) -> Title0
              pdf(FileName1, width=20, height=20 )
              plot.igraph(
                NetWork0, main=Title0, cex.main=2,
                edge.width=0.2, edge.weight=0.1, edge.color=E(NetWork0)$colorEdges,
                layout=LayOut, edge.curved=0.1, arrow.mode=3, edge.arrow.size=0.1, edge.arrow.width=0.5,
                vertex.label=V(NetWork0)$Symbol %>% as.character, vertex.size=(log2(Deg0)+4)*0.31, vertex.label.cex=V(NetWork0)$LabelCex,
                vertex.frame.color=V(NetWork0)$colorNode, vertex.label.font=2, vertex.label.color=V(NetWork0)$colorName, vertex.color=V(NetWork0)$colorNode )
              graphics.off()
            }
          } 
      }
  }
}
