#' @title stringDB interaction network
#' @param nodelist character list of Symbols to make network
#' @param foldchange data.frame Symbol list with fold-change values to display
#' @param species character for species hs mm rn dr
#' @param intmaxdh numeric maximum number of interaction to use for Davidson and Harel algorithm layout
#' @param confidence numeric for minus interaction confidence score
#' @param intmax numeric maximum number of interaction to use
#' @param nodesize numeric change Symbol size
#' @param title character plot title
#' @param path character file path
#' @param dirname character results directory name
#' @return result folder
#' @examples
#' # not run
#' # network(symbollist)
#' @author Florent Dumont <florent.dumont@universite-paris-saclay.fr>
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by left_join slice select slice_max
#' @importFrom rlang .data
#' @importFrom stats setNames
#' @importFrom igraph graph_from_data_frame degree V layout_as_tree layout_in_circle layout_on_grid layout_on_sphere layout_with_dh
#' @importFrom foreach foreach
#' @import Rgraphviz
#' @import moalstringdbhs
#' @import moalstringdbmm
#' @import moalstringdbrn
#' @import moalstringdbdr
#' @import moalstringdbss
#' @export
network <- function(
    nodelist = NULL , foldchange = NULL , species = NULL, confidence = 400,
    intmaxdh = 2500, intmax = 5000, nodesize = 0.39, title = "Network", path = ".", dirname = NULL )
{
  i=j=1
  moal:::orthoinfo -> orthoinfo
  # species selection
  ifelse(
    is.null(species),
    orthoinfo[[6]] -> Species0,
    orthoinfo %>% sapply("[[",1) %>% grep(species,.) %>% "[["(orthoinfo,.) -> Species0)
  # Nodes
  nodelist %>% as.character -> NodeList0
  NodeList0 %>% annot( species=Species0[1], dboutput="ebi", idtype = "SYMBOL" ) -> NodeList1
  # check
  if( nrow(NodeList1) > 0 &
      ( ( ( ( ( NodeList1$ENSGID %>% is.na ) | ( NodeList1$ENSGID %>% "=="("") ) | ( NodeList1$ENSGID %>% "=="("NA") ) ) %>% all ) |
      ( ( ( NodeList1$ENSPIDs %>% is.na ) | ( NodeList1$ENSPIDs %>% "=="("") ) | ( NodeList1$ENSPIDs %>% "=="("NA") ) ) %>% all ) ) %>% any %>% "!"(.) ) )
    { NWNODE1 <- TRUE }else{ NWNODE1 <- FALSE }
  #
  if( NWNODE1 )
  {
    # Nodes
    NodeList1$ENSPIDs %>% strsplit("\\|") %>% unlist %>% unique -> NodeList2
    NodeList2 %>% is.na %>% which -> sel ; if(length(sel) > 0){ NodeList2[-sel] -> NodeList2 }
    NodeList2 %>% "=="(.,"NA") %>% which -> sel ; if(length(sel) > 0){ NodeList2[-sel] -> NodeList2 }
    NodeList2 %>% as.character %>% annot(species=Species0[1], idtype="ENSP") %>% dplyr::select(.data$ENSPID,.data$Symbol) %>% data.frame -> NodeList3
    # check
    if( nrow(NodeList3) > 0 & NodeList3$Symbol %>% is.na %>% all %>% "!"(.) ){ NWNODE3 <- TRUE }else{ NWNODE3 <- FALSE }
    #
    if( NWNODE3 )
    {
      # stringdb loading
      if( Species0[1] == "hs" ){ moalstringdbhs::stringdbhs -> EdgeList0 }
      if( Species0[1] == "mm" ){ moalstringdbmm::stringdbmm -> EdgeList0 }
      if( Species0[1] == "rn" ){ moalstringdbrn::stringdbrn -> EdgeList0 }
      if( Species0[1] == "dr" ){ moalstringdbdr::stringdbdr -> EdgeList0 }
      if( Species0[1] == "ss" ){ moalstringdbss::stringdbss -> EdgeList0 }
      # Edges
      NodeList3 %>% dplyr::inner_join(EdgeList0, by =c("ENSPID"="NodeA")) %>% setNames(c("NodeA","SymbolA","NodeB","CombinedScore")) -> EdgeList1
      # Direct connection filtering
      NodeList3 %>% inner_join(EdgeList1, by=c("ENSPID"="NodeB")) %>% dplyr::select(c(1:3,5,4)) %>% setNames(c("NodeA","SymbolA","NodeB","CombinedScore","SymbolB")) -> EdgeList2
      # Confidence filtering
      EdgeList2 %>% dplyr::filter( .data$CombinedScore >= confidence ) -> EdgeList3
      # remove NA
      is.na(EdgeList3$SymbolA) %>% which -> sel ; if( length(sel) > 0 ){ EdgeList3[-sel,] -> EdgeList3  }
      EdgeList3$SymbolA %>% grep("^$|^NA$|^row",.) -> sel ; if( length(sel) > 0 ){ NodeList3[-sel,] -> NodeList3  }
      is.na(EdgeList3$SymbolB) %>% which -> sel ; if( length(sel) > 0 ){ EdgeList3[-sel,] -> EdgeList3  }
      EdgeList3$SymbolB %>% grep("^$|^NA$|^row",.) -> sel ; if( length(sel) > 0 ){ NodeList3[-sel,] -> NodeList3  }
      if( nrow(EdgeList3) > 0 & nrow(EdgeList3) < intmax )
      {
        # result directory
        ifelse( is.null(dirname), "nw" -> DirName, paste("nw_",dirname,sep="") -> DirName)
        path %>% file.path( DirName ) -> Path
        Path %>% dir.create
        #
        EdgeList3 %>% dplyr::select(.data$NodeA,.data$NodeB) -> EdgeList4
        EdgeList3 %>% dplyr::select(.data$NodeA,.data$SymbolA) %>% setNames(c("ENSPID","Symbol")) %>% data.frame ->  NodeListA
        EdgeList3 %>% dplyr::select(.data$NodeB,.data$SymbolB) %>% setNames(c("ENSPID","Symbol")) %>% data.frame ->  NodeListB
        rbind(NodeListA,NodeListB) %>% unique -> NodeList4
        #
        igraph::graph_from_data_frame(d=EdgeList4, vertices=NodeList4, directed=F) -> NetWork0
        igraph::degree(NetWork0, mode="all") -> Deg0
        rep("gray",length(V(NetWork0)$Symbol)) -> colorNode
        # layout
        list( c( "layout_as_tree","tree"),c("layout_in_circle","circle"),c("layout_on_grid","grid"),
              c("layout_on_sphere","sphere"),c("layout_with_fr","fr"),c("layout_with_dh","dh") ) -> NetworkLayout0
        if( nrow(EdgeList3) < intmaxdh ){ NetworkLayout0 -> NetworkLayout1 }else{ NetworkLayout0[-6] -> NetworkLayout1 }
        # Nodes
        # V(NetWork0)$Symbol %>% as.character %>% data.frame(Symbol=.,stringsAsFactors=F) %>% inner_join(NodeList1) %>% data.frame -> NodeList
        data.frame(Symbol=V(NetWork0)$Symbol %>% as.character,ENSPID=V(NetWork0) %>% attr("names") %>% as.character,stringsAsFactors=F) %>%
          inner_join(NodeList1) %>% data.frame -> NodeList
        # Nodes output
        c(EdgeList3$SymbolA,EdgeList3$SymbolA) %>% table %>% data.frame %>% setNames(c("Symbol","Connectivity")) -> Hub0
        NodeList %>% inner_join(Hub0,.) %>% dplyr::arrange(desc(.data$Connectivity)) -> NodeListoutput
        NodeListoutput %>% nrow -> NodeNumber
        paste("Nodes_",NodeNumber,".tsv",sep = "") -> FileName0
        file.path(Path,FileName0) -> FileName1
        NodeListoutput %>% output(FileName1)
        # Edges output
        paste("Edges_",nrow(EdgeList3),".tsv",sep ="") -> FileName0
        file.path(Path,FileName0) -> FileName1
        EdgeList3 %>% output(FileName1)
        # plot
        foreach( i=1:length(NetworkLayout1) ) %do%
          {
            #
            set.seed(1234567)
            eval( parse( text = paste( NetworkLayout1[[i]][1],"(NetWork0)",sep ="") ) ) -> LayOut
            # fold-change
            if( !is.null(foldchange) )
            {
              foreach( j=2:ncol(foldchange)) %do%
                {
                  networkcolor( foldchange=foldchange[,c(1,j)], colornode=colorNode, networkobj=NetWork0 ) -> NetworkColor0
                  NetworkColor0[[1]] -> colorNodeN
                  NodeList %>% left_join( NetworkColor0[[2]] ) %>% data.frame -> NodeList
                  #
                  paste( NetworkLayout1[[i]][2],"_", sub( "fc_(.*)", "\\1", colnames(foldchange)[j] ),
                         "_", nrow(EdgeList3), "_", NodeNumber, ".pdf", sep="" ) -> FileName0
                  file.path(Path,FileName0) -> FileName1
                  ifelse(is.null(title), dirname -> Title0, title -> Title0)
                  pdf(FileName1 , width = 10 , height = 10 )
                  plot( NetWork0, rescale=T, main=Title0, cex.main=2,
                        edge.width=0.2, edge.weight=0.1, edge.color="gray",
                        layout=LayOut, edge.curved=0.1, arrow.mode=3, edge.arrow.size=0.1, edge.arrow.width=0.5,
                        vertex.label=V(NetWork0)$Symbol %>% as.character, vertex.size=(log2(Deg0)+5)*0.5, vertex.label.cex=nodesize,
                        vertex.frame.color=colorNodeN, vertex.label.font=2, vertex.label.color="blue", vertex.color=colorNodeN )
                  graphics.off()
                }
            }else{
              paste(NetworkLayout1[[i]][2],"_",nrow(EdgeList3),"_",NodeNumber,".pdf",sep="") -> FileName0
              file.path(Path,FileName0) -> FileName1
              pdf(FileName1 , width = 10 , height = 10 )
              plot( NetWork0, rescale=T, main=dirname, cex.main=2,
                    edge.width=0.2, edge.weight=0.1, edge.color="gray",
                    layout=LayOut, edge.curved=0.1, arrow.mode=3, edge.arrow.size=0.1, edge.arrow.width=0.5,
                    vertex.label=V(NetWork0)$Symbol %>% as.character, vertex.size=(log2(Deg0)+5)*0.5, vertex.label.cex=nodesize,
                    vertex.frame.color=colorNode, vertex.label.font=2, vertex.label.color="blue", vertex.color=colorNode )
              graphics.off() }
          }
      }
    }# fi NWNODE3
  }# fi NWNODE1
}