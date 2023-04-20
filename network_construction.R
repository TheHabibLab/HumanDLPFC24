# Title     : TODO
# Objective : TODO
# Created by: OWNER
# Created on: 27/01/2021

# Red color
# Blue color
# Gray color
EDGE.STYLE <- c(NORMAL = 1, DASHED = 2, DOTTED = 3)
COR.EDGE.STYLE <- EDGE.STYLE["NORMAL"]
ANTI.EDGE.STYLE <- EDGE.STYLE["DASHED"]
RED.COLOR <- "#FF0000"
BLUE.COLOR <- "#0000FF"
GRAY_COLOR = "#CCCCCC"

pallete.bg <- c("#fea3a0", "#ffba99","#fed595", "#ffe696", "#fffb97", "#eaff98", "#b3eb92","#adf4ec", "#84d5fd", "#a8c5fc", "#d5adf3", "#ffadd1")[c(5, 11, 1, 8, 10,7, 4,9,  3, 6, 2, 12)]
pallete.border <- rev(c("#ff4954", "#fe794e", "#ffa84e", "#ffc54f", "#ffeb54", "#b8e74d", "#67d14b","#00cfc7", "#57c0fd", "#4d7ff2", "#9355da", "#ff58aa","#d9d9d9"))[c(5, 11, 1, 8, 10,7, 4,9,  3, 6, 2, 12)]


formate.matrix <- function(cor.matrix, thres.edge, thres.edge.neg){
  require(reshape2)
  thres.edge.neg <- if (thres.edge.neg == 0) -thres.edge else thres.edge.neg
  cor.matrix.m <- melt(cor.matrix)[melt(upper.tri(cor.matrix))$value,]
  colnames(cor.matrix.m) <- c("u", "v", "w")
  cor.matrix.m <- rbind(cor.matrix.m[cor.matrix.m$w >= thres.edge,],
                        cor.matrix.m[cor.matrix.m$w <= thres.edge.neg,])
  cor.matrix.m$sign <- plyr::mapvalues(sign(cor.matrix.m$w), c(-1, 1), c("neg", "pos"))
  return(cor.matrix.m)
}

check.nodes.infos <- function(cor.matrix, cell.types.info, features){
  if (length(setdiff(colnames(cor.matrix), rownames(cell.types.info))) > 0) {
    warning("Not all cell types present in infos dataframe")
  }
  if (!(all((features %in% c(colnames(cell.types.info), ""))))){
    warning("Feature not in cell types infos")
  }
}

calculate.layout <- function(cor.matrix, thres.edge, thres.edge.neg, cell.types.info, layout.file = '', save.layout = ''){
  require(igraph)
  matrix.plot <- formate.matrix(cor.matrix, thres.edge, thres.edge.neg)
  if (layout.file != "") return (read.csv(layout.file))
  matrix.plot.layout <-   matrix.plot[matrix.plot$sign == "pos", c("u", "v")]
  nodes.display <- rownames(cell.types.info)
  network.layout <- graph_from_data_frame(d=matrix.plot.layout, vertices=nodes.display, directed=F)
  layout <- layout_with_fr(network.layout)
  rownames(layout) <- nodes.display
  if (save.layout != "") write.csv(layout, save.layout)
  return (list(layout = layout, network = network.layout))
}

assign.nodes.color <- function(nodes.df.label){
  bg.cats <- unique(nodes.df.label$bg)
  nodes.df.label$bg <- plyr::mapvalues(nodes.df.label$bg, bg.cats, pallete.bg[1:length(bg.cats)])
  bd.cats <- unique(nodes.df.label$border)
  nodes.df.label$border <- plyr::mapvalues(nodes.df.label$border, bd.cats, pallete.border[1:length(bd.cats)])
  return(nodes.df.label)
}

calculate.nodes.df <- function(nodes.names, nodes.info, background, border, gray.nodes.by){
  if (length(gray.nodes.by) != 2){
    return()
  }
  nodes.df.label <- nodes.info[, c(background, border)]
  colnames(nodes.df.label) <- c("bg", "border")
  nodes.df <- assign.nodes.color(nodes.df.label)
  if (gray.nodes.by[1] != ""){
    to.grey <- rownames(nodes.info)[nodes.info[,gray.nodes.by[1]] == gray.nodes.by[2]]
    nodes.df[to.grey,"bg"] <- GRAY_COLOR
  }
  return (nodes.df)
}

calculate.edges.df <- function(matrix.plot, edge.colors, plot.by.weight){
  matrix.plot$style <- as.numeric(plyr::mapvalues(matrix.plot$sign, c("neg", "pos"), c(ANTI.EDGE.STYLE, 1)))
  if (!is.null(edge.colors)){
    matrix.plot$color <- plyr::mapvalues(matrix.plot$sign, c("neg", "pos"),
                                         c(edge.colors["neg"], edge.colors["pos"]))
  }
  matrix.plot$weight <- if (plot.by.weight) abs(matrix.plot$w) else 2
  return (matrix.plot)
}

plot.network <- function(edges.df, nodes.df, layout, save.fig){
  require(igraph);
  nodes.df <- cbind(rownames(nodes.df), nodes.df)
  colnames(nodes.df)[1] <- "subtype"
  edges.df$u <- as.vector(edges.df$u)
  edges.df$v <- as.vector(edges.df$v)
  network <- graph_from_data_frame(d=edges.df, vertices=nodes.df, directed=F)
  if (save.fig != "") while (!is.null(dev.list()))  dev.off()
  if (save.fig != "") pdf(save.fig)
  
  plot(network, layout = layout,
       edge.width=E(network)$weight, edge.lty = E(network)$style, edge.color=E(network)$color,
       vertex.color=V(network)$bg, vertex.frame.color=V(network)$border,
       label.cex = 0.5, size = 0.3, size2 = 0.3, vertex.size = 12)
  if (save.fig != "") while (!is.null(dev.list()))  dev.off()
  return(network)
}

#' @param cor.matrix A correlation matrix
#' @param thres.edge Threshold to build the layout
#' @param thres.edge.neg If the threshold to display negative weights is different than the positive weight
networks_graphs <- function(cor.matrix, layout.lst, cell.types.info = NULL, thres.edge = 0.3, thres.edge.neg = 0,
                             plot.fig = T, save.fig = '',
                            edge.colors = c(pos = RED.COLOR, neg = BLUE.COLOR), background.nodes.by = "", 
                            border.color.by = "", gray.nodes.by = c("", ""), edge.by.weight = F){
  if (border.color.by == "") {
    border.color.by = background.nodes.by
  } 
  check.nodes.infos(cor.matrix, cell.types.info, c(background.nodes.by, border.color.by, gray.nodes.by[1]))
  matrix.plot <- formate.matrix(cor.matrix, thres.edge, thres.edge.neg)
  if(!plot.fig) return()
  nodes.df <- calculate.nodes.df(rownames(layout.lst$layout), cell.types.info, background.nodes.by, border.color.by, gray.nodes.by)
  edges.df <- calculate.edges.df(matrix.plot, edge.colors, edge.by.weight)
  if (save.fig != "") {
    if(gray.nodes.by[1] == ""){
      save.fig <- paste(save.fig, "network_bg_", background.nodes.by, "_border_", border.color.by, ".pdf", sep = "")   
    } else {
      save.fig <- paste(save.fig, "network_bg_", background.nodes.by, "_border_", border.color.by,"_gray_by_",  gray.nodes.by[1], "_", gray.nodes.by[2], ".pdf", sep = "")   
    }
    
  }
  plot.network(edges.df, nodes.df, layout.lst$layout, save.fig)
}

calculate.communities <- function(network, cor.matrix, thres.edge = 0.3, thres.edge.neg = 0){
  matrix.plot <- formate.matrix(cor.matrix, thres.edge, thres.edge.neg)
  membership.obj <- leading.eigenvector.community(network)
  names(membership.obj$membership) <- membership.obj$names
  return(membership.obj$membership)
}

