source("network_construction.R")
library(pheatmap)

cor.m.colors <- function(){
  colors <- c("#003BC0", "#FFFFFF", "#D80D29")
  paletteLength <- 100
  colors_step <- colorRampPalette(colors)(paletteLength)
  breaks <- c(seq(-1, 0, length.out=50), 
                   seq(0.01, 1, length.out=50))
  return(list(colors = colors_step, breaks = breaks))
}

props.m.colors <- function(props.m){
  colors <- c("#0A6265","#5BC3AE", "#E7E1BD" , "#ECA24D", "#DB4325")
  paletteLength <- 100
  seqs <- seq(0, max(props.m), length.out = 5)
  colors_step <- colorRampPalette(colors)(paletteLength)
  breaks <-   c(seq(seqs[1], seqs[2], length.out=25), 
                seq(seqs[2]+0.01, seqs[3], length.out=25),
                seq(seqs[3]+0.01, seqs[4], length.out=25),
                seq(seqs[4]+0.01, seqs[5], length.out=25))
  return(list(colors = colors_step, breaks = breaks))
}

calc.prop.by.feature <- function(dfs, group.by, ident.by){
  all.groups <- unique(Reduce(c, lapply(dfs, function(df) unique(df[,group.by]))))
  freqs.by.type <- lapply(names(dfs), function(df.name){
    df.tmp <- dfs[[df.name]][,c(group.by, ident.by)]
    #if(is.numeric(df.tmp[,group.by]))  df.tmp[,group.by] <- paste("g.", df.tmp[,group.by], sep = "")
    return(as.data.frame.matrix(prop.table(table(df.tmp[,c(group.by, ident.by)]), margin = 1))[all.groups,])
  })
  names(freqs.by.type) <- names(dfs)
  return(freqs.by.type)
}

plot.freq.by.type <- function(dfs, group.by, ident.by, path){
  dfs.freq <- calc.prop.by.feature(dfs, group.by, ident.by)
  dfs.freq.env <<- dfs.freq
  path.freqs <- paste(path, "/by.type/", sep = "")
  dir.create(path.freqs, recursive = T)
  for (df.name in names(dfs.freq)){
    df <- dfs.freq[[df.name]]
    colors.heatmap <- props.m.colors(df)
    while (!is.null(dev.list()))  dev.off()
    pdf(paste(path.freqs, df.name, ".pdf", sep = ""), height = 0.1*nrow(df) + 2, width = 0.5*ncol(df)+2)
    
    pheatmap(df, color = colors.heatmap$colors, breaks = colors.heatmap$breaks,
             main = paste(df.name, " - No scaling", sep = ""),
             )
    pheatmap(df, color = colors.heatmap$colors, breaks = colors.heatmap$breaks,
             main = paste(df.name, " - Scaling by row", sep = ""), scale= "row")
    pheatmap(df, color = colors.heatmap$colors, breaks = colors.heatmap$breaks,
             main = paste(df.name, " - Scaling by column", sep = ""), scale = "column")
    while (!is.null(dev.list()))  dev.off()
  }
}

plot.all.freqs <- function(dfs, group.by, ident.by, path){
  df.prop.by.type <- do.call(cbind, unname(calc.prop.by.feature(dfs, group.by, ident.by)))
  path.all.types <- paste(path, "/all.types/", sep = "")
  dir.create(path.all.types, recursive = T)
  scale.colors <- props.m.colors(df.prop.by.type)
  while (!is.null(dev.list()))  dev.off()
  pdf(paste(path.all.types, "frequencies.pdf", sep = ""), height = 0.1*nrow(df.prop.by.type) + 2, width = 0.5*ncol(df.prop.by.type)+2)
  pheatmap(df.prop.by.type, main = paste("No scaling", sep = ""), color = scale.colors$colors, breaks = scale.colors$breaks)
  pheatmap(df.prop.by.type, main = paste("Scaling by row", sep = ""), scale= "row", color = scale.colors$colors, breaks = scale.colors$breaks)
  pheatmap(df.prop.by.type, main = paste("Scaling by column", sep = ""), scale= "column", color = scale.colors$colors, breaks = scale.colors$breaks)
  while (!is.null(dev.list()))  dev.off()
}

plot.frequencies <- function(dfs, group.by, ident.by, path){
  plot.freq.by.type(dfs, group.by, ident.by, path)
  plot.all.freqs(dfs, group.by, ident.by, path)
}

plot_networks.of.subtypes <- function(dfs,infos.by.subtype, path, group.by, ident.by, 
                                      bg.by = "", frame.by = "", gray.by = c("", ""), plot.community = F,
                                      thres.edge = 0.3, thres.edge.neg = 0){
  cor.freqs <- cor(do.call(cbind, unname(calc.prop.by.feature(dfs, group.by, ident.by))),
                  method = "spearman")
  layout.lst <- calculate.layout(cor.freqs, thres.edge, thres.edge.neg, infos.by.subtype)
  
  networks_graphs(cor.freqs, layout.lst, infos.by.subtype, thres.edge, thres.edge.neg, 
                            save.fig = path,
                            background.nodes.by = bg.by,
                            border.color.by = frame.by,
                            gray.nodes.by = gray.by)
  
  if (plot.community){
    coms <- calculate.communities(layout.lst$network, cor.freqs, thres.edge, thres.edge.neg)
    infos.by.subtype$community <- "out" 
    infos.by.subtype[names(coms),"community"] <- coms
    networks_graphs(cor.freqs, layout.lst, infos.by.subtype, thres.edge, thres.edge.neg, 
                    plot.fig = T, save.fig = path,
                    background.nodes.by = "community",
                    border.color.by = frame.by,
                    gray.nodes.by = gray.by,
                    edge.by.weight = FALSE)
  }
}

plot.corr.matrix <- function(cor.matrix, path){
  path.all.types <- paste(path, "/all.types/", sep = "")
  dir.create(path.all.types, recursive = T)
  colors.heatmap <- cor.m.colors()
  while (!is.null(dev.list()))  dev.off()
  pdf(paste(path.all.types, "cor.pdf", sep = ""), height = 0.2*ncol(cor.matrix)+2, width = 0.2*ncol(cor.matrix)+2)
  pheatmap(cor.matrix, main = "Correlation Matrix", color = colors.heatmap$colors, breaks = colors.heatmap$breaks)
  while (!is.null(dev.list()))  dev.off()
}

plot_frequencies_and_cor <- function(dfs, group.by, ident.by, path){
  plot.frequencies(dfs, group.by, ident.by, path)
  cor.freqs <- cor(do.call(cbind, calc.prop.by.feature(dfs, group.by, ident.by)), method = "spearman")
  plot.corr.matrix(cor.freqs, path)
}

generate_fake_data <- function(){
  cell.types <- list("Astrocytes" = list(),
                     "Microglia" = list(),
                     "Neurons" = list())
  g1.n <- 10
  g2.n <- 10
  g3.n <- 10
  projids <- paste("projid.", seq(1, 30), sep = "")
  cell.types[["Astrocytes"]][["df"]] <- data.frame(
  astrocytes.0 = c(sample(seq(2000, 3000), g1.n), sample(seq(8000, 9000), g2.n), sample(seq(5000, 6000), g3.n)),
  astrocytes.1 = c(sample(seq(8000, 9000), g1.n), sample(seq(2000, 3000), g2.n), sample(seq(5000, 6000), g3.n)),
  astrocytes.2 = c(sample(seq(4000, 6000), g1.n), sample(seq(12000, 15000), g2.n), sample(seq(7000, 9000), g3.n)),
  astrocytes.3 = c(sample(seq(12000, 15000), g1.n), sample(seq(4000, 6000), g2.n), sample(seq(7000, 9000), g3.n)),
  astrocytes.4 = sample(seq(2500, 8000), g1.n + g2.n + g3.n),
  row.names = projids)
  cell.types[["Microglia"]][["df"]] <- data.frame(
  microglia.0 = c(sample(seq(2000, 3000), g1.n), sample(seq(8000, 9000), g2.n), sample(seq(5000, 6000), g3.n)),
  microglia.1 = c(sample(seq(8000, 9000), g1.n), sample(seq(2000, 3000), g2.n), sample(seq(5000, 6000), g3.n)),
  microglia.2 = c(sample(seq(4000, 6000), g1.n), sample(seq(12000, 15000), g2.n), sample(seq(7000, 9000), g3.n)),
  microglia.3 = sample(seq(2500, 10000), g1.n + g2.n + g3.n),
  row.names = projids)
  cell.types[["Neurons"]][["df"]] <- data.frame(
  neurons.0 = c(sample(seq(2000, 4000), g1.n), sample(seq(2000, 4000), g2.n), sample(seq(10000, 15000), g3.n)),
  neurons.1 = c(sample(seq(10000, 15000), g1.n), sample(seq(10000, 15000), g2.n), sample(seq(2000, 2500), g3.n)),
  neurons.2 = c(sample(seq(6000, 10000), g1.n), sample(seq(8000, 10000), g2.n), sample(seq(8000, 9000), g3.n)),
  neurons.3 = sample(seq(3000, 5000), g1.n + g2.n + g3.n),
  row.names = projids)
  n.subtypes <- sum(Reduce(c, lapply(cell.types, function(x){ncol(x$df)})))
  
  dfs <- list()
  for (cell.type.name in names(cell.types)){
    df <- cell.types[[cell.type.name]][["df"]]
    cell.subtypes <- Reduce(c, sapply(colnames(df), function(subtype){rep(subtype, sum(df[,subtype]))}))
    cell.projids <- Reduce(c, apply(df, 2, function(subtype){return(Reduce(c, sapply(names(subtype), function(projid){rep(projid, subtype[projid])})))}))

    dfs[[cell.type.name]] <- data.frame(cell.name = paste("astr.cell.",
                                                          1:sum(df), sep = ""),
                                        ident.is = cell.subtypes,
                                        projid.is = cell.projids)
  }
  infos.by.type <- data.frame(ident.subtype =        Reduce(c, sapply(cell.types, function(cell.type){colnames(cell.type$df)})),
                              n.cells =              Reduce(c, sapply(cell.types, function(cell.type){colSums(cell.type$df)})),
                              association.to.trait = sample(c(0, 1, 2), n.subtypes, replace = T),
                              cell.type = substr(Reduce(c, sapply(cell.types, function(cell.type){colnames(cell.type$df)})), 1, 3)
  )
  
  for (df.n in names(dfs)){
    df <- dfs[[df.n]]
    write.table(df, file = paste("input_for_frequencies/props_of_", df.n, ".txt", sep = ""))
  }
  write.table(infos.by.type, "input_for_frequencies/infos.txt")
  return(list(dfs = dfs, infos.by.subtype = infos.by.type, infos.by.projid = NULL))
}

main <- function(dfs.of.props, infos.by.subtype){
  plot_frequencies_and_cor(dfs.of.props, group.by = "projid.is", ident.by = "ident.is", path = "results_frequencies/")
  plot_networks.of.subtypes(dfs.of.props,  infos.by.subtype, "results_frequencies/all.types/", 
                            group.by = "projid.is", 
                            ident.by = "ident.is",
                            thres.edge = 0.3,
                            thres.edge.neg = -0.3,
                            bg.by = "cell.type", frame.by = "association.to.trait")

  plot_networks.of.subtypes(dfs.of.props,  infos.by.subtype, "results_frequencies/all.types/", 
                            group.by = "projid.is", 
                            ident.by = "ident.is",
                            thres.edge = 0.3,
                            thres.edge.neg = -0.3,
                            bg.by = "association.to.trait")
  
  

  plot_networks.of.subtypes(dfs.of.props,  infos.by.subtype, "results_frequencies/all.types/", 
                            group.by = "projid.is", 
                            ident.by = "ident.is",
                            thres.edge = 0.3,
                            thres.edge.neg = -0.3,
                            bg.by = "cell.type", plot.community = T)
  
  plot_networks.of.subtypes(dfs.of.props,  infos.by.subtype, "results_frequencies/all.types/", 
                            group.by = "projid.is", 
                            ident.by = "ident.is",
                            thres.edge = 0.3,
                            thres.edge.neg = -0.3,
                            bg.by = "association.to.trait", 
                            gray.by = c("association.to.trait", 0))
  }

#dfs.example <- list(astrocytes = read.table("input_for_frequencies/props_of_Astrocytes.txt"),
            #microglia = read.table("input_for_frequencies/props_of_Microglia.txt"),
            #neurons = read.table("input_for_frequencies/props_of_Neurons.txt"))

#infos.by.type <- read.table("input_for_frequencies/infos.txt")
#main(dfs.example, infos.by.type)

