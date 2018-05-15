## These are plotting functions to complement the dpt R package in order to make specific types of plots
## Many are simply lightly modified versions of the functions in the dpt R package (Haghverdi et al, 2016, Nat Methods, 13:845)

plot_dpt_ACR_all_col <- function(
  ts, pt, branches, shape_cat, col_vec, shape_name, col_name, 
  x = 1L, y = 2L,
  w_width = .1,
  path_col = c('red', 'darkgreen', 'purple'),
  ...,
  dcs = eig_decomp(ts@transitions, max(x, y))$vectors[, -1]
) {
  require(viridis)
  require(destiny)
  require(scran)
  require(smoother)
  stopifnot(is(ts, 'Transitions'))
  stopifnot(length(x) == 1L, length(y) == 1L)
  stopifnot(is.integer(branches))
  
  stopifnot(all(c('Branch', 'DPT', 'DPT.1', 'DPT.2')[1:(length(branches+1))] %in% names(pt)))
  idx <- vector('list', length(branches))
  idx <- lapply(branches, function(x) {as.integer(pt$Branch) %in% c(1L, 2L, x + 2L)})
  ## so this grabbed all uncertain, unassigned, and correct branch
  
  if (is.null(colnames(dcs))) {
    colnames(dcs) <- paste0('DC', seq_len(ncol(dcs)))
  }
  evs <- as.data.frame(as.matrix(dcs))[, c(x, y)]
  nms <- names(evs)
  
  minDPT <- which.min(pt$DPT)
  print(minDPT)
  if(evs[minDPT, 1] >0) {
    evs[,1] <- -evs[,1]
  }
  if(evs[minDPT, 2] >0) {
    evs[,2] <- -evs[,2]
  }
  
  DPT <- lapply(branches, function(x) {switch(x, pt$DPT, pt$DPT.1, pt$DPT.2)})
  ## so if branch is 1, this grabs DPT; if branch is 2, this grabs DPT.1; branch 3 DPT.2
  ## even if DPT.2 doesn't exist, this is fine unless branch 3 is called, which should have been stopped above
  
  path <- vector('list', length(DPT))
  for(i in 1:length(DPT)){
    path[[i]] <- average_path(DPT[[i]][idx[[i]]], evs[idx[[i]], ], w_width)
  }
  ## so now path includes all uncertain and unassigned
  
  if(length(branches)==3){
    return(ggplot(cbind(evs, DPT = DPT[[1]], shape_cat=shape_cat, col_vec=col_vec), aes_string(nms[[1L]], nms[[2L]], colour = 'col_vec'))
           + geom_point(aes_string(shape='shape_cat', colour='col_vec'), show.legend=TRUE)
           + labs(shape=shape_name, colour=col_name)
           + scale_shape_manual(values=c(1, 16, 17, 15, 3, 7, 8, 3))
           + scale_colour_gradientn(colours = viridis(256, option = "D"))
           + geom_path(data = path[[1]], colour = path_col[1])
           + geom_path(data = path[[2]], colour = path_col[2])
           + geom_path(data = path[[3]], colour = path_col[3])
           + theme_bw() 
           + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
    )
  } 
  if(length(branches)==2){
    return(ggplot(cbind(evs, DPT = DPT[[1]], shape_cat=shape_cat, col_vec=col_vec), aes_string(nms[[1L]], nms[[2L]], colour = 'col_vec'))
           + scale_colour_gradientn(colours = viridis(256, option = "D"))
           + geom_point(aes_string(shape='shape_cat', colour='col_vec'), show.legend=TRUE)
           + labs(shape=shape_name, colour=col_name)
           + scale_shape_manual(values=c(1, 16, 17, 15, 3, 7, 8, 3))
           + geom_path(data = path[[1]], colour = path_col[1])
           + geom_path(data = path[[2]], colour = path_col[2])
           + theme_bw() 
           + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
    )
  }
  if(length(branches)==1){
    return(ggplot(cbind(evs, DPT = DPT[[1]], shape_cat=shape_cat, col_vec=col_vec), aes_string(nms[[1L]], nms[[2L]], colour = 'col_vec'))
           + scale_colour_gradientn(colours = viridis(256, option = "D"))
           + geom_point(aes_string(shape='shape_cat', colour='col_vec'), show.legend=TRUE)
           + labs(shape=shape_name, colour=col_name)
           + scale_shape_manual(values=c(1, 16, 17, 15, 3, 7, 8, 3))
           + geom_path(data = path[[1]], colour = path_col[1])
           + theme_bw() 
           + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
    )
  }
  
}




plot_dpt_ACR_all_colcat <- function(
  ts, pt, branches, col_cat, cat_cols, col_name,
  x = 1L, y = 2L,
  w_width = .1,
  path_col = c('red', 'darkgreen', 'purple'),
  ...,
  dcs = eig_decomp(ts@transitions, max(x, y))$vectors[, -1]
) {
  require(destiny)
  require(scran)
  require(smoother)
  stopifnot(is(ts, 'Transitions'))
  stopifnot(length(x) == 1L, length(y) == 1L)
  stopifnot(is.integer(branches))
  
  stopifnot(all(c('Branch', 'DPT', 'DPT.1', 'DPT.2')[1:(length(branches+1))] %in% names(pt)))
  idx <- vector('list', length(branches))
  idx <- lapply(branches, function(x) {as.integer(pt$Branch) %in% c(1L, 2L, x + 2L)})
  ## so this grabbed all uncertain, unassigned, and correct branch
  
  if (is.null(colnames(dcs))) {
    colnames(dcs) <- paste0('DC', seq_len(ncol(dcs)))
  }
  evs <- as.data.frame(as.matrix(dcs))[, c(x, y)]
  nms <- names(evs)
  minDPT <- which.min(pt$DPT)
  print(minDPT)
  if(evs[minDPT, 1] >0) {
    evs[,1] <- -evs[,1]
  }
  if(evs[minDPT, 2] >0) {
    evs[,2] <- -evs[,2]
  }

  DPT <- lapply(branches, function(x) {switch(x, pt$DPT, pt$DPT.1, pt$DPT.2)})
  ## so if branch is 1, this grabs DPT; if branch is 2, this grabs DPT.1; branch 3 DPT.2
  ## even if DPT.2 doesn't exist, this is fine unless branch 3 is called, which should have been stopped above
  
  path <- vector('list', length(DPT))
  for(i in 1:length(DPT)){
    path[[i]] <- average_path(DPT[[i]][idx[[i]]], evs[idx[[i]], ], w_width)
  }
  ## so now path includes all uncertain and unassigned
  
  ## ideally this plotting would be adaptable to the number of branches
  if(length(branches)==3){
    return(ggplot(cbind(evs, DPT = DPT[[1]], col_vec=col_cat), aes_string(nms[[1L]], nms[[2L]], colour = 'col_vec'))
           + geom_point(aes_string(colour='col_vec'), show.legend=TRUE)
           + labs(colour=col_name)
           + scale_color_manual(values=cat_cols)
           + geom_path(data = path[[1]], colour = path_col[1])
           + geom_path(data = path[[2]], colour = path_col[2])
           + geom_path(data = path[[3]], colour = path_col[3])
           + theme_bw() 
           + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
    )
  } 
  if(length(branches)==2){
    return(ggplot(cbind(evs, DPT = DPT[[1]], col_vec=col_cat), aes_string(nms[[1L]], nms[[2L]], colour = 'col_vec'))
           + scale_color_manual(values=cat_cols)
           + geom_point(aes_string(colour='col_vec'), show.legend=TRUE)
           + labs(colour=col_name)
           + geom_path(data = path[[1]], colour = path_col[1])
           + geom_path(data = path[[2]], colour = path_col[2])
           + theme_bw() 
           + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
    )
  }
  if(length(branches)==1){
    return(ggplot(cbind(evs, DPT = DPT[[1]], col_vec=col_cat), aes_string(nms[[1L]], nms[[2L]], colour = 'col_vec'))
           + scale_color_manual(values=cat_cols)
           + geom_point(aes_string(colour='col_vec'), show.legend=TRUE)
           + labs(colour=col_name)
           + geom_path(data = path[[1]], colour = path_col[1])
           + theme_bw() 
           + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
    )
  }
  
}


plot_protein_DPT <- function(sce, pt, namepdf='null', scale=FALSE){
  require(destiny)
  require(scran)
  require(RColorBrewer)
  
  ylim <- c(0, 4)
  
  if(scale){
    sce$log10CD25 <- scale(sce$log10CD25)
    sce$log10CD69 <- scale(sce$log10CD69)
    sce$log10CD62L <- scale(sce$log10CD62L)
    sce$log10CD44 <- scale(sce$log10CD44)
    ylim <- c(-2, 2)
  }
  
  loCD25 <- loess(sce$log10CD25[order(pt$DPT)]~pt$DPT[order(pt$DPT)])
  loCD69 <- loess(sce$log10CD69[order(pt$DPT)]~pt$DPT[order(pt$DPT)])
  loCD62L <- loess(sce$log10CD62L[order(pt$DPT)]~pt$DPT[order(pt$DPT)])
  loCD44 <- loess(sce$log10CD44[order(pt$DPT)]~pt$DPT[order(pt$DPT)])
  
  cols2 <-  brewer.pal(4, 'Dark2')
  
  if(namepdf != 'null')
  {
    pdf(namepdf, height=6, width=3.5)
  }
  par(mgp=c(2,1,0), bty='l')
  plot(pt$DPT[order(pt$DPT)], predict(loCD62L), col=cols2[1], lwd=2, type='l', 
       ylim=ylim, xlab='pseudotime', ylab='protein expression')
  lines(pt$DPT[order(pt$DPT)], predict(loCD69), col=cols2[2], lwd=2)
  lines(pt$DPT[order(pt$DPT)], predict(loCD25), col=cols2[3], lwd=2)
  lines(pt$DPT[order(pt$DPT)], predict(loCD44), col=cols2[4], lwd=2)
  legend('topright', legend=c('CD62L', 'CD69', 'CD25', 'CD44'), fill=cols2, bty='n')
  if(namepdf != 'null')
  {
    dev.off()
  }
  
}

plot_enrich <- function(){
  require(edgeR)
  require(goseq)
  
}

dpt_distrib_plot <- function(pt, sce,
                             conditions=c('N4 6h', 'T4 6h', 'G4 6h', 'NP68 6h'), 
                             col=c('red', 'green3', 'blue', 'grey40'), 
                             trajectory_name){
  temp <- data.frame(DPT=pt$DPT, condition=sce$Condition)
  temp <- temp[temp$condition %in% conditions,]
  temp$condition <- factor(as.character(temp$condition), levels=conditions)
  ggplot(temp, aes(condition, DPT, fill=factor(condition))) + 
    geom_violin() + 
    scale_fill_manual(values=col) + 
    guides(fill=FALSE) + 
    theme_bw() + labs(y=trajectory_name) +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
}


quick_dpt <- function(sce, conditions, col, base_condition, stim_condition, pdf_name='null', trajectory_name='pseudotime'){
  ## It can happen that the root is assigned to a cell not from the least stimulated condition. 
  ## An extra set of arguments have been added to the function to use just the most null and most stimulated conditions 
  ## to find the root cell if a root is not found in the most unstimulated condition using all the data.
  
  ## Because these analyses were originally performed in 2016 before dpt was incorportated into the 
  ## destiny package, we implemented the old dpt function and wrote our own plotting functions.
  ## Since the updates to the destiny package, Transitions estimation of sigma now is different by a 
  ## factor of sqrt(2). To maintain compatibility with our old results, we have adjusted the sigmas by 
  ## this factor in our function below. 
  
  library(destiny)
  pcs <- prcomp(t(as.matrix(exprs(sce))), scale=TRUE)
  require(Matrix)
  dm <- DiffusionMap(data=pcs$x[,c(1:50)], distance='euclidean')
  sigs <- optimal_sigma(dm)*2^0.5
  ts <- Transitions(data=pcs$x[,c(1:50)], distance='euclidean', sigma=sigs)
  pt <- dpt(ts, branching=FALSE)
  rc_min <- which.min(pt$DPT)
  rc_max <- which.max(pt$DPT)
  if(sce$Condition[rc_min] %in% base_condition) {
    rc <- rc_min
    print('Root cell assigned.')
  } else if (sce$Condition[rc_max] %in% base_condition) {
    rc <- rc_max
    print('Root cell assigned.')
  } else {
    print('Unstimulated cell is not at an extreme. Using DPT fit of only the extreme conditions to identify root cell before applying to whole data.')
    sce_temp <- sce[,sce$Condition %in% c(base_condition,stim_condition)]
    pcs <- prcomp(t(as.matrix(exprs(sce_temp))), scale=TRUE)
    require(Matrix)
    dm <- DiffusionMap(data=pcs$x[,c(1:50)], distance='euclidean')
    sigs <- optimal_sigma(dm)*2^0.5
    ts <- Transitions(data=pcs$x[,c(1:50)], distance='euclidean', sigma=sigs)
    pt <- dpt(ts, branching=FALSE)
    rc_min <- which.min(pt$DPT)
    rc_max <- which.max(pt$DPT)
    if(sce_temp$Condition[rc_min] %in% base_condition) {
      rc <- rc_min
      print('Root cell assigned.')
    } else if (sce_temp$Condition[rc_max] %in% base_condition) {
      rc <- rc_max
      print('Root cell assigned.')
    }
    rc_name <- colnames(sce_temp)[rc]
    rc <- which(colnames(sce) %in% rc_name)
  }
  
  ## now back to main analysis
  pcs <- prcomp(t(as.matrix(exprs(sce))), scale=TRUE)
  require(Matrix)
  dm <- DiffusionMap(data=pcs$x[,c(1:50)], distance='euclidean')
  sigs <- optimal_sigma(dm)*2^0.5
  ts <- Transitions(data=pcs$x[,c(1:50)], distance='euclidean', sigma=sigs)
  pt <- dpt(ts, branching=FALSE, root=rc)
  
  col_cat=factor(sce$Condition,levels=conditions)
  y <- plot_dpt_ACR_all_colcat(ts=ts, pt=pt, branches=c(1L), path_col=c('black'),
                               col_cat=col_cat,
                               cat_cols=col,
                               col_name='Condition',
                               x = 1L, y = 2L, w_width = 0.4)
  a <- plot_dpt_ACR_all_col(ts=ts, pt=pt, branches=c(1L), path_col=c('black'),
                            shape_cat=col_cat,
                            col_vec=pt$DPT,
                            shape_name='Condition',
                            col_name='Pseudotime',
                            x = 1L, y = 2L, w_width = 0.4)
  z <- dpt_distrib_plot(pt, sce, conditions=conditions, col=col, trajectory_name=trajectory_name)
  return(list(p1=y, p2=z, p3=a, pcs=pcs, ts=ts, pt=pt))
}



