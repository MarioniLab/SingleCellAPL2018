## these are functions from the original dpt package so that I can use my wrappers after updating R

#' Return average path over pseudotime
#' 
#' @param pt       The pseudotime for the cells of interest
#' @param x        Coordinates (e.g. from PCA or DM). coordinates \eqn{\times} dimensions
#' @param w_width  Window width for smoothing the path (see \code{\link[smoother]{smth.gaussian}})
#' 
#' @return \code{\link[base]{data.frame}} of smoothed coordinates
#' 
#' @importFrom graphics plot
#' @importFrom smoother smth.gaussian
#' @export
average_path <- function(pt, x, w_width = .1) {
  stopifnot(identical(nrow(x), length(pt)))
  as.data.frame(apply(x[order(pt), ], 2, function(col) smth.gaussian(col, w_width, tails = TRUE)))
}


#' Split a branch into three
#' 
#' Assign cells to a branch b by maximizing finite Kendall correlation of dpts of the
#' other two branches on b + finite Kendall anticorrelation of dpts on the other branches.
#' 
#' Both dpt and indb contain an entry for each cell, as this function performs the cut.
#'
#' @param dpt  Pseudo time distances for each of the three branches     (\eqn{n \times 3} matrix of double)
#' @param bid  Cell indices in ascending pseudotime distance per branch (\eqn{n \times 3} matrix of integer)
#' @param b    Index of branch to cut (1, 2, or 3)
#' 
#' @importFrom smoother smth.gaussian
#' @export
branchcut <- function(dpt, bid, b) {
  n <- nrow(bid)
  all_branches <- seq_len(3L)
  
  # sanity checks
  stopifnot(b %in% all_branches)
  stopifnot(ncol(dpt) == 3L, ncol(bid) == 3L)
  stopifnot(nrow(dpt) == n)
  stopifnot(is.double(dpt), is.integer(bid))
  
  # find cell indexes per branch 
  other <- all_branches[all_branches != b]
  b1 <- other[[1L]]
  b2 <- other[[2L]]
  
  # DPT for other branches, sorted by b3
  b3_idxs <- bid[, b]
  dpt1 <- dpt[b3_idxs, b1]
  dpt2 <- dpt[b3_idxs, b2]
  
  kcor <- vapply(seq_len(n - 1L), function(s1) {
    s2 <- s1 + 1L
    l <- seq_len(s1)
    r <- seq(s2, n)
    
    k_l <- kendall_finite_cor(dpt1[l], dpt2[l], dpt1[[s2]], dpt2[[s2]])
    k_r <- kendall_finite_cor(dpt1[r], dpt2[r], dpt1[[s1]], dpt2[[s1]])
    
    k_l/s1 - k_r/(n - s1)
  }, double(1))
  
  kcor <- smth.gaussian(kcor, 5L)
  cut <- which.max(kcor)
  
  b3_idxs[seq_len(cut)]
}


#' Forget all cached data
#' 
#' For speedups, this package uses \code{\link[memoise]{memoise}}.
#' Since this is a memory leak, we provide an easy way to free this memory.
#' 
#' @importFrom memoise forget
#' @export
dpt_clear_cache <- function() {
  forget(dpt_to_cell)
  forget(propagation_matrix)
}


#' Calculate DPTs
#' 
#' Given a cell index, returns the DPT of this cell to all cells
#' 
#' @param ts    A \code{\link{Transitions}} object
#' @param cell  Index of a cell for which all DPTs should be calculated
#' 
#' @export
dpt_to_cell <- memoise::memoise(function(ts, cell) {
  propagations <- propagation_matrix(ts)
  cell_propagations <- propagations[cell, ]
  
  apply(propagations, 1, function(row)
    sqrt(sum((cell_propagations - row) ^ 2)))
})


#' Diffusion Pseudo Time
#'
#' Create pseudotime ordering and assigns cell to one of three branches
#' 
#' @param ts           A \code{\link{Transitions}} object
#' @param branching    Detect a branching? (\code{TRUE} or \code{FALSE})
#' @param tips         Tip cell indices for each branch (integer vector of length 1 or 3)
#' @param root         The root index from which to calculate the DPTs (integer of length 1)
#' 
#' @details
#' All parameters are optional, but at least \code{branching} or \code{tips} has to be there.
#' If unspecified: \describe{
#'  \item{\code{root}}{will be the furthest cell from a random one}
#' 	\item{\code{branching}}{will be \code{TRUE} if multiple \code{tips} are specified}
#'  \item{\code{tips}}{will be \code{root} and (\code{if (branching)}) the most distant cells from it}
#' }
#' 
#' @return A \code{\link[base]{data.frame}} with the rows: \describe{
#' 	\item{\code{Branch}}{Branch labels for each cell, 1,2,3 or NA for undeceided}
#'  \item{\code{DPT}}{Diffusion pseudotime in respect to the root cell}
#'  \item{(\code{if (branching)}) \code{DPT.1}, \code{DPT.2}}{Diffusion pseudotime in respect to the other tips}
#' }
#' 
#' @export
dpt <- function(ts,
                branching = length(tips) > 1L,
                tips = if (branching) find_tips(ts, root) else root,
                root = random_root(ts)) {
  
  if (missing(branching) && missing(tips))
    stop('you need to specify at least `branching` or `tips`')
  
  n <- length(ts@phi0)
  stopifnot(is.logical(branching), length(branching) == 1L)
  stopifnot(is.integer(tips), length(tips) %in% c(1L, 3L))
  
  tip_labels <- rep(NA_integer_, n)
  tip_labels[tips] <- tips
  
  if (branching) {
    dpt <- vapply(tips, function(cell) dpt_to_cell(ts, cell), double(n))
    bid <- apply(dpt, 2, order)
    
    # cut it into three branches
    branch <- lapply(seq_len(3), function(b) branchcut(dpt, bid, b))
    
    unassigned <- setdiff(seq_len(3L), Reduce(union, branch, integer()))
    data.frame(
      Branch = organize_branches(branch, unassigned),
      Tips = tip_labels,
      DPT = dpt[, 1],
      DPT = dpt[, 2:3])  # DPT.1 & DPT.2
  } else {
    data.frame(
      Branch = rep(factor('branch 1'), n),
      Tips = tip_labels,
      DPT = dpt_to_cell(ts, tips[[1]]))
  }
}


#' Find tips in a Transitions object
#' 
#' @param ts         A \code{\link{Transitions}} object
#' @param root       Root cell index from which to find tips. (default: random)
#' 
#' @return An integer vector of length 1 or 3, depending if it was called using \code{branching} or not.
#' 
#' @export
find_tips <- function(ts, root = random_root(ts)) {
  x <- root
  dx <- dpt_to_cell(ts, x)
  y <- which.max(dx)
  dy <- dpt_to_cell(ts, y)
  z <- which.max(dx + dy)
  
  c(x, y, z)
}


# given two orderings b1 and b2, compute the delta (e.i. finite) Kendall
# correlation of adding a new cell with bnew1, bnew2 to the orderings.
kendall_finite_cor <- function(b1, b2, new1, new2) {
  b11 <- numeric(length(b1))
  b11[b1 >= new1] <- 1
  b11[b1 <  new1] <- -1
  
  b22 <- numeric(length(b2))
  b22[b2 >= new2] <- 1
  b22[b2 <  new2] <- -1
  
  b11 %*% b22
}


# This function organizes the cell arrays branch[[i]] and unassigned (created in dpt)
# to build Branch labels (length <- number of cells) indicating the branch each cell belongs to.
# Cells which are assigned to more than one branch in dpt as well
# as cells which are not assigned to any branch are defined as undeceided (label NA)
#' @importFrom utils combn
organize_branches <- function(branch, unassigned) {
  n <- do.call(max, branch)
  
  intersect_branches <- function(bs) intersect(branch[[bs[[1]]]], branch[[bs[[2]]]])
  branch_intersections <- lapply(combn(3L, 2L, simplify = FALSE), intersect_branches)
  inters <- Reduce(union, branch_intersections, integer())
  
  branch <- lapply(branch, function(b) setdiff(b, inters))
  branch_nums <- seq_along(branch)  # TODO: change
  
  branch_labels <- paste('branch', branch_nums)
  branches_label <- paste(branch_nums, collapse = ',')
  unassigned_label <- paste('unassigned', branches_label)
  uncertain_label <- paste('uncertain', branches_label)
  
  levels <- c(uncertain_label, unassigned_label, branch_labels)
  
  labels <- factor(rep(uncertain_label, n), levels)
  for (b in seq_along(branch)) {
    labels[branch[[b]]] <- paste('branch', branch_nums[[b]])
  }
  labels[unassigned] <- unassigned_label
  labels
}


#' Plot a DPT data.frame
#' 
#' Plots diffusion components from a Diffusion Map and the accompanying Diffusion pseudo time (\code{\link{dpt}})
#' 
#' @param ts        A \code{\link{Transitions}} object from which the dpt object was created
#' @param branches  Two numeric Branch IDs to use for the path. The first one will be used as the start of the DPT
#' @param pt        A \code{\link[base]{data.frame}} as returned by \code{\link{dpt}} (must contain the columns Branch and DPT)
#' @param x,y       The dimensions to use from the DiffusionMap
#' @param w_width   Window width for smoothing the path (see \code{\link[smoother]{smth.gaussian}})
#' @param path_col  Color for the path
#' @param ...       All parameters after this have to be specified by (full) name
#' @param dcs       Diffusion components created by eigen decomposition of the ts object
#' 
#' @importFrom methods is
#' @importFrom destiny eig_decomp
#' @importFrom ggplot2 ggplot geom_point geom_path aes_string
#' @export
plot_dpt <- function(
  ts, pt, branches,
  x = 1L, y = 2L,
  w_width = .1,
  path_col = 'red',
  ...,
  dcs = eig_decomp(ts@transitions, max(x, y))$vectors[, -1]
) {
  stopifnot(all(c('Branch', 'DPT', 'DPT.1', 'DPT.2') %in% names(pt)))
  stopifnot(is(ts, 'Transitions'))
  stopifnot(length(x) == 1L, length(y) == 1L)
  stopifnot(is.integer(branches))
  
  idx <- as.integer(pt$Branch) %in% c(1L, 2L, branches + 2L)
  
  if (is.null(colnames(dcs))) {
    colnames(dcs) <- paste0('DC', seq_len(ncol(dcs)))
  }
  evs <- as.data.frame(as.matrix(dcs))[, c(x, y)]
  nms <- names(evs)
  
  DPT <- switch(branches[[1L]], pt$DPT, pt$DPT.1, pt$DPT.2)
  
  path <- average_path(DPT[idx], evs[idx, ], w_width)
  
  (ggplot(cbind(evs, DPT = DPT), aes_string(nms[[1L]], nms[[2L]], colour = 'DPT'))
    + geom_point()
    + geom_path(data = path, colour = path_col))
}


#' @importFrom Matrix Diagonal
#' @importMethodsFrom Matrix solve
propagation_matrix <- memoise::memoise(function(ts) {
  n <- nrow(ts@transitions)
  inv <- solve(Diagonal(n) - ts@transitions + ts@phi0 %*% t(ts@phi0))
  inv - Diagonal(n)
})


#' Find a random root cell index
#' 
#' Finds a cell index whose cell has the maximum DPT distance from a randomly selected one.
#' 
#' @param ts  A \code{\link{Transitions}} object
#' 
#' @export
random_root <- function(ts) which.max(dpt_to_cell(ts, sample.int(length(ts@phi0), 1)))


#' Transitions class
#' 
#' @slot transitions  Transition matrix
#' @slot phi0         First eigenvector of transitions matrix
#' @slot d_rot        Diagonal normalization matrix
#' 
#' @name Transitions class
#' @aliases Transitions-class
#' @export
setClass(
  'Transitions',
  slots = c(
    transitions = 'dsCMatrix',
    phi0 = 'numeric',
    d_rot = 'ddiMatrix'
  ))

#' @param data      \code{\link{matrix}}, \code{\link{data.frame}}, or \code{\link[Biobase]{ExpressionSet}}
#' @param sigma     Parameter for gaussian kernel. Number or \code{\link[destiny]{Sigmas}} object
#' @param k         Number of nearest neighbors to use. NULL/Default: Whole matrix
#' @param distance  Distance metric to use
#' 
#' @importFrom methods new
#' @importFrom Matrix Diagonal
#' @importFrom destiny DiffusionMap
#' @name Transitions class
#' @export
Transitions <- function(
  data,
  sigma = NULL,
  k = NULL,
  distance = c('euclidean', 'cosine', 'rankcor')
) {
  dm <- DiffusionMap(data, sigma, k, distance = distance)
  phi0 <- dm@d_norm / sqrt(sum(dm@d_norm ^ 2))
  d_rot <- Diagonal(x = dm@d_norm ^ -.5)
  
  new('Transitions', transitions = dm@transitions, phi0 = phi0, d_rot = d_rot)
}
