# Used by:

# add_corners
#######################################################
# Get the corner coordinates
# 
# Gets the coordinates of the corners when the tree is plotted.
#
# Because this function needs to use a modified version of the APE plot.phylo
# function, and for complex reasons APE's .C functions cannot be used 
# elsewhere without causing problems with R CMD check, this function is
# left up to user specification.  Basically, the user puts in
# the name of the function, which is available in the extension data
# (extdata) directory of the package.  The defaults work on the developer's 
# machine, other users may have to e.g. change "manual" to "auto".
#
# @param tr A tree object in \code{\link[ape]{phylo}} format.
# @param coords_fun The name of the function to use to get node coordinates. Default: 
# "plot_phylo3_nodecoords".
# @param location Default is "manual", which looks for whatever is hard-coded. Most users
# should probably use "auto", which will look in /extdata.




#######################################################
# plot_phylo3_nodecoords
#######################################################
#' Modified plot.phylo to return the node coordinates
#' 
#' Just the standard plot.phylo (for phylo3 APE objects)
#' In addition to the standard outputs, this function returns xx and yy,
#' the plotting coordinates for those nodes.
#' 
#' @param x An ape phylo object
#' @param type Default "phylogram"
#' @param use.edge.length Default TRUE
#' @param node.pos Default NULL
#' @param show.tip.label Default TRUE
#' @param show.node.label Default FALSE
#' @param edge.color Default "black"
#' @param edge.width Default 1
#' @param edge.lty Default 1
#' @param font Default 3
#' @param cex Default par("cex")
#' @param adj Default NULL
#' @param srt Default 0
#' @param no.margin Default FALSE
#' @param root.edge Default FALSE
#' @param label.offset Default 0
#' @param underscore Default FALSE
#' @param x.lim Default NULL
#' @param y.lim Default NULL
#' @param direction Default "rightwards"
#' @param lab4ut Default "horizontal"
#' @param tip.color Default "black"
#' @param plot Default FALSE
#' @param rotate.tree Default 0
#' @param ... Additional arguments to standard functions
#' @return A typical plot object, but with node coordinates added
#' @export
#' @seealso \code{\link[stats]{optim}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' testval=1
#'
plot_phylo3_nodecoords <- function (x, type = "phylogram", use.edge.length = TRUE, node.pos = NULL, 
    show.tip.label = TRUE, show.node.label = FALSE, edge.color = "black", 
    edge.width = 1, edge.lty = 1, font = 3, cex = par("cex"), 
    adj = NULL, srt = 0, no.margin = FALSE, root.edge = FALSE, 
    label.offset = 0, underscore = FALSE, x.lim = NULL, y.lim = NULL, 
    direction = "rightwards", lab4ut = "horizontal", tip.color = "black", 
    plot = FALSE, rotate.tree = 0, ...) 
{
    require("ape")
    
    # Try registering foreign functions as in here:
	# http://cran.r-project.org/doc/manuals/R-exts.html#System-and-foreign-language-interfaces
    # p_myCfun = R_GetCCallable("packA", "myCfun");
    
    Ntip <- length(x$tip.label)
    if (Ntip == 1) {
        warning("found only one tip in the tree")
        return(NULL)
    }
    if (any(tabulate(x$edge[, 1]) == 1)) 
        stop("there are single (non-splitting) nodes in your tree; you may need to use collapse.singles()")
    .nodeHeight <- function(Ntip, Nnode, edge, Nedge, yy) .C("node_height", 
        as.integer(Ntip), as.integer(Nnode), as.integer(edge[, 
            1]), as.integer(edge[, 2]), as.integer(Nedge), as.double(yy), 
        DUP = FALSE, PACKAGE = "ape")[[6]]
    .nodeDepth <- function(Ntip, Nnode, edge, Nedge) .C("node_depth", 
        as.integer(Ntip), as.integer(Nnode), as.integer(edge[, 
            1]), as.integer(edge[, 2]), as.integer(Nedge), double(Ntip + 
            Nnode), DUP = FALSE, PACKAGE = "ape")[[6]]
    .nodeDepthEdgelength <- function(Ntip, Nnode, edge, Nedge, 
        edge.length) .C("node_depth_edgelength", as.integer(Ntip), 
        as.integer(Nnode), as.integer(edge[, 1]), as.integer(edge[, 
            2]), as.integer(Nedge), as.double(edge.length), double(Ntip + 
            Nnode), DUP = FALSE, PACKAGE = "ape")[[7]]
    Nedge <- dim(x$edge)[1]
    Nnode <- x$Nnode
    ROOT <- Ntip + 1
    type <- match.arg(type, c("phylogram", "cladogram", "fan", 
        "unrooted", "radial"))
    direction <- match.arg(direction, c("rightwards", "leftwards", 
        "upwards", "downwards"))
    if (is.null(x$edge.length)) 
        use.edge.length <- FALSE
    if (type %in% c("unrooted", "radial") || !use.edge.length || 
        is.null(x$root.edge) || !x$root.edge) 
        root.edge <- FALSE
    if (type == "fan" && root.edge) {
        warning("drawing root edge with type = 'fan' is not yet supported")
        root.edge <- FALSE
    }
    phyloORclado <- type %in% c("phylogram", "cladogram")
    horizontal <- direction %in% c("rightwards", "leftwards")
    xe <- x$edge
    if (phyloORclado) {
        phyOrder <- attr(x, "order")
        if (is.null(phyOrder) || phyOrder != "cladewise") {
            x <- reorder(x)
            if (!identical(x$edge, xe)) {
                ereorder <- match(x$edge[, 2], xe[, 2])
                if (length(edge.color) > 1) {
                  edge.color <- rep(edge.color, length.out = Nedge)
                  edge.color <- edge.color[ereorder]
                }
                if (length(edge.width) > 1) {
                  edge.width <- rep(edge.width, length.out = Nedge)
                  edge.width <- edge.width[ereorder]
                }
                if (length(edge.lty) > 1) {
                  edge.lty <- rep(edge.lty, length.out = Nedge)
                  edge.lty <- edge.lty[ereorder]
                }
            }
        }
        yy <- numeric(Ntip + Nnode)
        TIPS <- x$edge[x$edge[, 2] <= Ntip, 2]
        yy[TIPS] <- 1:Ntip
    }
    z <- reorder(x, order = "pruningwise")
    if (phyloORclado) {
        if (is.null(node.pos)) {
            node.pos <- 1
            if (type == "cladogram" && !use.edge.length) 
                node.pos <- 2
        }
        if (node.pos == 1) 
            yy <- .nodeHeight(Ntip, Nnode, z$edge, Nedge, yy)
        else {
            ans <- .C("node_height_clado", as.integer(Ntip), 
                as.integer(Nnode), as.integer(z$edge[, 1]), as.integer(z$edge[, 
                  2]), as.integer(Nedge), double(Ntip + Nnode), 
                as.double(yy), DUP = FALSE, PACKAGE = "ape")
            xx <- ans[[6]] - 1
            yy <- ans[[7]]
        }
        if (!use.edge.length) {
            if (node.pos != 2) 
                xx <- .nodeDepth(Ntip, Nnode, z$edge, Nedge) - 
                  1
            xx <- max(xx) - xx
        }
        else {
            xx <- .nodeDepthEdgelength(Ntip, Nnode, z$edge, Nedge, 
                z$edge.length)
        }
    }
    else {
        rotate.tree <- 2 * pi * rotate.tree/360
        switch(type, fan = {
            TIPS <- x$edge[which(x$edge[, 2] <= Ntip), 2]
            xx <- seq(0, 2 * pi * (1 - 1/Ntip), 2 * pi/Ntip)
            theta <- double(Ntip)
            theta[TIPS] <- xx
            theta <- c(theta, numeric(Nnode))
            theta <- .nodeHeight(Ntip, Nnode, z$edge, Nedge, 
                theta)
            if (use.edge.length) {
                r <- .nodeDepthEdgelength(Ntip, Nnode, z$edge, 
                  Nedge, z$edge.length)
            } else {
                r <- .nodeDepth(Ntip, Nnode, z$edge, Nedge)
                r <- 1/r
            }
            theta <- theta + rotate.tree
            xx <- r * cos(theta)
            yy <- r * sin(theta)
        }, unrooted = {
            nb.sp <- .nodeDepth(Ntip, Nnode, z$edge, Nedge)
            XY <- if (use.edge.length) unrooted.xy(Ntip, Nnode, 
                z$edge, z$edge.length, nb.sp, rotate.tree) else unrooted.xy(Ntip, 
                Nnode, z$edge, rep(1, Nedge), nb.sp, rotate.tree)
            xx <- XY$M[, 1] - min(XY$M[, 1])
            yy <- XY$M[, 2] - min(XY$M[, 2])
        }, radial = {
            X <- .nodeDepth(Ntip, Nnode, z$edge, Nedge)
            X[X == 1] <- 0
            X <- 1 - X/Ntip
            yy <- c((1:Ntip) * 2 * pi/Ntip, rep(0, Nnode))
            Y <- .nodeHeight(Ntip, Nnode, z$edge, Nedge, yy)
            xx <- X * cos(Y + rotate.tree)
            yy <- X * sin(Y + rotate.tree)
        })
    }
    if (phyloORclado) {
        if (!horizontal) {
            tmp <- yy
            yy <- xx
            xx <- tmp - min(tmp) + 1
        }
        if (root.edge) {
            if (direction == "rightwards") 
                xx <- xx + x$root.edge
            if (direction == "upwards") 
                yy <- yy + x$root.edge
        }
    }
    if (no.margin) 
        par(mai = rep(0, 4))
    if (is.null(x.lim)) {
        if (phyloORclado) {
            if (horizontal) {
                x.lim <- c(0, NA)
                pin1 <- par("pin")[1]
                strWi <- strwidth(x$tip.label, "inches")
                xx.tips <- xx[1:Ntip] * 1.04
                alp <- try(uniroot(function(a) max(a * xx.tips + 
                  strWi) - pin1, c(0, 1e+06))$root, silent = TRUE)
                if (is.character(alp)) 
                  tmp <- max(xx.tips) * 1.5
                else {
                  tmp <- if (show.tip.label) 
                    max(xx.tips + strWi/alp)
                  else max(xx.tips)
                }
                x.lim[2] <- tmp
            }
            else x.lim <- c(1, Ntip)
        }
        else switch(type, fan = {
            if (show.tip.label) {
                offset <- max(nchar(x$tip.label) * 0.018 * max(yy) * 
                  cex)
                x.lim <- c(min(xx) - offset, max(xx) + offset)
            } else x.lim <- c(min(xx), max(xx))
        }, unrooted = {
            if (show.tip.label) {
                offset <- max(nchar(x$tip.label) * 0.018 * max(yy) * 
                  cex)
                x.lim <- c(0 - offset, max(xx) + offset)
            } else x.lim <- c(0, max(xx))
        }, radial = {
            if (show.tip.label) {
                offset <- max(nchar(x$tip.label) * 0.03 * cex)
                x.lim <- c(-1 - offset, 1 + offset)
            } else x.lim <- c(-1, 1)
        })
    }
    else if (length(x.lim) == 1) {
        x.lim <- c(0, x.lim)
        if (phyloORclado && !horizontal) 
            x.lim[1] <- 1
        if (type %in% c("fan", "unrooted") && show.tip.label) 
            x.lim[1] <- -max(nchar(x$tip.label) * 0.018 * max(yy) * 
                cex)
        if (type == "radial") 
            x.lim[1] <- if (show.tip.label) 
                -1 - max(nchar(x$tip.label) * 0.03 * cex)
            else -1
    }
    if (phyloORclado && direction == "leftwards") 
        xx <- x.lim[2] - xx
    if (is.null(y.lim)) {
        if (phyloORclado) {
            if (horizontal) 
                y.lim <- c(1, Ntip)
            else {
                y.lim <- c(0, NA)
                pin2 <- par("pin")[2]
                strWi <- strwidth(x$tip.label, "inches")
                yy.tips <- yy[1:Ntip] * 1.04
                alp <- try(uniroot(function(a) max(a * yy.tips + 
                  strWi) - pin2, c(0, 1e+06))$root, silent = TRUE)
                if (is.character(alp)) 
                  tmp <- max(yy.tips) * 1.5
                else {
                  tmp <- if (show.tip.label) 
                    max(yy.tips + strWi/alp)
                  else max(yy.tips)
                }
                y.lim[2] <- tmp
            }
        }
        else switch(type, fan = {
            if (show.tip.label) {
                offset <- max(nchar(x$tip.label) * 0.018 * max(yy) * 
                  cex)
                y.lim <- c(min(yy) - offset, max(yy) + offset)
            } else y.lim <- c(min(yy), max(yy))
        }, unrooted = {
            if (show.tip.label) {
                offset <- max(nchar(x$tip.label) * 0.018 * max(yy) * 
                  cex)
                y.lim <- c(0 - offset, max(yy) + offset)
            } else y.lim <- c(0, max(yy))
        }, radial = {
            if (show.tip.label) {
                offset <- max(nchar(x$tip.label) * 0.03 * cex)
                y.lim <- c(-1 - offset, 1 + offset)
            } else y.lim <- c(-1, 1)
        })
    }
    else if (length(y.lim) == 1) {
        y.lim <- c(0, y.lim)
        if (phyloORclado && horizontal) 
            y.lim[1] <- 1
        if (type %in% c("fan", "unrooted") && show.tip.label) 
            y.lim[1] <- -max(nchar(x$tip.label) * 0.018 * max(yy) * 
                cex)
        if (type == "radial") 
            y.lim[1] <- if (show.tip.label) 
                -1 - max(nchar(x$tip.label) * 0.018 * max(yy) * 
                  cex)
            else -1
    }
    if (phyloORclado && direction == "downwards") 
        yy <- max(yy) - yy
    if (phyloORclado && root.edge) {
        if (direction == "leftwards") 
            x.lim[2] <- x.lim[2] + x$root.edge
        if (direction == "downwards") 
            y.lim[2] <- y.lim[2] + x$root.edge
    }
    asp <- if (type %in% c("fan", "radial", "unrooted")) 
        1
    else NA
    if (plot)
    	{  	plot(0, type = "n", xlim = x.lim, ylim = y.lim, ann = FALSE, 
        	axes = FALSE, asp = asp, ...) 	}
    if (plot) {
        if (is.null(adj)) 
            adj <- if (phyloORclado && direction == "leftwards") 
                1
            else 0
        if (phyloORclado && show.tip.label) {
            MAXSTRING <- max(strwidth(x$tip.label, cex = cex))
            loy <- 0
            if (direction == "rightwards") {
                lox <- label.offset + MAXSTRING * 1.05 * adj
            }
            if (direction == "leftwards") {
                lox <- -label.offset - MAXSTRING * 1.05 * (1 - 
                  adj)
            }
            if (!horizontal) {
                psr <- par("usr")
                MAXSTRING <- MAXSTRING * 1.09 * (psr[4] - psr[3])/(psr[2] - 
                  psr[1])
                loy <- label.offset + MAXSTRING * 1.05 * adj
                lox <- 0
                srt <- 90 + srt
                if (direction == "downwards") {
                  loy <- -loy
                  srt <- 180 + srt
                }
            }
        }
        if (type == "phylogram") {
            phylogram.plot(x$edge, Ntip, Nnode, xx, yy, horizontal, 
                edge.color, edge.width, edge.lty)
        }
        else {
            if (type == "fan") {
                ereorder <- match(z$edge[, 2], x$edge[, 2])
                if (length(edge.color) > 1) {
                  edge.color <- rep(edge.color, length.out = Nedge)
                  edge.color <- edge.color[ereorder]
                }
                if (length(edge.width) > 1) {
                  edge.width <- rep(edge.width, length.out = Nedge)
                  edge.width <- edge.width[ereorder]
                }
                if (length(edge.lty) > 1) {
                  edge.lty <- rep(edge.lty, length.out = Nedge)
                  edge.lty <- edge.lty[ereorder]
                }
                circular.plot(z$edge, Ntip, Nnode, xx, yy, theta, 
                  r, edge.color, edge.width, edge.lty)
            }
            else cladogram.plot(x$edge, xx, yy, edge.color, edge.width, 
                edge.lty)
        }
        if (root.edge) 
            switch(direction, rightwards = segments(0, yy[ROOT], 
                x$root.edge, yy[ROOT]), leftwards = segments(xx[ROOT], 
                yy[ROOT], xx[ROOT] + x$root.edge, yy[ROOT]), 
                upwards = segments(xx[ROOT], 0, xx[ROOT], x$root.edge), 
                downwards = segments(xx[ROOT], yy[ROOT], xx[ROOT], 
                  yy[ROOT] + x$root.edge))
        if (show.tip.label) {
            if (is.expression(x$tip.label)) 
                underscore <- TRUE
            if (!underscore) 
                x$tip.label <- gsub("_", " ", x$tip.label)
            if (phyloORclado) 
                text(xx[1:Ntip] + lox, yy[1:Ntip] + loy, x$tip.label, 
                  adj = adj, font = font, srt = srt, cex = cex, 
                  col = tip.color)
            if (type == "unrooted") {
                if (lab4ut == "horizontal") {
                  y.adj <- x.adj <- numeric(Ntip)
                  sel <- abs(XY$axe) > 0.75 * pi
                  x.adj[sel] <- -strwidth(x$tip.label)[sel] * 
                    1.05
                  sel <- abs(XY$axe) > pi/4 & abs(XY$axe) < 0.75 * 
                    pi
                  x.adj[sel] <- -strwidth(x$tip.label)[sel] * 
                    (2 * abs(XY$axe)[sel]/pi - 0.5)
                  sel <- XY$axe > pi/4 & XY$axe < 0.75 * pi
                  y.adj[sel] <- strheight(x$tip.label)[sel]/2
                  sel <- XY$axe < -pi/4 & XY$axe > -0.75 * pi
                  y.adj[sel] <- -strheight(x$tip.label)[sel] * 
                    0.75
                  text(xx[1:Ntip] + x.adj * cex, yy[1:Ntip] + 
                    y.adj * cex, x$tip.label, adj = c(adj, 0), 
                    font = font, srt = srt, cex = cex, col = tip.color)
                }
                else {
                  adj <- abs(XY$axe) > pi/2
                  srt <- 180 * XY$axe/pi
                  srt[adj] <- srt[adj] - 180
                  adj <- as.numeric(adj)
                  xx.tips <- xx[1:Ntip]
                  yy.tips <- yy[1:Ntip]
                  if (label.offset) {
                    xx.tips <- xx.tips + label.offset * cos(XY$axe)
                    yy.tips <- yy.tips + label.offset * sin(XY$axe)
                  }
                  font <- rep(font, length.out = Ntip)
                  tip.color <- rep(tip.color, length.out = Ntip)
                  cex <- rep(cex, length.out = Ntip)
                  for (i in 1:Ntip) text(xx.tips[i], yy.tips[i], 
                    cex = cex[i], x$tip.label[i], adj = adj[i], 
                    font = font[i], srt = srt[i], col = tip.color[i])
                }
            }
            if (type %in% c("fan", "radial")) {
                xx.tips <- xx[1:Ntip]
                yy.tips <- yy[1:Ntip]
                angle <- atan2(yy.tips, xx.tips)
                if (label.offset) {
                  xx.tips <- xx.tips + label.offset * cos(angle)
                  yy.tips <- yy.tips + label.offset * sin(angle)
                }
                s <- xx.tips < 0
                angle <- angle * 180/pi
                angle[s] <- angle[s] + 180
                adj <- as.numeric(s)
                font <- rep(font, length.out = Ntip)
                tip.color <- rep(tip.color, length.out = Ntip)
                cex <- rep(cex, length.out = Ntip)
                for (i in 1:Ntip) text(xx.tips[i], yy.tips[i], 
                  x$tip.label[i], font = font[i], cex = cex[i], 
                  srt = angle[i], adj = adj[i], col = tip.color[i])
            }
        }
        if (show.node.label) 
            text(xx[ROOT:length(xx)] + label.offset, yy[ROOT:length(yy)], 
                x$node.label, adj = adj, font = font, srt = srt, 
                cex = cex)
    }
    L <- list(type = type, use.edge.length = use.edge.length, 
        node.pos = node.pos, show.tip.label = show.tip.label, 
        show.node.label = show.node.label, font = font, cex = cex, 
        adj = adj, srt = srt, no.margin = no.margin, label.offset = label.offset, 
        x.lim = x.lim, y.lim = y.lim, direction = direction, 
        tip.color = tip.color, Ntip = Ntip, Nnode = Nnode, edge = xe, xx = xx, yy = yy)
    assign("last_plot.phylo", c(L, list(edge = xe, xx = xx, yy = yy)), 
        envir = .PlotPhyloEnv)
    invisible(L)
}

