#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use File::Copy;
use File::Basename;
use FindBin qw($Bin);
my $r_path = "$Bin/../../bin/R/bin";
$ENV{'PATH'} = "$r_path;$ENV{'PATH'}";

my %opts;
GetOptions (\%opts,"i=s","d=s","w=i","h=i","cex=f","eg=f","subtre=s","tretype=s","coltype=s","o=s","sl=s","ma=s","laboff=f","rule=f", "lcex=f", "html=s");

my $usage = <<"USAGE";
       Discription:   
       Usage :perl $0 [options]
				-i	newwick format tree file.
				-o	output pdf file name.
				-d	lable design file.
				-w	width  defalt: 5
				-h	heigth defalt:6
				-cex	label cex default:1
                                -lcex       legend cex default:1
				-eg	edge width default:1
				-tretype	tree type:["phylogram", "fan", "radial","cladogram"]. default:phylogram.
				-coltype	tree branch color type.[clu/tip], default: clu.
				-subtre	[T/F] output subtre file if lables in design file less than input treefile. default :F.
				 -sl [T/F]   show labels or not defalt:T.
				 -ma [T/F]   make branch length the same.defalt:F.
				 -laboff	label offset.default:0.0002.
				 -rule	ruler length.default:0.01

USAGE

die $usage if ( !($opts{i}&&$opts{o}));


$opts{h}=defined $opts{h}?$opts{h}:6;
$opts{w}=defined $opts{w}?$opts{w}:5;
$opts{cex}=defined $opts{cex}?$opts{cex}:1;
$opts{lcex}=defined $opts{lcex}?$opts{lcex}:1;
$opts{eg}=defined $opts{eg}?$opts{eg}:1;
$opts{d}=defined $opts{d}?$opts{d}:"none";
$opts{tretype}=defined $opts{tretype}?$opts{tretype}:"phylogram";
$opts{coltype}=defined $opts{coltype}?$opts{coltype}:"clu";
$opts{subtre}=defined $opts{subtre}?$opts{subtre}:"F";
$opts{gl}=defined $opts{gl}?$opts{gl}:"F";
$opts{sl}=defined $opts{sl}?$opts{sl}:"T";
$opts{ma}=defined $opts{ma}?$opts{ma}:"F";
$opts{laboff}=defined $opts{laboff}?$opts{laboff}:0.0002;
$opts{rule}=defined $opts{rule}?$opts{rule}:0.01;

copy( $opts{i}, '__input__' );
$opts{i} = '__input__';

if( $opts{d} ne "none" ){
	copy( $opts{d}, '__design__' );
	$opts{d} = '__design__';
}

my $cmd="
treefile=\"$opts{i}\"
design =\"$opts{d}\"
tretype=\"$opts{tretype}\"
coltype=\"$opts{coltype}\"
subtre=\"$opts{subtre}\"

#mycol <-c(\"#000000\",\"#3A89CC\",\"#769C30\",\"#CD0000\",\"#D99536\",\"#7B0078\",\"#BFBC3B\",\"#6E8B3D\",\"#00688B\",\"#C10077\",\"#CAAA76\",\"#EEEE00\",\"#458B00\",\"#8B4513\",\"#008B8B\",\"#6E8B3D\",\"#8B7D6B\",\"#7FFF00\",\"#CDBA96\",\"#ADFF2F\")
mycol <- c(34, 51, 27, 31,142, 23, 50, 75, 525, 62, 119, 46, 475, 554, 622, 483, 657, 545, 402, 477, 503, 40, 115, 5, 376,473,546,482)
mycol <-colors()[rep(mycol,20)]
library(ape)

plot.phylo.yg <-function (x, type = \"phylogram\", use.edge.length = TRUE, node.pos = NULL, 
    show.tip.label = TRUE, show.node.label = FALSE, edge.color = \"black\",gline=TRUE, 
    edge.width = 1, edge.lty = 1, font = 1, cex = par(\"cex\"), magic =FALSE,
    adj = NULL, srt = 0, no.margin = FALSE, root.edge = FALSE, legend =FALSE,
    label.offset = 0, underscore = TRUE, x.lim = NULL, y.lim = NULL, 
    direction = \"rightwards\", lab4ut = \"horizontal\", tip.color = \"black\", 
    plot = TRUE, ...) 
{



####phylogram.yg.plot##################

phylogram.yg.plot <-function (edge, Ntip, Nnode, xx, yy, horizontal, edge.color, 
    edge.width, edge.lty,gline) 
{#print(gline)
    nodes <- (Ntip + 1):(Ntip + Nnode)
    if (!horizontal) {
        tmp <- yy
        yy <- xx
        xx <- tmp
    }
#print(xx)
#print(yy)
#points(xx,yy,col=\"darkblue\")
    x0v <- xx[nodes]
    y0v <- y1v <- numeric(Nnode)
    NodeInEdge1 <- vector(\"list\", Nnode)
    for (i in nodes) {
        ii <- i - Ntip
        j <- NodeInEdge1[[ii]] <- which(edge[, 1] == i)
        tmp <- range(yy[edge[j, 2]])
        tmpx <- max(xx[edge[j, 2]])
        y0v[ii]  <- tmp[1]
        y1v[ii]  <- tmp[2] 

    }

### gline #######
    lab.width <- rep(0,Ntip)
    if(show.tip.label){
         lab.width <- strwidth(x\$tip.label,\"inches\",adj = adj, font = font, srt = srt, cex = cex)*1.01        
    }else{ lab.width<-0.1 }

if(gline){  
    gnod <- vector(\"list\")
    i <- Ntip+1
    g <- 0
    gno <-c()
    while (i <= length(tip.color)) {
        if(tip.color[i]>1){ g <-g+1; #print(i);gno <-c(gno,i)
                 gnod[[g]] <- i          
             while(any(gnod[[g]]> Ntip)){
                 gs <- which(gnod[[g]] <=Ntip)
                 gl <- which(gnod[[g]] >Ntip)
                 xgs <-gnod[[g]][gs];
                 xgl <-gnod[[g]][gl]; 
                 gend <-max(xgl)
                 j <-c() 
                 for(m in 1:length(xgl)){
                     j <- c(j,which(edge[, 1] == xgl[m]))
                 }                
                 xj <- edge[j,2];
                 gnod[[g]] <- c(xgs,xj)
             }
             i <- gend+1                              
        }else { i <- i+1 }                 
    }
#print(gnod)
#print(gno)
#print(g)
color.l <- tip.color[gno]
#print(color.l)
x0l <-y0l <-y1l <- numeric(g)
os <- 0.4*(strwidth(\"a\",\"inches\"))
  if(g >0){
    for(i in 1:g){
         x0l[i] <-max(xx[gnod[[i]]])+max(lab.width[gnod[[i]]])+3*os
         y0l[i] <-yy[min(gnod[[i]])]-os
         y1l[i] <-yy[max(gnod[[i]])]+os
    }
  }
}else {color.l <- \"white\"
y0l<-x0l<- y1l <-x0l <-0
}

############

    x0h <- xx[edge[, 1]]
    x1h <- xx[edge[, 2]]
    y0h <- yy[edge[, 2]]

    nc <- length(edge.color)
    nw <- length(edge.width)
    nl <- length(edge.lty)
    if (nc + nw + nl == 3) {
        color.v  <-edge.color
        width.v <- edge.width
        lty.v <- edge.lty
    }
    else {
        Nedge <- dim(edge)[1]
        edge.color <- rep(edge.color, length.out = Nedge)
        edge.width <- rep(edge.width, length.out = Nedge)
        edge.lty <- rep(edge.lty, length.out = Nedge)
        DF <- data.frame(edge.color, edge.width, edge.lty, stringsAsFactors = FALSE)
        #print(DF)
        color.v <- rep(1, Nnode)
        width.v <- rep(1, Nnode)
        lty.v <- rep(1, Nnode)
        for (i in 1:Nnode) {
            br <- NodeInEdge1[[i]]
            #print (NodeInEdge1[[i]])

            if (length(br) > 2) {
                x <- unique(DF[br, 1])
                if (length(x) == 1) { 
                  color.v[i] <-x 
                
                 
                }
                x <- unique(DF[br, 2])
                if (length(x) == 1) 
                  width.v[i] <- x
                x <- unique(DF[br, 3])
                if (length(x) == 1) 
                  lty.v[i] <- x
            }
            else {
                A <- br[1]
                B <- br[2]
                if (any(DF[A, ] != DF[B, ])) {
                  #color.v[i] <- edge.color[B]
                  color.v[i] <- 1  
                  width.v[i] <- edge.width[B]
                  lty.v[i] <- edge.lty[B]
                  y0v <- c(y0v, y0v[i])
                  y1v <- c(y1v, yy[i + Ntip])
                  x0v <- c(x0v, x0v[i])
                  #color.v <- c(color.v, edge.color[A])
                  color.v <- c(color.v,color.v[i] )
                  width.v <- c(width.v, edge.width[A])
                  lty.v <- c(lty.v, edge.lty[A])
                  y0v[i] <- yy[i + Ntip]
                }
                else {
                  color.v[i] <- edge.color[A]
                  width.v[i] <- edge.width[A]
                  lty.v[i] <- edge.lty[A]
                }
            }
        }
    }

#print(x0h)
#print(x1h)

    if (horizontal) {
        segments(x0h, y0h, x1h, y0h, col = edge.color, lwd = edge.width, 
            lty = edge.lty)
        segments(x0v, y0v, x0v, y1v, col = color.v, lwd = width.v, 
            lty = lty.v)
        segments(x0l, y0l, x0l, y1l, col = color.l, lwd = 2, 
            lty = lty.v)
    }
    else {
        segments(y0h, x0h, y0h, x1h, col = edge.color, lwd = edge.width, 
            lty = edge.lty)
        segments(y0v, x0v, y1v, x0v, col = color.v, lwd = width.v, 
            lty = lty.v)
        segments(y0l, x0l, y1l, x0l, col = color.l, lwd = 2, 
            lty = lty.v)
    }
}
#<environment: namespace:ape>

########phylogram.yg.plot##################end



    Ntip <- length(x\$tip.label)
    if (Ntip == 1) {
        warning(\"found only one tip in the tree\")
        return(NULL)
    }
    if (any(tabulate(x\$edge[, 1]) == 1)) 
        stop(\"there are single (non-splitting) nodes in your tree; you may need to use collapse.singles()\")
    .nodeHeight <- function(Ntip, Nnode, edge, Nedge, yy) .C(\"node_height\", 
        as.integer(Ntip), as.integer(Nnode), as.integer(edge[, 
            1]), as.integer(edge[, 2]), as.integer(Nedge), as.double(yy), 
        DUP = FALSE, PACKAGE = \"ape\")[[6]]
    .nodeDepth <- function(Ntip, Nnode, edge, Nedge) .C(\"node_depth\", 
        as.integer(Ntip), as.integer(Nnode), as.integer(edge[, 
            1]), as.integer(edge[, 2]), as.integer(Nedge), double(Ntip + 
            Nnode), DUP = FALSE, PACKAGE = \"ape\")[[6]]
    .nodeDepthEdgelength <- function(Ntip, Nnode, edge, Nedge, 
        edge.length) .C(\"node_depth_edgelength\", as.integer(Ntip), 
        as.integer(Nnode), as.integer(edge[, 1]), as.integer(edge[, 
            2]), as.integer(Nedge), as.double(edge.length), double(Ntip + 
            Nnode), DUP = FALSE, PACKAGE = \"ape\")[[7]]
    Nedge <- dim(x\$edge)[1]
    Nnode <- x\$Nnode
    ROOT <- Ntip + 1
    type <- match.arg(type, c(\"phylogram\", \"cladogram\", \"fan\", 
        \"unrooted\", \"radial\"))
    direction <- match.arg(direction, c(\"rightwards\", \"leftwards\", 
        \"upwards\", \"downwards\"))
    if (is.null(x\$edge.length)) 
        use.edge.length <- FALSE
    if (type \%in\% c(\"unrooted\", \"radial\") || !use.edge.length || 
        is.null(x\$root.edge) || !x\$root.edge) 
        root.edge <- FALSE
    if (type == \"fan\" && root.edge) {
        warning(\"drawing root edge with type = 'fan' is not yet supported\")
        root.edge <- FALSE
    }
    phyloORclado <- type \%in\% c(\"phylogram\", \"cladogram\")
    horizontal <- direction \%in\% c(\"rightwards\", \"leftwards\")
    xe <- x\$edge
    if (phyloORclado) {
        phyOrder <- attr(x, \"order\")
        if (is.null(phyOrder) || phyOrder != \"cladewise\") {
            x <- reorder(x)
            if (!identical(x\$edge, xe)) {
                ereorder <- match(x\$edge[, 2], xe[, 2])
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
        TIPS <- x\$edge[x\$edge[, 2] <= Ntip, 2]
        yy[TIPS] <- 1:Ntip
    }
    z <- reorder(x, order = \"pruningwise\")
    if (phyloORclado) {
        if (is.null(node.pos)) {
            node.pos <- 1
            if (type == \"cladogram\" && !use.edge.length) 
                node.pos <- 2
        }
        if (node.pos == 1) 
            yy <- .nodeHeight(Ntip, Nnode, z\$edge, Nedge, yy)
        else {
            ans <- .C(\"node_height_clado\", as.integer(Ntip), 
                as.integer(Nnode), as.integer(z\$edge[, 1]), as.integer(z\$edge[, 
                  2]), as.integer(Nedge), double(Ntip + Nnode), 
                as.double(yy), DUP = FALSE, PACKAGE = \"ape\")
            xx <- ans[[6]] - 1
            yy <- ans[[7]]
        }
        if (!use.edge.length) {
            if (node.pos != 2) 
                xx <- .nodeDepth(Ntip, Nnode, z\$edge, Nedge) - 
                  1
            xx <- max(xx) - xx
        }
        else {
            xx <- .nodeDepthEdgelength(Ntip, Nnode, z\$edge, Nedge, 
                z\$edge.length)
        }
    }
    else switch(type, fan = {
        TIPS <- x\$edge[which(x\$edge[, 2] <= Ntip), 2]
        xx <- seq(0, 2 * pi * (1 - 1/Ntip), 2 * pi/Ntip)
        theta <- double(Ntip)
        theta[TIPS] <- xx
        theta <- c(theta, numeric(Nnode))
        theta <- .nodeHeight(Ntip, Nnode, z\$edge, Nedge, theta)
        if (use.edge.length) {
            r <- .nodeDepthEdgelength(Ntip, Nnode, z\$edge, Nedge, 
                z\$edge.length)
        } else {
            r <- .nodeDepth(Ntip, Nnode, z\$edge, Nedge)
            r <- 1/r
        }
				if(magic){	
					r[1:Ntip] <-max(r)				
				}		
        xx <- r * cos(theta)
        yy <- r * sin(theta)		
    }, unrooted = {
        nb.sp <- .nodeDepth(Ntip, Nnode, z\$edge, Nedge)
        XY <- if (use.edge.length) unrooted.xy(Ntip, Nnode, z\$edge, 
            z\$edge.length, nb.sp) else unrooted.xy(Ntip, Nnode, 
            z\$edge, rep(1, Nedge), nb.sp)
        xx <- XY\$M[, 1] - min(XY\$M[, 1])
        yy <- XY\$M[, 2] - min(XY\$M[, 2])
    }, radial = {
        X <- .nodeDepth(Ntip, Nnode, z\$edge, Nedge)
        X[X == 1] <- 0
        X <- 1 - X/Ntip
        yy <- c((1:Ntip) * 2 * pi/Ntip, rep(0, Nnode))
        Y <- .nodeHeight(Ntip, Nnode, z\$edge, Nedge, yy)
        xx <- X * cos(Y)
        yy <- X * sin(Y)
    })
    if (phyloORclado) {
        if (!horizontal) {
            tmp <- yy
            yy <- xx
            xx <- tmp - min(tmp) + 1
        }
        if (root.edge) {
            if (direction == \"rightwards\") 
                xx <- xx + x\$root.edge
            if (direction == \"upwards\") 
                yy <- yy + x\$root.edge
        }
    }
    if (no.margin) 
        par(mai = rep(0, 4))
    if (is.null(x.lim)) {
        if (phyloORclado) {
            if (horizontal) {
                x.lim <- c(0, NA)
                pin1 <- par(\"pin\")[1]
                strWi <- strwidth(x\$tip.label, \"inches\")
                xx.tips <- xx[1:Ntip] * 1.05
                alp <- try(uniroot(function(a) max(a * xx.tips + 
                  strWi) - pin1, c(0, 1e+06))\$root, silent = TRUE)
                if (is.character(alp)) 
                  tmp <- max(xx.tips) * 1.5
                else {
                  tmp <- if (show.tip.label) 
                    max(xx.tips + strWi/alp*1.05)
                  else max(xx.tips)
                }
                x.lim[2] <- tmp
				if(legend) x.lim[2] <-x.lim[2]*1.25
            }
            else x.lim <- c(1, Ntip)
        }
        else switch(type, fan = {
            if (show.tip.label) {
                offset <- max(nchar(x\$tip.label) * 0.025 * max(yy) * 
                  cex)
                x.lim <- c(min(xx) - offset, max(xx) + offset)
				if(legend){ x.lim <- c(min(xx) - offset, max(xx) *1.25+ offset)}
            } else x.lim <- c(min(xx), max(xx))
        }, unrooted = {
            if (show.tip.label) {
                offset <- max(nchar(x\$tip.label) * 0.025 * max(yy) * 
                  cex)
                x.lim <- c(0 - offset, max(xx) + offset*3)
            } else x.lim <- c(0, max(xx))
        }, radial = {
            if (show.tip.label) {
                offset <- max(nchar(x\$tip.label) * 0.03 * cex)
                x.lim <- c(-1 - offset, 1 + offset*3)
            } else x.lim <- c(-1, 1)
        })
    }
    else if (length(x.lim) == 1) {
        x.lim <- c(0, x.lim)
        if (phyloORclado && !horizontal) 
            x.lim[1] <- 1
        if (type \%in\% c(\"fan\", \"unrooted\") && show.tip.label) 
            x.lim[1] <- -max(nchar(x\$tip.label) * 0.025 * max(yy) * 
                cex)
        if (type == \"radial\") 
            x.lim[1] <- if (show.tip.label) 
                -1 - max(nchar(x\$tip.label) * 0.03 * cex)
            else -1
    }
    if (phyloORclado && direction == \"leftwards\") 
        xx <- x.lim[2] - xx
    if (is.null(y.lim)) {
        if (phyloORclado) {
            if (horizontal) 
                y.lim <- c(0.8, Ntip*1.05)
            else {
                y.lim <- c(0.8, NA)
                pin2 <- par(\"pin\")[2]
                strWi <- strwidth(x\$tip.label, \"inches\")
                yy.tips <- yy[1:Ntip] * 1.05
                alp <- try(uniroot(function(a) max(a * yy.tips + 
                  strWi) - pin2, c(0, 1e+06))\$root, silent = TRUE)
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
                offset <- max(nchar(x\$tip.label) * 0.10 * max(yy) * 
                  cex)
                y.lim <- c(min(yy) - offset, max(yy) + offset)
            } else y.lim <- c(min(yy), max(yy))
        }, unrooted = {
            if (show.tip.label) {
                offset <- max(nchar(x\$tip.label) * 0.025 * max(yy) * 
                  cex)
                y.lim <- c(0 - offset, max(yy) + offset)
            } else y.lim <- c(0, max(yy))
        }, radial = {
            if (show.tip.label) {
                offset <- max(nchar(x\$tip.label) * 0.03 * cex)
                y.lim <- c(-1 - offset, 1 + offset)
            } else y.lim <- c(-1, 1)
        })
    }
    else if (length(y.lim) == 1) {
        y.lim <- c(0, y.lim)
        if (phyloORclado && horizontal) 
            y.lim[1] <- 1
        if (type \%in\% c(\"fan\", \"unrooted\") && show.tip.label) 
            y.lim[1] <- -max(nchar(x\$tip.label) * 0.05 * max(yy) * 
                cex)
        if (type == \"radial\") 
            y.lim[1] <- if (show.tip.label) 
                -1 - max(nchar(x\$tip.label) * 0.05 * max(yy) * 
                  cex)
            else -1
    }
    if (phyloORclado && direction == \"downwards\") 
        yy <- y.lim[2] - yy
    if (phyloORclado && root.edge) {
        if (direction == \"leftwards\") 
            x.lim[2] <- x.lim[2] + x\$root.edge
        if (direction == \"downwards\") 
            y.lim[2] <- y.lim[2] + x\$root.edge
    }
    asp <- if (type \%in\% c(\"fan\", \"radial\", \"unrooted\")) 
        1
    else NA
    plot(0, type = \"n\", xlim = x.lim, ylim = y.lim, ann = FALSE, 
        axes = FALSE, asp = asp, ...)
    if (plot) {
        if (is.null(adj)) 
            adj <- if (phyloORclado && direction == \"leftwards\") 
                1
            else 0
        if (phyloORclado && show.tip.label) {
            MAXSTRING <- max(strwidth(x\$tip.label, cex = cex))
            loy <- 0
            if (direction == \"rightwards\") {
                lox <- label.offset + MAXSTRING * 1.05 * adj
            }
            if (direction == \"leftwards\") {
                lox <- -label.offset - MAXSTRING * 1.05 * (1 - 
                  adj)
            }
            if (!horizontal) {
                psr <- par(\"usr\")
                MAXSTRING <- MAXSTRING * 1.09 * (psr[4] - psr[3])/(psr[2] - 
                  psr[1])
                loy <- label.offset + MAXSTRING * 1.05 * adj
                lox <- 0
                srt <- 90 + srt
                if (direction == \"downwards\") {
                  loy <- -loy
                  srt <- 180 + srt
                }
            }
        }

################################################
        

#####################
        
        if (type == \"phylogram\") { 
				gline <-gline
				if(magic ){
						if(direction == \"rightwards\") {xx[1:Ntip] <- max(xx[1:Ntip])}
						if(direction == \"downwards\") {yy[1:Ntip] <- min(yy[1:Ntip])} 
				}					
            phylogram.yg.plot(x\$edge, Ntip, Nnode, xx, yy, horizontal, 
                edge.color, edge.width, edge.lty,gline)
        }
        else {
            if (type == \"fan\") {
                ereorder <- match(z\$edge[, 2], x\$edge[, 2])
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
                circular.plot(z\$edge, Ntip, Nnode, xx, yy, theta, 
                  r, edge.color, edge.width, edge.lty)
            }
            else cladogram.plot(x\$edge, xx, yy, edge.color, edge.width, 
                edge.lty)
        }
        if (root.edge) 
            switch(direction, rightwards = segments(0, yy[ROOT], 
                x\$root.edge, yy[ROOT]), leftwards = segments(xx[ROOT], 
                yy[ROOT], xx[ROOT] + x\$root.edge, yy[ROOT]), 
                upwards = segments(xx[ROOT], 0, xx[ROOT], x\$root.edge), 
                downwards = segments(xx[ROOT], yy[ROOT], xx[ROOT], 
                  yy[ROOT] + x\$root.edge))
        if (show.tip.label) {
            if (is.expression(x\$tip.label)) 
                underscore <- TRUE
            if (!underscore) 
                x\$tip.label <- gsub(\"_\", \" \", x\$tip.label)
            if (phyloORclado) 
                text(xx[1:Ntip] + lox, yy[1:Ntip] + loy, x\$tip.label, 
                  adj = adj, font = font, srt = srt, cex = cex, 
                  col = tip.color)
#print(xx)
#print(yy)
#points(xx[1:Ntip] + lox, yy[1:Ntip] + loy,col=\"green\",pch=2)

            if (type == \"unrooted\") {
                if (lab4ut == \"horizontal\") {
                  y.adj <- x.adj <- numeric(Ntip)
                  sel <- abs(XY\$axe) > 0.75 * pi
                  x.adj[sel] <- -strwidth(x\$tip.label)[sel] * 
                    1.05
                  sel <- abs(XY\$axe) > pi/4 & abs(XY\$axe) < 0.75 * 
                    pi
                  x.adj[sel] <- -strwidth(x\$tip.label)[sel] * 
                    (2 * abs(XY\$axe)[sel]/pi - 0.5)
                  sel <- XY\$axe > pi/4 & XY\$axe < 0.75 * pi
                  y.adj[sel] <- strheight(x\$tip.label)[sel]/2
                  sel <- XY\$axe < -pi/4 & XY\$axe > -0.75 * pi
                  y.adj[sel] <- -strheight(x\$tip.label)[sel] * 
                    0.75
                  text(xx[1:Ntip] + x.adj * cex, yy[1:Ntip] + 
                    y.adj * cex, x\$tip.label, adj = c(adj, 0), 
                    font = font, srt = srt, cex = cex, col = tip.color)
                }
                else {
                  adj <- abs(XY\$axe) > pi/2
                  srt <- 180 * XY\$axe/pi
                  srt[adj] <- srt[adj] - 180
                  adj <- as.numeric(adj)
                  xx.tips <- xx[1:Ntip]
                  yy.tips <- yy[1:Ntip]
                  if (label.offset) {
                    xx.tips <- xx.tips + label.offset * cos(XY\$axe)
                    yy.tips <- yy.tips + label.offset * sin(XY\$axe)
                  }
                  font <- rep(font, length.out = Ntip)
                  tip.color <- rep(tip.color, length.out = Ntip)
                  cex <- rep(cex, length.out = Ntip)
                  for (i in 1:Ntip) text(xx.tips[i], yy.tips[i], 
                    cex = cex[i], x\$tip.label[i], adj = adj[i], 
                    font = font[i], srt = srt[i], col = tip.color[i])
                }
            }
            if (type \%in\% c(\"fan\", \"radial\")) {
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
                  x\$tip.label[i], font = font[i], cex = cex[i], 
                  srt = angle[i], adj = adj[i], col = tip.color[i])
            }
        }
        if (show.node.label) 
            text(xx[ROOT:length(xx)] + label.offset, yy[ROOT:length(yy)], 
                x\$node.label, adj = adj, font = font, srt = srt, 
                cex = cex)
    }

                
    L <- list(type = type, use.edge.length = use.edge.length, 
        node.pos = node.pos, show.tip.label = show.tip.label, 
        show.node.label = show.node.label, font = font, cex = cex, 
        adj = adj, srt = srt, no.margin = no.margin, label.offset = label.offset, 
        x.lim = x.lim, y.lim = y.lim, direction = direction, 
        tip.color = tip.color, Ntip = Ntip, Nnode = Nnode)
    assign(\"last_plot.phylo\", c(L, list(edge = xe, xx = xx, yy = yy)), 
        envir = .PlotPhyloEnv)
    invisible(L)
}
#<environment: namespace:ape>



tree <-read.tree(file=treefile)
egcol=\"#000000\"
labcol=\"#000000\"

####### read in a Newick  Tree ########
design_tre <-function(tree,design=\"none\",tiplen=0,coltype=\"clu\"){
		newtree <- tree
		if(design!=\"none\"){
						de <-read.table(design,sep=\"\\t\",fill=T, colClasses=c(\"character\",\"character\"))
						labs=de[,1]						
						de <-de[,2]	
						names(de)<-labs
						class_count <-as.matrix(table(de))
						class_color  <-mycol[1:(length(class_count))]
						class <-data.frame(count=class_count,color=as.character(class_color))
									
						
						## select labs ##########
						nlab <- which(tree\$tip.label \%in\% labs )
						tnlab <-which(tree\$tip.label \%in\% setdiff(tree\$tip.label,labs))
						if(length(tnlab)>0){
							   tnedge <-which(tree\$edge[,2] \%in\% tnlab)
							   eddge <- tree\$edge
							   tchop <-tnedge
							   s1 <-0
							   s2 <-length(tchop)
							   while(s2-s1>0){
							   s1 <-length(tchop)
							   chop <- c(tnlab,names(which(table(eddge[tchop,1])==2)))
							   tchop <- which(eddge[,2] \%in\% chop)
							   s2 <-length(tchop)
							   }
							   newtree\$edge <-newtree\$edge[-tchop,]
							   newtree\$edge.length <-newtree\$edge.length[-tchop] 
							   lastnlabs <-which(newtree\$edge[,2]<=length(newtree\$tip.label))
							   newtree\$edge[lastnlabs,2] <- sapply(newtree\$edge[lastnlabs,2],function(x) which(nlab \%in\% x)) 
							   newtree\$tip.label <-newtree\$tip.label[nlab]
							   newtree <-collapse.singles(newtree)
							   eg <-c(1:length(unique(newtree\$edge[,1])))+length(nlab)
							   names(eg) <-unique(sort(newtree\$edge[,1]))
							   newtree\$edge[,1] <-as.vector(sapply(as.character(newtree\$edge[,1]),function(x) eg[[x]]))
							   ed2 <- which(newtree\$edge[,2]>length(newtree\$tip.label))
							   newtree\$edge[ed2,2] <-as.vector(sapply(as.character(newtree\$edge[ed2,2]),function(x) eg[[x]]))
							   newtree\$Nnode <-max(newtree\$edge[,1])-min(newtree\$edge[,1])+1
						}
						##################
						
						## set branch colors ##########											
						labcol <-vector(\"list\")
						for (i in 1:length(newtree\$tip.label)){
								labcol[[newtree\$tip.label[i]]] <- class\$color[which(rownames(class) \%in\% de[newtree\$tip.label[i]])]			
						}
						
						if(coltype==\"tip\"){
								for(i in 1:(nrow(newtree\$edge))){
									egcol[i]=\"#000000\"
									  n <-newtree\$edge[i,2]
									  if(n <=length(newtree\$tip.label)){
											egcol[i]=as.character(labcol[[newtree\$tip.label[n]]])			 
									  }
								}
								labcol <-as.vector(unlist(labcol))
						}else if(coltype==\"clu\"){
								labcol <-as.vector(unlist(labcol))
								xt <-newtree
								Ntip <- length(xt\$tip.label)
								Nnode <- xt\$Nnode
								xedge <- cbind(1:nrow(xt\$edge),xt\$edge)
								xeo <- xedge[order(xedge[,2],decreasing = T),]
								xlb <- xedge[order(xedge[,3],decreasing = F),]
								m <-j <-max(xeo[,2])
								while(j>Ntip+1){
									   i=2*(m-j)+1
									   if(labcol[xeo[i,3]]==labcol[xeo[i+1,3]]){
											labcol[xeo[i,2]]=labcol[xeo[i,3]]
									   }else{
											labcol[xeo[i,2]]=\"#000000\"
									   }
									   j <-j-1
								}
								i <- i+2
								root <-xeo[i:nrow(xeo),]
								if( any(labcol[root[,3]] != labcol[root[1,3]]) ){
									 labcol[root[1,2]]=1 
								}else{
											labcol[root[1,2]] <-labcol[root[1,3]] 
								}
								for(i in 1:(nrow(xt\$edge))){
									  n <-xt\$edge[i,2]
									  egcol[i]= labcol[n] 
								}
						}
						#labcol <-as.vector(unlist(labcol))
						####################
		}
		
		#make the branch longger or shorter
		for(i in 1:(nrow(newtree\$edge))){
					mal <-max(newtree\$edge.length)
					if(tiplen>=1){
						 newtree\$edge.length[i] <-newtree\$edge.length[i]+(tiplen-1)*mal 
					}else if(tiplen>0){ 
						 newtree\$edge.length[i] <-newtree\$edge.length[i]*tiplen
					} 
		}
		##########################		
		ntr <- list(newtree=newtree,labcol=labcol,egcol=egcol,class=class)
}

	

ntre <-design_tre(tre=tree,design=design,coltype=coltype,)
tree <- ntre\$newtree
egcol <- ntre\$egcol
labcol <- ntre\$labcol

if(subtre==\"T\"){
write.tree(tree,\"__select.tre__\")
}

pdf(\"__output__\",width=$opts{w},height=$opts{h})
legend=FALSE

plot.phylo.yg(tree,edge.color=egcol,edge.width=$opts{eg},cex=$opts{cex},magic=$opts{ma},show.tip.label = $opts{sl},tip.color=labcol,label.offset=$opts{laboff},type=tretype,legend=legend)
axis(side=1,at=c(0,$opts{rule}),labels=F,tcl=0.2)
axis(side=1,at=$opts{rule}/2,labels=\"$opts{rule}\",tick=F,cex.axis=$opts{cex},mgp=c(0, 0, 0))

if(design!=\"none\"){
if(length(rownames(ntre\$class))>1){    
#plot.new()
#par(mar=c(0,0,0,0))
legend(\"topleft\",legend=rownames(ntre\$class),col=as.character(ntre\$class\$color),fill=as.character(ntre\$class\$color),ncol=1,cex=$opts{lcex},bty=\"n\")
}
}
dev.off()
";
open CMD ,">trecmd.r";
print CMD $cmd;
close CMD;


my $ret = system("R --restore --no-save < trecmd.r");
if ( $ret ) {
        die "Error, died with $ret";
}

unlink 'trecmd.r';

unlink $opts{i};
if ( $opts{d} ne 'none'){
	unlink $opts{d};
}

move( '__output__', $opts{o} );


&get_html_link_v2($opts{html}, "Plot Tree Results", $opts{o}, "pdf");


use File::Basename;
use File::Copy;
#tabH, tabN, txt
sub get_html_link_v2{
	my @param = @_;
	my $html = shift @param;
	my $title = shift @param;
	my @par = @param;

    my $out_dir = "$html.files";
    mkdir $out_dir;


	open HTML, ">$html";
	print HTML "<html><head><title>$title</title></head><body><h3>Output Files:</h3><p><ul>\n";

    for(my $i = 0; $i <= $#par; $i+=2){
            my $file = $par[$i];
            my $type = $par[$i+1];
            my $name = basename $file;
            if( $type eq "tabH" ){
                get_html_table($file, "$out_dir/$name.html", "y");
            }elsif(  $type eq "tabN" ){
                get_html_table($file, "$out_dir/$name.html", "n");
            }elsif(  $type eq "txt" ){
                txt2html($file, "$out_dir/$name.html");
            }else{
                copy($file, "$out_dir/$name");
            }
    }

    my $basename = basename $out_dir;
	for (my $i = 0; $i <= $#par; $i+=2){
            my $file = $par[$i];
            my $type = $par[$i+1];
            my $name = basename $file;
            if( $type eq "tabH" ){
                print HTML "<li><a href=\"$basename/$name.html\" target=\"_parent\">$name</a></li>\n";
            }elsif(  $type eq "tabN" ){
                print HTML "<li><a href=\"$basename/$name.html\" target=\"_parent\">$name</a></li>\n";
            }elsif(  $type eq "txt" ){
                print HTML "<li><a href=\"$basename/$name.html\" target=\"_parent\">$name</a></li>\n";
            }else{
                print HTML "<li><a href=\"$basename/$name\" target=\"_parent\">$name</a></li>\n";
            }
	        
	}
	print HTML "</ul></p>\n";


}



sub get_html_table{
	my $in = shift;
	my $out = shift;
	my $isHeader = shift;
	open IN,"$in";
	open OUT,">$out";
    print OUT "<!DOCTYPE html>\n";
    print OUT "<html>\n";
	print OUT "<table border=\"1\">\n";
	if ( $isHeader eq "y") {
		my $h = <IN>;
		chomp $h;
		my @s = split /\t/, $h;
		print OUT "<tr align=\"center\">";
		for(@s){
			print OUT "<th>$_</th>";
		}
		print OUT "</tr>"
	}

	while (<IN>) {
		chomp;
		my @s = split /\t/, $_;
		print OUT "<tr align=\"center\">";
		for(@s){
			print OUT "<td>$_</td>";
		}
		print OUT "</tr>"
	}

	print OUT "</table>";
    print OUT "</html>";
}




sub txt2html{
	my $txt = shift;
	my $html = shift;
	open IN,"$txt";
	open OUT,">$html";
    print OUT "<!DOCTYPE html>\n";
    print OUT "<html>";
	print OUT "<body>\n";
	while (<IN>) {
		chomp;
		print OUT "$_<br />\n";
	}

	print OUT "</body>\n";
    print OUT "</html>\n";
	close IN;
	close OUT;
}


sub getRTable{
    my $in = shift;
    my $out = shift;
    open IN,"$in";
    open OUT,">$out";
    my $h = <IN>;
    print OUT "\t$h";
    while (<IN>){
          print OUT "$_";
    }
    close IN;
    close OUT;
}