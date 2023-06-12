#' Draw 2-d probability simplex, with model lines for T3 or T1 model
#'
#' Outline the 2-d probability simplex, and draw the T1 or T3 model points for quartet frequencies.
#' The models "T1" and "T3" are described more fully by \insertCite{MAR19;textual}{MSCquartets}.
#'
#' @references
#' \insertRef{MAR19}{MSCquartets}
#' 
#' @param model \code{"T1"} or \code{"T3"}, for 1-tree or 3-tree model
#' @param maintitle main title for plot
#' @param titletext additional text for title
#' @return NULL
#'
#' @examples
#'    simplexPrepare("T3",maintitle="Main title",titletext="further text")
#'
#' @seealso \code{\link{simplexLabels}},
#'          \code{\link{simplexPoint}},
#'          \code{\link{simplexSegment}},
#'          \code{\link{simplexText}},
#'          \code{\link{simplexCoords}}
#'
#' @export
simplexPrepare <- function(model = "T3", 
                           maintitle = NULL, 
                           titletext = NULL) {
  if (!(model %in% c("T1", "T3", "cut")))
    stop("Invalid model name; use 'T1','T3', or 'cut'.") #check parameters
  lineWidth = 2
  
  oldpar=par(mar = c(0, 0, 4, 0) + 0.1)# set margin
  on.exit(par(oldpar))# and restore when done         
  
  top = c(1, 0, 0) # Simplex vertices
  left = c(0, 1, 0)
  right = c(0, 0, 1)
  mid = rep(1, 3) #Simplex barycenter
  lim <- c(-1, 1)
  yylim <- c(-.1, .6)
  
  plot(
    0,
    0,
    xlim = lim,
    ylim = yylim,
    type = 'n',
    asp = 1,
    axes = FALSE,
    xlab = '',
    ylab = '',
    main = maintitle
  ) # create plot
  mtext(eval(bquote(.(titletext))),
        side = 3,
        line = 0,
        cex = 1.2) # add rest of title text
  simplexSegment(top, left, lty = 'dotted', lwd = lineWidth) # outline Simplex
  simplexSegment(top, right, lty = 'dotted', lwd = lineWidth)
  simplexSegment(left, right, lty = 'dotted', lwd = lineWidth)
  simplexSegment(top, mid, lwd = lineWidth) # Plot lines giving model for quartet probabilities under the coalsecent on a species tree
  if (model == "T3") {
    # for 3-tree model use 3 lines stopping at centroid
    simplexSegment(left, mid, lwd = lineWidth)
    simplexSegment(right, mid, lwd = lineWidth)
  }
  if (model =="cut") {
    # for cut model use 3 lines through centroid
    botmid=c(0,.5,.5)
    rightmid=c(.5,0,.5)
    leftmid=c(.5,.5,0)
    simplexSegment(mid, botmid, lwd = lineWidth)
    simplexSegment(left, rightmid, lwd = lineWidth)
    simplexSegment(right, leftmid, lwd = lineWidth)
  }
}

###############################################

#' Plot point in 2-d probability simplex
#'
#' Normalizes a point given in 3-d non-normalized coordinates, then plots it in 
#' the 2-d probability simplex.
#' 
#' @param v a 3-d point in non-negative orthant, coordinates not summing to 0
#' 
#' @param ... other options to pass to graphics::points function
#' 
#' @return NULL
#'
#' @examples
#'    simplexPrepare("T3","Example Plot")
#'    simplexPoint(c(15,65,20),pch=3,col="blue")
#'
#' @seealso \code{\link{simplexLabels}},
#'          \code{\link{simplexPrepare}},
#'          \code{\link{simplexSegment}},
#'          \code{\link{simplexText}},
#'          \code{\link{simplexCoords}}
#' @export
simplexPoint <- function(v, 
                         ...) {
  coords <- simplexCoords(v)
  points(coords$x, coords$y, ...)
}

##################################################

#' Plot line segment in 2-d probability simplex
#'
#' Normalizes two points in 3-d, and draws line segment between them in 2-d probability simplex.
#' 
#' @param v,w  3-d endpoints of line segment in non-negative orthant, coords not summing to 0
#' 
#' @param ... other options to pass to graphics::segments function
#' 
#' @return NULL
#'
#' @examples
#'    simplexPrepare("T3","Example Plot")
#'    simplexSegment(c(15,65,20),c(15,70, 15),col="green")
#'
#' @seealso \code{\link{simplexLabels}},
#'          \code{\link{simplexPoint}},
#'          \code{\link{simplexPrepare}},
#'          \code{\link{simplexText}},
#'          \code{\link{simplexCoords}}
#' @export
simplexSegment <- function(v, 
                           w, 
                           ...) {
  # Draw line segment in planar simplex
  # Args:
  #      v,w = vectors giving endpoints in R^3
  #
  coords0 <- simplexCoords(v)
  coords1 <- simplexCoords(w)
  segments(coords0$x, coords0$y, coords1$x, coords1$y, ...)
}

########################################################

#' Add text at a point in 2-d probability simplex
#' 
#' Add text to a 2-d probability simplex plot, at specified location.
#'
#' @param v a 3-d point in non-negative orthant, coordinates not summing to 0
#' 
#' @param  label text to add to plot
#' 
#' @param ... other options to pass to graphics::text function
#' 
#' @return NULL
#'
#' @examples
#'    simplexPrepare("T3","Example Plot")
#'    simplexText(c(15,65,20),"tree ac|bd")
#'
#' @seealso \code{\link{simplexLabels}},
#'          \code{\link{simplexPoint}},
#'          \code{\link{simplexPrepare}},
#'          \code{\link{simplexSegment}},
#'          \code{\link{simplexCoords}}
#' @export
simplexText <- function(v, 
                        label = '', 
                        ...) {
  coords <- simplexCoords(v)
  text(coords$x, coords$y, label = label, ...)
}

#######################################################

#' Label vertices of 2-d probability simplex  
#' 
#' Add labels to vertices of the probability simplex.
#'
#' @param top label for top
#' @param left label for left bottom
#' @param right label for right bottom
#' @return NULL
#'
#' @examples
#'    simplexPrepare("T3","Example Plot")
#'    simplexLabels("ab|cd","ac|bd","ad|bc")
#'
#' @seealso
#'          \code{\link{simplexPoint}},
#'          \code{\link{simplexPrepare}},
#'          \code{\link{simplexSegment}},
#'          \code{\link{simplexText}},
#'          \code{\link{simplexCoords}}
#' @export
simplexLabels <- function(top = '',
                    left = '',
                    right = '') {
  topCoord   <- simplexCoords(c(1, 0, 0))
  leftCoord  <- simplexCoords(c(0, 1, 0))
  rightCoord <- simplexCoords(c(0, 0, 1))
  text(topCoord$x, topCoord$y, top, pos = 3)
  text(leftCoord$x, leftCoord$y, left, pos = 1)
  text(rightCoord$x, rightCoord$y, right, pos = 1)
}

###############################################

#' Convert 3-d coordinates to 2-d probability simplex coordinates
#' 
#' Convert from 3-d Cartesian coordinates to 2-d coordinates suitable for plotting in the probability simplex.
#'
#' @details Applies an affine coordinate trandformation that maps the centroid (1/3,1/3,1/3) to the origin (0,0), and 
#' rescales so that the line segments between (1,0,0), (0,1,0), and (0,0,1) are mapped to segments of length 1.
#' 
#' An input vector \code{v} is first normalized so its component sum to 1 before the map is applied.
#' 
#' @param v vector of 3 non-negative numbers, not summing to 0
#' 
#' @return 2-d coordinates to plot normalized point in simplex
#'
#' @examples
#'      simplexCoords(c(15,65,20))
#'
#' @seealso \code{\link{simplexLabels}},
#'          \code{\link{simplexPoint}},
#'          \code{\link{simplexPrepare}},
#'          \code{\link{simplexSegment}},
#'          \code{\link{simplexText}}
#' @export
simplexCoords <- function(v) {
  sum <- sum(v)
  
  if ((length(v) != 3) ||
      (sum(v < 0) > 0) || (sum == 0))
    stop("Invalid coordinates")
  
  v <- v / sum  #normalize so ternary coord sum to 1
  x <-
    (.5 * v[1] + v[3] - .5)  #map to planar triangle, with sides 1, barycenter at origin
  y <- (.8660254 * v[1] - 0.2886751)
  xyCoords <- data.frame(cbind(x = x, y = y))
  return(xyCoords)
}
