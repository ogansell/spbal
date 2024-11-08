# findFirstStudyRegionPoint.R

#' @name findFirstStudyRegionPoint
#'
#' @title Get a randomly chosen Halton point from within the study area and the associated seeds.
#'
#' @description This function repeatedly calls function spbal::getBASSample
#' to generate the Halton frame sample. This function selects the first point at random from those
#' points in the study area. This point and the seeds used to generate the sample are returned to
#' the caller.
#'
#' @author This function was written by Phil Davies.
#'
#' @param shapefile Shape file as a polygon (sp or sf) of the study area(s).
#' @param seeds A vector of 2 seeds, u1 and u2. If not specified, the default is NULL and will
#' be defined randomly using function \code{generateUVector}.
#' @param bb Bounding box which defines the Master Sample. A bounding box must be
#' supplied.
#' @param verbose Boolean if you want to see any output printed to screen. Helpful if taking a
#' long time. Default is FALSE i.e. no informational messages are displayed.
#'
#' @return A list containing three variables:
#'
#' \itemize{
#' \item \code{seeds} The u1 and u2 seeds used to generate the first point.
#' \item \code{k} The index of the first point in the initial sample.
#' }
#'
# 1. Set J1 = 4 and J2 = 3.
# 2. Generate B = 2^J1 x 3^J2 points from a random-start Halton sequence H
#    with a random seed (u1, u2).
# 3. Find points from H in the study area. Call this set S. If S is empty,
#    increment J1 and J2 and go to step 2.
# 4. Randomly choose a point from S. Let xk be this point where k is the 'site index'
#    (I think that's what we call it).
# 5. Set the seed to (u1 + k - 1, u2 + k - 1).
# 6. Re-number ID by subtracting k (re-generate sample using seeds for 5 - first point must also be ID=1)

# For example, let (u1, u2) = (1, 5) and S = {x2, x6, x7}.
# If x6 is randomly chosen, then the new seed is (1 + 6 - 1, 5 + 6 - 1) = (6, 10)
# (the sixth point in H).

# The only difference is that the random-start Halton sequence must be length B.
#' @keywords internal
findFirstStudyRegionPoint <- function(shapefile, bb, seeds, verbose = FALSE){
  # must not be called without seeds! (also checked in getBASSample).
  if(base::is.null(seeds)){
    msg <- "spbal(findFirstStudyRegionPoint) The seeds parameter must not be NULL."
    msgs <- base::sprintf(msg)
    base::stop(msgs)
  }

  # Initialise variables.
  J <- base::c(4, 3)
  bases <- base::c(2, 3)
  crs <- sf::st_crs(shapefile)

  # default number of sample points to find.
  n <- (bases[1]^J[1]) * (bases[2]^J[2])

  pts_in_intersection <- 0
  call.getBASSample.cnt <- 0

  while(pts_in_intersection < 1){
    # shapefile, bb, n, seeds
    call.getBASSample.cnt <- call.getBASSample.cnt + 1
    result <- getBASSample(shapefile = shapefile, bb = bb, n = n, seeds = seeds)
    diff_ <- result$sample
    seeds <- result$seed

    # find number of points within our study area/bb intersection.
    pts_in_intersection <- base::length(diff_$SiteID)
    n <- n * 2
  }

  if(verbose){
    msg <- "spbal(findFirstStudyRegionPoint) Needed %s call(s) to locate first study area point."
    msgs <- base::sprintf(msg, call.getBASSample.cnt)
    base::message(msgs)
  }

  # select a point in the study area at random.
  base::set.seed(seeds[1] + seeds[2])
  k <- base::sample(pts_in_intersection, 1)
  k <- diff_$SiteID[k]

  # select our random first point. # return SiteID = 1, know what k is.
  #first.pt <- diff_ #[1,]

  if(verbose){
    msg <- "spbal(findFirstStudyRegionPoint) Random point selected: %s."
    msgs <- base::sprintf(msg, k)
    base::message(msgs)
  }

  result <- base::list(seeds = seeds,
                       k     = k)
  return(result)
}

## Potentially faster than sf st_sample simple random sample of polygon.
## Not to be exported for now. Likely slower in some situations.
SRSPoly <- function(n = 1, shapefile, bb, verbose = FALSE){
  bb.bounds <- sf::st_bbox(bb)
  n_found <- 0
  ndraw <- n + 10
  while(n_found < n){
    xy <- cbind(runif(ndraw, bb.bounds["xmin"], bb.bounds["xmax"]), runif(ndraw, bb.bounds["ymin"], bb.bounds["ymax"]))
    pts.coord <- sf::st_as_sf(base::data.frame(SiteID = 1:ndraw, xy), coords = c(2, 3))
    
    sf::st_crs(pts.coord) <- sf::st_crs(bb)
    # find the intersection. Generates the same as sf::st_intersection(pts.coord, shapefile)
    if(n_found == 0) {
      pts.intersect <- pts.coord[shapefile,]
    }else{ 
      pts.intersect <- rbind( pts.intersect, pts.coord[shapefile,] )
    }
    n_found <- nrow(pts.intersect)
    ndraw <- ndraw*2
  }
  return(sf::st_coordinates(pts.intersect[1:n,]))
}

## Function to find a random point in the polygon and to then find where in the Halton Sequence it is.
## Blair to provide better algebra than my current method.
setBASSeed <- function(shapefile, bb, verbose = FALSE){
  bases <- base::c(2, 3)
  crs <- sf::st_crs(shapefile)
  
  bb.bounds <- sf::st_bbox(bb)
  scale.bas <- bb.bounds[3:4] - bb.bounds[1:2]
  shift.bas <- bb.bounds[1:2]
  
  bases <- c(2,3)
  J <- c(7,5)
  Bxy <- bases^J
  seeds <- c(0,0)

  ## Get a single random sample from the polygon.
  pts.unif <- SRSPoly(n = 1, shapefile, bb, verbose)

  upts <- (pts.unif - shift.bas)/scale.bas

  Axy <- cbind(floor((upts[,1] + 2*.Machine$double.eps)*bases[1]^J[1]), 
               floor((upts[,2] + 2*.Machine$double.eps)*bases[2]^J[2]))
  ix <- round(spbal:::cppBASpts(n = Bxy[1], seeds = 0, bases = bases[1], FALSE)$pts*Bxy[1],0)  ## Make sure it's an integer
  iy <- round(spbal:::cppBASpts(n = Bxy[2], seeds = 0, bases = bases[2], FALSE)$pts*Bxy[2],0)
  seeds[1] <- which(ix == Axy[1]) - 1  ## Starts at 0.
  seeds[2] <- which(iy == Axy[2]) - 1
  
  U <- sample(10000, 2, replace = TRUE)
  seeds <- seeds + U*Bxy

  result <- base::list(seeds = seeds, k = 1) ## Keep k for backward compatability.
  return(result)
}