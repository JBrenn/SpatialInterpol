\name{OrdKrig}
\alias{OrdKrig}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Ordinary local kriging and inverse distance weighting
}
\description{
Ordinary local kriging and inverse distance weighting. Using global variogram estimation (see \code{\link[gstat]{variogram}}) and krige::krige.locations, krige::idw.locations (see \code{\link[gstat]{krige}}), providing GEOtif raster output according to a given raster mask.  
}
\usage{
OrdKrig(wpath = "/home/jbre/R/OrdKrig", 
        datafolder = "raw", rastermask = "mask/Mask_master.tif",
        inverseDistWeigths = FALSE, local = TRUE,
        variable = "Humus____",  npix = 100,
        cutoff = c("AdigeVenosta"=400, "Adige"=400, "Venosta"=450), 
        anis_deg = c("AdigeVenosta"=0, "Adige"=0, "Venosta"=90), 
        anis_ax = c("AdigeVenosta"=.5, "Adige"=.5, "Venosta"=.5),
        psill = c("AdigeVenosta"=1, "Adige"=1, "Venosta"=1), 
        nugget = c("AdigeVenosta"=1, "Adige"=1, "Venosta"=1),
        nmax = c("AdigeVenosta"=12, "Adige"=12, "Venosta"=12), 
        nmin = c("AdigeVenosta"=1, "Adige"=1, "Venosta"=1),
        omax = c("AdigeVenosta"=3, "Adige"=3, "Venosta"=3),
        idp = c("AdigeVenosta"=2.0, "Adige"=2.0, "Venosta"=2.0),
        var_model="Sph",
        validation = FALSE, kfold=5,
        coordsys = "+proj=utm +zone=32 ellps=WGS84")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{wpath}{
working path containing data files in \emph{datafolder}, as well as a raster mask in \emph{rastermask}.  
}
  \item{datafolder}{
folder within \emph{wpath} containing raw data files
}
  \item{rastermask}{
folder within \emph{wpath} containing raster mask
}
  \item{inverseDistWeigths}{
boolean, if TRUE method is invers distance weighting (IDW), else method is local kriging
}
  \item{local}{
boolean, if TRUE local kriging, else global kriging, see \code{\link[gstat]{krige}}
}
  \item{variable}{
character, variable name to interpolate, choose according to the column name of the input data files
}
  \item{npix}{
grid resolution of the input map. If the resolution of the raster mask coarser then \emph{npix} the raster mask is resampled to \emph{npix} in order to define estimation locations 
}
  \item{cutoff}{
see \code{\link[gstat]{variogram}}, also used for the \emph{range} argument in \code{\link[gstat]{vgm}}
}
  \item{anis_deg}{
see \emph{anis} argument in \code{\link[gstat]{vgm}}, main axis direction
}
  \item{psill}{
see \code{\link[gstat]{vgm}}
}
  \item{nugget}{
see \code{\link[gstat]{vgm}}
}
  \item{nmax}{
see \code{\link[gstat]{krige}}
}
  \item{nmin}{
see \code{\link[gstat]{krige}}
}
  \item{omax}{
see \code{\link[gstat]{krige}}
}
  \item{idp}{
see \code{\link[gstat]{idw}}
}
  \item{model}{
see \code{\link[gstat]{fit.variogram}}
}
  \item{validation}{
boolean, if TRUE fold validation is run, no map output
}
  \item{kfold}{
see \code{\link[caret]{createFolds}}
}
  \item{coordsys}{
characeter, definition of coordinate system, see \code{\link[raster]{crs}}
}

}
\details{
Each input file in datafolder has to contain geographic coordinates with the column names \emph{x_Coord} and \emph{y_Coord} and variable values to interpolate.
Parameters for the functions \link[gstat]{vgm}, \link[gstat]{fit.variagram}, \link[gstat]{variogam}, \link[gstat]{krige} and \link[gstat]{idw} are gven by named vectors. Names refer to the zone names (\emph{namezone}) interpolationis done for, data input files have to named accordingly: \emph{Masterfile_namezone.txt}.
An example data folder can be found here \href{https://cloud.scientificnet.org/index.php/s/950PW7aOFd49ScD}{here}.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Johannes Brenner, \email{johannesbrenner@gmail.com}
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
