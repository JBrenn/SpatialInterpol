\name{OrdKrig_optim_idw}
\alias{OrdKrig_optim_idw}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Optimization function inverse distance weighting
}
\description{
Optimization function performing 5-fold validation and returning RMSE of estimates for inverse distance weighting.
}
\usage{
OrdKrig_optim_idw(par = c(idp = 2.0, maxdist=300, nmax=12, omax=3),
                  wpath = "/home/jbre/R/OrdKrig", 
                  datafile = "master/Masterfile_AdigeVenosta.txt",
                  variable = "Humus____")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{par}{
see \code{\link[stats]{optim}} or \code{\link[hydroPSO]{hydroPSO}}. For names of vector \emph{par} see  \code{\link[gstat]{idw}}
}
  \item{wpath}{
see \code{\link{OrdKrig}}
}
  \item{datafile}{
path to and name of the data to optimize on
}
  \item{variable}{
see \code{\link{OrdKrig}}
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
RMSE of local inverse distance weights estimations vs. observed data
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
# do not run!

# hydroPSO::hydroPSO(fn = OrdKrig_optim_idw, method="spso2011",
#                    lower = c(1,100,8,1), upper = c(16,1000,100,25),
#                    control=list(drty.out = "/home/jbre/R/OrdKrig/PSO_idw", npart=40, 
#                                 parallel="none", par.pkgs = c("gstat","caret","hydroGOF","sp")))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
