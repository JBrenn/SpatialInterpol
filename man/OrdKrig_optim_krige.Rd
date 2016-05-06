\name{OrdKrig_optim_krige}
\alias{OrdKrig_optim_krige}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
OrdKrig_optim_krige(par = c(cutoff=300, anis_deg=0, anis_ax=.5, nmax=12, omax=3, nugget=1),
                            wpath = "/home/jbre/R/OrdKrig", 
                            datafile = "raw/Masterfile_Adige.txt",
                            variable = "Humus____",
                            model="Sph", kfold=5)
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
  \item{model}{
see \code{\link[gstat]{fit.variogram}}
}
  \item{kfold}{
see \code{\link[caret]{createFolds}}
}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
RMSE of local ordinary krigging estimations vs. observed data
}
\references{
Johannes Brenner, \email{johannesbrenner@gmail.com}
}
\author{
%%  ~~who you are~~
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

# hydroPSO::hydroPSO(fn = OrdKrig_optim_krige, method="spso2011",
#                   lower = c(0,0,0.01,8,1,0), upper = c(1000,359,1,100,25,10),
#                   control=list(drty.out = "/home/jbre/R/OrdKrig/PSO_krige", npart=40, 
#                                parallel="none", par.pkgs = c("gstat","caret","hydroGOF","sp")))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
