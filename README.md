# SpatialInterpol

A wrapper of krige::idw (inverse distance weighting) and krige::krige (local ordinary krigging). Function \emph{krige} is preparing simple .txt files for the use with \emph{krige} and \emph{idw}, performing interpolation and writing GEOtif maps. Optimisation functions are available for krige and idw and can be se with \emph{optim} or \emph{hydrPSO}.   

# How to start

First install the package with:

```R
install.packages("devtools")
library(devtools)
install_github("JBrenn/SpatialInterpol")
```

and then import the library with:

```R
library(SpatialInterpol)
```
