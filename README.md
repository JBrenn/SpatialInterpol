# SpatialInterpol

A wrapper of gstat::idw (inverse distance weighting) and gstat::krige (local ordinary krigging). Function **OrdKrig** is preparing simple .txt files for the use with _krige_ and _idw_, performing interpolation and writing GEOtif maps. Optimisation functions are available for krige and idw and can be used with _optim_ or _hydroPSO_.   

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
