# choose variables to use from colnames of namefile (see below)
# "Bodenart__" "Humus____"  "pH_Wert__i" "Karbonate_" "Kalkbedarf" "Phosphat__" "Kali__K_O_" "Magnesium_" "Bor__B__im" "Mangan__Mn" "Kupfer__Cu" "Zink__Zn__" "NUM"       
# "ID_string"  "ID_suolo"   "Num_Soil"   "Soil_newCl" "Schluff"    "Tonig"      "Sand"

# library(gstat)
# library(raster)
# library(caret)
# library(SpatialPosition)

OrdKrig <- function ( wpath = "/home/jbre/R/OrdKrig", 
                      datafolder = "raw", rastermask = "mask/Mask_master.tif",
                      inverseDistWeigths = FALSE, local=TRUE,
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
                      model="Sph",
                      validation = FALSE, kfold=5,
                      coordsys = "+proj=utm +zone=32 ellps=WGS84"
                    )
{
  
  filesIN <- dir(file.path(wpath, datafolder))
  
  if (is.null(names(cutoff)))
  {
    # get namezones according to input files 
    namezones <- c()
    for (namefile in filesIN)
    {
      nr_char <- nchar(filesIN)
      names(nr_char) <- filesIN
      namezone <- strsplit(substr(namefile,1,nr_char[namefile]-4), "_")[[1]][2]
      namezones <- c(namezones, namezone)
    }
    
    # list of arguments
    args <- list("cutoff"=cutoff, "anis_deg"=anis_deg,"anis_ax"=anis_ax,"psill"=psill,
                 "nugget"=nugget,"nmax"=nmax,"nmin"=nmin,"idp"=idp)
    
    # assing argument names to parameter value (named vector)
    for (arg in 1:length(args))
    {
      x <- args[[arg]]
      if (length(x)==length(namezones)) {
        names(x) <- namezones
        assign(x = names(args)[arg], value = x) 
      } else {
        print(paste("argument",  names(args)[arg], "is not given with name of investigated zone and differs in length from input .csv files in the folder", datafolder, ". Please check!", sep=" " ))
      }
    }
  }
    
  val_list <- list()
  
  for (namefile in filesIN)
  {
    
    nr_char <- nchar(filesIN)
    names(nr_char) <- filesIN
    namezone <- strsplit(substr(namefile,1,nr_char[namefile]-4), "_")[[1]][2]
    
    print(paste("processing zone", namezone, sep=" "))
    
    # read table 
    worktab <- read.table(file = file.path(wpath, datafolder, namefile), header = TRUE, sep = ",",dec = ".")
    worktab <- cbind(worktab$x_Coord,worktab$y_Coord,worktab[,variable])
    
    # matrix 2 data.frame
    worktab <- as.data.frame(worktab)
    # rename cols
    names(worktab) <- c("X","Y","VARIABLE")
    
    # zeros
    worktab[worktab$VARIABLE <= 0,"VARIABLE"] <- 0.001
    
    if (validation)
    {
      
      flds <- createFolds(y = worktab$VARIABLE, k = kfold, list = TRUE, returnTrain = FALSE)
      
      val_out_var <- list()
      val_out_df  <- data.frame()
      
      for (i in 1:length(flds))
      {
        # training
        train_set <- worktab[-flds[[i]],]
        
        # train
        
        if (!inverseDistWeigths)
        {
          # gstatVariogram - Calculate Sample variogram 
          my_var <- variogram(log(VARIABLE)~1, data=train_set, locations = ~X+Y, cutoff = cutoff[namezone])
          # Fit a Variogram Model to a Sample Variogram  
          m <- vgm(psill = psill[namezone], model = model, range = cutoff[namezone], nugget = nugget[namezone], 
                   anis = c(anis_deg[namezone], anis_ax[namezone]))
          my_var_fit <- fit.variogram(my_var, m)
        }

        # validation set
        valid_set <- worktab[flds[[i]],]
        
        Xnew <- valid_set[,c("X", "Y")]
        Xnew <- SpatialPoints(Xnew)
        
        myloc <- data.frame("X" = train_set$X,"Y" = train_set$Y)
        myloc <- SpatialPoints(myloc)
        
        if (inverseDistWeigths) {
          
          # Inverse Distance Weighting (local only)
          ord_krig <- gstat::idw(formula = train_set$VARIABLE~1, locations = myloc, newdata = Xnew, idp = idp[namezone],
                                   nmax = nmax[namezone], nmin = nmin[namezone], omax = omax[namezone], maxdist = cutoff[namezone])
          
          names(ord_krig) <- c("predict")
        
          val_out_df <- rbind(val_out_df, data.frame(fold=i ,valid_set, ord_krig$predict))
          
        } else {
          
          # Ordinary Kriging (local only)
          ord_krig <- gstat::krige(formula = train_set$VARIABLE~1, locations = myloc, newdata = Xnew, model = m,
                                   nmax = nmax[namezone], nmin = nmin[namezone], omax = omax[namezone], maxdist = m$range[2])
          
          names(ord_krig) <- c("predict", "variance")
          
          val_out_var[[paste("fold", i, sep="")]] <- list(vario = my_var, vario_fit = my_var_fit)
          val_out_df <- rbind(val_out_df, data.frame(fold=i ,valid_set, ord_krig$predict, ord_krig$variance))
        }
                                                    
      }
      
      val_list[[namezone]] <- list(variogram = val_out_var, df = val_out_df)
     
    } else {
      
    if (!inverseDistWeigths)
    {
      # Choose parameters as function arguments
      
      # gstatVariogram - Calculate Sample variogram 
      my_var <- variogram(log(VARIABLE)~1, data=worktab, locations = ~X+Y, cutoff = cutoff[namezone])
      # Fit a Variogram Model to a Sample Variogram  
      m <- vgm(psill = psill[namezone], model = model, range = cutoff[namezone], nugget = nugget[namezone], 
               anis = c(anis_deg[namezone], anis_ax[namezone]))
      my_var_fit <- fit.variogram(my_var, m)
    }
 
      coordinates(worktab) <- ~X+Y
      crs(worktab) <- coordsys
      
        if (!is.na(rastermask)) {
          
          # get raster mask
          # 1 read in raster
          mask <- raster(file.path(wpath, rastermask))
          # crop mask to data extent
          mask_A <- crop(mask,extent(worktab))
          
          fac <- round(res(mask_A)[2]/npix,0)
          
          # resample to npix
          if (res(mask_A)[2] != npix)
          {
            mask_A <- disaggregate(mask_A, fact=fac, method='')
          }
          
        } else {
          
          # 2 calculate convex hull for worktab data
          ch <- chull(cbind(worktab@coords[,1],worktab@coords[,2]))
          coords <- worktab@coords[c(ch, ch[1]), ] 
          # convert to Polygon
          sp_poly <- SpatialPolygons(list(Polygons(list(Polygon(coords)), ID=1)))
          # convert to raster
          rast <- raster()
          extent(rast) <- c(min(coords[,1]), max(coords[,1]), min(coords[,2]), max(coords[,2]))
          res(rast) <- npix
          mask_A <- rasterize(sp_poly, rast)
          
        }

        # raster 2 SpatialPixelsDataFrame
        mask_sppxdf <- as(mask_A, "SpatialPixelsDataFrame")
        crs(mask_sppxdf) <- coordsys
        
        # get locations to estimate VARIABLE in
        #Xnew <- SpatialPoints(coordinates(mask_A)[!is.na(values(mask_A)),])
        #crs(Xnew) <- coordsys
        
        # Ordinary krigging (gstat)
        if (inverseDistWeigths) {
          
          if (local) {
            ord_krig <- gstat::idw(formula = worktab$VARIABLE~1, worktab, mask_sppxdf, idp = idp[namezone],
                                   nmax = nmax[namezone], nmin = nmin[namezone], omax = omax[namezone], maxdist = cutoff[namezone])
            locglob <- "local"
          } else {
            ord_krig <- gstat::idw(formula = worktab$VARIABLE~1, worktab, mask_sppxdf, idp = idp[namezone])
            locglob <- "global"
          }
          
          names(ord_krig) <- c("predict", "variance")
          
          r_pred <- raster(ord_krig["predict"])
          
          # export raster as GeoTiff
          # different possible formats see ?writeRaster
          dir.create(file.path(wpath, variable, "maps"), recursive = T)
          print("write .tif map files")
          writeRaster(x = r_pred, filename = file.path(wpath, variable, "maps", paste(namezone, "_", variable, "_", npix, "_", locglob, "_predict_sp_idw.tif",sep="")),
                      overwrite=TRUE, format="GTiff")
          
          val_list[[namezone]] <- list(krig = ord_krig, map_pred = r_pred)
          
        } else {
          
          if (local) {
            ord_krig <- gstat::krige(formula = worktab$VARIABLE~1, locations = worktab, newdata = mask_sppxdf, model = m,
                                     nmax = nmax[namezone], nmin = nmin[namezone], omax = omax[namezone], maxdist = m$range[2])
            locglob <- "local"
          } else {
            ord_krig <- gstat::krige(formula = worktab$VARIABLE~1, locations = worktab, newdata = mask_sppxdf, model = m)
            locglob <- "global"
          }
          
          arg_spec <- paste(model, npix ,psill[namezone], cutoff[namezone], nugget[namezone], anis_deg[namezone], anis_ax[namezone],
                            nmax[namezone], nmin[namezone], omax[namezone], sep="_")
          
          names(ord_krig) <- c("predict", "variance")
          
          r_pred <- raster(ord_krig["predict"])
          r_vari <- raster(ord_krig["variance"])
          
          # export raster as GeoTiff
          # different possible formats see ?writeRaster
          dir.create(file.path(wpath, variable, "maps"), recursive = T)
          print("write .tif map files")
          writeRaster(x = r_pred, filename = file.path(wpath, variable, "maps", paste(namezone, "_", variable, "_", locglob, "_predict_sp_krige_",arg_spec,".tif",sep="")),
                      overwrite=TRUE, format="GTiff")
          writeRaster(x = r_vari, filename = file.path(wpath, variable, "maps", paste(namezone, "_", variable, "_", locglob, "_variance_sp_krige_",arg_spec,".tif",sep="")),
                      overwrite=TRUE, format="GTiff")
          
          val_list[[namezone]] <- list(vario = my_var, vario_fit = my_var_fit, krig = ord_krig, map_pred = r_pred, map_var = r_vari)
        }
      
    }
      
    }
    
  if (validation) {
    if (inverseDistWeigths) save(list = "val_list", file = file.path(wpath, variable, "validation_idw.RData"))
    if (!inverseDistWeigths) save(list = "val_list", file = file.path(wpath, variable, "validation_krige.RData"))
    return(val_list)
  } else {
    if (inverseDistWeigths) save(list = "val_list", file = file.path(wpath, variable, "NOvalidation_idw.RData"))
    if (!inverseDistWeigths) save(list = "val_list", file = file.path(wpath, variable, "NOvalidation_krige.RData"))
    return(val_list)
  }
  
}