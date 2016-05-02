# Optimistation function for ordinary kriging

# choose variables to use from colnames of namefile (see below)
# "Bodenart__" "Humus____"  "pH_Wert__i" "Karbonate_" "Kalkbedarf" "Phosphat__" "Kali__K_O_" "Magnesium_" "Bor__B__im" "Mangan__Mn" "Kupfer__Cu" "Zink__Zn__" "NUM"       
# "ID_string"  "ID_suolo"   "Num_Soil"   "Soil_newCl" "Schluff"    "Tonig"      "Sand"

# library(gstat)
# library(caret)
# library(hydroGOF)
# library(sp)

OrdKrig_optim_idw <- function(par = c(idp = 2.0),
                                wpath = "/home/lv70864/jbrenner/R/OrdKrig", 
                                datafile = "master/Masterfile_AdigeVenosta.txt",
                                rastermask = "mask/Mask_master.tif",
                                variable = "Humus____",
                                var_model = "Sph", kfold = 5,
                                coordsys = "+proj=utm +zone=32 ellps=WGS84"
                              )
{
  # read table 
  worktab <- read.table(file = file.path(wpath, datafile), header = TRUE, sep = ",",dec = ".")
  worktab <- cbind(worktab$x_Coord, worktab$y_Coord, worktab[,variable])
  
  # matrix 2 data.frame
  worktab <- as.data.frame(worktab)
  # rename cols
  names(worktab) <- c("X","Y","VARIABLE")
  worktab_save <- worktab
  
  coordinates(worktab) <- ~X+Y
  crs(worktab) <- coordsys
  
  # zeros
  if (variable == "Humus____") worktab[worktab$VARIABLE <= 0,"VARIABLE"] <- 0.001
  
  # get 5 folds
  flds <- createFolds(y = worktab$VARIABLE, k = kfold, list = TRUE, returnTrain = FALSE)
  
  val_out_df  <- data.frame()
  
  for (i in 1:length(flds))
  {
    # training
    train_set <- worktab_save[-flds[[i]],]
    
    coordinates(train_set) <- ~X+Y
    crs(train_set) <- coordsys
    
    # validation set
    valid_set <- worktab_save[flds[[i]],]
    
    coordinates(valid_set) <- ~X+Y
    crs(valid_set) <- coordsys
    
    if (!is.na(rastermask)) {
      
      # get raster mask
      # 1 read in raster
      mask <- raster(file.path(wpath, rastermask))
      # crop mask to data extent
      mask_A <- crop(mask,extent(worktab))
      
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
    
    # Inverse Distance Weights
    ord_krig <- gstat::idw(formula = train_set$VARIABLE~1, train_set, mask_sppxdf, idp = par[1])
    
    names(ord_krig) <- c("predict", "variance")
    
    # create raster 
    r_pred <- raster(ord_krig["predict"])
    
    # extract estimations for validation points
    estimates <- extract(r_pred, valid_set)
    
    val_out_df <- rbind(val_out_df, data.frame(coordinates(valid_set) ,valid_set@data, ESTIM = estimates))
    
  }
  
  return(RMSE(pred = val_out_df$ESTIM, obs = val_out_df$VARIABLE, na.rm = T))
  
}

# windows machine
# hydroPSO::hydroPSO(fn = OrdKrig_optim_idw, method="ipso",
#                    lower = c(0.01), upper = c(10),
#                    control=list(npart=40, parallel="parallelWin", par.pkgs = c("gstat","caret","hydroGOF","sp","raster")))

# linux server
