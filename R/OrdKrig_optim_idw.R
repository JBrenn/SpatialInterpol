# Optimistation function for ordinary kriging

# choose variables to use from colnames of namefile (see below)
# "Bodenart__" "Humus____"  "pH_Wert__i" "Karbonate_" "Kalkbedarf" "Phosphat__" "Kali__K_O_" "Magnesium_" "Bor__B__im" "Mangan__Mn" "Kupfer__Cu" "Zink__Zn__" "NUM"       
# "ID_string"  "ID_suolo"   "Num_Soil"   "Soil_newCl" "Schluff"    "Tonig"      "Sand"

# library(gstat)
# library(caret)
# library(hydroPSO)
# library(sp)

OrdKrig_optim_idw <- function(par = c(idp = 2.0, maxdist=300, nmax=12, omax=3),
                                wpath = "/home/jbre/R/OrdKrig", 
                                datafile = "master/Masterfile_AdigeVenosta.txt",
                                variable = "Humus____",
                                var_model = "Sph", kfold = 5
                              )
{
  # read table 
  worktab <- read.table(file = file.path(wpath, datafile), header = TRUE, sep = ",",dec = ".")
  worktab <- cbind(worktab$x_Coord, worktab$y_Coord, worktab[,variable])
  
  # matrix 2 data.frame
  worktab <- as.data.frame(worktab)
  # rename cols
  names(worktab) <- c("X","Y","VARIABLE")
  
  # zeros
  worktab[worktab$VARIABLE <= 0,"VARIABLE"] <- 0.001
  
  # get 5 folds
  flds <- createFolds(y = worktab$VARIABLE, k = kfold, list = TRUE, returnTrain = FALSE)
  
  val_out_df  <- data.frame()
  
  for (i in 1:length(flds))
    {
      # training
      train_set <- worktab[-flds[[i]],]
      
      # validation set
      valid_set <- worktab[flds[[i]],]
      
      # IDW
      
      Xnew <- valid_set[,c("X", "Y")]
      Xnew <- SpatialPoints(Xnew)
      
      myloc <- data.frame("X" = train_set$X,"Y" = train_set$Y)
      myloc <- SpatialPoints(myloc)
      
      ord_krig <- idw(formula = train_set$VARIABLE~1, locations = myloc, newdata = Xnew, 
                      idp = par[1], nmax = par[3], nmin = 1, omax = par[4], maxdist = par[2])
      
      names(ord_krig) <- c("predict")
      
      val_out_df <- rbind(val_out_df, data.frame(fold=i ,valid_set, ord_krig$predict))
      
    }
  
  return(RMSE(pred = val_out_df$ord_krig.predict, obs = val_out_df$VARIABLE, na.rm = T))
  
}

# # keep care: trade of between search distance and number of NA estimations
# # the smaller the search radius, the better the estimation - but lot of NAs
# # How to solve?

# hydroPSO::hydroPSO(fn = OrdKrig_optim_idw, method="spso2011",
#                    lower = c(1,100,8,1), upper = c(16,1000,100,25),
#                    control=list(drty.out = "/home/jbre/R/OrdKrig/PSO_idw", npart=40, 
#                                 parallel="none", par.pkgs = c("gstat","caret","hydroGOF","sp")))

