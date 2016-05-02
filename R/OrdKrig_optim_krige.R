# Optimistation function for ordinary kriging

# choose variables to use from colnames of namefile (see below)
# "Bodenart__" "Humus____"  "pH_Wert__i" "Karbonate_" "Kalkbedarf" "Phosphat__" "Kali__K_O_" "Magnesium_" "Bor__B__im" "Mangan__Mn" "Kupfer__Cu" "Zink__Zn__" "NUM"       
# "ID_string"  "ID_suolo"   "Num_Soil"   "Soil_newCl" "Schluff"    "Tonig"      "Sand"

# library(gstat)
# library(caret)
# library(hydroGOF)
# library(sp)

OrdKrig_optim_krige <- function(par = c(c_off=300, anis_deg=0, anis_ax=.5),
                                wpath = "C:/Users/JBrenner/Desktop/OrdKrig", 
                                datafile = "raw/Masterfile_Adige.txt",
                                variable = "Humus____",
                                var_model="Sph", kfold=5
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
  if (variable == "Humus____") worktab[worktab$VARIABLE <= 0,"VARIABLE"] <- 0.001
  
  # get 5 folds
  flds <- createFolds(y = worktab$VARIABLE, k = kfold, list = TRUE, returnTrain = FALSE)
  
  val_out_df  <- data.frame()
  
  for (i in 1:length(flds))
  {
    # training
    train_set <- worktab[-flds[[i]],]
    
    # train
    # gstatVariogram - Calculate Sample variogram 
    my_var <- variogram(log(VARIABLE)~1, data=train_set, locations = ~X+Y, cutoff = par[1])
    # Fit a Variogram Model to a Sample Variogram  
    my_var_fit <- fit.variogram(my_var, vgm(1, var_model, par[1], 1, 
                                            anis = c(par[2], par[3])))
    
    # validation set
    valid_set <- worktab[flds[[i]],]
    
    # Ordinary Kriging
    
    Xnew <- valid_set[,c("X", "Y")]
    Xnew <- SpatialPoints(Xnew)
    
    myloc <- data.frame("X" = train_set$X,"Y" = train_set$Y)
    myloc <- SpatialPoints(myloc)
    
    ord_krig <- krige(formula = train_set$VARIABLE~1, locations = myloc,
                      newdata = Xnew,
                      model = my_var_fit, nmax = 10, nmin = 1, maxdist = par[1])
    
    names(ord_krig) <- c("predict", "variance")
    
    val_out_df <- rbind(val_out_df, data.frame(fold=i ,valid_set, ord_krig$predict, ord_krig$variance))
    
  }
  
  return(RMSE(pred = val_out_df$ord_krig.predict, obs = val_out_df$VARIABLE, na.rm = T))
  
}

# # 
hydroPSO::hydroPSO(fn = OrdKrig_optim, method="ipso",
                  lower = c(0,0,0), upper = c(1000,360,1),  
                  control=list(npart=40, parallel="none", par.pkgs = c("gstat","caret","hydroGOF","sp")))
