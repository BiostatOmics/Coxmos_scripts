#Script for Garnatxa cluster
args <- commandArgs(trailingOnly = TRUE)

#### ### ### ### ### ### #
# UPDATING GLOBALS SIZE #
#### ### ### ### ### ### #
MB = 6000
bytes = MB*1024^2
options(future.globals.maxSize = bytes)

#libraries
library(Coxmos)
library(RColorConesa) #from GitHub #devtools #usethis and get install.packages('hrbrthemes', repos='http://cran.us.r-project.org')
library(openxlsx)

#Load ggplot theme
loadGgplotTheme <- function(path){
  file <- paste0(path,"ggplot_theme.R")
  source(file, echo = F)
}

path = "/home/salguero/Proyectos/" #edit to your current path
loadGgplotTheme(path)

#### ### ## #
# TO MODIFY #
#### ### ## #

load(args[1])
NAME = args[2]
data_type = args[3]
data_norm = args[4]

if(!data_type %in% c("1", "2", "3","CLINICAL", "OMIC", "MO")){
  stop("Third argument must be one of: '1', '2', '3', 'CLINICAL', 'OMIC' or 'MO'")
}

#### ### ### ### #
# CONCAT VERSION #
#### ### ### ### #

# Perform only for MO data

if(data_type %in% c("1","CLINICAL")){
  # Classical
  FLAG_COX = T
  # HD
  FLAG_COXSW = T
  FLAG_COXEN = T
  FLAG_sPLSICOX = T
  FLAG_sPLSDRCOX_PENALTY = T
  FLAG_sPLSDRCOX_DYNAMIC = T
  FLAG_sPLSDACOX_DYNAMIC = T
  # MO
  FLAG_SB.sPLSICOX = F
  FLAG_iSB.sPLSICOX = F

  FLAG_SB.sPLSDRCOX_PENALTY = F
  FLAG_iSB.sPLSDRCOX_PENALTY = F

  FLAG_SB.sPLSDRCOX_DYNAMIC = F
  FLAG_iSB.sPLSDRCOX_DYNAMIC = F

  FLAG_SB.sPLSICOX = F
  FLAG_iSB.sPLSICOX = F

  FLAG_SB.sPLSDACOX = F
  FLAG_iSB.sPLSDACOX = F

  FLAG_MB.sPLSDRCOX = F
  FLAG_MB.sPLSDACOX = F
}else if(data_type %in% c("3","OMIC")){ #### THIS HAS BEEN CHANGED !!!
  # Classical
  FLAG_COX = F
  # HD
  FLAG_COXSW = T
  FLAG_COXEN = T
  FLAG_sPLSICOX = T
  FLAG_sPLSDRCOX_PENALTY = T
  FLAG_sPLSDRCOX_DYNAMIC = T
  FLAG_sPLSDACOX_DYNAMIC = T
  # MO
  FLAG_SB.sPLSICOX = F
  FLAG_iSB.sPLSICOX = F

  FLAG_SB.sPLSDRCOX_PENALTY = F
  FLAG_iSB.sPLSDRCOX_PENALTY = F

  FLAG_SB.sPLSDRCOX_DYNAMIC = F
  FLAG_iSB.sPLSDRCOX_DYNAMIC = F

  FLAG_SB.sPLSICOX = F
  FLAG_iSB.sPLSICOX = F

  FLAG_SB.sPLSDACOX = F
  FLAG_iSB.sPLSDACOX = F

  FLAG_MB.sPLSDRCOX = F
  FLAG_MB.sPLSDACOX = F
}else if(data_type %in% c("2","MO")){ # THIS HAS BEEN CHANGED !!!
  # Classical
  FLAG_COX = F
  # HD
  FLAG_COXSW = F
  FLAG_COXEN = F
  FLAG_sPLSICOX = F
  FLAG_sPLSDRCOX_PENALTY = F
  FLAG_sPLSDRCOX_DYNAMIC = F
  FLAG_sPLSDACOX_DYNAMIC = F
  # MO
  FLAG_SB.sPLSICOX = T
  FLAG_iSB.sPLSICOX = T

  FLAG_SB.sPLSDRCOX_PENALTY = T
  FLAG_iSB.sPLSDRCOX_PENALTY = T

  FLAG_SB.sPLSDRCOX_DYNAMIC = T
  FLAG_iSB.sPLSDRCOX_DYNAMIC = T

  FLAG_SB.sPLSICOX = T
  FLAG_iSB.sPLSICOX = T

  FLAG_SB.sPLSDACOX = T
  FLAG_iSB.sPLSDACOX = T

  FLAG_MB.sPLSDRCOX = T
  FLAG_MB.sPLSDACOX = T
}

# METHODS
lst_evaluations <- c("survivalROC", "cenROC", "nsROC", "smoothROCtime_C", "smoothROCtime_I", "risksetROC")
names(lst_evaluations) <- lst_evaluations

# Scale Parameters
if(data_type %in% c("1","CLINICAL")){
  x.center = T
  x.scale = T
}else if(data_type %in% c("2","OMIC")){
  x.center = T
  x.scale = F
}else if(data_type %in% c("3","MO")){
  x_omics <- ls(pattern = "X\\.")
  y_omics <- ls(pattern = "Y\\.")
  lst_omics_names <- NULL
  x.center = list()
  x.scale = list()
  for(o in x_omics){
    omic_name <- strsplit(o, "\\.")[[1]][2]
    lst_omics_names <- c(lst_omics_names, omic_name)
    if(endsWith(x = o, "clinical")){
      x.scale[[omic_name]] = T
    }else{
      x.scale[[omic_name]] = F
    }
    x.center[[omic_name]] = T
  }
}

#### ### ### ### ### ### ### ### ### ##
# Set BLOCKS - INTERSECT OBSERVATIONS #
#### ### ### ### ### ### ### ### ### ##
if(data_type %in% c("3","MO")){
  x_omics <- ls(pattern = "X\\.")
  y_omics <- ls(pattern = "Y\\.")

  Y.full <- mget(y_omics[[1]])[[1]]
  for(idx in 2:length(y_omics)){
    aux <- mget(y_omics[[idx]])[[1]]
    Y.full <- rbind(Y.full, aux[!rownames(aux) %in% rownames(Y.full),])
  }

  X = list()
  for(o in x_omics){
    omic_name <- strsplit(o, "\\.")[[1]][2]
    X[[omic_name]] = mget(o)[[1]]
  }

  all_X <- purrr::map(X, ~rownames(.))

  #rownames 1st omic and in the others, including Y
  X1_NONAs = X[[1]]
  final_names <- rownames(X1_NONAs)
  for(i in 2:length(X)){
    final_names <- intersect(all_X[[i]], final_names)
  }
  final_names <- intersect(rownames(Y.full), final_names)

  #update X and Y
  for(i in 1:length(X)){
    X[[i]] = X[[i]][final_names,,drop=F]
  }

  Y <- Y.full[final_names,,drop=F]
}

# Survival Parameters
MIN_EPV = 5
pred.attr = "mean"

# Algorithm Parameters
remove_non_significant = T
remove_non_significant_models = F

remove_near_zero_variance = T
remove_zero_variance = F
toKeep.zv = NULL
remove_variance_at_fold_level = F
alpha = 0.05

returnData = T
verbose = T
PARALLEL = T #!!!

# cox
FORCE = T

# SW
BACKWARDS = T
initialModel = "NULL"
toKeep.sw = NULL
alpha_ENT = 0.1
alpha_OUT = 0.15
alpha_PH  = 0.05
check_PH = F

# coxEN
EN.alpha.list = seq(0.1,0.9,0.2)

# PLS
max.ncomp = 8
max.iter = 500
tol = 1e-15

# Cross Validation Parameters
times = NULL
return_models = F
seed = 123

# Weights Parameters
w_AIC = 0
w_C.Index = 0
w_I.BRIER = 0
w_AUC = 1

# sPLS-ICOX # sPLS-DRCOX
penalty.list = seq(0.1,0.9,0.2)

# mixOmics
design = NULL
vector = NULL
MIN_NVAR = 1 #in order to create a vector and not only get the maximum number of variables if lesser than 10
MAX_NVAR = NULL
n.cut_points = 10

# Eval stop detection
MIN_AUC_INCREASE = 0.01 # 1%
MIN_AUC = 0.85 # 75%
MIN_COMP_TO_CHECK = 3
EVAL_METHOD = "AUC"

# Model Lists
lst_models_full <- NULL
lst_models <- NULL
lst_models_dummy <- NULL
lst_models_pls <- NULL

# Evaluation multiple models
max_time_points = 15

# Kaplan-Meiers
minProp = 0.2

#### ### ### #### ### ### ###
#Cross Validation Parameters #
#### ### ### ### ### #### ###

# Function to optimize % TRAIN and % TEST and # Folds.
# queremos 50 obs por fold y menor que el test
# o bien, que el test final tenga más obs que un fold
PTRAIN = seq(0.7, 0.8, 0.05)
f = c(5:10)

EVENT_RATIO <- sum(Y$event) / nrow(Y)
#MIN 10 evento por fold
MIN_EVENTS_PER_FOLD <- 10
min_n_fold <- floor(MIN_EVENTS_PER_FOLD/EVENT_RATIO)

aux <- expand.grid(PTRAIN, f)
colnames(aux) <- c("PTRAIN", "folds")

n_fold <- NULL
n_test <- NULL
for(i in 1:nrow(aux)){
  n_fold <- c(n_fold, floor(floor(nrow(Y) * aux[i,1]) / aux[i,2]))
  n_test <- c(n_test, floor(nrow(Y) * (1-aux[i,1])))
}

aux <- cbind(aux, n_fold)
aux <- cbind(aux, n_test)

aux$diff <- aux$n_test-aux$n_fold

idx <- which(aux$n_fold > 50 & aux$n_fold < n_test & aux$n_fold >= min_n_fold)
if(length(idx)==0){
  idx <- which(aux$n_test > aux$n_fold & aux$n_fold >= min_n_fold)
  if(length(idx)!=0){
    subaux <- aux[idx,,drop=F]
    idx <- which.min(subaux$diff)
    best_fila <- subaux[idx,,drop=F]
  }else{
    idx <- which(aux$n_test > aux$n_fold)
    idx <- which.max(aux[idx,]$n_fold)
    best_fila <- aux[idx,,drop=F]
  }
}else{
  subaux <- aux[idx,,drop=F]
  idx <- which.min(subaux$diff)
  best_fila <- subaux[idx,,drop=F]
}

n_run = 3
k_folds = best_fila$folds
PTRAIN = best_fila$PTRAIN

####

fast_mode = F
pred.method = "cenROC"

todaydate <- format(Sys.time(), '%Y-%m-%d')
txt_folder <- paste0(NAME, ifelse(fast_mode, "FAST_", "CONCATENATE_"), data_norm, "_", pred.method, "_runs_", n_run, "_folds_", k_folds)
folder <- paste0(txt_folder,"_",todaydate,"/")

dir.create(folder)

#### ### ### ### ### ### ### ##
# IF CLINICAL GET BINARY DATA #
#### ### ### ### ### ### ### ##

FLAG_CLINIC_TRANSFORMATION <- FALSE
if(data_type %in% c("1","CLINICAL")){
  message("CLINICAL does have cox, coxSW. X matrix updated to k-1 variables for factors.")
  if(any(unlist(purrr::map(colnames(X), ~class(X[,.]))) %in% "factor")){
    FLAG_CLINIC_TRANSFORMATION = TRUE
    X.original <- X
    X.dummy <- factorToBinary(X = X, all = F, sep = "_") # create dummy-1 binary variables
    X.pls <- factorToBinary(X = X, all = T, sep = "_") # create dummy binary variables
  }else{
    X.dummy <- X
    X.pls <- X
  }
}else if(data_type %in% c("3","MO")){
  message("MO could have factors in clinical matrix if present.")
  X.dummy <- X
  X.pls <- X
  X.original <- X
  for(block in names(X)){
    if(any(unlist(purrr::map(colnames(X[[block]]), ~class(X[[block]][,.]))) %in% "factor")){
      FLAG_CLINIC_TRANSFORMATION = TRUE
      X.dummy[[block]] <- factorToBinary(X = X[[block]], all = F, sep = "_") # create dummy-1 binary variables
      X.pls[[block]] <- factorToBinary(X = X[[block]], all = T, sep = "_") # create dummy binary variables
    }
  }
}else{
  if(any(unlist(purrr::map(colnames(X), ~class(X[,.]))) %in% "factor")){
    message("OMIC should not have factors in X matrix.")
    FLAG_CLINIC_TRANSFORMATION = TRUE
    X.original <- X
    X.dummy <- factorToBinary(X = X, all = F, sep = "_") # create dummy-1 binary variables
    X.pls <- factorToBinary(X = X, all = T, sep = "_") # create dummy binary variables
  }else{
    X.dummy <- X
    X.pls <- X
  }
}

#### ### ### ### ###
#Set Train and Test#
#### ### ### ### ###
set.seed(123)
index_train <- caret::createDataPartition(Y$event,
                                          p = PTRAIN, #dynamic
                                          list = FALSE,
                                          times = 1)

#same test because complete test have all the variables
if(data_type %in% c("3","MO")){
  X_train <- list()
  X_train.dummy <- list()
  X_test <- list()
  for(omic in names(X.pls)){
    X_train[[omic]] <- X.pls[[omic]][index_train,]
    X_train.dummy[[omic]] <- X.dummy[[omic]][index_train,]
    X_test[[omic]] <- X.pls[[omic]][-index_train,]
  }

  Y_train <- Y[index_train,]
  Y_test <- Y[-index_train,]

}else{
  X_train <- X.pls[index_train,]
  X_train.dummy <- X.dummy[index_train,]
  X_test <- X.pls[-index_train,]

  Y_train <- Y[index_train,]
  Y_test <- Y[-index_train,]
}

#### ### ### ### ### ###
# Function CONCATENATE #
#### ### ### ### ### ###

# HACER LOS CAMBIOS EN TRAIN Y APLICAR ESO A TEST
Scaling.type <- function(reguVal, scaletype){
  res.mat = NULL
  for (ov in names(reguVal)){
    des.mat2 = reguVal[[ov]]
    colnames(des.mat2) <- paste0(colnames(des.mat2), "_", ov)
    if(scaletype=='softBlock'){
      #Escalado tipo pareto
      des.mat2 = des.mat2 / (ncol(des.mat2)^(1/4))
    }
    if(scaletype=='hardBlock'){
      #Cada bloque de variables tiene el mismo peso total en el modelo
      des.mat2 = des.mat2 / (ncol(des.mat2)^(1/2))
    }
    if (is.null(res.mat)) {
      res.mat = des.mat2
    } else {
      res.mat = cbind(res.mat, des.mat2)
    }
  }
  rownames(res.mat) = rownames(reguVal[[ov]])
  return(res.mat)
}

scaling.ISGL <- function(OmValues,des.mat, scaleType = 'auto'){
  ## Specify the scaling type
  scale <- ifelse(scaleType == 'none', FALSE, TRUE)
  center <- ifelse(scaleType == 'none', FALSE, TRUE)
  for (i in 1:length(OmValues)) {
    #Separar data.train y data.validate
    data.train = OmValues[[i]]$data.train
    data.validate = OmValues[[i]]$data.validate
    if (names(OmValues)[i] != 'response'){
      dataT.scaled = scale(data.train, scale = scale, center = center)
      # Centering and scaling factors
      center.at = attr(dataT.scaled, "scaled:center")
      scale.at = attr(dataT.scaled, "scaled:scale")
      # Scale validation data acording to scaling parameters of training data
      centered_data = sweep(data.validate, 2, center.at, FUN = "-")
      dataV.scaled = sweep(centered_data, 2, scale.at, FUN = "/")
      #Remove regulators lost by low variability and save as new values
      lowVar =  names(which(apply(dataT.scaled,2,function(x) all(is.nan(x)))))
      OmValues[[i]]$data.train = dataT.scaled[,!colnames(dataT.scaled) %in% lowVar]
      OmValues[[i]]$data.validate = dataV.scaled[,!colnames(dataV.scaled) %in% lowVar]
    } else{
      if(scaleType != 'none'){
        intercept = mean(data.train)
        OmValues[[i]]$data.train = data.train - intercept
        OmValues[[i]]$data.validate = data.validate - intercept
      } else{ intercept = 0 }
    }
  }
  # Once variables are scaled apply the group scaling if necessary
  data.train = lapply(OmValues, function(x) x$data.train)
  data.validate = lapply(OmValues, function(x) x$data.validate)
  dataT.scaled = Scaling.type(data.train[!grepl('condition|response', names(data.train))],scaletype = scaleType)
  dataV.scaled = Scaling.type(data.validate[!grepl('condition|response', names(data.validate))],scaletype = scaleType)
  data.train = list( x = as.matrix(if (is.null(des.mat)) dataT.scaled else data.frame(data.train$condition, dataT.scaled, check.names = FALSE)), y = data.train$response )
  data.validate = list( x = as.matrix(if (is.null(des.mat)) dataV.scaled else data.frame(data.validate$condition, dataV.scaled, check.names = FALSE)), y = data.validate$response )
  if(!is.null(des.mat)){
    colnames(data.train$x)[1:ncol(des.mat)]=colnames(des.mat)
    colnames(data.validate$x)[1:ncol(des.mat)]=colnames(des.mat)
  }
  return(list('data.train'=data.train, 'data.validate'=data.validate, 'intercept'= intercept))
}

# requiere que los datos estén centrados y escalados antes de la corrección
lst_center = list()
lst_scale = list()
for(b in names(X_train)){
  X_train.dummy[[b]] <- scale(X_train.dummy[[b]], center = TRUE, scale = TRUE)
  X_train[[b]] <- scale(X_train[[b]], center = TRUE, scale = TRUE)
  lst_center[[b]] <- attr(X_train[[b]], "scaled:center")
  lst_scale[[b]] <- attr(X_train[[b]], "scaled:scale")
  X_test[[b]] <- scale(X_test[[b]], center = lst_center[[b]], scale = lst_scale[[b]])
}

#data_norm: hardBlock - softBlock - auto
X_train <- Scaling.type(reguVal = X_train, scaletype = data_norm)
X_train.dummy <- Scaling.type(reguVal = X_train.dummy, scaletype = data_norm)
X_test <- Scaling.type(reguVal = X_test, scaletype = data_norm)

# ALL NORM DONE !!!
if(data_norm %in% c("hardBlock", "softBlock")){
  x.scale <- FALSE
  x.center <- FALSE
}else{
  x.scale <- TRUE
  x.center <- FALSE
}

save(list = c("X_train", "X_train.dummy", "X_test", "Y_train", "Y_test"), file = paste0(folder, "train_test_data.RData")) # eval_results

### ###
# COX #
### ###
if(FLAG_COX){
  aux_folder = paste0(folder, "cox_plot/")
  dir.create(aux_folder)

  best_cox <- cox(X = X_train.dummy, Y = Y_train,
                  x.center = x.center, x.scale = x.scale,
                  remove_near_zero_variance = remove_near_zero_variance,
                  remove_zero_variance = remove_zero_variance,
                  toKeep.zv = toKeep.zv,
                  remove_non_significant = remove_non_significant,
                  alpha = alpha,
                  MIN_EPV = MIN_EPV,
                  FORCE = FORCE,
                  returnData = returnData,
                  verbose = verbose)

  message(paste0("\nTIME FOR COX: ", best_cox$time, "\n"))

  save(list = c("best_cox"), file = paste0(aux_folder, "cox.RData"))
  gc()

  lst_models_full[[best_cox$class]] = best_cox
  lst_models[[best_cox$class]] = best_cox
  lst_models_dummy[[best_cox$class]] = best_cox
}

### ### ###
# COX - SW #
### ### ###

# Add LIMIT of 22k variables or too much time to be computed.
# for HNSC cnv (55k variables)
if(FLAG_COXSW){
  if(ncol(X_train.dummy)<=22000){
    aux_folder = paste0(folder, "coxSW_plot/")
    dir.create(aux_folder)

    best_coxSW <- coxSW(X = X_train.dummy, Y = Y_train,
                        x.center = x.center, x.scale = x.scale,
                        initialModel = initialModel, toKeep.sw = toKeep.sw,
                        remove_near_zero_variance = remove_near_zero_variance,
                        remove_zero_variance = remove_zero_variance,
                        toKeep.zv = toKeep.zv,
                        remove_non_significant = remove_non_significant,
                        max.variables = ncol(X_train.dummy),
                        alpha = alpha, alpha_ENT = alpha_ENT, alpha_OUT = alpha_OUT,
                        MIN_EPV = MIN_EPV,
                        BACKWARDS = BACKWARDS, verbose = verbose)

    message(paste0("\nTIME FOR COXSW: ", best_coxSW$time, "\n"))

    save(list = c("best_coxSW"), file = paste0(aux_folder, "coxSW.RData"))
    gc()

    lst_models_full[[best_coxSW$class]] = best_coxSW
    lst_models[[best_coxSW$class]] = best_coxSW
    lst_models_dummy[[best_coxSW$class]] = best_coxSW
  }
}

### ### ##
# COX EN #
### ### ##
if(FLAG_COXEN){
  aux_folder = paste0(folder, "coxEN_plot/")
  dir.create(aux_folder)

  cv.coxEN_res <- cv.coxEN(X = X_train.dummy, Y = Y_train,
                           EN.alpha.list = EN.alpha.list, #EN penalization
                           max.variables = ncol(X_train.dummy),
                           n_run = n_run, k_folds = k_folds,
                           alpha = alpha, remove_non_significant = remove_non_significant, times = times, max_time_points = max_time_points,
                           remove_near_zero_variance = remove_near_zero_variance,
                           remove_zero_variance = remove_zero_variance,
                           toKeep.zv = toKeep.zv,
                           w_AIC = w_AIC, w_C.Index = w_C.Index, w_I.BRIER = w_I.BRIER, w_AUC = w_AUC,
                           MIN_AUC_INCREASE = MIN_AUC_INCREASE, MIN_AUC = MIN_AUC, MIN_COMP_TO_CHECK = MIN_COMP_TO_CHECK,
                           x.center = x.center, x.scale = x.scale,
                           fast_mode = fast_mode, return_models = return_models, MIN_EPV = MIN_EPV,
                           pred.attr = pred.attr, pred.method = pred.method, seed = seed, PARALLEL = PARALLEL,
                           verbose = verbose)

  message(paste0("\nTIME FOR cv.COXEN: ", cv.coxEN_res$time, "\n"))

  save_ggplot(cv.coxEN_res$plot_AUC, folder = aux_folder, wide = T, name = "cv.coxEN_res_AUC", format = "svg")
  save_ggplot(cv.coxEN_res$plot_BRIER, folder = aux_folder, wide = T, name = "cv.coxEN_res_BRIER", format = "svg")
  save_ggplot(cv.coxEN_res$plot_c_index, folder = aux_folder, wide = T, name = "cv.coxEN_res_c_index", format = "svg")
  save_ggplot(cv.coxEN_res$plot_AIC, folder = aux_folder, wide = T, name = "cv.coxEN_res_AIC", format = "svg")

  best_coxEN <- coxEN(X = X_train.dummy, Y = Y_train,
                      EN.alpha = cv.coxEN_res$opt.EN.alpha,
                      max.variables = cv.coxEN_res$opt.nvar,
                      x.center = x.center, x.scale = x.scale,
                      remove_near_zero_variance = remove_near_zero_variance,
                      remove_zero_variance = remove_zero_variance,
                      toKeep.zv = toKeep.zv,
                      remove_non_significant = remove_non_significant,
                      alpha = alpha,
                      MIN_EPV = MIN_EPV,
                      returnData = returnData,
                      verbose = verbose)

  save(list = c("cv.coxEN_res", "best_coxEN"), file = paste0(aux_folder, "coxEN.RData"))
  gc()

  lst_models_full[[cv.coxEN_res$class]] = cv.coxEN_res
  lst_models_full[[best_coxEN$class]] = best_coxEN
  lst_models[[best_coxEN$class]] = best_coxEN
  lst_models_dummy[[best_coxEN$class]] = best_coxEN

}

### ### ### #
# sPLS-ICOX #
### ### ### #
if(FLAG_sPLSICOX){
  aux_folder = paste0(folder, "splsicox_plot/")
  dir.create(aux_folder)

  cv.splsicox_res <- cv.splsicox(X = X_train, Y = Y_train,
                                 max.ncomp = max.ncomp, penalty.list = penalty.list,
                                 n_run = n_run, k_folds = k_folds,
                                 alpha = alpha, remove_non_significant = remove_non_significant, times = times, max_time_points = max_time_points,
                                 remove_near_zero_variance = remove_near_zero_variance,
                                 remove_zero_variance = remove_zero_variance,
                                 toKeep.zv = toKeep.zv,
                                 w_AIC = w_AIC, w_C.Index = w_C.Index, w_I.BRIER = w_I.BRIER, w_AUC = w_AUC,
                                 MIN_AUC_INCREASE = MIN_AUC_INCREASE, MIN_AUC = MIN_AUC, MIN_COMP_TO_CHECK = MIN_COMP_TO_CHECK,
                                 x.center = x.center, x.scale = x.scale,
                                 fast_mode = fast_mode, return_models = return_models, MIN_EPV = MIN_EPV,
                                 pred.attr = pred.attr, pred.method = pred.method, seed = seed, PARALLEL = PARALLEL,
                                 verbose = verbose)

  message(paste0("\nTIME FOR cv.sPLS-ICOX: ", cv.splsicox_res$time, "\n"))

  if(!is.null(cv.splsicox_res$best_model_info)){
    save_ggplot(cv.splsicox_res$plot_AUC, folder = aux_folder, wide = T, name = "cv.splsicox_res_AUC", format = "svg")
    save_ggplot(cv.splsicox_res$plot_BRIER, folder = aux_folder, wide = T, name = "cv.splsicox_res_BRIER", format = "svg")
    save_ggplot(cv.splsicox_res$plot_c_index, folder = aux_folder, wide = T, name = "cv.splsicox_res_c_index", format = "svg")
    save_ggplot(cv.splsicox_res$plot_AIC, folder = aux_folder, wide = T, name = "cv.splsicox_res_AIC", format = "svg")

    best_splsicox <- splsicox(X = X_train, Y = Y_train,
                              n.comp = cv.splsicox_res$opt.comp,
                              penalty = cv.splsicox_res$opt.penalty,
                              x.center = x.center, x.scale = x.scale,
                              remove_near_zero_variance = remove_near_zero_variance,
                              remove_zero_variance = remove_zero_variance,
                              toKeep.zv = toKeep.zv,
                              remove_non_significant = remove_non_significant,
                              alpha = alpha,
                              MIN_EPV = MIN_EPV,
                              returnData = returnData,
                              verbose = verbose)

    save(list = c("cv.splsicox_res", "best_splsicox"), file = paste0(aux_folder, "splsicox.RData"))
    gc()

    lst_models_full[[cv.splsicox_res$class]] = cv.splsicox_res
    lst_models_full[[best_splsicox$class]] = best_splsicox
    lst_models[[best_splsicox$class]] = best_splsicox
    lst_models_pls[[best_splsicox$class]] = best_splsicox
  }

}

### ### ### ##
# sPLS-DRCOX #
### ### ### ##
if(FLAG_sPLSDRCOX_PENALTY){
  aux_folder = paste0(folder, "splsdrcox_penalty_plot/")
  dir.create(aux_folder)

  cv.splsdrcox_penalty_res <- cv.splsdrcox_penalty(X = X_train, Y = Y_train,
                                   max.ncomp = max.ncomp, penalty.list = penalty.list,
                                   n_run = n_run, k_folds = k_folds,
                                   alpha = alpha, remove_non_significant = remove_non_significant, times = times, max_time_points = max_time_points,
                                   remove_near_zero_variance = remove_near_zero_variance,
                                   remove_zero_variance = remove_zero_variance,
                                   toKeep.zv = toKeep.zv,
                                   w_AIC = w_AIC, w_C.Index = w_C.Index, w_I.BRIER = w_I.BRIER, w_AUC = w_AUC,
                                   MIN_AUC_INCREASE = MIN_AUC_INCREASE, MIN_AUC = MIN_AUC, MIN_COMP_TO_CHECK = MIN_COMP_TO_CHECK,
                                   x.center = x.center, x.scale = x.scale,
                                   fast_mode = fast_mode, return_models = return_models, MIN_EPV = MIN_EPV,
                                   pred.attr = pred.attr, pred.method = pred.method, seed = seed, PARALLEL = PARALLEL,
                                   verbose = verbose)

  message(paste0("\nTIME FOR cv.sPLSDRCOX: ", cv.splsdrcox_penalty_res$time, "\n"))

  save_ggplot(cv.splsdrcox_penalty_res$plot_AUC, folder = aux_folder, wide = T, name = "cv.splsdrcox_penalty_res_AUC", format = "svg")
  save_ggplot(cv.splsdrcox_penalty_res$plot_BRIER, folder = aux_folder, wide = T, name = "cv.splsdrcox_penalty_res_BRIER", format = "svg")
  save_ggplot(cv.splsdrcox_penalty_res$plot_c_index, folder = aux_folder, wide = T, name = "cv.splsdrcox_penalty_res_c_index", format = "svg")
  save_ggplot(cv.splsdrcox_penalty_res$plot_AIC, folder = aux_folder, wide = T, name = "cv.splsdrcox_penalty_res_AIC", format = "svg")

  best_splsdrcox_penalty <- splsdrcox_penalty(X = X_train, Y = Y_train,
                              n.comp = cv.splsdrcox_penalty_res$opt.comp,
                              penalty = cv.splsdrcox_penalty_res$opt.penalty,
                              x.center = x.center, x.scale = x.scale,
                              remove_near_zero_variance = remove_near_zero_variance,
                              remove_zero_variance = remove_zero_variance,
                              toKeep.zv = toKeep.zv,
                              remove_non_significant = remove_non_significant,
                              alpha = alpha,
                              MIN_EPV = MIN_EPV,
                              returnData = returnData,
                              verbose = verbose)

  save(list = c("cv.splsdrcox_penalty_res", "best_splsdrcox_penalty"), file = paste0(aux_folder, "splsdrcox.RData"))
  gc()

  lst_models_full[[cv.splsdrcox_penalty_res$class]] = cv.splsdrcox_penalty_res
  lst_models_full[[best_splsdrcox_penalty$class]] = best_splsdrcox_penalty
  lst_models[[best_splsdrcox_penalty$class]] = best_splsdrcox_penalty
  lst_models_pls[[best_splsdrcox_penalty$class]] = best_splsdrcox_penalty

}

### ### ### ### ### ##
# sPLS-DRCOX_DYNAMIC #
### ### ### ### ### ##
if(FLAG_sPLSDRCOX_DYNAMIC){
  aux_folder = paste0(folder, "splsdrcox_dynamic_plot/")
  dir.create(aux_folder)

  cv.splsdrcox_dynamic_res <- cv.splsdrcox(X = X_train, Y = Y_train,
                                           max.ncomp = max.ncomp,
                                           vector = vector, MIN_NVAR = MIN_NVAR, MAX_NVAR = MAX_NVAR,
                                           n.cut_points = n.cut_points, EVAL_METHOD = EVAL_METHOD,
                                           n_run = n_run, k_folds = k_folds,
                                           alpha = alpha, remove_non_significant = remove_non_significant,
                                           times = times, max_time_points = max_time_points,
                                           remove_near_zero_variance = remove_near_zero_variance,
                                           remove_zero_variance = remove_zero_variance,
                                           remove_non_significant_models = remove_non_significant_models,
                                           toKeep.zv = toKeep.zv,
                                           w_AIC = w_AIC, w_C.Index = w_C.Index, w_I.BRIER = w_I.BRIER, w_AUC = w_AUC,
                                           MIN_AUC_INCREASE = MIN_AUC_INCREASE, MIN_AUC = MIN_AUC, MIN_COMP_TO_CHECK = MIN_COMP_TO_CHECK,
                                           x.center = x.center, x.scale = x.scale,
                                           fast_mode = fast_mode, return_models = return_models, MIN_EPV = MIN_EPV,
                                           pred.attr = pred.attr, pred.method = pred.method, seed = seed, PARALLEL = PARALLEL,
                                           verbose = verbose)

  message(paste0("\nTIME FOR cv.sPLSDRCOX_DYNAMIC: ", cv.splsdrcox_dynamic_res$time, "\n"))

  save_ggplot(cv.splsdrcox_dynamic_res$plot_AUC, folder = aux_folder, wide = T, name = "cv.splsdrcox_dynamic_res_AUC", format = "svg")
  save_ggplot(cv.splsdrcox_dynamic_res$plot_BRIER, folder = aux_folder, wide = T, name = "cv.splsdrcox_dynamic_res_BRIER", format = "svg")
  save_ggplot(cv.splsdrcox_dynamic_res$plot_c_index, folder = aux_folder, wide = T, name = "cv.splsdrcox_dynamic_res_c_index", format = "svg")
  save_ggplot(cv.splsdrcox_dynamic_res$plot_AIC, folder = aux_folder, wide = T, name = "cv.splsdrcox_dynamic_res_AIC", format = "svg")

  best_splsdrcox_dynamic <- splsdrcox(X = X_train, Y = Y_train,
                                      n.comp = cv.splsdrcox_dynamic_res$opt.comp,
                                      vector = cv.splsdrcox_dynamic_res$opt.nvar,
                                      x.center = x.center, x.scale = x.scale,
                                      remove_near_zero_variance = remove_near_zero_variance,
                                      remove_zero_variance = remove_zero_variance,
                                      toKeep.zv = toKeep.zv,
                                      remove_non_significant = remove_non_significant,
                                      alpha = alpha,
                                      MIN_EPV = MIN_EPV,
                                      returnData = returnData,
                                      verbose = verbose)

  save(list = c("cv.splsdrcox_dynamic_res", "best_splsdrcox_dynamic"), file = paste0(aux_folder, "splsdrcox_dynamic.RData"))
  gc()

  lst_models_full[[cv.splsdrcox_dynamic_res$class]] = cv.splsdrcox_dynamic_res
  lst_models_full[[best_splsdrcox_dynamic$class]] = best_splsdrcox_dynamic
  lst_models[[best_splsdrcox_dynamic$class]] = best_splsdrcox_dynamic
  lst_models_pls[[best_splsdrcox_dynamic$class]] = best_splsdrcox_dynamic

}

### ### ### ### ### ##
# sPLS-DACOX_DYNAMIC #
### ### ### ### ### ##
if(FLAG_sPLSDACOX_DYNAMIC){
  aux_folder = paste0(folder, "splsdacox_dynamic_plot/")
  dir.create(aux_folder)

  cv.splsdacox_dynamic_res <- cv.splsdacox(X = X_train, Y = Y_train,
                                           max.ncomp = max.ncomp, vector = vector,
                                           MIN_NVAR = MIN_NVAR, MAX_NVAR = MAX_NVAR, n.cut_points = n.cut_points, EVAL_METHOD = EVAL_METHOD,
                                           n_run = n_run, k_folds = k_folds,
                                           alpha = alpha, remove_non_significant = remove_non_significant, times = times, max_time_points = max_time_points,
                                           remove_near_zero_variance = remove_near_zero_variance,
                                           remove_zero_variance = remove_zero_variance,
                                           remove_non_significant_models = remove_non_significant_models,
                                           toKeep.zv = toKeep.zv,
                                           w_AIC = w_AIC, w_C.Index = w_C.Index, w_I.BRIER = w_I.BRIER, w_AUC = w_AUC,
                                           MIN_AUC_INCREASE = MIN_AUC_INCREASE, MIN_AUC = MIN_AUC, MIN_COMP_TO_CHECK = MIN_COMP_TO_CHECK,
                                           x.center = x.center, x.scale = x.scale,
                                           fast_mode = fast_mode, return_models = return_models, MIN_EPV = MIN_EPV,
                                           pred.attr = pred.attr, pred.method = pred.method, seed = seed, PARALLEL = PARALLEL,
                                           verbose = verbose)

  message(paste0("\nTIME FOR cv.sPLSDACOX_DYNAMIC: ", cv.splsdacox_dynamic_res$time, "\n"))

  save_ggplot(cv.splsdacox_dynamic_res$plot_AUC, folder = aux_folder, wide = T, name = "cv.splsdacox_res_AUC", format = "svg")
  save_ggplot(cv.splsdacox_dynamic_res$plot_BRIER, folder = aux_folder, wide = T, name = "cv.splsdacox_res_BRIER", format = "svg")
  save_ggplot(cv.splsdacox_dynamic_res$plot_c_index, folder = aux_folder, wide = T, name = "cv.splsdacox_res_c_index", format = "svg")
  save_ggplot(cv.splsdacox_dynamic_res$plot_AIC, folder = aux_folder, wide = T, name = "cv.splsdacox_res_AIC", format = "svg")

  best_splsdacox_dynamic <- splsdacox(X = X_train, Y = Y_train,
                                      n.comp = cv.splsdacox_dynamic_res$opt.comp,
                                      vector = cv.splsdacox_dynamic_res$opt.nvar,
                                      x.center = x.center, x.scale = x.scale,
                                      remove_near_zero_variance = remove_near_zero_variance,
                                      remove_zero_variance = remove_zero_variance,
                                      toKeep.zv = toKeep.zv,
                                      remove_non_significant = remove_non_significant,
                                      alpha = alpha,
                                      MIN_EPV = MIN_EPV,
                                      returnData = returnData,
                                      verbose = verbose)

  save(list = c("cv.splsdacox_dynamic_res", "best_splsdacox_dynamic"), file = paste0(aux_folder, "splsdacox_dynamic.RData"))
  gc()

  lst_models_full[[cv.splsdacox_dynamic_res$class]] = cv.splsdacox_dynamic_res
  if(!all(is.na(best_splsdacox_dynamic))){
    lst_models_full[[best_splsdacox_dynamic$class]] = best_splsdacox_dynamic
    lst_models[[best_splsdacox_dynamic$class]] = best_splsdacox_dynamic
    lst_models_pls[[best_splsdacox_dynamic$class]] = best_splsdacox_dynamic
  }
}

### ### ### ## #
# SB.sPLS-ICOX #
### ### ### ## #
if(FLAG_SB.sPLSICOX){
  aux_folder = paste0(folder, "sb.splsicox_plot/")
  dir.create(aux_folder)

  cv.sb.splsicox_res <- cv.sb.splsicox(X = X_train, Y = data.matrix(Y_train),
                                       max.ncomp = max.ncomp, penalty.list = penalty.list,
                                       n_run = n_run, k_folds = k_folds, alpha = alpha,
                                       remove_near_zero_variance = remove_near_zero_variance, remove_zero_variance = remove_zero_variance, toKeep.zv = toKeep.zv,
                                       remove_non_significant_models = remove_non_significant_models,
                                       remove_non_significant = remove_non_significant,
                                       w_AIC = w_AIC, w_C.Index = w_C.Index, w_I.BRIER = w_I.BRIER, w_AUC = w_AUC, times = times, max_time_points = max_time_points,
                                       MIN_AUC_INCREASE = MIN_AUC_INCREASE, MIN_AUC = MIN_AUC, MIN_COMP_TO_CHECK = MIN_COMP_TO_CHECK,
                                       x.center = x.center, x.scale = x.scale,
                                       fast_mode = fast_mode, return_models = return_models, MIN_EPV = MIN_EPV,
                                       pred.attr = pred.attr, pred.method = pred.method, seed = seed, PARALLEL = PARALLEL, verbose = verbose)

  message(paste0("\nTIME FOR cv.SB.sPLS-ICOX: ", cv.sb.splsicox_res$time, "\n"))

  save_ggplot(cv.sb.splsicox_res$plot_AUC, folder = aux_folder, wide = T, name = "cv.sb.splsicox_res_AUC", format = "svg")
  save_ggplot(cv.sb.splsicox_res$plot_BRIER, folder = aux_folder, wide = T, name = "cv.sb.splsicox_res_BRIER", format = "svg")
  save_ggplot(cv.sb.splsicox_res$plot_c_index, folder = aux_folder, wide = T, name = "cv.sb.splsicox_res_c_index", format = "svg")
  save_ggplot(cv.sb.splsicox_res$plot_AIC, folder = aux_folder, wide = T, name = "cv.sb.splsicox_res_AIC", format = "svg")

  best_sb.splsicox <- sb.splsicox(X = X_train,
                                Y = data.matrix(Y_train),
                                n.comp = cv.sb.splsicox_res$opt.comp, #1
                                penalty = cv.sb.splsicox_res$opt.penalty,
                                x.center = x.center, x.scale = x.scale,
                                remove_non_significant = remove_non_significant, #cox variables
                                remove_near_zero_variance = remove_near_zero_variance, remove_zero_variance = remove_zero_variance, toKeep.zv = toKeep.zv,
                                returnData = returnData, verbose = verbose, MIN_EPV = MIN_EPV)

  lst_models_full[[cv.sb.splsicox_res$class]] = cv.sb.splsicox_res
  lst_models_full[[best_sb.splsicox$class]] = best_sb.splsicox
  lst_models[[best_sb.splsicox$class]] = best_sb.splsicox
  lst_models_pls[[best_sb.splsicox$class]] = best_sb.splsicox

  save(list = c("cv.sb.splsicox_res", "best_sb.splsicox"), file = paste0(aux_folder, "sb.splsicox.RData"))
  gc()
}


### ### ### ### #
# iSB.sPLS-ICOX #
### ### ### ### #
if(FLAG_iSB.sPLSICOX){
  aux_folder = paste0(folder, "isb.splsicox_plot/")
  dir.create(aux_folder)

  cv.isb.splsicox_res <- cv.isb.splsicox(X = X_train, Y = data.matrix(Y_train),
                                  max.ncomp = max.ncomp, penalty.list = penalty.list,
                                  n_run = n_run, k_folds = k_folds, alpha = alpha, remove_non_significant_models = remove_non_significant_models,
                                  w_AIC = w_AIC, w_C.Index = w_C.Index, w_I.BRIER = w_I.BRIER, w_AUC = w_AUC, times = times, max_time_points = max_time_points,
                                  remove_near_zero_variance = remove_near_zero_variance, remove_zero_variance = remove_zero_variance, toKeep.zv = toKeep.zv,
                                  remove_non_significant = remove_non_significant,
                                  MIN_AUC_INCREASE = MIN_AUC_INCREASE, MIN_AUC = MIN_AUC, MIN_COMP_TO_CHECK = MIN_COMP_TO_CHECK,
                                  x.center = x.center, x.scale = x.scale,
                                  fast_mode = fast_mode, return_models = return_models, MIN_EPV = MIN_EPV,
                                  pred.attr = pred.attr, pred.method = pred.method, seed = seed,
                                  returnData = returnData, verbose = verbose, PARALLEL = PARALLEL)

  message(paste0("\nTIME FOR cv.iSB.sPLS-ICOX: ", cv.isb.splsicox_res$time, "\n"))

  best_isb.splsicox <- isb.splsicox(X = X_train,
                                    Y = data.matrix(Y_train),
                                    cv.isb = cv.isb.splsicox_res,
                                    x.center = x.center, x.scale = x.scale,
                                    remove_non_significant = remove_non_significant, #cox variables
                                    remove_near_zero_variance = remove_near_zero_variance, remove_zero_variance = remove_zero_variance, toKeep.zv = toKeep.zv,
                                    returnData = returnData, verbose = verbose, MIN_EPV = MIN_EPV)

  lst_models_full[[cv.isb.splsicox_res$class]] = cv.isb.splsicox_res
  lst_models_full[[best_isb.splsicox$class]] = best_isb.splsicox
  lst_models[[best_isb.splsicox$class]] = best_isb.splsicox
  lst_models_pls[[best_isb.splsicox$class]] = best_isb.splsicox

  save(list = c("cv.isb.splsicox_res","best_isb.splsicox"), file = paste0(aux_folder, "isb.splsicox.RData"))
  gc()
}

### ### ### ### ### ### #
# SB.sPLS-DRCOX-Penalty #
### ### ### ### ### ### #
if(FLAG_SB.sPLSDRCOX_PENALTY){
  aux_folder = paste0(folder, "sb.splsdrcox_penalty_plot/")
  dir.create(aux_folder)

  cv.sb.splsdrcox_penalty_res <- cv.sb.splsdrcox_penalty(X = X_train, Y = data.matrix(Y_train),
                                         max.ncomp = max.ncomp, penalty.list = penalty.list,
                                         n_run = n_run, k_folds = k_folds, alpha = alpha,
                                         remove_near_zero_variance = remove_near_zero_variance, remove_zero_variance = remove_zero_variance, toKeep.zv = toKeep.zv,
                                         remove_non_significant_models = remove_non_significant_models,
                                         remove_non_significant = remove_non_significant,
                                         w_AIC = w_AIC, w_C.Index = w_C.Index, w_I.BRIER = w_I.BRIER, w_AUC = w_AUC, times = times, max_time_points = max_time_points,
                                         MIN_AUC_INCREASE = MIN_AUC_INCREASE, MIN_AUC = MIN_AUC, MIN_COMP_TO_CHECK = MIN_COMP_TO_CHECK,
                                         x.center = x.center, x.scale = x.scale,
                                         fast_mode = fast_mode, return_models = return_models, MIN_EPV = MIN_EPV,
                                         pred.attr = pred.attr, pred.method = pred.method, seed = seed, PARALLEL = PARALLEL)

  message(paste0("\nTIME FOR cv.SB.sPLSDRCOX: ", cv.sb.splsdrcox_penalty_res$time, "\n"))

  save_ggplot(cv.sb.splsdrcox_penalty_res$plot_AUC, folder = aux_folder, wide = T, name = "cv.sb.splsdrcox_penalty_res_AUC", format = "svg")
  save_ggplot(cv.sb.splsdrcox_penalty_res$plot_BRIER, folder = aux_folder, wide = T, name = "cv.sb.splsdrcox_penalty_res_BRIER", format = "svg")
  save_ggplot(cv.sb.splsdrcox_penalty_res$plot_c_index, folder = aux_folder, wide = T, name = "cv.sb.splsdrcox_penalty_res_c_index", format = "svg")
  save_ggplot(cv.sb.splsdrcox_penalty_res$plot_AIC, folder = aux_folder, wide = T, name = "cv.sb.splsdrcox_penalty_res_AIC", format = "svg")

  best_sb.splsdrcox_penalty <- sb.splsdrcox_penalty(X = X_train,
                                    Y = data.matrix(Y_train),
                                    n.comp = cv.sb.splsdrcox_penalty_res$opt.comp,
                                    penalty = cv.sb.splsdrcox_penalty_res$opt.penalty,
                                    x.center = x.center, x.scale = x.scale, alpha = alpha,
                                    remove_non_significant = remove_non_significant,
                                    remove_near_zero_variance = remove_near_zero_variance, remove_zero_variance = remove_zero_variance, toKeep.zv = toKeep.zv,
                                    returnData = returnData, MIN_EPV = MIN_EPV, verbose = verbose)

  lst_models_full[[cv.sb.splsdrcox_penalty_res$class]] = cv.sb.splsdrcox_penalty_res
  lst_models_full[[best_sb.splsdrcox_penalty$class]] = best_sb.splsdrcox_penalty
  lst_models[[best_sb.splsdrcox_penalty$class]] = best_sb.splsdrcox_penalty
  lst_models_pls[[best_sb.splsdrcox_penalty$class]] = best_sb.splsdrcox_penalty

  save(list = c("cv.sb.splsdrcox_penalty_res", "best_sb.splsdrcox_penalty"), file = paste0(aux_folder, "sb.splsdrcox.penalty.RData"))
  gc()
}

### ### ### ### ### ### #
# SB.sPLS-DRCOX-Dynamic #
### ### ### ### ### ### #
if(FLAG_SB.sPLSDRCOX_DYNAMIC){
  aux_folder = paste0(folder, "sb.splsdrcox_dynamic_plot/")
  dir.create(aux_folder)

  cv.sb.splsdrcox_dynamic_res <- cv.sb.splsdrcox(X = X_train, Y = data.matrix(Y_train),
                                         max.ncomp = max.ncomp,
                                         vector = vector, MIN_NVAR = MIN_NVAR, MAX_NVAR = MAX_NVAR,
                                         n.cut_points = n.cut_points, EVAL_METHOD = EVAL_METHOD,
                                         n_run = n_run, k_folds = k_folds,
                                         remove_near_zero_variance = remove_near_zero_variance, remove_zero_variance = remove_zero_variance, toKeep.zv = toKeep.zv,
                                         remove_non_significant = remove_non_significant,
                                         alpha = alpha, remove_non_significant_models = remove_non_significant_models,
                                         w_AIC = w_AIC, w_C.Index = w_C.Index, w_I.BRIER = w_I.BRIER, w_AUC = w_AUC, times = times, max_time_points = max_time_points,
                                         MIN_AUC_INCREASE = MIN_AUC_INCREASE, MIN_AUC = MIN_AUC, MIN_COMP_TO_CHECK = MIN_COMP_TO_CHECK,
                                         x.center = x.center, x.scale = x.scale,
                                         fast_mode = fast_mode, return_models = return_models, MIN_EPV = MIN_EPV,
                                         pred.attr = pred.attr, pred.method = pred.method, seed = seed, PARALLEL = PARALLEL)

  message(paste0("\nTIME FOR cv.SB.sPLSDRCOX: ", cv.sb.splsdrcox_dynamic_res$time, "\n"))

  save_ggplot(cv.sb.splsdrcox_dynamic_res$plot_AUC, folder = aux_folder, wide = T, name = "cv.sb.splsdrcox_dynamic_res_AUC", format = "svg")
  save_ggplot(cv.sb.splsdrcox_dynamic_res$plot_BRIER, folder = aux_folder, wide = T, name = "cv.sb.splsdrcox_dynamic_res_BRIER", format = "svg")
  save_ggplot(cv.sb.splsdrcox_dynamic_res$plot_c_index, folder = aux_folder, wide = T, name = "cv.sb.splsdrcox_dynamic_res_c_index", format = "svg")
  save_ggplot(cv.sb.splsdrcox_dynamic_res$plot_AIC, folder = aux_folder, wide = T, name = "cv.sb.splsdrcox_dynamic_res_AIC", format = "svg")

  best_sb.splsdrcox_dynamic <- sb.splsdrcox(X = X_train,
                                    Y = data.matrix(Y_train),
                                    n.comp = cv.sb.splsdrcox_dynamic_res$opt.comp,
                                    vector = cv.sb.splsdrcox_dynamic_res$opt.nvar,
                                    x.center = x.center, x.scale = x.scale,
                                    remove_non_significant = remove_non_significant,
                                    remove_near_zero_variance = remove_near_zero_variance, remove_zero_variance = remove_zero_variance, toKeep.zv = toKeep.zv,
                                    returnData = returnData, MIN_EPV = MIN_EPV, verbose = verbose)

  lst_models_full[[cv.sb.splsdrcox_dynamic_res$class]] = cv.sb.splsdrcox_dynamic_res
  lst_models_full[[best_sb.splsdrcox_dynamic$class]] = best_sb.splsdrcox_dynamic
  lst_models[[best_sb.splsdrcox_dynamic$class]] = best_sb.splsdrcox_dynamic
  lst_models_pls[[best_sb.splsdrcox_dynamic$class]] = best_sb.splsdrcox_dynamic

  save(list = c("cv.sb.splsdrcox_dynamic_res", "best_sb.splsdrcox_dynamic"), file = paste0(aux_folder, "sb.splsdrcox.dynamic.RData"))
  gc()
}


### ### ### ### ### ### ##
# iSB.sPLS-DRCOX-Penalty #
### ### ### ### ### ### ##
if(FLAG_iSB.sPLSDRCOX_PENALTY){
  aux_folder = paste0(folder, "isb.splsdrcox_penalty_plot/")
  dir.create(aux_folder)

  cv.isb.splsdrcox_penalty_res <- cv.isb.splsdrcox_penalty(X = X_train, Y = data.matrix(Y_train),
                                                           max.ncomp = max.ncomp, penalty.list = penalty.list,
                                                           n_run = n_run, k_folds = k_folds, alpha = alpha, remove_non_significant_models = remove_non_significant_models,
                                                           w_AIC = w_AIC, w_C.Index = w_C.Index, w_I.BRIER = w_I.BRIER, w_AUC = w_AUC, times = times, max_time_points = max_time_points,
                                                           remove_near_zero_variance = remove_near_zero_variance, remove_zero_variance = remove_zero_variance, toKeep.zv = toKeep.zv,
                                                           remove_non_significant = remove_non_significant,
                                                           MIN_AUC_INCREASE = MIN_AUC_INCREASE, MIN_AUC = MIN_AUC, MIN_COMP_TO_CHECK = MIN_COMP_TO_CHECK,
                                                           x.center = x.center, x.scale = x.scale,
                                                           fast_mode = fast_mode, return_models = return_models, MIN_EPV = MIN_EPV,
                                                           pred.attr = pred.attr, pred.method = pred.method, seed = seed,
                                                           returnData = T, verbose = T, PARALLEL = PARALLEL)

  message(paste0("\nTIME FOR cv.iSB.sPLSDRCOX: ", cv.isb.splsdrcox_penalty_res$time, "\n"))

  best_isb.splsdrcox_penalty <- isb.splsdrcox_penalty(X = X_train,
                                                      Y = data.matrix(Y_train),
                                                      cv.isb = cv.isb.splsdrcox_penalty_res, alpha = alpha,
                                                      x.center = x.center, x.scale = x.scale,
                                                      remove_non_significant = remove_non_significant, #cox variables
                                                      remove_near_zero_variance = remove_near_zero_variance, remove_zero_variance = remove_zero_variance, toKeep.zv = toKeep.zv,
                                                      returnData = returnData, verbose = verbose, MIN_EPV = MIN_EPV)

  lst_models_full[[cv.isb.splsdrcox_penalty_res$class]] = cv.isb.splsdrcox_penalty_res
  lst_models_full[[best_isb.splsdrcox_penalty$class]] = best_isb.splsdrcox_penalty
  lst_models[[best_isb.splsdrcox_penalty$class]] = best_isb.splsdrcox_penalty
  lst_models_pls[[best_isb.splsdrcox_penalty$class]] = best_isb.splsdrcox_penalty

  save(list = c("cv.isb.splsdrcox_penalty_res","best_isb.splsdrcox_penalty"), file = paste0(aux_folder, "isb.splsdrcox.penalty.RData"))
  gc()
}

### ### ### ### ### ### ##
# iSB.sPLS-DRCOX-Dynamic #
### ### ### ### ### ### ##
if(FLAG_iSB.sPLSDRCOX_DYNAMIC){
  aux_folder = paste0(folder, "isb.splsdrcox_dynamic_plot/")
  dir.create(aux_folder)

  cv.isb.splsdrcox_dynamic_res <- cv.isb.splsdrcox(X = X_train, Y = data.matrix(Y_train),
                                                   max.ncomp = max.ncomp,
                                                   vector = vector, MIN_NVAR = MIN_NVAR, MAX_NVAR = MAX_NVAR,
                                                   n.cut_points = n.cut_points, EVAL_METHOD = EVAL_METHOD,
                                                   n_run = n_run, k_folds = k_folds,
                                                   alpha = alpha, remove_non_significant_models = remove_non_significant_models,
                                                   w_AIC = w_AIC, w_C.Index = w_C.Index, w_I.BRIER = w_I.BRIER, w_AUC = w_AUC, times = times, max_time_points = max_time_points,
                                                   remove_near_zero_variance = remove_near_zero_variance, remove_zero_variance = remove_zero_variance, toKeep.zv = toKeep.zv,
                                                   remove_non_significant = remove_non_significant,
                                                   MIN_AUC_INCREASE = MIN_AUC_INCREASE, MIN_AUC = MIN_AUC, MIN_COMP_TO_CHECK = MIN_COMP_TO_CHECK,
                                                   x.center = x.center, x.scale = x.scale,
                                                   fast_mode = fast_mode, return_models = return_models, MIN_EPV = MIN_EPV,
                                                   pred.attr = pred.attr, pred.method = pred.method, seed = seed,
                                                   returnData = T, verbose = T, PARALLEL = PARALLEL)

  message(paste0("\nTIME FOR cv.iSB.sPLSDRCOX: ", cv.isb.splsdrcox_dynamic_res$time, "\n"))

  best_isb.splsdrcox_dynamic <- isb.splsdrcox(X = X_train,
                                              Y = data.matrix(Y_train),
                                              cv.isb = cv.isb.splsdrcox_dynamic_res, alpha = alpha,
                                              x.center = x.center, x.scale = x.scale,
                                              remove_non_significant = remove_non_significant, #cox variables
                                              remove_near_zero_variance = remove_near_zero_variance, remove_zero_variance = remove_zero_variance, toKeep.zv = toKeep.zv,
                                              returnData = returnData, verbose = verbose, MIN_EPV = MIN_EPV)

  lst_models_full[[cv.isb.splsdrcox_dynamic_res$class]] = cv.isb.splsdrcox_dynamic_res
  lst_models_full[[best_isb.splsdrcox_dynamic$class]] = best_isb.splsdrcox_dynamic
  lst_models[[best_isb.splsdrcox_dynamic$class]] = best_isb.splsdrcox_dynamic
  lst_models_pls[[best_isb.splsdrcox_dynamic$class]] = best_isb.splsdrcox_dynamic

  save(list = c("cv.isb.splsdrcox_dynamic_res","best_isb.splsdrcox_dynamic"), file = paste0(aux_folder, "isb.splsdrcox.dynamic.RData"))
  gc()
}

### ### ### ### ### ### ##
# SB.sPLS-DACOX-Dynamic #
### ### ### ### ### ### ##
if(FLAG_SB.sPLSDACOX){
  aux_folder = paste0(folder, "sb.splsdacox_dynamic_plot/")
  dir.create(aux_folder)

  cv.sb.splsdacox_dynamic_res <- cv.sb.splsdacox(X = X_train, Y = data.matrix(Y_train),
                                                 max.ncomp = max.ncomp, vector = vector, MIN_NVAR = MIN_NVAR, MAX_NVAR = MAX_NVAR,
                                                 n.cut_points = n.cut_points, EVAL_METHOD = EVAL_METHOD,
                                                 n_run = n_run, k_folds = k_folds, alpha = alpha, remove_non_significant_models = remove_non_significant_models,
                                                 w_AIC = w_AIC, w_C.Index = w_C.Index, w_I.BRIER = w_I.BRIER, w_AUC = w_AUC, times = times, max_time_points = max_time_points,
                                                 remove_near_zero_variance = remove_near_zero_variance, remove_zero_variance = remove_zero_variance, toKeep.zv = toKeep.zv,
                                                 remove_non_significant = remove_non_significant,
                                                 MIN_AUC_INCREASE = MIN_AUC_INCREASE, MIN_AUC = MIN_AUC, MIN_COMP_TO_CHECK = MIN_COMP_TO_CHECK,
                                                 x.center = x.center, x.scale = x.scale,
                                                 fast_mode = fast_mode, return_models = return_models, MIN_EPV = MIN_EPV,
                                                 pred.attr = pred.attr, pred.method = pred.method, seed = seed,
                                                 returnData = T, verbose = T, PARALLEL = PARALLEL)

  message(paste0("\nTIME FOR cv.SB.sPLSDACOX: ", cv.sb.splsdacox_dynamic_res$time, "\n"))

  best_sb.splsdacox_dynamic <- sb.splsdacox(X = X_train,
                                            Y = data.matrix(Y_train),
                                            n.comp = cv.sb.splsdacox_dynamic_res$opt.comp, #2
                                            vector = cv.sb.splsdacox_dynamic_res$opt.nvar,
                                            alpha = alpha,
                                            x.center = x.center, x.scale = x.scale,
                                            remove_non_significant = remove_non_significant, #cox variables
                                            remove_near_zero_variance = remove_near_zero_variance, remove_zero_variance = remove_zero_variance, toKeep.zv = toKeep.zv,
                                            returnData = returnData, verbose = verbose, MIN_EPV = MIN_EPV)

  lst_models_full[[cv.sb.splsdacox_dynamic_res$class]] = cv.sb.splsdacox_dynamic_res
  lst_models_full[[best_sb.splsdacox_dynamic$class]] = best_sb.splsdacox_dynamic
  lst_models[[best_sb.splsdacox_dynamic$class]] = best_sb.splsdacox_dynamic
  lst_models_pls[[best_sb.splsdacox_dynamic$class]] = best_sb.splsdacox_dynamic

  save(list = c("cv.sb.splsdacox_dynamic_res","best_sb.splsdacox_dynamic"), file = paste0(aux_folder, "sb.splsdacox.dynamic.RData"))
  gc()
}

### ### ### ### ### ### ##
# iSB.sPLS-DACOX-Dynamic #
### ### ### ### ### ### ##
if(FLAG_iSB.sPLSDACOX){
  aux_folder = paste0(folder, "isb.splsdacox_dynamic_plot/")
  dir.create(aux_folder)

  cv.isb.splsdacox_dynamic_res <- cv.isb.splsdacox(X = X_train, Y = data.matrix(Y_train),
                                                   max.ncomp = max.ncomp, vector = vector, MIN_NVAR = MIN_NVAR, MAX_NVAR = MAX_NVAR,
                                                   n.cut_points = n.cut_points, EVAL_METHOD = EVAL_METHOD,
                                                   n_run = n_run, k_folds = k_folds, alpha = alpha, remove_non_significant_models = remove_non_significant_models,
                                                   w_AIC = w_AIC, w_C.Index = w_C.Index, w_I.BRIER = w_I.BRIER, w_AUC = w_AUC, times = times, max_time_points = max_time_points,
                                                   remove_near_zero_variance = remove_near_zero_variance, remove_zero_variance = remove_zero_variance, toKeep.zv = toKeep.zv,
                                                   remove_non_significant = remove_non_significant,
                                                   MIN_AUC_INCREASE = MIN_AUC_INCREASE, MIN_AUC = MIN_AUC, MIN_COMP_TO_CHECK = MIN_COMP_TO_CHECK,
                                                   x.center = x.center, x.scale = x.scale,
                                                   fast_mode = fast_mode, return_models = return_models, MIN_EPV = MIN_EPV,
                                                   pred.attr = pred.attr, pred.method = pred.method, seed = seed,
                                                   returnData = T, verbose = T, PARALLEL = PARALLEL)

  message(paste0("\nTIME FOR cv.iSB.sPLSDACOX: ", cv.isb.splsdacox_dynamic_res$time, "\n"))

  best_isb.splsdacox_dynamic <- isb.splsdacox(X = X_train,
                                              Y = data.matrix(Y_train),
                                              cv.isb = cv.isb.splsdacox_dynamic_res, alpha = alpha,
                                              x.center = x.center, x.scale = x.scale,
                                              remove_non_significant = remove_non_significant, #cox variables
                                              remove_near_zero_variance = remove_near_zero_variance, remove_zero_variance = remove_zero_variance, toKeep.zv = toKeep.zv,
                                              returnData = returnData, verbose = verbose, MIN_EPV = MIN_EPV)

  lst_models_full[[cv.isb.splsdacox_dynamic_res$class]] = cv.isb.splsdacox_dynamic_res
  lst_models_full[[best_isb.splsdacox_dynamic$class]] = best_isb.splsdacox_dynamic
  lst_models[[best_isb.splsdacox_dynamic$class]] = best_isb.splsdacox_dynamic
  lst_models_pls[[best_isb.splsdacox_dynamic$class]] = best_isb.splsdacox_dynamic

  save(list = c("cv.isb.splsdacox_dynamic_res","best_isb.splsdacox_dynamic"), file = paste0(aux_folder, "isb.splsdacox.dynamic.RData"))
  gc()
}

### ### ### ### ## ### ##
# MO.sPLSDRCOX_MixOmics #
### ### ### ### ## ### ##
if(FLAG_MB.sPLSDRCOX){
  aux_folder = paste0(folder, "mb.splsdrcox_dynamic_plot/")
  dir.create(aux_folder)

  cv.mb.splsdrcox_dynamic_res <- cv.mb.splsdrcox(X = X_train, Y = data.matrix(Y_train),
                                                 max.ncomp = max.ncomp, n_run = n_run, k_folds = k_folds,
                                                 vector = vector, design = design, MIN_NVAR = MIN_NVAR, MAX_NVAR = MAX_NVAR, n.cut_points = n.cut_points, EVAL_METHOD = EVAL_METHOD,
                                                 remove_near_zero_variance = remove_near_zero_variance, remove_zero_variance = remove_zero_variance, toKeep.zv = toKeep.zv,
                                                 remove_non_significant = remove_non_significant,
                                                 alpha = alpha, remove_non_significant_models = remove_non_significant_models,
                                                 w_AIC = w_AIC, w_C.Index = w_C.Index, w_I.BRIER = w_I.BRIER, w_AUC = w_AUC, times = times, max_time_points = max_time_points,
                                                 MIN_AUC_INCREASE = MIN_AUC_INCREASE, MIN_AUC = MIN_AUC, MIN_COMP_TO_CHECK = MIN_COMP_TO_CHECK,
                                                 x.center = x.center, x.scale = x.scale,
                                                 fast_mode = fast_mode, return_models = return_models, MIN_EPV = MIN_EPV,
                                                 pred.attr = pred.attr, pred.method = pred.method, seed = seed, PARALLEL = PARALLEL, verbose = verbose)

  message(paste0("\nTIME FOR cv.MB.sPLSDRCOX_DYNAMIC: ", cv.mb.splsdrcox_dynamic_res$time, "\n"))

  save_ggplot(cv.mb.splsdrcox_dynamic_res$plot_AUC, folder = aux_folder, wide = T, name = "cv.mb.splsdrcox_dynamic_res_AUC", format = "svg")
  save_ggplot(cv.mb.splsdrcox_dynamic_res$plot_BRIER, folder = aux_folder, wide = T, name = "cv.mb.splsdrcox_dynamic_res_BRIER", format = "svg")
  save_ggplot(cv.mb.splsdrcox_dynamic_res$plot_c_index, folder = aux_folder, wide = T, name = "cv.mb.splsdrcox_dynamic_res_c_index", format = "svg")
  save_ggplot(cv.mb.splsdrcox_dynamic_res$plot_AIC, folder = aux_folder, wide = T, name = "cv.mb.splsdrcox_dynamic_res_AIC", format = "svg")

  best_mb.splsdrcox_dynamic <- mb.splsdrcox(X = X_train,
                                            Y = Y_train,
                                            n.comp = cv.mb.splsdrcox_dynamic_res$opt.comp,
                                            vector = cv.mb.splsdrcox_dynamic_res$opt.nvar,
                                            design = cv.mb.splsdrcox_dynamic_res$design,
                                            x.center = x.center, x.scale = x.scale,
                                            remove_non_significant = remove_non_significant,
                                            remove_near_zero_variance = remove_near_zero_variance, remove_zero_variance = remove_zero_variance, toKeep.zv = toKeep.zv,
                                            pred.method = pred.method, MIN_EPV = MIN_EPV,
                                            returnData = returnData, verbose = verbose)

  lst_models_full[[cv.mb.splsdrcox_dynamic_res$class]] = cv.mb.splsdrcox_dynamic_res
  lst_models_full[[best_mb.splsdrcox_dynamic$class]] = best_mb.splsdrcox_dynamic
  lst_models[[best_mb.splsdrcox_dynamic$class]] = best_mb.splsdrcox_dynamic
  lst_models_pls[[best_mb.splsdrcox_dynamic$class]] = best_mb.splsdrcox_dynamic

  save(list = c("cv.mb.splsdrcox_dynamic_res", "best_mb.splsdrcox_dynamic"), file = paste0(aux_folder, "mb.splsdrcox.RData"))
  gc()

}

### ### ### ### ## ### ##
# MO.sPLSDACOX_MixOmics #
### ### ### ### ## ### ##

if(FLAG_MB.sPLSDACOX){
  aux_folder = paste0(folder, "mb.splsdacox_dynamic_plot/")
  dir.create(aux_folder)

  cv.mb.splsdacox_dynamic_res <- cv.mb.splsdacox(X = X_train, Y = data.matrix(Y_train),
                                                 max.ncomp = max.ncomp, n_run = n_run, k_folds = k_folds,
                                                 vector = vector, design = design, MIN_NVAR = MIN_NVAR, MAX_NVAR = MAX_NVAR, n.cut_points = n.cut_points, EVAL_METHOD = EVAL_METHOD,
                                                 remove_near_zero_variance = remove_near_zero_variance, remove_zero_variance = remove_zero_variance, toKeep.zv = toKeep.zv,
                                                 remove_non_significant = remove_non_significant,
                                                 alpha = alpha, remove_non_significant_models = remove_non_significant_models,
                                                 w_AIC = w_AIC, w_C.Index = w_C.Index, w_I.BRIER = w_I.BRIER, w_AUC = w_AUC, times = times, max_time_points = max_time_points,
                                                 MIN_AUC_INCREASE = MIN_AUC_INCREASE, MIN_AUC = MIN_AUC, MIN_COMP_TO_CHECK = MIN_COMP_TO_CHECK,
                                                 x.center = x.center, x.scale = x.scale,
                                                 #y.center = y.center, y.scale = y.scale,
                                                 fast_mode = fast_mode, return_models = return_models, MIN_EPV = MIN_EPV,
                                                 pred.attr = pred.attr, pred.method = pred.method, seed = seed, PARALLEL = PARALLEL, verbose = verbose)

  message(paste0("\nTIME FOR cv.MB.sPLSDACOX_DYNAMIC: ", cv.mb.splsdacox_dynamic_res$time, "\n"))

  save_ggplot(cv.mb.splsdacox_dynamic_res$plot_AUC, folder = aux_folder, wide = T, name = "cv.mb.splsdacox_dynamic_res_AUC", format = "svg")
  save_ggplot(cv.mb.splsdacox_dynamic_res$plot_BRIER, folder = aux_folder, wide = T, name = "cv.mb.splsdacox_dynamic_res_BRIER", format = "svg")
  save_ggplot(cv.mb.splsdacox_dynamic_res$plot_c_index, folder = aux_folder, wide = T, name = "cv.mb.splsdacox_dynamic_res_c_index", format = "svg")
  save_ggplot(cv.mb.splsdacox_dynamic_res$plot_AIC, folder = aux_folder, wide = T, name = "cv.mb.splsdacox_dynamic_res_AIC", format = "svg")

  best_mb.splsdacox_dynamic <- mb.splsdacox(X = X_train,
                                            Y = Y_train,
                                            n.comp = cv.mb.splsdacox_dynamic_res$opt.comp,
                                            vector = cv.mb.splsdacox_dynamic_res$opt.nvar,
                                            design = cv.mb.splsdacox_dynamic_res$design,
                                            x.center = x.center, x.scale = x.scale,
                                            #y.center = y.center, y.scale = y.scale,
                                            remove_non_significant = remove_non_significant,
                                            remove_near_zero_variance = remove_near_zero_variance, remove_zero_variance = remove_zero_variance, toKeep.zv = toKeep.zv,
                                            pred.method = pred.method, MIN_EPV = MIN_EPV,
                                            returnData = returnData, verbose = verbose)

  lst_models_full[[cv.mb.splsdacox_dynamic_res$class]] = cv.mb.splsdacox_dynamic_res
  lst_models_full[[best_mb.splsdacox_dynamic$class]] = best_mb.splsdacox_dynamic
  lst_models[[best_mb.splsdacox_dynamic$class]] = best_mb.splsdacox_dynamic
  lst_models_pls[[best_mb.splsdacox_dynamic$class]] = best_mb.splsdacox_dynamic

  save(list = c("cv.mb.splsdacox_dynamic_res", "best_mb.splsdacox_dynamic"), file = paste0(aux_folder, "mb.splsdacox.RData"))
  gc()

}

## ## ## ## ##
# SAVE MODELS #
## ## ## ## ##

message("\n## COMPUTING SAVING MODELS ##\n")

save(list = c("lst_models_full"), file = paste0(folder, "lst_models_full.RData")) # CV, dummy and pls
save(list = c("lst_models"), file = paste0(folder, "lst_models.RData")) # dummy and pls
save(list = c("lst_models_dummy"), file = paste0(folder, "lst_models_dummy.RData")) # dummy
save(list = c("lst_models_pls"), file = paste0(folder, "lst_models_pls.RData")) # pls


## ## ## ## ##
# EVALUATION #
## ## ## ## ##

message("\n## COMPUTING EVALUATION ##\n")

eval_results <- purrr::map(lst_evaluations, ~eval_Coxmos_models(lst_models = lst_models,
                                                                X_test = X_test, Y_test = Y_test, pred.method = .,
                                                                pred.attr = pred.attr, verbose = verbose,
                                                                times = times, max_time_points = max_time_points, PARALLEL = PARALLEL))

save(list = c("eval_results"), file = paste0(folder, "eval_results.RData")) # eval_results
gc()

#### ### ## #
# EVALPLOTS #
#### ### ## #

message("\n## COMPUTING EVALPLOTS AUC ##\n")

evaluation_folder = paste0(folder, "evaluation_plot_AUC/")
dir.create(evaluation_folder)

lst_evaluation_plot <- plot_evaluation.list(lst_eval_results = eval_results, evaluation = "AUC")

for(eval_name in "cenROC"){
  if(all(is.na(lst_evaluation_plot[[eval_name]]))){
    next
  }

  write.csv(x = lst_evaluation_plot[[eval_name]]$df, file = paste0(evaluation_folder,"results_", eval_name, ".csv"), row.names = FALSE, quote = F)
  openxlsx::write.xlsx(x = lst_evaluation_plot[[eval_name]]$df, file = paste0(evaluation_folder,"results_", eval_name, ".xlsx"))
  #save evaluation plots
  save_ggplot(plot = lst_evaluation_plot[[eval_name]]$lst_plots$lineplot.mean, folder = evaluation_folder, name = paste0(eval_name, "_", "cenROC"), wide = T, format = "svg")

}

save(list = c("eval_results", "lst_evaluation_plot"), file = paste0(evaluation_folder, "eval_results.RData"))
gc()

#### ### ## #
# EVALPLOTS #
#### ### ## #

message("\n## COMPUTING EVALPLOTS BRIER ##\n")

evaluation_folder = paste0(folder, "evaluation_plot_Brier/")
dir.create(evaluation_folder)

lst_evaluation_plot <- plot_evaluation.list(eval_results, evaluation = "IBS")

for(eval_name in names(lst_evaluation_plot)){
  if(all(is.na(lst_evaluation_plot[[eval_name]]))){
    next
  }
  write.csv(x = lst_evaluation_plot[[eval_name]]$df, file = paste0(evaluation_folder,"results_", eval_name, ".csv"), row.names = FALSE, quote = F)
  openxlsx::write.xlsx(x = lst_evaluation_plot[[eval_name]]$df, file = paste0(evaluation_folder,"results_", eval_name, ".xlsx"))
  #save evaluation plots
  for(comp_name in names(lst_evaluation_plot[[eval_name]]$lst_plots)){
    save_ggplot(plot = lst_evaluation_plot[[eval_name]]$lst_plots[[comp_name]], folder = evaluation_folder, name = paste0(eval_name, "_", comp_name), wide = T, format = "svg")
  }
  #save test plots
  for(comp_name in names(lst_evaluation_plot[[eval_name]]$lst_plot_comparisons)){
    save_ggplot(plot = lst_evaluation_plot[[eval_name]]$lst_plot_comparisons[[comp_name]], folder = evaluation_folder, name = paste0("comparison_", eval_name, "_", comp_name), wide = T, format = "svg")
  }
}

save(list = c("eval_results", "lst_evaluation_plot"), file = paste0(evaluation_folder, "eval_results.RData"))
gc()

## ## ##
# TIME #
## ## ##

message("\n## COMPUTING TIME ##\n")

time_folder = paste0(folder, "time_plot/")
dir.create(time_folder)

ggp_time <- plot_time.list(lst_models_full)

save_ggplot(ggp_time, folder = time_folder, name = "time", format = "svg", wide = T)

save(list = c("ggp_time"), file = paste0(time_folder, "time_results.RData"))

gc()

## ## ## ## ## ##
## FOREST PLOTS #
## ## ## ## ## ##

message("\n## COMPUTING FOREST PLOTS ##\n")

forest_folder = paste0(folder, "forest_plot/")
dir.create(forest_folder)

lst_forest_plot <- plot_forest.list(lst_models = lst_models)
save_ggplot_lst(lst_plots = lst_forest_plot, folder = forest_folder, wide = T, format = "svg", prefix = "forest_")

save(list = c("lst_forest_plot"), file = paste0(forest_folder, "forest_results.RData"))

gc()

## ## ## ## ## ## ## ## #
# PH - DIAGNOSTIC PLOTS #
## ## ## ## ## ## ## ## #

message("\n## COMPUTING PH - DIAGNOSTIC PLOTS ##\n")

ph_folder = paste0(folder, "ph_plot/")
dir.create(ph_folder)

lst_ph_ggplot <- plot_proportionalHazard.list(lst_models = lst_models)

save_ggplot_lst(lst_plots = lst_ph_ggplot, folder = ph_folder, wide = T, format = "svg", prefix = "ph_")

save(list = c("lst_ph_ggplot"), file = paste0(ph_folder, "ph_results.RData"))

gc()

## ## ## ##
## BIPLOT #
## ## ## ##

message("\n## COMPUTING BIPLOTS - PLS ##\n")

biplot_folder = paste0(folder, "biplot_plot/")
dir.create(biplot_folder)

biplot.list <- purrr::map(.x = lst_models_pls, ~plot_sPLS_Coxmos(model = ., comp = c(1,2), mode = "biplot",
                                                                top = 20, overlaps = 20, only_top = TRUE, names = TRUE))
for(m in names(biplot.list)){

  if(is.null(biplot.list[[m]])){
    next
  }else{
    for(b in names(biplot.list[[m]]$plot_block)){
      if("plot" %in% names(biplot.list[[m]]$plot_block[[b]])){
        save_ggplot(plot = biplot.list[[m]]$plot_block[[b]]$plot, format = "svg",
                    folder = biplot_folder, name = paste0(m, "_", b, "_biplot_1vs2"))
      }else{
        save_ggplot(plot = biplot.list[[m]]$plot_block[[b]], format = "svg",
                    folder = biplot_folder, name = paste0(m, "_", b, "_biplot_1vs2"))
      }
    }
  }
}

save(list = c("biplot.list"), file = paste0(biplot_folder, "biplot_results.RData"))

## ## ## ## ## ##
# DENSITY PLOTS #
## ## ## ## ## ##

message("\n## COMPUTING DENSITY PLOTS ##\n")

density_folder = paste0(folder, "density_plot/")
dir.create(density_folder)

density.plots.lp <- plot_cox.event.list(lst_models, type = "lp")
density.plots.risk <- plot_cox.event.list(lst_models, type = "risk")
density.plots.expected <- plot_cox.event.list(lst_models, type = "expected")
density.plots.survival <- plot_cox.event.list(lst_models, type = "survival")

for(cn in names(density.plots.lp)){

  if(is.null(density.plots.lp[[cn]])){
    next
  }

  if(cn %in% names(density.plots.lp) && !is.na(density.plots.lp[[cn]])){
    save_ggplot(plot = density.plots.lp[[cn]]$plot.density, format = "svg", folder = density_folder, name = paste0("pred.lp_density_",cn))
    save_ggplot(plot = density.plots.lp[[cn]]$plot.histogram, format = "svg", folder = density_folder, name = paste0("pred.lp_histogram_",cn))
  }

  if(cn %in% names(density.plots.risk) && !is.na(density.plots.risk[[cn]])){
    save_ggplot(plot = density.plots.risk[[cn]]$plot.density, format = "svg", folder = density_folder, name = paste0("pred.risk_density_",cn))
    save_ggplot(plot = density.plots.risk[[cn]]$plot.histogram, format = "svg", folder = density_folder, name = paste0("pred.risk_histogram_",cn))
  }

  if(cn %in% names(density.plots.expected) && !is.na(density.plots.expected[[cn]])){
    save_ggplot(plot = density.plots.expected[[cn]]$plot.density, format = "svg", folder = density_folder, name = paste0("pred.expected_density_",cn))
    save_ggplot(plot = density.plots.expected[[cn]]$plot.histogram, format = "svg", folder = density_folder, name = paste0("pred.expected_histogram_",cn))
  }

  if(cn %in% names(density.plots.survival) && !is.na(density.plots.survival[[cn]])){
    save_ggplot(plot = density.plots.survival[[cn]]$plot.density, format = "svg", folder = density_folder, name = paste0("pred.expected_density_",cn))
    save_ggplot(plot = density.plots.survival[[cn]]$plot.histogram, format = "svg", folder = density_folder, name = paste0("pred.expected_histogram_",cn))
  }
}

save(list = c("density.plots.lp", "density.plots.risk", "density.plots.expected", "density.plots.survival"), file = paste0(density_folder, "density_results.RData"))

gc()

## ## ## ## ## #
## EVENT PLOTS #
## ## ## ## ## #

message("\n## COMPUTING EVENT PLOTS ##\n")

event_folder = paste0(folder, "event_plot/")
dir.create(event_folder)

ggp_density.event <- plot_events(Y = rbind(Y_train, Y_test),
                                 categories = c("Censored","Death"), #name for FALSE/0 (Censored) and TRUE/1 (Event)
                                 y.text = "Number of observations",
                                 roundTo = 0.5,
                                 max.breaks = 15)

save_ggplot(plot = ggp_density.event$plot, folder = event_folder, name = "density_events")

save(list = c("ggp_density.event"), file = paste0(event_folder, "event_results.RData"))

gc()

## ## ## ## ## ## #
## PSEUDO-BETA COX #
## ## ## ## ## ## #

message("\n## COMPUTING PSEUDO-BETA COX ##\n")

psbeta_folder = paste0(folder, "pbetacox_plot/")
dir.create(psbeta_folder)

ggp.simulated_beta <- plot_pseudobeta.list(lst_models = lst_models_pls,
                                           error.bar = T, onlySig = T, alpha = 0.05,
                                           zero.rm = T, auto.limits = T, top = 20,
                                           show_percentage = T, size_percentage = 3, verbose = F)

if(!data_type %in% c("3","MO")){
  save_ggplot_lst(lst_plots = ggp.simulated_beta, object_name = "plot", folder = psbeta_folder, wide = T, format = "svg", prefix = "pbetacox_")
}else{
  for(m in names(ggp.simulated_beta)){
    for(b in names(X_train)){
      if(b %in% names(ggp.simulated_beta[[m]]$plot)){
        save_ggplot(plot = ggp.simulated_beta[[m]]$plot[[b]], folder = psbeta_folder, wide = T, name = paste0("pbetacox_",m,"_",b), format = "svg")
      }
    }
  }
}

save(list = c("ggp.simulated_beta"), file = paste0(psbeta_folder, "pseudobeta_results.RData"))

gc()

## ## ## ## ## #
# KAPLAN-MEIER #
## ## ## ## ## #

message("\n## COMPUTING KAPLAN-MEIER ##\n")

km_folder = paste0(folder, "km_train_plot/")
dir.create(km_folder)

LST_KM_RES_LP <- getAutoKM.list(type = "LP",
                                lst_models = lst_models,
                                comp = max.ncomp,
                                top = 20,
                                ori_data = F,
                                BREAKTIME = NULL, minProp = minProp,
                                only_sig = T, alpha = 0.05)

LST_KM_RES_COMP <- getAutoKM.list(type = "COMP",
                                  lst_models = lst_models,
                                  comp = max.ncomp,
                                  top = 20,
                                  ori_data = F,
                                  BREAKTIME = NULL, minProp = minProp,
                                  only_sig = T, alpha = 0.05)

LST_KM_RES_VAR <- getAutoKM.list(type = "VAR",
                                 lst_models = lst_models,
                                 comp = max.ncomp,
                                 top = 20,
                                 ori_data = TRUE,
                                 BREAKTIME = NULL, minProp = minProp,
                                 only_sig = T, alpha = 0.05)

for(m in names(LST_KM_RES_LP)){
  if(!all(is.na(LST_KM_RES_LP[[m]]))){
    save_ggplot_lst(lst_plots = LST_KM_RES_LP[[m]]$LST_PLOTS, object_name = NULL, folder = km_folder, wide = T, format = "svg", prefix = paste0("km_", m, "_"))
  }
  if(m %in% names(LST_KM_RES_COMP) && !all(is.na(LST_KM_RES_COMP[[m]]))){
    if(!data_type %in% c("3","MO")){
      save_ggplot_lst(lst_plots = LST_KM_RES_COMP[[m]]$LST_PLOTS, object_name = NULL, folder = km_folder, wide = T, format = "svg", prefix = paste0("km_", m, "_"))
    }else{
      for(b in names(LST_KM_RES_COMP[[m]]$LST_PLOTS)){
        save_ggplot_lst(lst_plots = LST_KM_RES_COMP[[m]]$LST_PLOTS[[b]], object_name = NULL, folder = km_folder, wide = T, format = "svg", prefix = paste0("km_", m, "_"))
      }
    }
  }
  if(m %in% names(LST_KM_RES_VAR) && !all(is.na(LST_KM_RES_VAR[[m]]))){
    if(!data_type %in% c("3","MO")){
      save_ggplot_lst(lst_plots = LST_KM_RES_VAR[[m]]$LST_PLOTS, object_name = NULL, folder = km_folder, wide = T, format = "svg", prefix = paste0("km_", m, "_"))
    }else{
      for(b in names(LST_KM_RES_COMP[[m]]$LST_PLOTS)){
        save_ggplot_lst(lst_plots = LST_KM_RES_VAR[[m]]$LST_PLOTS[[b]], object_name = NULL, folder = km_folder, wide = T, format = "svg", prefix = paste0("km_", m, "_"))
      }
    }
  }
}

save(list = c("LST_KM_RES_LP", "LST_KM_RES_COMP", "LST_KM_RES_VAR"), file = paste0(km_folder, "KM_results.RData"))

## ## ## ## ## ## ## #
# Kaplan-Meier TEST  #
## ## ## ## ## ## ## #

message("\n## COMPUTING KAPLAN-MEIER TESTING ##\n")

km_folder = paste0(folder, "km_test_plot/")
dir.create(km_folder)

lst_cutoff <- getCutoffAutoKM.list(LST_KM_RES_LP)
LST_KM_TEST_LP <- getTestKM.list(lst_models = lst_models,
                                 X_test = X_test, Y_test = Y_test,
                                 type = "LP",
                                 ori_data = FALSE,
                                 BREAKTIME = NULL, n.breaks = 20, title = "Test Data",
                                 lst_cutoff = lst_cutoff)

lst_cutoff <- getCutoffAutoKM.list(LST_KM_RES_COMP)
LST_KM_TEST_COMP <- getTestKM.list(lst_models = lst_models,
                                   X_test = X_test, Y_test = Y_test,
                                   type = "COMP",
                                   ori_data = FALSE,
                                   BREAKTIME = NULL, n.breaks = 20, title = "Test Data",
                                   lst_cutoff = lst_cutoff)

lst_cutoff <- getCutoffAutoKM.list(LST_KM_RES_VAR)
LST_KM_TEST_VAR <- getTestKM.list(lst_models = lst_models,
                                  X_test = X_test, Y_test = Y_test,
                                  type = "VAR",
                                  ori_data = TRUE,
                                  BREAKTIME = NULL, n.breaks = 20, title = "Test Data",
                                  lst_cutoff = lst_cutoff)

for(m in names(LST_KM_TEST_LP)){

  if(!all(is.na(LST_KM_TEST_LP[[m]]))){
    save_ggplot(plot = LST_KM_TEST_LP[[m]], folder = km_folder, wide = T, name = paste0("km_test_", m))
  }

  if(m %in% names(LST_KM_TEST_COMP) && !all(is.na(LST_KM_TEST_COMP[[m]]))){
    save_ggplot_lst(lst_plots = LST_KM_TEST_COMP[[m]], object_name = NULL, folder = km_folder, wide = T, format = "svg", prefix = paste0("km_", m, "_"))
  }

  if(m %in% names(LST_KM_TEST_VAR) && !all(is.na(LST_KM_TEST_VAR[[m]]))){
    if(!data_type %in% c("3","MO")){
      save_ggplot_lst(lst_plots = LST_KM_TEST_VAR[[m]], object_name = NULL, folder = km_folder, wide = T, format = "svg", prefix = paste0("km_", m, "_"))
    }else{
      for(b in names(LST_KM_RES_COMP[[m]]$LST_PLOTS)){
        save_ggplot_lst(lst_plots = LST_KM_TEST_VAR[[m]][[b]], object_name = NULL, folder = km_folder, wide = T, format = "svg", prefix = paste0("km_", m, "_"))
      }
    }
  }
}

save(list = c("LST_KM_TEST_LP", "LST_KM_TEST_COMP", "LST_KM_TEST_VAR"), file = paste0(km_folder, "KM_test_results.RData"))

## ## ## ## ## ## ## ## ## ## #
# AUC relevance per Variable  #
## ## ## ## ## ## ## ## ## ## #

message("\n## COMPUTING AUC relevance per Variable ##\n")

aucVariable_folder = paste0(folder, "AUC_per_Variable/")
dir.create(aucVariable_folder)

LST_AUC_VARIABLE_RES_VAR_TRAIN <- eval_Coxmos_model_per_variable.list(lst_models = lst_models,
                                                                X_train, Y_train,
                                                                pred.method = "cenROC", pred.attr = "mean",
                                                                times = NULL, max_time_points = 15,
                                                                PARALLEL = FALSE, verbose = FALSE)

LST_PLOT_AUC_VARIABLE_RES_VAR_TRAIN <- plot_evaluation.list(LST_AUC_VARIABLE_RES_VAR_TRAIN, evaluation = "AUC")

LST_AUC_VARIABLE_RES_VAR_TEST <- eval_Coxmos_model_per_variable.list(lst_models = lst_models,
                                                                 X_test, Y_test,
                                                                 pred.method = "cenROC", pred.attr = "mean",
                                                                 times = NULL, max_time_points = 15,
                                                                 PARALLEL = FALSE, verbose = FALSE)


LST_PLOT_AUC_VARIABLE_RES_VAR_TEST <- plot_evaluation.list(LST_AUC_VARIABLE_RES_VAR_TEST, evaluation = "AUC")

for(m in names(LST_PLOT_AUC_VARIABLE_RES_VAR_TEST)){

  if(!all(is.na(LST_PLOT_AUC_VARIABLE_RES_VAR_TRAIN[[m]]$lst_plots$lineplot.mean))){
    save_ggplot(plot = LST_PLOT_AUC_VARIABLE_RES_VAR_TRAIN[[m]]$lst_plots$lineplot.mean,
                folder = aucVariable_folder, wide = T, name = paste0("AUC_Variable_train_", m))
  }

  if(!all(is.na(LST_PLOT_AUC_VARIABLE_RES_VAR_TEST[[m]]$lst_plots$lineplot.mean))){
    save_ggplot(plot = LST_PLOT_AUC_VARIABLE_RES_VAR_TEST[[m]]$lst_plots$lineplot.mean,
                folder = aucVariable_folder, wide = T, name = paste0("AUC_Variable_test_", m))
  }
}

save(list = c("LST_AUC_VARIABLE_RES_VAR_TRAIN", "LST_PLOT_AUC_VARIABLE_RES_VAR_TRAIN",
              "LST_AUC_VARIABLE_RES_VAR_TEST", "LST_PLOT_AUC_VARIABLE_RES_VAR_TEST"), file = paste0(aucVariable_folder, "AUC_variable_results.RData"))

## ## ## ## ## ##
## SAVE RESULTS #
## ## ## ## ## ##

message("\n## COMPUTING SAVING ALL RESULTS ##\n")

##Save all results
save.image(paste0(folder,"output.performance_",todaydate,".RData"))

message("The script has successfully finished!")
