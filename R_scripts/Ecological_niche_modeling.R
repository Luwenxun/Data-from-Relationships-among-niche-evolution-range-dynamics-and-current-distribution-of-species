#Ecological niche modeling
#Load the required packages for the analysis
library(raster)
library(biomod2)
library(SDMtune)
# Acquire environmental variables
bio01 = raster("D:/climate_data/bio01.tif")
bio07 = raster("D:/climate_data/bio07.tif")
bio12 = raster("D:/climate_data/bio12.tif")
bio15 = raster("D:/climate_data/bio15.tif")
Clim_Cur <- stack(bio01, bio07, bio12, bio15)
#Generate pseudo-absences for the whole C. zawadskii complex
occ.cza.com <- read.delim("cza_complex.csv", h=T, sep=",")
cza.pa.com <- BIOMOD_FormatingData(resp.var = rep(1, nrow(occ.cza.com)),
                     expl.var = Clim_Cur,
                     resp.xy = occ.cza.com[,2:3],
                     resp.name = 'cza.com',
                     PA.nb.rep = 1,
                     PA.nb.absences = 10000,
                     PA.strategy = 'random',
                     PA.table = NULL,
                     na.rm = TRUE)
write.csv(cza.pa.com@coord, file="cza-com-pa.csv")
#Tune hyperparameters of the BRT model (using C. zawadskii as an example)
#Prepare presence and background locations
cza.p <- read.delim("D:/occ_point/obs-cza-thin.csv", h=T, sep=",")#presence
cza.a <- read.delim("D:/occ_point/pa.csv", h=T, sep=",")#background
#Create an SWD object
cza.tune <- prepareSWD(species = "cza", p = cza.p, a = cza.a, env = Clim_Cur)
cza.datasets <- trainValTest(cza.tune, test = 0.2, seed = 25)
cza.train <- cza.datasets[[1]]
cza.test <- cza.datasets[[2]]
cza.folds <- randomFolds(cza.train, k = 4, seed = 25)
#Tune model hyperparameters
BRT.h <- list(distribution='bernoulli', bag.fraction = 0.5, shrinkage = c(0.01, 0.005, 0.001), n.trees = seq(5000, 10000, 500),  interaction.depth = c(1, 3, 5, 7, 9))
cza.BRT.om <- optimizeModel(cza.BRT.model, hypers = BRT.h, metric = "auc", seed = 25)
cza.BRT.best_model <- cza.BRT.om@models[[1]]
cza.BRT.om@results[1, ]
write.csv(cza.BRT.om@results[1, ], file= "cza.BRT-tuning.csv")
#Evaluate the final model
set.seed(25)
cza.BRT.final_model <- train("BRT", data = cza.train, distribution = cza.BRT.om@results[1, 1], n.trees = cza.BRT.om@results[1, 2], interaction.depth = cza.BRT.om@results[1, 3], shrinkage = cza.BRT.om@results[1, 4], bag.fraction=cza.BRT.om@results[1, 5])
plotROC(cza.BRT.final_model, test = cza.test)
#Generate BRT model for  C. zawadskii 
# Load species occurrences and pseudo-absences
cza.pa <- read.delim("D:/evolutionaryrescue/multi_EM_model/occ/pa-cza.csv", h=T, sep=",")
#Format Data with true absences and pseudo-absences
cza.Resp <- as.numeric(cza.pa[, 'Occ'])
cza.XY<- cza.pa[, c('Longitude','Latitude')]
cza.data <- BIOMOD_FormatingData(resp.var = cza.Resp,
	expl.var = Clim_Cur,
	resp.xy = cza.XY,
	resp.name = 'cza')
#Set modeling options
cza.opt <- BIOMOD_ModelingOptions(GBM=list(shrinkage = 0.001, n.trees = 7000,  interaction.depth = 9, n.cores = 3))
# Model single models
cza.model <- BIOMOD_Modeling(bm.format = cza.data,
                                    modeling.id = 'cza_gbm',
                                    models = c('GBM'),
                                    bm.options = cza.opt,
                                    nb.rep = 10,
                                    data.split.perc = 80,
                                    metric.eval = c('TSS','ROC','KAPPA'),
                                    var.import = 3,
                                    do.full.models = FALSE,
                                    seed.val = 42,
                                    nb.cpu= 3)
#Model ensemble models
cza.EM <- BIOMOD_EnsembleModeling(bm.mod = cza.model,
                                      models.chosen = 'all',
                                      em.by = 'all',
                                      metric.select = c('TSS'),
                                      metric.select.thresh = c(0.7),
                                      metric.eval = c('KAPPA', 'TSS', 'ROC'),
                                      var.import = 3,
                                      prob.mean = TRUE,
                                      prob.ci.alpha = 0.05,
                                      committee.averaging = TRUE,
                                      prob.mean.weight.decay = 'proportional',
                                      seed.val = 42)
# Get evaluation scores
cza.EM.evaluations <- get_evaluations(cza.EM, as.data.frame = TRUE)
write.csv(cza.EM.evaluations, file= "cza.EM.evaluations.csv")
#Project ensemble models
cza.EMProj <- BIOMOD_EnsembleForecasting(bm.em = cza.EM,
                                             proj.name = 'cza.CurrentEM',
                                             new.env = Clim_Cur,
                                             models.chosen = 'all',
                                             metric.binary = 'TSS',
                                             metric.filter = 'TSS')
