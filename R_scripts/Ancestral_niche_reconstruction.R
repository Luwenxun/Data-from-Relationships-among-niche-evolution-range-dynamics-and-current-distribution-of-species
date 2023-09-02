#Ancestral niche reconstruction
#Test phylosignal
#Load the required packages for the analysis
library(phytools)
library(phylosignal)
library(phylobase)
library(scales)
library(geiger)
#Load tree file and climatic data
cza_tree <- read.nexus("cza_tree.treefile")
bioclim <- read.delim("cza_com_clim.csv", h=T, sep=",")
#Test phylosignal
cza.data <- phylo4d(cza_tree, tip.data = bioclim[,2:17])
barplot.phylo4d(cza.data, tree.type = "phylo", tree.ladderize = FALSE)
cza.phySignal <- phyloSignal(cza.data = cza.data, method = "all")
write.csv(cza.phySignal, file="cza.phySignal.csv")
cza.phylosim <- phyloSim(tree = cza_tree, method = "all", nsim = 1000, reps = 99)
plot(cza.phylosim, stacked.methods = FALSE, quantiles = c(0.05, 0.95))
plot.phylosim(cza.phylosim, what = "pval", stacked.methods = TRUE)
#Test evolution models
fitCZA=function(trait=c("MAT","MWMT","MCMT","TD","MAP","AHM","DD_0","DD5","DD_18","DD18","NFFD","PAS","EMT","EXT","Eref","CMD")){

	trait=match.arg(trait, c("MAT","MWMT","MCMT","TD","MAP","AHM","DD_0","DD5","DD_18","DD18","NFFD","PAS","EMT","EXT","Eref","CMD"))

	# define set of models to compare
	models=c("BM", "OU", "EB", "white")
	summaries=c("diffusion", "Ornstein-Uhlenbeck", "early burst", "white noise")

	## ESTIMATING measurement error ##
	aic.se=numeric(length(models))
	lnl.se=numeric(length(models))

	for(m in 1:length(models)){
		cat("\n\n\n\n\t*** ", paste(toupper(summaries[m]),": fitting ", sep=""), models[m],
			" with SE *** \n", sep="")
		tmp=fitContinuous(cza_tree,cza_clim_mean[,trait],SE=NA, model=models[m],
                                    bounds=list(SE=c(0,0.5)), ncores=2)
		print(tmp)
		aic.se[m]=tmp$opt$aicc
		lnl.se[m]=tmp$opt$lnL
	}


	## ASSUMING no measurement error ##
	aic=numeric(length(models))
	lnl=numeric(length(models))

	for(m in 1:length(models)){
		cat("\n\n\n\n\t*** ", paste(toupper(summaries[m]),": fitting ", sep=""), models[m],
			 " *** \n", sep="")
		tmp=fitContinuous(cza_tree,cza_clim_mean[,trait],SE=0,model=models[m], ncores=2)
		print(tmp)
		aic[m]=tmp$opt$aicc
		lnl[m]=tmp$opt$lnL
	}

	## COMPARE AIC ##
	names(aic.se)<-names(lnl.se)<-names(aic)<-names(lnl)<-models
	delta_aic<-function(x) x-x[which(x==min(x))]

	# no measurement error
	daic=delta_aic(aic)
	cat("\n\n\n\t\t\t\t*** MODEL COMPARISON: ",trait," *** \n",sep="")
	cat("\tdelta-AIC values for models assuming no measurement error
    \t\t\t\t zero indicates the best model\n\n")
	print(daic, digits=2)

		# measurement error
	daic.se=delta_aic(aic.se)
	cat("\n\n\n\n\t\t\t\t*** MODEL COMPARISON: ",trait," ***\n",sep="")
	cat("\t\t   delta-AIC values for models estimating SE
    \t\t\t\t zero indicates the best model\n\n")
	print(daic.se, digits=2)
	cat("\n\n\n")

	res_aicc=rbind(aic, aic.se, daic, daic.se)
	rownames(res_aicc)=c("AICc","AICc_SE","dAICc", "dAICc_SE")

	return(res_aicc)
}
res_MAT=fitCZA("MAT")
print(res_MAT)
#Ancestral niche reconstruction (using Mean annual temperature as an example)
#Mean annual temperature
MAT<-as.matrix(cza_clim_mean)[,1]
fitBM_MAT <- anc.ML(cza_tree, MAT, model = "BM")
phenogram(cza_tree,c(MAT,fitBM_MAT$ace), col=alpha("red", 0.5), xlab="Mya", ylab="Mean annual temperature")

