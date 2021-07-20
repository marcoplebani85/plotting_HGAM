# code by Marco Plebani - July 19th, 2021



rm(list=ls()) # wipe R's memory clean
library(pacman) # load packages, installing them from CRAN if needed
p_load(RCurl) # allows accessing data from URL
p_load(mgcv) # for HGAMs
ss <- read.delim(text=getURL("https://raw.githubusercontent.com/marcoplebani85/datasets/master/flower_color_spectra.txt"))
head(ss)
ss$density <- ifelse(ss$density<0, 0, ss$density) # set spurious negative reflectance values to zero
ss$clr <- ifelse(ss$Taxon=="SpeciesB", "red", "black")
ss <- with(ss, ss[order(Locality, wl), ])






# GS: models with a global smoother for all observations, plus group-level smoothers, all with the same wiggliness 
start_time <- Sys.time()
gam_GS1 <- mgcv::bam(density ~ Taxon # main effect
		+ s(wl, by = Taxon, k = 20) # interaction
		+ s(wl, by = Locality, bs="fs", m=1),
		# "fs" is short for "factor-smoother [interaction]"
		family="quasipoisson",
		data = ss, method = 'REML'
		)
end_time <- Sys.time()
end_time - start_time
# gam.check(gam_GS1)
# k.check(gam_GS1)
# MuMIn::AICc(gam_GS1)


# Predict gam_GS1

# by Locality

nn <- unique(ss[, c("wl", "Taxon", "Locality", "clr")])
pred <- predict(object= gam_GS1, newdata=nn, type="response", se.fit=T)
nn$fit <- pred$fit
nn$se <- pred$se.fit
	
nn0 <- subset(nn, Locality==unique(nn$Locality)[1])

par(mfrow=c(1,2))
plot(fit ~ wl, xlab="Wavelength (nm)", ylab="Mean reflectance (GAM fit)",
	data = nn0,
	ylim=c(0,100),
	type="l", lwd=2, col=scales::alpha(unique(nn0$clr), 0.8)
	)
yysep<- nn0$fit+1.96*nn0$se
yysem<- nn0$fit-1.96*nn0$se
# Confidence interval as a shaded area:
polygon(c(nn0$wl,rev(nn0$wl)),c(yysep,rev(yysem)),col=scales::alpha(unique(nn0$clr), 0.3),
	border=NA)
	
for(i in 2:length(unique(nn$Locality))){
nn0 <- subset(nn, Locality==unique(nn$Locality)[i])
lines(fit ~ wl,
	data = nn0,
	type="l", lwd=2, col=scales::alpha(unique(nn0$clr), 0.8)
	)
yysep<- nn0$fit+1.96*nn0$se
yysem<- nn0$fit-1.96*nn0$se
# Confidence interval as a shaded area:
polygon(c(nn0$wl,rev(nn0$wl)),c(yysep,rev(yysem)),
	col=scales::alpha(unique(nn0$clr), 0.3),
	border=NA
	)
}

# gam_GS1 by Taxon
	
# gratia::draw(gam_GS1)
# plot(gam_GS1, pages=1)
nn <- unique(ss[, c("wl", 
				"Taxon", 
				"Locality",
				"clr")])
nn$Locality=0 # turns random effect off
# https://stats.stackexchange.com/questions/131106/predicting-with-random-effects-in-mgcv-gam
pred <- predict(object = gam_GS1, 
				type="response", 
				# exclude="c(Locality)", # turns random effect off, allegedly
				newdata=nn, 
				se.fit=T)
nn$fit <- pred$fit
nn$se <- pred$se.fit

nn0 <- subset(nn, Taxon==unique(nn$Taxon)[1])

plot(fit ~ wl, xlab="Wavelength (nm)", ylab="Mean reflectance (GAM fit)",
	data = nn0,
	ylim=c(0,100),
	type="l", lwd=2, col=scales::alpha(unique(nn0$clr), 0.8)
	)
yysep<- nn0$fit+1.96*nn0$se
yysem<- nn0$fit-1.96*nn0$se
# Confidence interval as a shaded area:
polygon(c(nn0$wl,rev(nn0$wl)),c(yysep,rev(yysem)),col=scales::alpha(unique(nn0$clr), 0.3),
	border=NA)
	
for(i in 2:length(unique(nn$Taxon))){
nn0 <- subset(nn, Taxon ==unique(nn$Taxon)[i])
lines(fit ~ wl,
	data = nn0,
	type="l", lwd=2, col=scales::alpha(unique(nn0$clr), 0.8)
	)
yysep<- nn0$fit+1.96*nn0$se
yysem<- nn0$fit-1.96*nn0$se
# Confidence interval as a shaded area:
polygon(c(nn0$wl,rev(nn0$wl)),c(yysep,rev(yysem)),
	col=scales::alpha(unique(nn0$clr), 0.3),
	border=NA
	)
}
