# code by Marco Plebani - July 19th, 2021



rm(list=ls()) # wipe R's memory clean
library(pacman) # load packages, installing them from CRAN if needed
p_load(RCurl) # allows accessing data from URL
ss <- read.delim(text=getURL("https://raw.githubusercontent.com/marcoplebani85/datasets/master/flower_color_spectra.txt"))
head(ss)
ss$density <- ifelse(ss$density<0, 0, ss$density) # set spurious negative reflectance values to zero
ss$clr <- ifelse(ss$Taxon=="SpeciesB", "red", "black")
ss <- with(ss, ss[order(Locality, wl), ])

# calculate population means

ss.smr <- doBy::summaryBy(density ~ wl + Locality + Taxon + clr,
						FUN=c(mean),
						data=ss
						)

names(ss.smr) <- c("wl", "Locality", "Taxon", "clr", "dens_mean")
ss.smr$LT <- paste(ss.smr$Locality, ss.smr$Taxon, sep="_")






# plot rolling means per population

ddd <- subset(ss.smr, LT ==unique(ss.smr$LT)[1])
plot(zoo::rollmean(ddd$dens_mean, k=5), 
		ylim=c(0, 80),
		data=ddd, xlab="Wavelength [nm]", ylab="Reflectance [%]",
		col=scales::alpha(ddd$clr, 0.5),
		type="l", lwd=3
		)	
for(i in 2:length(unique(ss.smr$LT))){
	ddd <- subset(ss.smr, LT ==unique(ss.smr$LT)[i])
	points(zoo::rollmean(ddd$dens_mean, k=5),
		data=ddd,
		col=scales::alpha(ddd$clr, 0.5),
		type="l", lwd=3
		)
}






# Plot HGAMs:

gam_G1 <- mgcv::bam(density ~ Taxon # main effect
		+ s(wl, by = Taxon, k = 20) # interaction
		+ s(Locality, bs="re"), # "re" is short for "random effect"
		data = ss, method = 'REML',
		family="quasipoisson"
		)
# gam.check(gam_G1)
# k.check(gam_G1)	
MuMIn::AICc(gam_G1)

# Predict gam_G1

# by Locality

nn <- unique(ss[, c("wl", "Taxon", "Locality", "clr")])
pred <- predict(object= gam_G1, newdata=nn, type="response", se.fit=T)
nn$fit <- pred$fit
nn$se <- pred$se.fit
	
nn0 <- subset(nn, Locality==unique(nn$Locality)[1])
par(mfrow=c(1,2))
plot(fit ~ wl, xlab="Wavelength (nm)", ylab="Mean reflectance (GAM fit)",
	data = nn0,
	ylim=c(-20,100),
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

# G1 by Taxon
	
# gratia::draw(gam_G1)
# plot(gam_G1, pages=1)
nn <- unique(ss[, c("wl", 
				"Taxon", 
				"Locality",
				"clr")])
nn$Locality=0 # turn the random effects "off"
# https://stats.stackexchange.com/questions/131106/predicting-with-random-effects-in-mgcv-gam
pred <- predict(object = gam_G1, 
				type="response", 
				newdata=nn, 
				se.fit=T)
nn$fit <- pred$fit
nn$se <- pred$se.fit
nn <- with(nn, nn[order(Taxon, wl), ])
nn0 <- subset(nn, Taxon==unique(nn$Taxon)[1])
plot(fit ~ wl, xlab="Wavelength (nm)", ylab="Mean reflectance (GAM fit)",
	data = nn0,
	ylim=c(-20,100),
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