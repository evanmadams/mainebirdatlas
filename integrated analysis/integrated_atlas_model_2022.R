require(R2jags)
require(coda)
require(reshape2)
#require(car)
require(stringr)
require(plyr)
require(dplyr)
require(ggplot2)
require(snow)
require(rgdal)
require(rgeos)
require(sp)
require(mgcv)
require(matrixcalc)
require(spdep)
require(igraph)
require(parallel)
require(nimble)
require(bayesplot)
require(raster)
#require(readr)

#select the species you want to model
spp <- 'Black-throated Green Warbler'

#lets add atlas data as a covariate to the previous pcar model

#atlas data (not needed currently, just importing data from the Data vizualzation section)

# bba2 <- readOGR(dsn = getwd(), layer = 'atlas_block_breeding_code_join')
# 
# bba2

#####BEGIN MODEL code#########

#build a distance matrix for all Blocks

amap <- readOGR(dsn = 'data/blockupdate2019', layer = 'BBA_Blocks_20190418_OldManIsland')

#convert to utm

amap <- spTransform(amap,  '+proj=utm +zone=19 +datum=NAD83 +units=m +no_defs +ellps=GRS80
                    +towgs84=0,0,0')

#dissolve to quad level and make a new polygon
#block scale kriging is not working quickly

amap@data$QuadCode <- factor(gsub('[0-9]+', '', amap@data$QuadNameCo))

amap@data$QuadName <- factor(str_sub(amap@data$CLOBlockNa, end = -4))

quadmap <- gUnaryUnion(amap, id = amap@data$QuadName)

quad.dat <- data.frame(quadID = 1:length(levels(amap@data$QuadName)), QuadCode = levels(amap@data$QuadName))
rownames(quad.dat) <- levels(amap@data$QuadName)

quadmap <- SpatialPolygonsDataFrame(quadmap, quad.dat)
cents <- coordinates(quadmap)
quadmap.dat <- quadmap@data

plot(quadmap)

#remove two quads that don't have any adjacent quads

quadmap <- quadmap[c(-717, -718), ]
quadmap@data$quadID <- 1:nrow(quadmap@data)

#create adjancy list from polygon data

#create adjacency list from polygon data

adj <- poly2nb(quadmap)

adjm <- nb2mat(adj, style = 'B')

adj.sparse <- unlist(adj)

num <- unlist(lapply(adj, length))

pcdat.eff <- read.csv("data/pc_effort_2018_22.csv")
pcdat.eff$Date <- as.POSIXct(pcdat.eff$Bird.Survey.Date..MM.DD.YYYY, format = '%m/%d/%Y')

pcdat.obs <- read.csv("data/pc_obs_2018_22.csv")

eff <- pcdat.eff
obs <- pcdat.obs

obs <- obs[which(is.na(obs$First.Minute.Detected) == FALSE),]

levels(obs$Distance..meters.)

obs$dist <- recode(obs$Distance..meters., '<5' = '5', '<10' = '10', '`25'= '25', '>100'= '100', '>150' = '150', '>250' = '250', '>300' = '300', '>350' = '350', '>400' = '400', '>500' = '500')
obs$dist <- as.numeric(as.character(obs$dist))

hist(obs$dist)

#recode observer id's

eff$Observer <- recode(eff$Observer, 'MATTHEW DICKINSON' = 'Matthew Dickinson', 'C WEST' = 'Christopher West')

#create a site index and covariates for effort data

eff$PC_ID <- as.factor(paste0(eff$Atlas.Block.Name, eff$Point.Count..))

#create a site index for observations and link to effort index

obs$PC_ID <- as.factor(paste0(obs$Atlas.Block.Name, obs$Point.Count..))

#obs <- droplevels(obs)

sites <- levels(eff$PC_ID)
sites <- data.frame(sites, siteID = 1:length(sites))

write.csv(sites, 'site_key.csv')

mean(levels(eff$PC_ID) %in% levels(sites$sites))

#not all obs IDs are found in the sites -- need to check this
mean(levels(obs$PC_ID) %in% levels(sites$sites))

eff <- merge(eff, sites, by.x = 'PC_ID', by.y = 'sites', all.x = TRUE, sort = FALSE)
obs <- merge(obs, sites, by.x = 'PC_ID', by.y = 'sites', all.x = TRUE, sort = FALSE)

obs.m <- merge(eff, obs, by = 'Sampling.Event.UNIQUE.ID', all.x = TRUE)

write.csv(obs.m, 'effort_linked_obs.csv')

eff.dups <- dcast(eff, PC_ID ~ 'N', fun.aggregate = length)

dups <- eff.dups[which(eff.dups$N > 1),]

#View(obs[which(obs$PC_ID == as.character(dups$PC_ID[1])), ])
#View(obs[which(obs$PC_ID == as.character(dups$PC_ID[2])), ])

# #remove the duplicate surveys from eff and obs if needed
# rem <- c('LE1164', 'LE1163', 'LE1162', 'LE1158')
# 
# eff <- eff[which(!(eff$Sampling.Event.UNIQUE.ID %in% rem)), ]
# obs <- obs[which(!(obs$Sampling.Event.UNIQUE.ID %in% rem)), ]

# site.obsID <- unique(data.frame(siteID = eff$siteID))
# site.obsID <- site.obsID[order(site.obsID$siteID),]

#pull the count data for the species in question


obs <- obs[which(obs$Species == spp),]
eff <- droplevels(eff)
obs <- droplevels(obs)

#old code/not needed anymore???
#remove the sites visited multiple times -- just picking the ones completed by all the same person
#removing the later visits in eff
#eff <- eff[c(-1176, -1172, -1173, -1919, -1920),]

#removing the later visits in obs


#####now process the data and summarize for JAGS###############

#turn them into distance bins
db1 <- cut(obs$dist, c(0, 20, 40, 60, 80, 100))
plot(db1)
db1 <- as.numeric(db1)

#assign observers to each obs

obsList <- levels(droplevels(eff[which(eff$PC_ID %in% obs$PC_ID), 'Observer']))
obsList <- data.frame(Observer = obsList, obsID = 1:length(obsList))

eff <- merge(eff, obsList, by = 'Observer', all.x = TRUE, sort = FALSE)

eff <- eff[order(eff$siteID),]

eff$QuadCode <- factor(str_sub(eff$Atlas.Block.Name, start = 0, end = -4))
eff$Atlas.Block.Name <- str_replace_all(eff$Atlas.Block.Name, ' ', '_')

#old code for linking, maybe remove when we are sure we don't need it...
#amap.dat$BlockCode <- str_replace_all(amap.dat$BlockCode, ' ', '_')

#eff.blockID <- left_join(eff, amap.dat, by = c('Atlas.Block.Name' = 'BlockCode'))
#eff$BlockID <- eff.blockID$blockID

eff$QuadCode <- str_replace_all(eff$QuadCode, ' ', '_')
quad.dat$QuadCode <- str_replace_all(quad.dat$QuadCode, ' ', '_')

eff.quadID <- left_join(eff, quad.dat, by = c('QuadCode'))

eff$QuadID <- eff.quadID$quadID

#add time of day and time of year variables to the effort data

eff$date <- as.POSIXlt(eff$Bird.Survey.Date..MM.DD.YYYY, format = "%m/%d/%Y")
eff$ydate <-eff$date$yday

eff$time <- as.POSIXlt(eff$Bird.Survey.start.time..in.HH.MM.24.hour.format, format = "%H:%M")
eff$ytime <- eff$time - as.POSIXlt('4:00', format = '%H:%M')

#pull in habitat variables

covs <- read.csv('mebba_envcovs.csv')

eff$PC_ID2 <- str_replace_all(eff$PC_ID, ' ', '_')
covs$PC_ID2 <- str_replace_all(covs$PC_ID, ' ', '_')

eff.covs <- left_join(eff, covs, by = c('PC_ID2'))

#we are missing some sites, for now we are going to replace with mean values, but we need to add all sites eventually
#looks like a few names are missing. In the future I'll be using the actual lat/lons to pull in data but we won't do that until all points have been surveyed

eff.covs$X11 <- ifelse(is.na(eff.covs$X11), round(mean(eff.covs$X11, na.rm = TRUE)), eff.covs$X11)
eff.covs$X21 <- ifelse(is.na(eff.covs$X21), round(mean(eff.covs$X21, na.rm = TRUE)), eff.covs$X21)
eff.covs$X22 <- ifelse(is.na(eff.covs$X22), round(mean(eff.covs$X22, na.rm = TRUE)), eff.covs$X22)
eff.covs$X23 <- ifelse(is.na(eff.covs$X23), round(mean(eff.covs$X23, na.rm = TRUE)), eff.covs$X23)
eff.covs$X24 <- ifelse(is.na(eff.covs$X24), round(mean(eff.covs$X24, na.rm = TRUE)), eff.covs$X24)
eff.covs$X31 <- ifelse(is.na(eff.covs$X31), round(mean(eff.covs$X31, na.rm = TRUE)), eff.covs$X31)
eff.covs$X41 <- ifelse(is.na(eff.covs$X41), round(mean(eff.covs$X41, na.rm = TRUE)), eff.covs$X41)
eff.covs$X42 <- ifelse(is.na(eff.covs$X42), round(mean(eff.covs$X42, na.rm = TRUE)), eff.covs$X42)
eff.covs$X43 <- ifelse(is.na(eff.covs$X43), round(mean(eff.covs$X43, na.rm = TRUE)), eff.covs$X43)
eff.covs$X52 <- ifelse(is.na(eff.covs$X52), round(mean(eff.covs$X52, na.rm = TRUE)), eff.covs$X52)
eff.covs$X71 <- ifelse(is.na(eff.covs$X71), round(mean(eff.covs$X71, na.rm = TRUE)), eff.covs$X71)
eff.covs$X81 <- ifelse(is.na(eff.covs$X81), round(mean(eff.covs$X81, na.rm = TRUE)), eff.covs$X81)
eff.covs$X82 <- ifelse(is.na(eff.covs$X82), round(mean(eff.covs$X82, na.rm = TRUE)), eff.covs$X82)
eff.covs$X90 <- ifelse(is.na(eff.covs$X90), round(mean(eff.covs$X90, na.rm = TRUE)), eff.covs$X90)
eff.covs$X95 <- ifelse(is.na(eff.covs$X95), round(mean(eff.covs$X95, na.rm = TRUE)), eff.covs$X95)
eff.covs$ch <- ifelse(is.na(eff.covs$ch), mean(eff.covs$ch, na.rm = TRUE), eff.covs$ch)
eff.covs$imp <- ifelse(is.na(eff.covs$imp), mean(eff.covs$imp, na.rm = TRUE), eff.covs$imp)
eff.covs$ele <- ifelse(is.na(eff.covs$ele), mean(eff.covs$ele, na.rm = TRUE), eff.covs$ele)

#derive some variables
eff.covs$cult <- eff.covs$X81 + eff.covs$X82
eff.covs$wet <- eff.covs$X90 + eff.covs$X95

envcovs <- as.matrix(eff.covs[, c('X21', 'X41', 'X42', 'X43', 'X52', 'X71', 'cult', 'wet', 'ch', 'imp', 'ele')])

envcovs <- apply(envcovs, 2, function(x) scale(x))

#quad-level covs

quad.covs <- read.csv('data/mebba_envcovs_quad.csv')

quad.covs$X11 <- ifelse(is.na(quad.covs$X11), round(mean(quad.covs$X11, na.rm = TRUE)), quad.covs$X11)
quad.covs$X21 <- ifelse(is.na(quad.covs$X21), round(mean(quad.covs$X21, na.rm = TRUE)), quad.covs$X21)
quad.covs$X22 <- ifelse(is.na(quad.covs$X22), round(mean(quad.covs$X22, na.rm = TRUE)), quad.covs$X22)
quad.covs$X23 <- ifelse(is.na(quad.covs$X23), round(mean(quad.covs$X23, na.rm = TRUE)), quad.covs$X23)
quad.covs$X24 <- ifelse(is.na(quad.covs$X24), round(mean(quad.covs$X24, na.rm = TRUE)), quad.covs$X24)
quad.covs$X31 <- ifelse(is.na(quad.covs$X31), round(mean(quad.covs$X31, na.rm = TRUE)), quad.covs$X31)
quad.covs$X41 <- ifelse(is.na(quad.covs$X41), round(mean(quad.covs$X41, na.rm = TRUE)), quad.covs$X41)
quad.covs$X42 <- ifelse(is.na(quad.covs$X42), round(mean(quad.covs$X42, na.rm = TRUE)), quad.covs$X42)
quad.covs$X43 <- ifelse(is.na(quad.covs$X43), round(mean(quad.covs$X43, na.rm = TRUE)), quad.covs$X43)
quad.covs$X52 <- ifelse(is.na(quad.covs$X52), round(mean(quad.covs$X52, na.rm = TRUE)), quad.covs$X52)
quad.covs$X71 <- ifelse(is.na(quad.covs$X71), round(mean(quad.covs$X71, na.rm = TRUE)), quad.covs$X71)
quad.covs$X81 <- ifelse(is.na(quad.covs$X81), round(mean(quad.covs$X81, na.rm = TRUE)), quad.covs$X81)
quad.covs$X82 <- ifelse(is.na(quad.covs$X82), round(mean(quad.covs$X82, na.rm = TRUE)), quad.covs$X82)
quad.covs$X90 <- ifelse(is.na(quad.covs$X90), round(mean(quad.covs$X90, na.rm = TRUE)), quad.covs$X90)
quad.covs$X95 <- ifelse(is.na(quad.covs$X95), round(mean(quad.covs$X95, na.rm = TRUE)), quad.covs$X95)
quad.covs$ch <- ifelse(is.na(quad.covs$ch), mean(quad.covs$ch, na.rm = TRUE), quad.covs$ch)
quad.covs$imp <- ifelse(is.na(quad.covs$imp), mean(quad.covs$imp, na.rm = TRUE), quad.covs$imp)
quad.covs$ele <- ifelse(is.na(quad.covs$ele), mean(quad.covs$ele, na.rm = TRUE), quad.covs$ele)

quad.covs$cult <- quad.covs$X81 + quad.covs$X82
quad.covs$wet <- quad.covs$X90 + quad.covs$X95

quad.covs <- as.matrix(quad.covs[, c('X21', 'X41', 'X42', 'X43', 'X52', 'X71', 'cult', 'wet', 'ch', 'imp', 'ele')])

quad.covs <- apply(quad.covs, 2, function(x) scale(x))


#make the observation covariates
obscovs <- eff[, c('Wind.Speed.Code', 'Sky.Condition.Code', 'Background.Noise.Code', 'ydate', 'ytime')]

#convert NAs into median values

obscovs[which(is.na(obscovs$Wind.Speed.Code) == TRUE), 'Wind.Speed.Code'] <- 0
obscovs[which(is.na(obscovs$Sky.Condition.Code) == TRUE), 'Sky.Condition.Code'] <- 1
obscovs[which(is.na(obscovs$Background.Noise.Code) == TRUE), 'Background.Noise.Code'] <- 0
obscovs$ydate <- scale(obscovs$ydate)
obscovs$ytime <- scale(obscovs$ytime)

#define the distance observation bands for each survey
maxd <- 100

#define the number of distance breaks for each survey
db <- 5

#define the midpoints for each distance break for each survey
mdpts <- c(10, 30, 50, 70, 90)

#define the difference between distance breaks for each survey (delta)
delta <- c(20, 20, 20, 20, 20)

#create an idealized dataset
# f<-(2*mdpts*delta)/(maxd*maxd)
# g<-c(1/4.5, 1/4.5, 1/4.5, 1/4.5, .5/4.5)
# pi.sim<-f * g
# pi.sim/sum(pi.sim)
# simdat<-rmultinom(1000, 1, pi.sim)
# simsum<-rowSums(simdat)
# db1.sim<-c(rep(1, simsum[1]), rep(2, simsum[2]), rep(3, simsum[3]), rep(4, simsum[4]), rep(5, simsum[5]) )
# 

#histograms of each observer

obs.m <- merge(obs, eff, 'siteID', all.x = TRUE, sort = FALSE)

p <- ggplot(obs.m, aes(x = dist))
p <- p + geom_histogram(binwidth = 20) + facet_grid(~Observer)
print(p)

p <- ggplot(obs.m, aes(y = dist, x = Observer))
p <- p + geom_violin() #+ facet_grid(~Observer)
print(p)

#time to detection data

summary(obs$First.Minute.Detected)

time <- recode(obs$First.Minute.Detected, '10' = 9L)

hist(time)

#counts per site a quick count just for testing -- i need to double check vs. effort data later

y <- dcast(obs, siteID ~ 'num', fun.aggregate = length)
y <- merge(sites, y, 'siteID', all.x = TRUE, sort = FALSE)
y$y <- ifelse(is.na(y$num)== TRUE, 0, y$num)
y <- y[order(y$siteID),]

write.csv(eff, 'mebba_effort.csv')

#pull in atlas data to add as a covariate and for integrated modeling

bba <- read.csv('data/quad_spp_atlas_summary.csv')
bba$QuadCode <- factor(gsub('[0-9]+', '', bba$QuadNameCo))
bba$QuadCode <- str_replace_all(bba$QuadCode, ' ', '_')
bba <- left_join(quad.dat, bba, 'QuadCode')

bba.det <- rep(NA, times = nrow(bba))     #make a vector to house the results
bba.det <- ifelse(bba$Black.throated.Green.Warbler >= 1, 1, 0) #summarize the data as detections ignoring NAs for now
bba.det <- data.frame(QuadCode = bba$QuadCode, bba.det = bba.det)

bba.br <- rep(NA, times = nrow(bba))
bba.br <- ifelse(bba$Black.throated.Green.Warbler > 2, 1, 0)   #up the threshold for inclusion, only probable and greater
bba.br <- data.frame(QuadCode = bba$QuadCode, bba.br = bba.br)

eff.bbadet <- left_join(eff[, 1:34], bba.det, by = 'QuadCode')
eff.bbabr <- left_join(eff[, 1:34], bba.br, by = 'QuadCode')

#create be data for each quad in the random effect

quadmap.dat <- quadmap@data
quadmap.dat$QuadName <- str_replace_all(quadmap.dat$QuadCode, ' ', '_')

quad.bba <- left_join(quadmap.dat, bba, by = c('QuadName' = 'QuadCode'))

be <- ifelse(quad.bba$Black.throated.Green.Warbler > 2, 1, 0) 

quadmap.be <- quadmap
quadmap.be@data <- quad.bba

spplot(quadmap.be, 'Black.throated.Green.Warbler')

#creating CAR data

CM2<-as.carCM(adj.sparse, rep(1, length(adj.sparse)), num)
C<- CM2$C
M <- CM2$M

bou<-carBounds(C, adj.sparse, num, M)

mumod<-rep(0,length(num))

R = structure(.Data = c(1, 0, 0, 1), .Dim = c(2, 2))

#compile the data for nimble

dat <- list(maxd=maxd,midpts=mdpts, nbreaks=db, nobs=length(db1), delta=delta, nperiods = 10,
            site.idx = obs$siteID, nsites = nrow(eff), L = length(adj.sparse), wei = rep(1, length(adj.sparse)), adj = adj.sparse, num = num, bou = bou, C= C, M = M,
            distclass=db1, tinterval = time + 1, distcovs = as.matrix(obscovs[,1:3]), avcovs = as.matrix(obscovs[,4:5]), ndistcovs = ncol(obscovs[,1:3]), navcovs = ncol(obscovs[,4:5]),
            quadID = eff$QuadID, nquads = length(num),  envcovs = envcovs, nenvcovs = ncol(envcovs), quadcovs = quad.covs, be = be,
            R = R,
            y = y$y)

hist(dat$y)
hist(rpois(length(dat$y), mean(dat$y)))

nimbrun <- function(seed, ndata){
  
  require(nimble)
  
  ########BUILD A MODEL IN NIMBLE############################################################
  mod <- nimbleCode(
    {
      #priors
      #distance model
      sigma.mu ~ dunif(0,500)
      theta.mu ~ dunif(0, 10)
      
      for(i in 1:ndistcovs){
        sigma.b0[i] ~ dnorm(0, 0.01)
        #theta.b0[i] ~ dnorm(0, 0.01)
      }
      
      #availability model
      beta.a0 ~ dnorm(0, 0.1)
      
      for(i in 1:navcovs){
        pa.b0[i] ~ dnorm(0, 0.01)
      }
      
      #count model
      for(i in 1:nenvcovs){
        lambda.b[i] ~ dnorm(0, 0.01)
      }
      
      beta11 ~ dnorm(0, 0.01)
      beta12 ~ dnorm(0, 0.01)
      beta21 ~ dnorm(0, 0.01)
      beta22 ~ dnorm(0, 0.01)
      
      #mcar setup 
      for(k in 1:2){
      u[1:nquads, k] ~ dcar_proper(muu[1:nquads], C[1:L], adj[1:L], num[1:nquads], M[1:nquads], tauu[k], gammau[k])
      gammau[k] ~ dunif(bou[1], bou[2])
      tauu[k] <- 1
      mu0[k] ~ dnorm(0, 0.1) #this is currently not being used, double-check the setup
      a0[k]~dnorm(0, tau0[k])
      tauv[k]~dgamma(2, 0.5)
      tau0[k]~dgamma(2, 0.5)
      
      for(n in 1:nquads){

      v[n, k] ~ dnorm(0, tauv[k])
        
      }#n
      }#k
      
      Prec[1:2, 1:2] ~ dwish(R[1:2, 1:2], 2)
      Cov[1:2, 1:2] <- inverse(Prec[1:2, 1:2])
      Achol[1:2, 1:2] <- t(chol(Cov[1:2, 1:2]))

      
      for(n in 1:nquads){
        
      phi[n, 1:2] <- Achol[1:2,1:2] %*% u[n, 1:2]
      }#n
      
      sig[1] <- sqrt(Cov[1,1])
      sig[2] <- sqrt(Cov[2,2])
      cor12 <- Cov[1,2]/(sig[1]*sig[2])    # between point count and atlas data
      

      theta.psi ~ dnorm(0, 0.01)
      r ~ dunif(0, 500)
      #r ~ dgamma(0.01, 0.01) 
      
      #confirming that the mvcar effect is zero
      
      u1 <- mean(u[1:nquads, 1])
      u2 <- mean(u[1:nquads, 2])
      
      #distance detection model
      
      for(n in 1:nsites){
        log(sigma[n]) <- log(sigma.mu) + inprod(distcovs[n,], sigma.b0[1:ndistcovs])
        log(theta[n]) <- log(theta.mu) #+ inprod(obscovs[n,], theta.b0[]) 
        
        for(b in 1:nbreaks){
          cloglog(g[n, b]) <- theta[n]*log(sigma[n]) - theta[n]*log(midpts[b])
          f[n, b] <- (2*midpts[b]*delta[b])/(maxd*maxd)
          pi.pd[n, b] <- g[n, b]*f[n, b]
          pi.pd.c[n, b] <- pi.pd[n, b]/pdet[n]
        }#b
        pdet[n] <- sum(pi.pd[n, 1:nbreaks])
        
        #availability model
        
        logit(p.a[n]) <- beta.a0 + inprod(avcovs[n,], pa.b0[1:navcovs])
        
        for(j in 1:nperiods){
          pi.pa[n, j] <- p.a[n]*pow(1 - p.a[n], (j - 1))
          pi.pa.c[n, j] <- pi.pa[n, j]/pavail[n]
        }
        pavail[n] <- sum(pi.pa[n, 1:nperiods])
        
      }#n
      
      
      #linking the distance observations to this detection framework
      for(i in 1:nobs){
        distclass[i] ~ dcat(pi.pd.c[site.idx[i], 1:nbreaks])
        distclassnew[i] ~ dcat(pi.pd.c[site.idx[i], 1:nbreaks])
        
        E.dc[i] <- pow(1 - pow(pi.pd.c[site.idx[i], distclass[i]], 0.5), 2)
        E.dcnew[i] <- pow(1 - pow(pi.pd.c[site.idx[i], distclassnew[i]], 0.5), 2)
        
        tinterval[i] ~ dcat(pi.pa.c[site.idx[i], 1:nperiods])
        tintervalnew[i] ~ dcat(pi.pa.c[site.idx[i], 1:nperiods])
        
        E.ti[i] <- pow(1 - pow(pi.pa.c[site.idx[i], tinterval[i]], 0.5), 2)
        E.tinew[i] <- pow(1 - pow(pi.pa.c[site.idx[i], tintervalnew[i]], 0.5), 2)
        
      }#i
      
      Tobs <- sum(E.dc[1:nobs])
      Tobsnew <- sum(E.dcnew[1:nobs])
      # 
      Ttime <- sum(E.ti[1:nobs])
      Ttimenew <- sum(E.tinew[1:nobs])
      
      #point count abundance estimation
      
      for(n in 1:nsites){
        
        y[n] ~ dbin(pdet[n], navail[n])
        navail[n] ~ dbin(pavail[n], N[n])
        N[n] ~ dpois(lambda.star[n])
        lambda.star[n] <- lambda[n] #* x[n] # * rho[n]
        #rho[n] ~ dgamma(r, r)
        x[n] ~ dbern(psi[n])
        log(lambda[n]) <- a0[1] + phi[quadID[n], 1] #+ beta11*quad.covs[quadID[n], 9] # + beta12*quad.covs[quadID[n], 11] #+ inp+ v[quadID[n], 1] od(envcovs[n,], lambda.b[1:nenvcovs])
        logit(psi[n]) <- theta.psi
        
        # 
        ynew[n] ~ dbin(pdet[n], navail[n])
        # 
        e.y[n] <- pdet[n] * pavail[n] * N[n]
        E.y[n] <- pow(y[n] - e.y[n], 2)/(e.y[n] + 0.5)
        E.ynew[n] <- pow(ynew[n] - e.y[n], 2)/(e.y[n] + 0.5)
      }
      
      fit.y <- sum(E.y[1:nsites])
      fit.ynew <- sum(E.ynew[1:nsites])
      
      #atlas breeding bird spatial modeling
      
      for(m in 1:nquads){
      
      be[m] ~ dbern(phi.be[m])
      logit(phi.be[m]) <- a0[2] + phi[m, 2] #+ beta21*quad.covs[m, 9] # + beta22*quad.covs[m, 11] #+ v[m, 2]
      
      muu[m] <- 0
      
      }
      
      ######## Summary stats
      meanpavail<-mean(pavail[1:nsites]) # mean probability of availability
      meanpdet<-mean(pdet[1:nsites]) # mean probability of perceptibility
      # #bayesp.pd<-step(fit.new.pd-fit.pd) # Bayesian p-value for perceptibility model
      # #bayesp.pa<-step(fit.new.pa-fit.pa) # Bayesian p-value for availability model
      meanN<-mean(N[1:nsites]) # mean site-level abundance
      totN<-sum(N[1:nsites])  # population size of total area surveyed
      meansig<-mean(sigma[1:nsites]) # mean scale parameter across sites
      meanthet<-mean(theta[1:nsites])
      # dens<-meanN/(maxd*maxd*3.14159/10000) # density of birds per ha
      
      # for(b in 1:nbreaks){
      # meanpipdc[b] <- mean(pi.pd.c[1:nsites, b])
      # }
      # 
      # for(j in 1:nperiods){
      # meanpipac[j] <- mean(pi.pa.c[1:nsites, j])
      # }
      # 
    })
  ######GATHER DATA AND RUN THE JAGS MODEL####################################################
  dat <- ndata
  
  params <- (c('sigma.mu', 'sigma.b0', 'theta.mu', 'beta.a0', 'theta.psi', 'v', 'a0', 'phi', 'cor12', 'beta11', 'beta12', 'beta21', 'beta22',
               'meanpavail', 'meanpdet', 'meanN', 'totN', 'meansig', 'meanthet', 'ynew', 'r', 'sig', 'u1', 'u2',
               'Tobs', 'Tobsnew', 'Ttime', 'Ttimenew', 'fit.y', 'fit.ynew', 'gammau', 'tauu', 'mu0', 'lambda.b'))
  
  Nst <- dat$y + 1
  xst <-ifelse(dat$y >= 1, 1, 0)
  Nstu <- log(mean(Nst))
  
  inits <- list(N = Nst + 50,
                navail = Nst + 5,
                x = xst,
                sigma.mu = runif(1, 35, 50),
                theta.mu = runif(1, 3, 5),
                sigma.b0 = runif(3, -.1, .1),
                mu0 = Nstu,
                theta.psi = 0,
                beta.a0 = 0,
                pa.b0 = runif(2, -.1, .1),
                u = structure(.Data= runif(720 * 2, -.1, .1),.Dim=c(720, 2)),
                v = structure(.Data= runif(720 * 2, -.1, .1),.Dim=c(720, 2)),
                a0 = runif(2, -.1, .1),
                tauv = runif(2, .5, 5),
                tau0 = runif(2, .5, 5),
                Prec=structure(.Data=c(runif(1, 1, 3), 0, 0, runif(1, 1, 3)), .Dim=c(2,2)),
                lambda.b = runif(ncol(dat$envcovs), -.1, .1),
                beta11 = runif(1, -.1, .1),
                beta21 = runif(1, -.1, .1),
                beta12 = runif(1, -.1, .1),
                beta22 = runif(1, -.1, .1))
  
  #dput(dat, paste0('mebba_spp_', spp, '.dat'))
  
  ###nimble model fitting
  
  #quick attempt
  
  ni <- 90000
  na <- 100
  nb <- 60000
  nt <- 10
  nc <- 1
  
  
  fm <- nimbleModel(code = mod, constants = dat, inits = inits, name = 'MEBBA_IM', debug = TRUE)
  md <- modelDefClass$new(name = 'MEBBA_IM')
  test <- md$setupModel(code = mod, constants = dat, inits = inits, name = 'MEBBA_IM', debug = TRUE)
  fmc <- compileNimble(fm)
  fm.conf <- configureMCMC(fm, enableWAIC = TRUE)
  # fm.conf <- configureMCMC(fm, enableWAIC = TRUE)
  # fm.conf$calculateWAIC(nburnin = nb)
  fm.conf$addMonitors(params)
  fm.mcmc <- buildMCMC(fm.conf)
  fmc.mcmc <- compileNimble(fm.mcmc)
  #set.seed(seed)
  samples <- runMCMC(fmc.mcmc, niter = ni, nburnin = nb, nchains = nc, thin = nt, 
                     setSeed = seed,
                     summary = FALSE, samplesAsCodaMCMC = TRUE)
  
  
  return(samples)
  
}

start.time <- Sys.time()
this_cluster <- makeCluster(3)

chain_out <- parLapply(cl = this_cluster, X = 1:3, fun = nimbrun, ndata = dat)

stopCluster(this_cluster)

stop.time <- Sys.time()
stop.time - start.time

#combine together

out_mcmc <- as.mcmc.list(list(chain_out[[1]], chain_out[[2]], chain_out[[3]]))

saveRDS(out_mcmc, file= paste0('./mcar_results/', spp, 'pcar_im.rds'))


plot(out_mcmc[, 'totN'])
plot(out_mcmc[, 'a0[1]'])
plot(out_mcmc[, 'a0[2]'])
#plot(out_mcmc[, 'tauu[1]'])
#plot(out_mcmc[, 'tauu[2]'])
plot(out_mcmc[, 'tauv[2]'])
plot(out_mcmc[, 'tau0[1]'])
plot(out_mcmc[, 'tau0[2]'])
plot(out_mcmc[, 'sig[1]'])
plot(out_mcmc[, 'sig[2]'])
plot(out_mcmc[, 'cor12'])
plot(out_mcmc[, 'beta11'])
plot(out_mcmc[, 'beta21'])


# spp <- 'Black-throated Green Warbler'
#load(paste0('./mcar_results/mebbacar_spp_mpcar_zip_', spp, 'mcar_testing_pois_quadcov_nov_bc_prob_testing2.Rdata'))

#now let's make some quad-level predictions with the model

fit.sum <- summary(out_mcmc)
View(fit.sum$statistics)
View(fit.sum$quantiles)

#look at the fit model

fit <- out_mcmc

fit.rhat <- gelman.diag(fit, multivariate = FALSE)
View(fit.rhat[[1]])
summary(fit.rhat[[1]])

#dispersion
summary(fit[, 'fit.y'])$statistics[1]/summary(fit[, 'fit.ynew'])$statistics[1]

mean(unlist(fit[,'Tobsnew']) > unlist(fit[,'Tobs']))
plot(unlist(fit[,'Tobsnew']), unlist(fit[,'Tobs']))
abline(0,1,lwd=2)

mean(unlist(fit[,'Ttimenew']) > unlist(fit[,'Ttime']))
plot(unlist(fit[,'Ttimenew']), unlist(fit[,'Ttime']))
abline(0,1,lwd=2)

mean(unlist(fit[,'fit.ynew']) > unlist(fit[,'fit.y']))
plot(unlist(fit[,'fit.ynew']), unlist(fit[,'fit.y']))
abline(0,1,lwd=2)


#k-fold cross validation

## activate when need this tool
# start.time <- Sys.time()
# fm <- nimbleModel(code = mod, constants = dat, inits = inits)
# fm.conf <- configureMCMC(fm)
# fit.cv <- runCrossValidate(fm.conf, k = 5, nCores = 1)
# end.time <- Sys.time()
# end.time - start.time

#bayesplot model checking

#check ynews

#ynew <- matrix(cbind(fit[[1]][, 755:4430], fit[[2]][, 755:4430], fit[[3]][, 755:4430]), ncol = length(dat$y) * 3)
# ynew <- rbind(fit[[1]][, 755:4430], fit[[2]][, 755:4430], fit[[3]][, 755:4430])
# 
# bayesplot::pp_check(dat$y, ynew[sample(1:nrow(ynew), 100, replace = FALSE), ], ppc_dens_overlay)
# 
# #looking at individual
# 
# plot(fit[, 'meanpdet'])
# plot(fit[, 'meanpavail'])
# plot(fit[, 'dens'])
# plot(fit[, 'totN'])
# plot(fit[, 'sigma.mu'])
# plot(fit[, 'theta.mu'])
# plot(fit[, 'beta.a0'])
# plot(fit[, 'lambda.mu'])
# plot(fit[, 'lambda.sigma'])
# plot(fit[, 'theta.psi'])
# plot(fit[, 'r'])
# plot(fit[, 'mu0'])
# plot(fit[, 'pa.b0[1]'])
# plot(fit[, 'pa.b0[2]'])
# plot(fit[, 'sigma.b0[1]'])
# plot(fit[, 'sigma.b0[2]'])
# plot(fit[, 'sigma.b0[3]'])
# plot(fit[, 'u[3]'])
# plot(fit[, 'gammau'])
# plot(fit[, 'lambda.b[1]'])
# plot(fit[, 'lambda.b[2]'])
# plot(fit[, 'lambda.b[3]'])
# plot(fit[, 'lambda.b[4]'])
# plot(fit[, 'lambda.b[5]'])
# plot(fit[, 'lambda.b[6]'])
# plot(fit[, 'lambda.b[7]'])
# plot(fit[, 'lambda.b[8]'])
# plot(fit[, 'lambda.b[9]'])
# 
# plot(fit[, c('lambda.beta[1]', 'lambda.beta[2]', 'lambda.beta[3]', 'lambda.beta[4]')])
# 
# summary(fit[, c('meanpipdc[1]', 'meanpipdc[2]', 'meanpipdc[3]', 'meanpipdc[4]', 'meanpipdc[5]')])
# 
# summary(fit[, c('meanpipac[1]', 'meanpipac[2]', 'meanpipac[3]', 'meanpipac[4]', 'meanpipac[5]', 'meanpipac[6]', 'meanpipac[7]',
#                 'meanpipac[8]', 'meanpipac[9]', 'meanpipac[10]')])
# 
# plot(fit[, 'pa.b0[5]'])

#compare original y's to model developed y's

# ycomp <- data.frame(yobs = y$y, ynew = fit.sum$quantiles[226:2224, 3], ynew2 = fit.sum$statistics[226:2224, 1])
# ycomp <- melt(ycomp)
# 
# p <- ggplot(ycomp, aes(x = value))
# p <- p + geom_histogram() + facet_wrap(~variable)
# p <- p + theme_bw()
# 
# print(p)

#connect model predictions to atlas blocks

#spp <- 'Eastern Wood-Pewee'
#load('mebba_spp_Eastern Wood-Pewee')
#fit.sum <- summary(fit)

#pull in the prediction map

pgrid <- readOGR(getwd(), 'pgrid_2020_corrected')
pcovs <- read.csv('mebba_envcovs_pgrid.csv')

pgrid@data <- cbind(pgrid@data, pcovs)

pgrid@data$QuadName <- factor(str_sub(pgrid@data$CLOBlcN, end = -4))
pgrid@data$QuadCode <- factor(gsub('[0-9]+', '', pgrid@data$QuadNmC))

pgrid_qc <- left_join(pgrid@data, quadmap@data, by = c('QuadName' = 'QuadCode'))

pgrid@data$quadID <- pgrid_qc$quadID
pgrid@data$wet <- pgrid@data$X90 + pgrid@data$X95
pgrid@data$cult <- pgrid@data$X81 + pgrid@data$X82

#where are the covariate NAs?

#spplot(pgrid, zcol = 'ele')

#key that connects blocks to JAGS model

#key <- read.csv('site_key.csv')
#eff <- read.csv('mebba_effort.csv')

#quads <- unique(eff[, c('Atlas.Quad.Name', 'QuadID')])

#original data set for bringing in original count data

#dat <- dget('mebba_test.dat')

#make density predictinos (# birds/ha) for each large atlas block

quad.ps <- data.frame(p = fit.sum$statistics[37:(37 + 719), 1], quadID = 1:length(fit.sum$statistics[37:(37 + 719), 1]))

#ps <- exp(fit.sum$statistics[35:(35 + 719), 1])/(pi * (.1^2))/100

#ps <- data.frame(ps, QuadID = 1:length(ps))

#link back to the siteID key

#psm <- merge(quads, ps, by = 'QuadID', all.x = TRUE)
p.quadcar <- left_join(pgrid@data, quad.ps, by = 'quadID')

cov.params <- fit.sum$statistics[9:(9 + ncol(dat$envcovs) - 1), 1]
pgrid.covs <- as.matrix(pgrid@data[, c('X21', 'X41', 'X42', 'X43', 'X52', 'X71', 'cult', 'wet', 'ch', 'imp', 'ele')])

pgrid.covs[, 1:7] <- pgrid.covs[, 1:7] * (pi * (.1)^2)
pgrid.covs[, 1] <- (pgrid.covs[, 1] - attr(scale(eff.covs$X21), 'scaled:center'))/attr(scale(eff.covs$X21), 'scaled:scale')
pgrid.covs[, 2] <- (pgrid.covs[, 2] - attr(scale(eff.covs$X41), 'scaled:center'))/attr(scale(eff.covs$X41), 'scaled:scale')
pgrid.covs[, 3] <- (pgrid.covs[, 3] - attr(scale(eff.covs$X42), 'scaled:center'))/attr(scale(eff.covs$X42), 'scaled:scale')
pgrid.covs[, 4] <- (pgrid.covs[, 4] - attr(scale(eff.covs$X43), 'scaled:center'))/attr(scale(eff.covs$X43), 'scaled:scale')
pgrid.covs[, 5] <- (pgrid.covs[, 5] - attr(scale(eff.covs$X52), 'scaled:center'))/attr(scale(eff.covs$X52), 'scaled:scale')
pgrid.covs[, 6] <- (pgrid.covs[, 6] - attr(scale(eff.covs$X71), 'scaled:center'))/attr(scale(eff.covs$X71), 'scaled:scale')
pgrid.covs[, 7] <- (pgrid.covs[, 7] - attr(scale(eff.covs$cult), 'scaled:center'))/attr(scale(eff.covs$cult), 'scaled:scale')
pgrid.covs[, 8] <- (pgrid.covs[, 8] - attr(scale(eff.covs$wet), 'scaled:center'))/attr(scale(eff.covs$wet), 'scaled:scale')
pgrid.covs[, 9] <- (pgrid.covs[, 9] - attr(scale(eff.covs$ch), 'scaled:center'))/attr(scale(eff.covs$ch), 'scaled:scale')
pgrid.covs[, 10] <- (pgrid.covs[, 10] - attr(scale(eff.covs$imp), 'scaled:center'))/attr(scale(eff.covs$imp), 'scaled:scale')
pgrid.covs[, 11] <- (pgrid.covs[, 11] - attr(scale(eff.covs$ele), 'scaled:center'))/attr(scale(eff.covs$ele), 'scaled:scale')
cov.ps <- t(apply(pgrid.covs, 1, function(x) cov.params*x))

p <- exp(p.quadcar$p + rowSums(cov.ps)) / (pi * 0.1^2) / 100

#now link this to the shapefile

pgrid@data$p <- p


## old prediction attempt using the atlas block shapefile
# amap.dat <- amap@data
# amap.dat$x <- cents[,2]
# amap.dat$y <- cents[,1]
# amap.dat <- dplyr::left_join(amap.dat, psm, by = c('QuadName' = 'QuadCode'))

#amap@data <- amap.dat

# amap <- SpatialPointsDataFrame(data.frame(amap@data$x, amap@data$y), amap@data,
#                                     proj4string = CRS(proj4string(amap)))
# 
# dat.merge@data$X1 <- coordinates(dat.merge)[, 1]
# dat.merge@data$Y1 <- coordinates(dat.merge)[, 2]

# spplot(pgrid, zcol = 'p', main=list(label= paste0(spp," density (N/ha) by Atlas Quad",cex=2)))
# png(file = paste0(spp,'dens_map.png'))
# spplot(amap, zcol = 'ps', main=list(label= paste0(spp," density (N/ha) by Atlas Quad",cex=2)))
# dev.off()
# 
#writeOGR(pgrid, dsn = getwd(), layer = paste0(spp,"_car_dens_1km2_2020"), driver = 'ESRI Shapefile', overwrite = TRUE)

#as a Raster
p.ras <- raster(extent(pgrid))
res(p.ras) <- 1000
proj4string(p.ras) <- proj4string(pgrid)
pred.ras <- rasterize(pgrid, p.ras, 'p')

plot(pred.ras)

writeRaster(pred.ras, paste0('raster_results_2020/', spp,"_car_dens_1km2_2020_test"), prj = TRUE, format = 'ascii',  overwrite = TRUE)

#sum(values(pred.ras), na.rm = TRUE) * 100
sum(pred.ras@data@values, na.rm = TRUE) * 100

#SE of predictions

quad.sds <- data.frame(sd = fit.sum$statistics[37:(37 + 719), 2], quadID = 1:length(fit.sum$statistics[37:(37 + 719), 2]))

#link back to the siteID key

#psm <- merge(quads, ps, by = 'QuadID', all.x = TRUE)
sd.quadcar <- left_join(pgrid@data, quad.sds, by = 'quadID')

cov.params <- fit.sum$statistics[9:(9 + ncol(dat$envcovs) - 1), 1]
pgrid.covs <- as.matrix(pgrid@data[, c('X21', 'X41', 'X42', 'X43', 'X52', 'X71', 'cult', 'wet', 'ch', 'imp', 'ele')])

pgrid.covs[, 1:7] <- pgrid.covs[, 1:7] * (pi * (.1)^2)
pgrid.covs[, 1] <- (pgrid.covs[, 1] - attr(scale(eff.covs$X21), 'scaled:center'))/attr(scale(eff.covs$X21), 'scaled:scale')
pgrid.covs[, 2] <- (pgrid.covs[, 2] - attr(scale(eff.covs$X41), 'scaled:center'))/attr(scale(eff.covs$X41), 'scaled:scale')
pgrid.covs[, 3] <- (pgrid.covs[, 3] - attr(scale(eff.covs$X42), 'scaled:center'))/attr(scale(eff.covs$X42), 'scaled:scale')
pgrid.covs[, 4] <- (pgrid.covs[, 4] - attr(scale(eff.covs$X43), 'scaled:center'))/attr(scale(eff.covs$X43), 'scaled:scale')
pgrid.covs[, 5] <- (pgrid.covs[, 5] - attr(scale(eff.covs$X52), 'scaled:center'))/attr(scale(eff.covs$X52), 'scaled:scale')
pgrid.covs[, 6] <- (pgrid.covs[, 6] - attr(scale(eff.covs$X71), 'scaled:center'))/attr(scale(eff.covs$X71), 'scaled:scale')
pgrid.covs[, 7] <- (pgrid.covs[, 7] - attr(scale(eff.covs$cult), 'scaled:center'))/attr(scale(eff.covs$cult), 'scaled:scale')
pgrid.covs[, 8] <- (pgrid.covs[, 8] - attr(scale(eff.covs$wet), 'scaled:center'))/attr(scale(eff.covs$wet), 'scaled:scale')
pgrid.covs[, 9] <- (pgrid.covs[, 9] - attr(scale(eff.covs$ch), 'scaled:center'))/attr(scale(eff.covs$ch), 'scaled:scale')
pgrid.covs[, 10] <- (pgrid.covs[, 10] - attr(scale(eff.covs$imp), 'scaled:center'))/attr(scale(eff.covs$imp), 'scaled:scale')
pgrid.covs[, 11] <- (pgrid.covs[, 11] - attr(scale(eff.covs$ele), 'scaled:center'))/attr(scale(eff.covs$ele), 'scaled:scale')
cov.sds <- t(apply(pgrid.covs, 1, function(x) abs(cov.params)*abs(x)))

sds <- sd.quadcar$sd + rowSums(cov.sds) #/ (pi * 0.1^2) / 100

sds <- (abs(p) * sds)/ (pi * 0.1^2) / 100

#now link this to the shapefile

pgrid@data$sd <- sds

pgrid@data$cv <- pgrid@data$sd/pgrid@data$p

cv.ras <- raster(extent(pgrid))
res(cv.ras) <- 1000
proj4string(cv.ras) <- proj4string(pgrid)
cv.ras <- rasterize(pgrid, cv.ras, 'cv')

plot(cv.ras)

writeRaster(cv.ras, paste0('raster_results_2020/', spp,"_car_dens_1km2_cv_test_2020"), prj = TRUE, format = 'ascii', overwrite = TRUE)


#####END#####

#TO-DO

#1 look at tau notation for dnorm distributions like a0 or v, we should  convert them to precision estimates?

