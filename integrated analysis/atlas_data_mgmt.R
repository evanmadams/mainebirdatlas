require(reshape2)
require(plyr)
require(dplyr)
require(ggplot2)
require(rgdal)
require(rgeos)
require(sp)

#summarize atlas data at the quad scale
dat <- readRDS('data/ebd_amap.rds')

#pull in first atlas data as a grid file

# bba1 <- readOGR('P:/MEBBA/Data Visualization/First Atlas Blocks', layer = 'MBBA1')
# 
# plot(bba1)
# 
# #pull in 2nd atlas grid
# 
# amap <- readOGR(dsn = 'data/blockupdate2019', layer = 'BBA_Blocks_20190418_OldManIsland')

dat$bcNum <- as.numeric(as.factor(dat$breeding_category))

bc.sum <- reshape2::dcast(dat, QuadNameCo + QuadCode ~ common_name, value.var = 'bcNum', fun.aggregate = max)

bc.sum[bc.sum < 1] <- 0
bc.sum[is.na(bc.sum)] <- 0

write.csv(bc.sum, 'data/quad_bc_sum.csv')




####SUMMARY STATS#############
#pull in the atlas data

at <- read.csv("data/Evan's_Atlas_eBird_Data.csv")
at <- at[which(at$Flag.Record != 'Exclude'), ]

at.comp <- read.csv('data/comp_blocks_0120.csv')

#make species x block summary table which provides the highest breeding code

atsum <- dcast(at, CLOBlockNa + NEW_COMMON_NAME + Updated.Final.Breeding.Category ~ 'Count', fun.aggregate = length)

at.mat <- dcast(at, CLOBlockNa ~  NEW_COMMON_NAME + Updated.Final.Breeding.Category, fun.aggregate = length)

#highest breeding code for each species/block combo

blocks <- unique(at$CLOBlockNa)
spp <- unique(at$NEW_COMMON_NAME)

#create a table of highest breeding codes with block completion data too

ab <- data.frame(BlockCode = rep(blocks, each = length(spp)), Species = spp)

ab$HighestBreedingCode <- NA
ab$BlockCompletion <- NA

for(i in 1:nrow(ab)){
  
  print(i)
  
  tmp <- atsum[which(atsum$CLOBlockNa == ab$BlockCode[i] & atsum$NEW_COMMON_NAME == ab$Species[i]), ]
  
  if(nrow(tmp) > 0){
    ab$HighestBreedingCode[i] <- paste0(tmp$Updated.Final.Breeding.Category[nrow(tmp)])
  }
  
  tmp <- at.comp[which(at.comp$Atlas.Block.Name == as.character(ab$BlockCode[i])), ]
  
  if(nrow(tmp) > 0){
    ab$BlockCompletion[i] <- paste0(tmp$Block.Completion)
  }
  
}


write.csv(ab, 'atlas_breeding_code_sum.csv')

ab <- read.csv('atlas_breeding_code_sum.csv')

ab.long <- dcast(ab, BlockCode + BlockCompletion ~ Species, value.var = 'HighestBreedingCode')

#pull in atlas block map

amap <- readOGR(dsn = 'data/newblocksshapefile', layer = 'BBA_Blocks_20180216')

amap.dat <- amap@data

amap.dat$BlockCode <- amap.dat$CLOBlockNa

ab.join <- join(amap.dat, ab.long, by = 'BlockCode')

amap@data <- ab.join

writeOGR(amap, dsn = getwd(), layer = 'atlas_block_breeding_code_join', driver = 'ESRI Shapefile', overwrite = TRUE)

#effort linked point count data point count data

pcdat.eff <- read.csv("data/pc_effort_2018_2019.csv")
pcdat.eff$Date <- as.POSIXct(pcdat.eff$Bird.Survey.Date..MM.DD.YYYY, format = '%m/%d/%Y')

pcdat.obs <- read.csv("data/pc_obs_2018_2019.csv")

eff <- pcdat.eff
obs <- pcdat.obs

eff$Observer <- recode(eff$Observer, 'MATTHEW DICKINSON' = 'Matthew Dickinson', 'C WEST' = 'Christopher West')

#create a site index and covariates for effort data

eff$PC_ID <- as.factor(paste0(eff$Atlas.Block.Name, eff$Point.Count..))

#create a site index for observations and link to effort index

obs.eff <- read.csv('effort_linked_obs.csv')

pc.eff <- data.frame(BlockCode = rep(sites$sites, each = length(spp)), Species = spp)
pc.eff$key <- paste0(pc.eff$BlockCode, pc.eff$Species)

y <- dcast(obs.eff, Atlas.Block.Name.x + Species ~ 'num', fun.aggregate = length)
y$key <- paste0(y$Atlas.Block.Name.x, y$Species)
y <- merge(pc.eff, y, by = 'key', all.x = TRUE, sort = FALSE)
y$y <- ifelse(is.na(y$num)== TRUE, 0, y$num)
y$pa <- ifelse(y$y > 0, 1, 0)


#combine the two data sets together

ab$key <- paste0(ab$BlockCode, ab$Species)

pc.ab <- merge(ab, y, by = 'key', all.x = TRUE)

summary(pc.ab)

pc.ab$HighestBreedingCode <- as.factor(pc.ab$HighestBreedingCode)
pc.ab$BlockCompletion <- as.factor(pc.ab$BlockCompletion)

pc.ab <- pc.ab[which(is.na(pc.ab$pa) == FALSE), ]

pc.ab$C2_pa <- ifelse((pc.ab$HighestBreedingCode == 'C2'| pc.ab$HighestBreedingCode == 'C3'| pc.ab$HighestBreedingCode == 'C4') & pc.ab$pa == 1, 1, 0)
pc.ab$C3_pa <- ifelse((pc.ab$HighestBreedingCode == 'C3'| pc.ab$HighestBreedingCode == 'C4') & pc.ab$pa == 1, 1, 0)
pc.ab$C4_pa <- ifelse(pc.ab$HighestBreedingCode == 'C4' & pc.ab$pa == 1, 1, 0)

mean(pc.ab$C2_pa, na.rm = TRUE)
mean(pc.ab$C3_pa, na.rm = TRUE)
mean(pc.ab$C4_pa, na.rm = TRUE)

mean(pc.ab$C2_pa[which(pc.ab$BlockCompletion == 'Complete')], na.rm = TRUE)
mean(pc.ab$C3_pa[which(pc.ab$BlockCompletion == 'Complete')], na.rm = TRUE)
mean(pc.ab$C4_pa[which(pc.ab$BlockCompletion == 'Complete')], na.rm = TRUE)


#bootstrap pc data per block to estimate spatial variation in dens or pa species by species
#look for correlations with this measurement and the relationship between atlas and pc data


