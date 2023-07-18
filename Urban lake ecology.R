# R script for manuscript
# Title:    Macrophyte communities affect the density and community composition of macroinvertebrates in artificial lakes, France.
# Authors:  Chaozhong Tan, Sabine Greulich, Valentin Medina, Xue Zheng, Pao Canu, Alan Fritsch and Karl Matthias Wantzen 

# set working direction
setwd('D:\\PhD\\Lake ecological study\\Statistical analysis') # Set working direction where you store your data in .csv format
library(vegan)
library(ecodist) # MRM
library(corrplot) # Correlation plot
library(gamlss)
library(rdacca.hp)
# Create useful functions ##### 
{
  # Calculate the mean and SE of a matrix
  Mean.SE <- function(data){
    # Mean
    nrow <- nrow(data)
    mean <- colSums(data)/nrow
    
    # SD
    ncol <- ncol(data)
    i=0;result.F <- NULL
    repeat{
      i=i+1
      result <- sd(data[,i])
      result.F <- c(result.F,result)
      if (i == ncol) break
    }
    
    # SE
    SE <- result.F/(nrow^(1/2))
    
    result <- cbind(mean,SE)
    return(result)
  }
  
  Mean.SE.grp <- function(matrix,grp){
    title <- Element.grp(grp)
    n <- length(title)
    i=0;result.F <- NULL
    repeat{
      i=i+1
      data <- matrix[grp==title[i],]
      result <- Mean.SE(data)
      result.F <- cbind(result.F,result)
      if(i==n)break
    }
    return(result.F)
    return('col names are grp variables')
  }

  # Conduct wilcox.test of vector variable among groups
  Wilcox.test.grp <- function(data,grp){
    
    # See how many levels of the grp
    level <- Element.grp(grp) 
    
    n <- length(level)
    # built two matrix to store the statistics of 'W' and 'P value'
    
    i=0;result.F.F <- NULL
    repeat{
      i=i+1
      
      # use Data[i] to compare with all the Data behind
      j=i;result.F <- NULL
      repeat{
        j=j+1
        
        test <- wilcox.test(data[grp==level[i]],data[grp==level[j]])
        t.statistic <- test$statistic
        p.value <- test$p.value
        result <- c(t.statistic,p.value)
        result.F <- c(result.F,result)
        if(j==n)break
      }
      result.F <- c(rep(0,(i-1)*2),result.F)
      result.F.F <- rbind(result.F.F,result.F)
      if(i==n-1)break
    }
    
    # give row name to the result
    rownames(result.F.F) <- level[-n]
    
    # give 'w' and 'p' to the result
    result.F.F <- rbind(
      (rep(c('w','p'),n-1))
      ,result.F.F)
    
    # give col name to the reesult
    i=0;colname <- NULL
    repeat{
      i=i+1
      result <- c(level[i+1],'')
      colname <- c(colname,result)
      if(i+1 == n) break
    }
    
    colnames(result.F.F) <- colname
    
    return(result.F.F)
  }
  
  # Define how many elements in a vector group
  Element.grp <- function(data){
    
    element.F <- NULL
    repeat{
      
      element <- data[1]
      
      # create a index from 1 to n
      n <- length(data);vector <- 1:n
      
      # get the index of the elements that are same as the first one
      index <- vector[data==element]
      
      # remove the same elements, and formed a new data list
      data <- data[-index]
      
      # record the elements
      element.F <- c(element.F,element)
      
      if (length(data)==0)break
    }
    return(element.F)
  }
  
  # Replace the strings within a vector group
  Replace <- function(original,replace){
    old <- Element.grp(original)
    new <- NULL
    # replace the grp var to color 
    i=0
    repeat{
      i=i+1
      
      # For each element in original
      # Check which they are
      # And then replace then accordingly
      j=0
      repeat{
        j=j+1
        if(original[i] == old[j]) new[i] <- replace[j]
        if(j==length(old))break
      }
      
      if(i == length(original) ) break
    }
    
    # record which replaced which
    i=0; alter.F <- NULL
    repeat{
      i=i+1
      alteration <- c(replace[i],'-->',old[i])
      alter.F <- rbind(alter.F,alteration)
      if(i==length(replace))break
    }
    alter.F
    
    result <- list(alter.F,new)
    return(result)
  }
  # To take the number n of letters from the left of the string spe.name
  Left <- function(spe.name,n){
    j=0;initial <- NULL
    repeat{
      j=j+1
      name <- spe.name[j]
      new.name <- strsplit(name,"")[[1]]
      new.name <- paste(new.name[1:n],collapse='')
      initial <- c(initial,new.name)
      if(j==length(spe.name))break
    }
    return(initial) # initial of the species
  }
  
  # Calculate abundance-based Jaccard dissimilarity
  # The Jaccard dissimilarity was calculated by Bray-curtis distance:
  # Jac.dis = 2B/(b+1)
  Jac.dis <- function(data){
    BC.dis <- vegdist(data,method = 'bray')
    Jaccard.distance = (2*BC.dis)/(BC.dis+1)
    return(Jaccard.distance)
  }
}

#1. Input data of environmental variables #####
env <- read.csv('Environmental data.csv')
env.title <- c('Site', 'Lake', 'Site.number', 'type', 'water transparency', 'macrophyte breath'
               ,rep('substrate',4),rep('water quality',7)
               ,rep('canopy',1))

# Head of the environmental variables (Ommited some variables because lack of space)
# Note: NA represent no values. Because, the 'substrate type' is only measured for both eulittoral and sublittoral zone but not for water surface.

#         Site Lake Site.number type Water.transparency Macrophyte.breadth Silt.mud Wooden.debris.leaf.litters Macrophytes Gravel DO... EC.us.cm. 
# 1 L1.S4.June   L1          S4  eu.                 53                4.0        0                         10          90      0  99.5       456
# 2 L1.S4.June   L1          S4 sub.                 NA                 NA        0                         10          90      0    NA        NA
# 3 L1.S1.June   L1          S1  eu.                 47                5.0       10                         20          65      5 106.9       328
# 4 L1.S1.June   L1          S1 sub.                 NA                 NA       20                         10          20     70    NA        NA
# 5 L1.S2.June   L1          S2  eu.                 60                1.5       25                         50          25      0 106.1       304
# 6 L1.S2.June   L1          S2 sub.                 NA                 NA       95                          5           0      0    NA        NA


#1.2 Environmental variables #####

# First 4 columns are site elements
env.substrate.title <- c(rep('yes',4),rep('no',2),rep('yes',4),rep('no',8))
env.other.env.title <- c(rep('yes',4),rep('yes',2),rep('no',4),rep('yes',8))

substrate <- env[,env.substrate.title == 'yes']
other.env <- env[,env.other.env.title == 'yes'] # other env that are not substrate

judge <- is.na(rowSums(other.env[,-c(1:5)]))
other.env.rm.na <- other.env[judge==FALSE,] # remove the empty value

# Mean and SE of variables of three lakes
title <- c('L1','L2','L3') # Three lakes in the study: Lake 1: Lake Bretonniere; Lake 2: Lake Bergeonnerie; lake 3: Lake Peupleraies.
i=0;result.F <- NULL
repeat{
  i=i+1
  data.substrate <- substrate[substrate$Lake == title[i],]
  data.other.env <- other.env.rm.na[other.env.rm.na$Lake == title[i],]
  
  data.substrate.onlydata <- data.substrate[,-c(1:4)]
  data.other.env.onlydata <- data.other.env[,-c(1:4)]
  
  # Substrate
    # substrate eulittoral zone
      data.eu <- data.substrate.onlydata[data.substrate$type=='eu.',]
      mean.se.substrate.eu<- Mean.SE(data.eu)
    
   # substrate sublittoral zone
      data.sub <- data.substrate.onlydata[data.substrate$type=='sub.',]
      mean.se.substrate.sub <- Mean.SE(data.sub)

  # Other env  
      mean.se.otherenv<- Mean.SE(data.other.env.onlydata)
  
  # Combine the data from eulittoral and sublittoral zone
  result <- rbind(mean.se.otherenv
                  ,mean.se.substrate.eu
                  ,mean.se.substrate.sub)
  
  result.F<-cbind(result.F,result)
  
  if(i==length(title))break
}
result.F <- round(result.F,2) # Take two decimal places
result.F # show the results

#1.3 Pair-wised comparison of environmental variables among three lakes #####

# 1.3.1 Variables of eulittoral zone among lakes (water quality, macrophyte and substrate type)
env.eu <- env[env$type == 'eu.',]

i=4 # the first 4 column are site labels, so the comparison starts from the 5th column
result.F <- NULL
repeat{
  i=i+1
  title.env <- colnames(env.eu)[i]
  data = env.eu[,i]
  
  result <- Wilcox.test.grp(data,env.eu$Lake)
  result <- cbind(title.env,result)
  result.F <- rbind(result.F,result)
  
  if(i== ncol(env.eu) )break
}
"The Mann-Whitney U test of environmental variables (water quality, macrophyte and substrate type) among lakes:"
result.F
#'title.env' indicates which variables used for comparison. The L1, L2, L3 indicate which two of lakes are used for comparison


# 1.3.2 Variables in the sublittoral zone among lakes (just substrate)
env.sub <- env[env$type == 'sub.',]
env.sub <- env.sub[,c(1:4,7:10)]

i=4;result.F <- NULL
repeat{
  i=i+1
  title.env <- colnames(env.sub)[i]
  data = env.sub[,i]
  result <- Wilcox.test.grp(data,env.sub$Lake)
  result <- cbind(title.env,result)
  result.F <- rbind(result.F,result)
  
  if(i== ncol(env.middle) )break
}
" The Mann-Whitney U test of substrate types among lakes:"
result.F
#'title.env' indicates which variables used for comparison. The L1, L2, L3 indicate which two of lakes are used for comparison

#2. MACROPHYTE SECTION#####
macrophyte <- read.csv('Macrophyte cover & weight.csv')
macrophyte.title <- c('site','lake','site.number','type'
                      ,'weight','weight'
                      ,rep('species',15))

# Head of the macrophyte matrix:

#         Site Lake Site.number   Type Wet_weight Dry_weight sp1 sp2 sp3 sp4 sp5 sp6 sp7 sp8 sp9 sp10 sp11 sp12 sp13 sp14 sp15
# 1 L1.S4.June   L1          S4    eu.        0.0        0.0   0   0   1   5   0   0   1  60   0    0    0    0   30    0    0
# 2 L1.S4.June   L1          S4 w.surf        0.0        0.0   0   0   1   0   0   0   1  60   0    0    1    0   30    5    0
# 3 L1.S1.June   L1          S1    eu.        0.0        0.0   0   0   0   0   0   0   2  60   0    0    0    0   20   10    0
# 4 L1.S1.June   L1          S1 w.surf       25.8        3.3   0   0   2   0   0   0   0   0   0    0    8    0    0    0    0
# 5 L1.S2.June   L1          S2    eu.        0.0        0.0   0   0   0  30   0   0   0  20   0    0    1    0    5    1    1
# 6 L1.S2.June   L1          S2 w.surf        0.0        0.0   0   0   1   0   0   0   1   0   0    0    5    0   10    0    0



macrophyte.data <- macrophyte[,macrophyte.title=='species']

# 2.1 SR at lake level #####
macrophyte.col <- colSums(macrophyte.data)
macrophyte.col[macrophyte.col>0]=1 # if value >0, then assign one
SR <- sum(macrophyte.col);SR
"Total species richness in three lakes:"
SR

# 2.2 SR at eulittoral zone and water surface #####

# Separate the data into different lakes
title <- Element.grp(macrophyte$Lake)
n <- length(title) # n = how many lakes
i=0;result.F <- NULL
repeat{
  i=i+1;data.title <- title[i]

# Calculate the SR at lake level
  # separate different lakes
  data = macrophyte.data[macrophyte$Lake == title[i],]
  data.col <- colSums(data);data.col[data.col>0]=1
  data.SR <- sum(data.col)
  
  # separate eulittoral and sublittoral zone in a lake
  title.eu.sublittoral <- (macrophyte$Type)[macrophyte$Lake == title[i]]
  data.eu <- data[title.eu.sublittoral=='eu.',]; data.w.surf <- data[title.eu.sublittoral=='w.surf',]
  
  data.col <- colSums(data.eu); data.col[data.col>0]=1; SR.eu <- sum(data.col)
  data.col <- colSums(data.w.surf);data.col[data.col>0]=1; SR.w.surf <- sum(data.col)
  
  result<- c(data.title,data.SR,SR.eu,SR.w.surf)
  result.F <- rbind(result.F,result)
  
  if(i==n)break
}
colnames(result.F) <- c('Lake','Total', 'Eulittoral','Sublittoral')
rownames(result.F) <- NULL
result.F


# 2.3 Mean SE of weight, cover and SR at site level #####
  # weight
wet.dry_weight <- macrophyte[,macrophyte.title=='weight']
title.weight <- macrophyte$Lake[macrophyte$Wet_weight>0] # remove the empty row for the title factor
wet.dry_weight <- wet.dry_weight[macrophyte$Wet_weight>0,] # remove the empty row for the data
  
  # cover and SR of each site
cover <- rowSums(macrophyte.data)
data <- macrophyte.data;data[data>0]=1
SR <- rowSums(data)

cover.SR <- cbind(cover,SR) # function Mean.SE only applies to matrix
'Vegetation Cover and SR at site level of three lakes:'
Mean.SE(cover.SR)

# cover and SR among littoral zones
title <- macrophyte$Lake;
title.level <- Element.grp(title) # L1, L2 and L3

i=0;result.cover.SR.F<-NULL;result.weight.F <- NULL # first define the variables that are about to use
repeat{
  i=i+1
  # separate SR and weight by lakes
  data.cover.SR <- cover.SR[title == title.level[i],]
  weight <- wet.dry_weight[title.weight == title.level[i],]
  title.eu.w.surf <- macrophyte$Type[title == title.level[i]] # titles for eulittoral zone and water surface
  
  # in each lake, the mean and SE of cover and biomass
  title.level[i]
  mean.se.lake <- Mean.SE(data.cover.SR)
  mea.se.weight <- Mean.SE(weight)
  # separate the eulittoral and sublittoral zone
  data.eu <- data.cover.SR[title.up.mid=='eu.',]; data.w.surf <- data.cover.SR[title.up.mid=='w.surf',]
  
  # the cover and SR for eulittoral and sublittoral zone
  mean.se.eu <- Mean.SE(data.eu)
  mean.se.w.surf <- Mean.SE(data.w.surf)
  
  # combine the results at lake level, eulittoral zone and water surface
  result.cover.SR <- cbind(mean.se.lake,mean.se.eu,mean.se.w.surf)
  
  result.weight.F <- rbind(result.weight.F,mea.se.weight)
  result.cover.SR.F <- rbind(result.cover.SR.F,result.cover.SR)
  
  if(i==3)break
}
result.cover.SR.F <- round(result.cover.SR.F,2) # Take two decimal places
result.cover.SR.F <- rbind(c('Lake level','','Eulittoral zone','','Water surface','')
                  ,result.cover.SR.F) # Assign names for each column
result.cover.SR.F
result.weight.F

#  2.4 Comparison of macrophyte cover among lakes #####
title <- c('eu.','w.surf')
i=0;result.F <- NULL
repeat{
  i=i+1
  data.grp <- macrophyte$Lake[macrophyte$Type == title[i]] # the group factor of eulittoral/water surface
  data.cover <- cover[macrophyte$Type == title[i]] # the macrophyte cover of eulittoral/water surface
  data.SR <- SR[macrophyte$Type == title[i]] # the macrophyte richness of eulittoral/water surface
  
  test.cover <- Wilcox.test.grp(data.cover,data.grp)
  test.SR <- Wilcox.test.grp(data.SR,data.grp)
  
  result <- cbind(test.cover,test.SR) # cover on the left, and richness on the right
  result.F <- rbind(result.F,result) # eulittoral zone at the upper part of the matrix, and water surface on the lower part
  
  if(i==2)break
}
result.F 
# left: cover; right: SR
# Upper: Eulittoral; Lower: Water surface 

#   2.4 Comparison of macrophyte weight among lakes
Wilcox.test.grp(wet.dry_weight$Wet_weight,title.weight)
Wilcox.test.grp(wet.dry_weight$Dry_weight,title.weight)


#3. BENTHOS SECTION#####
benthos <- read.csv('Site-species data.csv')  # without the Qualitative data
ben.title <- c('Site','Lake','Site.number','Type','Efforts',rep('Spe',66))
benthos.data <- benthos[,ben.title=='Spe']

# Example of the data:
#         Site Lake Site.number Type  Efforts.m2. z.Cli1 z.Cli2 Gast3 Gast4 Gast5
# 1 L1.S4.June   L1         S4    eu.       0.484      1      0     1     0    44
# 2 L1.S4.June   L1         S4   sub.       0.240      0      1     5     0     8
# 3 L1.S4.June   L1         S4 w.suf.       0.000      0      0     0     0     0
# 4 L1.S1.June   L1         S1    eu.       0.484      0      0     2     0     4
# 5 L1.S1.June   L1         S1   sub.       0.240      0      0     0     0     4

# 3.1 SR and density at site level #####
ben.lake.eu <- benthos[benthos$Type=='eu.',]
ben.lake.sub <- benthos[benthos$Type == 'sub.',]
ben.lake.w.surf <- benthos[benthos$Type == 'w.surf.',]

# combine the data as list
ben.list <- list(ben.lake.eu,ben.lake.sub,ben.lake.w.surf)

# Compare of SR and density among lakes
i=0;result.F <- NULL;mean.se.F <- NULL
repeat{
  i=i+1
  data.set <- ben.list[[i]]
  
  # remove the empty row
  data.set <- data.set[data.set$Efforts.m2.>0,]
  
  data <- data.set[,ben.title=='Spe']
 
  # Density
  density <- rowSums(data)/(data.set$Efforts.m2)
  # SR
  data[data>0]=1
  SR <- rowSums(data)
  grp <- data.set$Lake # L1, L2, L3
  
  density.SR <- cbind(density,SR)
  mean.se <- Mean.SE.grp(density.SR,grp)
  
  # Wilcox test among groups (i.e., L1, L2, L3)
  t1 <- Wilcox.test.grp(density,grp)
  t2 <- Wilcox.test.grp(SR,grp)
  result <- cbind(t1,t2)                # Left: Density; Right:SR
  result.F <- rbind(result.F,result)    # Row names: Eulittoral, Sublittoral, Water surface
  mean.se.F <- rbind(mean.se.F,mean.se) # Row names: Eulittoral, Sublittoral, Water surface
  if(i==3) break
}
result.F  # colume : Lake Eulittoral, Sublittoral, Water surface; Row : Density, SR
mean.se.F # colume : Lake (L1, L2, L3);  Row : Eulittoral, Sublittoral, Water surface

#4. NMDS and ANONSIM ####

# 4.1 NMDS calculation #####
ben.data <- benthos[,ben.title=="Spe"]
ben.rm.empty <- ben.data[rowSums(ben.data)>0,] 

ben.lg <- log(ben.rm.empty +1) # have to log transformed, otherwise the NMDS cannot display properly
Jac.dist <- Jac.dis(ben.lg)
nmds <- metaMDS(Jac.dist)
nmds$stress

# The coordinates of each sites
score <- scores(nmds)               
nmds.coord <- as.data.frame(score)  # convert coordinates to data frame

# 4.2 Plot NMDS #####

# Setting two groups

#   4.2.3 Setting color factors based on the Lake and sampling zone
{
  grp1 <- benthos$Type[rowSums(ben.data)>0]
  grp2 <- benthos$Lake[rowSums(ben.data)>0]

  main <- 'Lake 1 2 3 and eulittoral, sublittoral and water surface'
  
  pch <- c(15,16,17)
  Replace(grp1,pch)
  pch <- Replace(grp1,pch)[[2]]
  
  color <- c('#CC79A7','#0072B2','#E69F00')
  Replace(grp2,color)
  col <- Replace(grp2,color)[[2]]
}

#Plot NMDS with group factors of Lake and Littoral zone
plot(nmds.coord     # use the coordinate take out from nMDS 
     ,col= col      # separate the sites use different color based on grp
     ,pch=pch
     ,main=main
)
legend('bottomleft',c('Lake 1','Lake 2','Lake 3'
                      ,'Eulitoal zone','Sublittoral zone', 'Water surface')
       ,col = c('#CC79A7','#0072B2','#E69F00'
                ,'black','black','black')
       ,pch = c(1,1,1,
                15,16,17))

# Fit the species distribution on the nMDS plot
fit <- envfit(nmds,log(ben.rm.empty+1),perm = 999)

plot(fit, p.max = 0.05, col = "black")
nmds.spe <- fit$vectors
nmds.spe.p <- cbind(nmds.spe$arrows,nmds.spe$r,nmds.spe$pvals)
colnames(nmds.spe.p) <- c('NMDS1','NMDS2','R2','p-value')
nmds.spe.p

# 4.3 ANOSIM analysis #####

#Dissimilarity of eulittoral zone
Type.title <- benthos$Type
Lake.title <- benthos$Lake

Type.title <- Type.title[rowSums(ben.data)>0] # remove empty
Lake.title <- Lake.title[rowSums(ben.data)>0] # remove empty


ben.eu.lg <- ben.lg[Type.title=='eu.',]
lake.title.eu <- Lake.title[Type.title=='eu.']

#Dissimilarity of sublittoral zone
ben.sub.lg <- ben.lg[Type.title=='sub.',]
lake.title.sub <- Lake.title[Type.title=='sub.']

#Dissimilarity of water surface
ben.w.surf.lg <- ben.lg[Type.title=='w.surf.',]
lake.title.w.surf <- Lake.title[Type.title=='w.surf.']

#   4.3.2 ANOSIM of community dissimilarity among lakes #####
Janc.dis <- Jac.dis(ben.eu.lg)  # BC and Jaccard distance are very similiar in terms of nMDS.
mod1 <- anosim(Janc.dis, grouping = lake.title.eu)
mod1

Janc.dis <- Jac.dis(ben.sub.lg)
mod2 <- anosim(Janc.dis, grouping = lake.title.sub)
mod2

Janc.dis <- Jac.dis(ben.w.surf.lg)
mod3 <- anosim(Janc.dis, grouping = lake.title.w.surf)
mod3

par(mfrow = c(3,1))
plot(mod1);plot(mod2);plot(mod3)
par(mfrow = c(1,1))

#6. CCA #####

# 6.1 Prepare the variables
#   6.1.1 Other measured variables 
### Attention ###
"I did not put breadth and canopy cover into the analysis
because there are limited number of variables in the CCA
and the variables are not the directly impact the macroinvertebrates"

  w.trans <- env$Water.transparency; w.trans <- w.trans[is.na(w.trans)==FALSE]

  macrophyte.cover.eu <- rowSums(macrophyte[,-c(1:6)])[env$type=='eu.']
  macrophyte.cover.w.surf<- rowSums(macrophyte[,-c(1:6)])[env$type=='sub.']

#   6.1.2 Macrophyte cover  
  # Eulittoral zone
  macrophyte.eu <- cbind(macrophyte[,macrophyte.title=='species'],macrophyte.cover.eu)[macrophyte$Type=='eu.',]
  macrophyte.binary <- macrophyte.eu; macrophyte.binary[macrophyte.binary>0]=1
  
  # Remove rare species and add the total cover of macrophyte
  macrophyte.eu <- macrophyte.eu[,colSums(macrophyte.binary)>1]
  
  # Water surface
  macrophyte.w.surf <- cbind(macrophyte[,macrophyte.title=='species'],macrophyte.cover.w.surf)[macrophyte$Type=='w.surf',]
  macrophyte.binary <- macrophyte.w.surf; macrophyte.binary[macrophyte.binary>0]=1
  
  # Remove rare species and add the total cover of macrophyte
  macrophyte.w.surf <- macrophyte.w.surf[,colSums(macrophyte.binary)>1]
  
#   6.1.3 Water quality data
  water.quality <- (env[,env.title== 'water qaulity'])[env$type=='eu.',]

#   6.1.4 Substrate
  substrate.eu <- (env[,env.title=='substrate'])[env$type=='eu.',]; substrate.eu.scale <- scale(substrate.eu,center=TRUE,scale=TRUE)
  substrate.sub <- (env[,env.title=='substrate'])[env$type=='sub.',]; substrate.sub.scale <- scale(substrate.sub,center=TRUE,scale=TRUE)

#   6.1.5 Benthos
  ben.eu <- benthos[benthos$Type=='eu.',];ben.eu.data <- ben.eu[,-c(1:6)]
  ben.sub <- benthos[benthos$Type=='sub.',];ben.sub.data <- ben.sub[,-c(1:6)]
  ben.w.surf <- benthos[benthos$Type=='w.surf.',];ben.w.surf.data <- ben.w.surf[,-c(1:6)]


# 6.2 Gather all the variables in lists. Therefore, it would be easier to be included in a loop
  env.eu    <- cbind(water.quality,w.trans)
  env.sub <- cbind(water.quality,w.trans)
  env.w.surf <- cbind(water.quality,w.trans)
  env.eu.scale <- scale(env.eu)
  env.sub.scale <- scale(env.sub)
  env.w.surf.scale <- scale(env.w.surf)
  
  ben.list <- list(ben.eu,ben.sub,ben.w.surf)
  
  ben.data.list <- list(ben.eu.data,ben.sub.data,ben.w.surf.data)
  
  # Organize the independen variable based on:
  #     Eulittoral zone: water quality, substrate type and macrophyte
  #     Sublittoral zone: water quality and substrate type
  #     Water surface: water quality and macrophyte 
  explanatory.list <- list(
    list(env.eu.scale[,-7],substrate.eu.scale[,-4],macrophyte.eu[,-c(6,13)])# rm NaNO2, gravel stone, sp6 and sp15, becaus of low occurence and value
    ,list(env.sub.scale,substrate.sub)
    ,list(env.w.surf.scale,macrophyte.w.surf[,-c(4,7,8)]) #rm sp7, sp13, sp14 because of low occurence and percentage
  )


title <- c('Eulittoral zone','Sublittoral zone', 'Water surface')

# 6.3 A loop for CCA of three zones#####
cca.hp <- function(){
  i=0;result.F <- NULL
  par(mfrow=c(2,2))
  repeat{
    exp1 <- NULL
    exp2 <- NULL
    exp3 <- NULL
    
    i=i+1
    
    # input the response variables
    ben.data <- ben.data.list[[i]]
    
    # input the independent variable
    exp <- explanatory.list[[i]]
    
    n <- length(exp) # how many response variables. If n=3, it is for the Eulittoral zone
    
    # remove empty rows
    
    ###              ###
    ### Attention!!! ###   
    ##               ###
    # Also need to rm the empty rows in the lake type, which are used to group the site in the CCA plot 
    lake.type <- (ben.list[[i]])$Lake; lake.type <- lake.type[rowSums(ben.data)>0]
    
    exp1 <- (exp[[1]])[rowSums(ben.data)>0,]
    exp2 <- (exp[[2]])[rowSums(ben.data)>0,] 
    if (n==3) exp3 <- (exp[[3]])[rowSums(ben.data)>0,]
    
    ben.data <- ben.data[rowSums(ben.data)>0,]
    ben.data <- ben.data[,colSums(ben.data)>0]
    
    # re-organize the independent variables into one matrix for further analyses
    if(n==2) explanatory <- cbind(exp1,exp2)
    if(n==3) explanatory <- cbind(exp1,exp2,exp3)
  
    ben.lg <- log(ben.data+1)
    decorana(ben.lg) # DCA analysis (<3 use RDA; 4 use CCA)
  
#   6.3.1 CCA analysis ######
     explanatory <- as.data.frame(explanatory)
     spe.cca <- cca(ben.lg~.,data=explanatory) # ben.lg with all the column in explanatory vairables
     vif.cca(spe.cca) # co-linear coeffecience (the small the better, generally need < 10)
     
#   6.3.2 CCA plot #####
      
      #       6.4.1 default plot ####
      pl <- plot(spe.cca,scaling=3,main=title[i])

      #       6.4.2 plot from a blank ####
      pl <- plot(spe.cca, type="none", scaling=3, correlation=TRUE)
      # 'arg' should be one of “text”, “points”, “none”

      #           6.4.2.1 plot sites####
      pch <- Replace(lake.type,c(1,2,3))[[2]] #  replace the lake type with 1,2,3
      points(pl, "site", pch=pch
             ,lwd=1                       # line width
             ,col=1,cex=1.5)

      #           6.4.2.2 plot species ####
      fit <- envfit(spe.cca,ben.data,perm = 9999)
      
      plot(fit, p.max = 0.05, col = "black")
      
      # save the results from envfit
      cca.spe <- fit$vectors
      cca.spe.p <- cbind(cca.spe$arrows,cca.spe$r,cca.spe$pvals)

      #           6.4.2.3 plot variables ####
      text(spe.cca, display ="bp", scaling=3,col='blue') # plot arrows and variables
      # 'arg' should be one of “sites”, “species”, “wa”, “lc”, “bp”, “reg”, “cn”
      # bp=variables and arrow
      
#   6.3.3 CCA statistical analysis #####
  # Heriacical partition - HP

    # P.S: I waited four days to finish the analysis, which includes 12 variables. It took too much time...

    # 6.3.3.1 This is what others normally do
    # sig.var <- anova.cca(spe.cca, by="term", step=1000)   # The sig. of every each variable
    sig.model <- anova.cca(spe.cca, step=999)            # Permutation test of result of CCA/RDA
    sig.model
    
    # 6.3.3.2 other methods: RDACCA.hp to calculate the explained variance of each variable
    mod <- rdacca.hp(ben.lg,as.data.frame(explanatory), method="CCA", var.part=TRUE, type="adjR2")
    r.square <- mod$Total_explained_variation
    hier.part <- mod$Hier.part
    permu  <- permu.hp(ben.lg, as.data.frame(explanatory), method="CCA",type="adjR2", permutations=999)
    
    # the interaction among different group of variables (Macrophyte : Substrate : Water quality)
    if(n==2) iv <- list(as.data.frame(exp1),as.data.frame(exp2)) # water quality, susbtrate/macrophyte
    if(n==3) iv <- list(as.data.frame(exp1),as.data.frame(exp2)
                        ,as.data.frame(exp3))                        # water quality, substrate, macrophyte

    model.interaction <- rdacca.hp(ben.lg, iv, method="CCA",var.part=TRUE, type="adjR2")
    r.square.inter <- model.interaction$Total_explained_variation # R square = 1
    hier.part.inter <- model.interaction$Hier.part
    permu.inter <- permu.hp(ben.lg, iv, method="CCA", type="adjR2", permutations=999)

    result <- list(title[i]
                  #,sig.model,r.square,hier.part
                  ,permu
                  #,r.square.inter,hier.part.inter
                  ,permu.inter
                  )

    result
    result.F <- list(result.F,result)
      
     if(i==length(explanatory.list))break
  }
  par(mfrow=c(1,1))
  return(result.F)
}
cca.hp()

result.list <- cca.hp()
result.list

# 6.4 output results #####
capture.output(result.list,file='XX.text')

#7. MACROPHYTE cover and BENTHOS abundance #####
# 7.1 indi/sampling area ##### 

# ben
ben.w.surf <- benthos[benthos$Type == 'w.surf.',]
judge <- ben.w.surf$Efforts.m2

ben.eu <- benthos[benthos$Type == 'eu.',]
# rm empty
ben.rm <- ben.w.surf[judge>0,]
ben.data <- ben.rm[,ben.title == 'Spe']

# benthos density and SR
abundance.w.surf <- rowSums(ben.data)
density.w.surf <- abundance.w.surf/0.24  # indi/m2

ben.data.binary <- ben.data
ben.data.binary[ben.data.binary>0]=1
SR.w.surf <- rowSums(ben.data.binary)

# 7.2 macrophyte cover
cover <- rowSums(macrophyte[,macrophyte.title=='species'])
cover.w.surf <- cover[macrophyte$Type=='w.surf'];cover.w.surf <- cover.w.surf[judge>0]

lake <- ben.rm$Lake
pch.mid <- Replace(lake,c(1,2,3))[[2]]
plot(density.w.surf~cover.w.surf,pch=pch.mid
     #,ylim=c(0,15000)
     ,main='Mid P: Density(indi/m2)/Cover')
legend('topleft',legend=c('Lake 1','Lake 2','Lake 3'),pch=c(1,2,3))

# 7.3 Wilcox.test of density among lakes
Wilcox.test.grp(density.w.surf,lake)

matrix <- cbind(density.w.surf,rep(1,length(density.w.surf))) # to combine a column to perform the function Mean.SE.grp. Because this function can only deal with matrix

Mean.SE.grp(matrix,lake)


# 7.5 macroinvertebrate density and macrophyte in both Eulittoral zone and Water surface
#   Eulittoral zone
ben.eu.data <- ben.eu[,ben.title == 'Spe']
abundance.eu <- rowSums(ben.eu.data)
density.eu <- abundance.eu/ben.up$Efforts.m2
cover.eu <- cover[macrophyte$Type=='eu.']

lake <- ben.eu$Lake
pch.up <- Replace(lake,c(1,2,3))[[2]]
plot(density.eu ~ cover.eu,pch=pch.up
     #,ylim=c(0,3000)
     )

# Wilcox.test of density among lakes
Wilcox.test.grp(density.eu,lake)

lm.eu <- gamlss(density.eu ~ cover.eu)
lm.w.surf <- gamlss(density.w.surf~ cover.w.surf)

summary(lm.up);Rsq(lm.up)
summary(lm.midP);Rsq(lm.midP)

plot(lm.up)
plot(lm.midP)

par(mfrow=c(1,2))

plot(density.eu ~ cover.eu,pch=pch.up, main='Eullitoral zone: Density(indi m-2)/Cover')
abline(lm.eu,col='red')

plot(density.w.surf~cover.w.surf,pch=pch.mid, main='Water surface: Density(indi m-2)/Cover')
legend('topleft',legend=c('Lake 1','Lake 2','Lake 3'),pch=c(1,2,3))
abline(lm.w.surf,col='red')

par(mfrow=c(1,1))
