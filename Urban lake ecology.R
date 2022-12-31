# Urban lake ecology
setwd('D:\\PhD\\20211204 Lake ecological study\\Statistical analysis')
#library(eoffice)
#library(ggplotify)
library(vegan)
library(ecodist) # MRM
library(corrplot) # Correlation plot
library(gamlss)
library(betapart)
# Define functions #####
# Calculate the mean and SE of a matrix
# 'D:\PhD\Methods in PhD\Mean.SE.R'
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
# 'D:\PhD\Methods in PhD\Wilcox.test.grp.R'
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
# 'D:\PhD\Methods in PhD\Element.grp.R'
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
# 'D:\PhD\Methods in PhD\Replace.R'
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
# draw the correlation plot among environmental variables
Correlation.plot <-function(data){
  cor <-cor(data)
  test.p.value <- cor.mtest(data)$p # p value
  corrplot(cor
           #,method = 'number'   #just show the values of correlation
           ,method = 'color'     # fill the space with color
           ,order = 'AOE'        # order the values of correlation
           ,addCoef.col = 'black'# add the values of correlation
           ,type = c('lower')    # type: display the lower part of the matrix
           ,diag = TRUE          # display the correlation coefficients
           ,cl.pos = 'n'         # remove the correlation table
           ,tl.pos = 'd'         # show the variables in the middle 
           ,p.mat = test.p.value, sig.level = 0.05
           ,insig ='blank'       # values higher than 0.05 are not shown
  ) 
}
# Calculate abundance-based Jaccard dissimilarity
# The Jaccard dissimilarity was calculated by Bray-curtis distance:
# Jac.dis = 2B/(b+1)
Jac.dis <- function(data){
  BC.dis <- vegdist(data,method = 'bray')
  Jaccard.distance = (2*BC.dis)/(BC.dis+1)
  return(Jaccard.distance)
}
# Function Boxplot.grp() here you just need to input the data list and grp factor
Boxplot.grp <- function(list.of.data,grp,main,col,lty){
  n = length(list.of.data)
  
  i=0;group <- NULL;data.F <- NULL
  repeat{
    i=i+1
    data <- list.of.data[[i]]
    
    # transform the data into vactors
    data.F <- c(data.F,as.vector(data))
    
    # enlarge the grp factor
    group <- c(group,rep(grp[i],length(data)))
    if(i==n)break
  }
  
  # sort the group based on what you want
  group <- factor(group,levels = grp)
  boxplot(data.F~group,main=main,col=col,lty=lty)
  
}
# 1. ENVironmental variables #####
env <- read.csv('Environmental data.csv')
env.title <- c('Site', 'Lake', 'Site.number', 'Month', 'type', 'water transparency', 'macrophyte breath'
               ,rep('substrate',4),rep('water quality',7)
               ,rep('canopy',2))
#   1.1 Do water qualities are different from middle and upper zone? #####
{
  # water.q <- env[,env.title == 'water qaulity']
  # 
  # title <- c('Lake 1','Lake 2','Lake 3')
  # list <- c('L1','L2','L3')
  # i=0;Result.F <- NULL
  # repeat{
  #   i=i+1
  #   title <- title[i]; data <- water.q[env$Lake == list[i],]
  #   up.mid <- env$type[env$Lake == list[i]]
  #   data.up <- data[up.mid=='upper',]; data.mid <- data[up.mid=='middle',]
  #   
  #   j=0;rowname<-NULL;Result.lake<-NULL
  #   repeat{
  #     j=j+1
  #     test <- wilcox.test(data.up[,j],data.mid[,j]) # compare the j th column
  #     w.statistic <- test$statistic; p.value <- test$p.value
  #     rowname <- c(rowname,colnames(data)[j])
  #     test.result <- c( w.statistic, p.value )
  #     Result.lake <- rbind(Result.lake,test.result)
  #     if(j == ncol(data)) break
  #   }
  #   rownames(Result.lake) <- rowname
  #   
  #   Result.F <- list(Result.F,title,Result.lake)
  #   if(i==3) break
  # }
  # Result.F
}
#   Results of comparison
  {
  # Lake 1
  # DO...       36.0 0.71298949
  # EC..s.cm.   55.5 0.01533316 **
  # Water.temp  35.5 0.75236127
  # pH          28.5 0.75057571
  # NO3.N.mg.L. 32.0        NaN
  # NH3.N.mg.L. 26.5 0.58022394
  # P2O5.mg.L.  34.0 0.86917024
  # NaNO2..g.L. 25.0 0.43280248
  
  # Lake 2
  # DO...       37.0 0.6453768
  # EC..s.cm.   35.0 0.7927465
  # Water.temp  36.5 0.6730617
  # pH          33.0 0.9580911
  # NO3.N.mg.L. 32.0       NaN
  # NH3.N.mg.L. 22.0 0.3144789
  # P2O5.mg.L.  25.5 0.5150163
  # NaNO2..g.L. 32.5 1.0000000
  
  # Lake 3
  # DO...       46.0 0.15594939
  # EC..s.cm.   43.5 0.24729945
  # Water.temp  38.5 0.52524234
  # pH          41.0 0.37132498
  # NO3.N.mg.L. 32.0        NaN
  # NH3.N.mg.L. 31.0 0.95696844
  # P2O5.mg.L.  38.5 0.52647600
  # NaNO2..g.L. 12.0 0.01277699 **

}
# only EC in the L1, and NO2 in the Lake 3 have significant different
# We then use the env at upper only for the whole site.

#   1.2 Mean ± SE of each env among lakes and littoral zones #####
    # All env occurs in Upper, only substrate and macrophyte.cover occur in Middle.

# 1.2.1 At lake level

# we calculate the substrate differently

# first 5 are site elements
env.substrate.title <- c(rep('yes',5),rep('no',2),rep('yes',4),rep('no',10))
env.other.env.title <- c(rep('yes',5),rep('yes',2),rep('no',4),rep('yes',10))

substrate <- env[,env.substrate.title == 'yes']
other.env <- env[,env.other.env.title == 'yes'] # other env that are not substrate

judge <- is.na(rowSums(other.env[,-c(1:5)]))
other.env.rm.na <- other.env[judge==FALSE,]

# Mean and SE of variables of three lakes
title <- c('L1','L2','L3')
i=0;result.F <- NULL
repeat{
  
  i=i+1
  data.substrate <- substrate[substrate$Lake == title[i],]
  data.other.env <- other.env.rm.na[other.env.rm.na$Lake == title[i],]
  
  data.substrate.onlydata <- data.substrate[,-c(1:5)]
  data.other.env.onlydata <- data.other.env[,-c(1:5)]
  
  # Substrate
    # substrate upper
      data.upper <- data.substrate.onlydata[data.substrate$type=='upper',]
      mean.se.substrate.upper<- Mean.SE(data.upper)
    
   # substrate middle
      data.middle <- data.substrate.onlydata[data.substrate$type=='middle',]
      mean.se.substrate.middle<- Mean.SE(data.middle)

  # Other env  
      mean.se.otherenv<- Mean.SE(data.other.env.onlydata)
  
  result <- rbind(mean.se.substrate.upper
                  ,mean.se.substrate.middle
                  ,mean.se.otherenv)
  
  result.F<-cbind(result.F,result)
  
  if(i==length(title))break
}

result.F <- round(result.F,2)
result.F
# Substrate
#     uppper
#     middle
# Macrophyte
# Water quality
# Canopy cover

#write.csv(result.F,'221128 mean se of env variables.csv')

#   1.3 Comparisons of variables among three lakes #####

# 1.3.1 Variables in the upper among lakes
env.upper <- env[env$type == 'upper',]

i=5 # because the comparison starts from the 6th column
result.F <- NULL
repeat{
  i=i+1
  title.env <- colnames(env.upper)[i]
  data = env.upper[,i]
  
  result <- Wilcox.test.grp(data,env.upper$Lake)
  result.F <- rbind(result.F,title.env,result)
  
  if(i== ncol(env.upper) )break
}
result.F
#write.csv(result.F,'221128 comparison of env among lakes.upper.csv')

# 1.3.2 Variables in the middle among lakes (just substrate)
env.middle <- env[env$type == 'middle',]
env.middle <- env.middle[,-c(6,7,12:21)]

i=5;result.F <- NULL
repeat{
  i=i+1
  title.env <- colnames(env.middle)[i]
  data = env.middle[,i]
  # t.1.2 <- wilcox.test(data[env.middle$Lake=='L1'],data[env.middle$Lake=='L2']) # L1 vs L2
  # t.1.3 <- wilcox.test(data[env.middle$Lake=='L1'],data[env.middle$Lake=='L3']) # L1 vs L3
  # t.2.3 <- wilcox.test(data[env.middle$Lake=='L2'],data[env.middle$Lake=='L3']) # L2 vs L3
  # 
  # result <- c(title.env,t.1.2$statistic,t.1.2$p.value
  #             ,t.1.3$statistic,t.1.3$p.value
  #             ,t.2.3$statistic,t.2.3$p.value)
  
  result <- Wilcox.test.grp(data,env.middle$Lake)
  
  result.F <- rbind(result.F,title.env,result)
  
  if(i== ncol(env.middle) )break
}
result.F
#write.csv(result.F,'221128 comparison of env among lakes middle.csv')

# 2. MACROPHYTE SECTION#####
macrophyte <- read.csv('Macrophyte cover & weight.csv')
macrophyte.title <- c('site','lake','site.number','month','type'
                      ,'weight','weight'
                      ,rep('species',15))

macrophyte.data <- macrophyte[,macrophyte.title=='species']

#   2.1 SR at lake level #####
macrophyte.col <- colSums(macrophyte.data)
macrophyte.col[macrophyte.col>0]=1
SR <- sum(macrophyte.col);SR
"Total species richness in three lakes:"
SR

# Separate the data into different lakes
title <- Element.grp(macrophyte$Lake)
n <- length(title)
i=0;result.F <- NULL
repeat{
  i=i+1;data.title <- title[i]

# Calculate the SR at lake level
  # separate different lakes
  data = macrophyte.data[macrophyte$Lake == title[i],]
  data.col <- colSums(data);data.col[data.col>0]=1
  data.SR <- sum(data.col)
  
  # # species accumulation curve
  # rarefaction <- specaccum(data, method = "rarefaction"
  #                          , permutations = 1000,conditioned =TRUE
  #                          , gamma = "chao1",  w = NULL)
  # plot(rarefaction
  #      , add = FALSE, random = FALSE
  #      ,col = "Blue", lty = 1       # the color of the accmulation curve
  #      ,ci = 2, ci.type = c("line") # setting of the ci (confidence interval)
  #      ,ci.col = "red", ci.lty = 2  # setting of the ci
  #      , ci.length = 10             # setting of the ci
  #      #,xlim=c(0,300), ylim=c(0,35)# limits of x- and y- axis
  #      , ylim=c(0,13)
  #      ,xvar = c("site")            # based on "individual" or 'site'
  #      ,main = title[i]
  # ) # ci - confidence interval)
  
  # separate upper and middle
  title.up.mid <- (macrophyte$Type)[macrophyte$Lake == title[i]]
  data.up <- data[title.up.mid=='upper',]; data.mid <- data[title.up.mid=='middle',]
  
  data.col <- colSums(data.up); data.col[data.col>0]=1; SR.up <- sum(data.col)
  data.col <- colSums(data.mid);data.col[data.col>0]=1; SR.mid <- sum(data.col)
  
  result<- c(data.title,data.SR,SR.up,SR.mid)
  result.F <- rbind(result.F,result)
  colnames(result.F) <- c('Lake','Total', 'Upper','Middle')
  
  if(i==n)break
}
result.F
#write.csv(result.F,'221128 Species richness at lake level.csv')

# Results of macrophyte SR: tota; Upper; Middle
#   Lake 1: 12;11;8
#   Lake 2: 9;9;4
#   Lake 3: 4;4;2


#   2.2 Mean SE of weight, cover and SR at site level #####
  # weight
wet.dry_weight <- macrophyte[,macrophyte.title=='weight']
title.weight <- macrophyte$Lake[macrophyte$Wet_weight>0]
wet.dry_weight <- wet.dry_weight[macrophyte$Wet_weight>0,]
  
  # cover and SR
cover <- rowSums(macrophyte.data)
data <- macrophyte.data;data[data>0]=1
SR <- rowSums(data)

'Vegetation Cover and SR at site level of three lakes:'
cover.SR <- cbind(cover,SR)
Mean.SE(cover.SR)

# cover and SR among littoral zone and lakes
title <- macrophyte$Lake;
title.level <- Element.grp(title)

i=0;result.cover.SR.F<-NULL;result.weight.F <- NULL
repeat{
  i=i+1
  # separate the lake
  data.cover.SR <- cover.SR[title == title.level[i],]
  weight <- wet.dry_weight[title.weight == title.level[i],]
  title.up.mid <- (macrophyte$Type)[title == title.level[i]]
  
  mean.se.lake <- Mean.SE(data.cover.SR)
  mea.se.weight <- Mean.SE(weight)
  # separate the upper and middle
  data.up <- data.cover.SR[title.up.mid=='upper',]; data.mid <- data.cover.SR[title.up.mid=='middle',]
  
  mean.se.up <- Mean.SE(data.up)
  mean.se.mid <- Mean.SE(data.mid)
  
  result.cover.SR <- cbind(mean.se.lake,mean.se.up,mean.se.mid)
  result.weight.F <- rbind(result.weight.F,mea.se.weight)
  
  result.cover.SR.F <- rbind(result.cover.SR.F,result.cover.SR)
  
  if(i==3)break
}
result.cover.SR.F <- round(result.cover.SR.F,2)
result.cover.SR.F <- rbind(c('Lake level','','Upper','','Middle','')
                  ,result.cover.SR.F)
result.cover.SR.F
result.weight.F
#write.csv(result.cover.SR.F,'221128 macrophyte cover and SR at site level.csv')
#write.csv(result.weight.F,'221130 macrophyte dry and wet weight at site level.csv')

# Results of cover at the finest level (littoral zone)
#   Lake 1: 44.625 ± 9.16464
#   Lake 2: 58.875 ± 7.19657
#   Lake 3: 76.9375 ± 4.596619

#   2.3 Boxplot of Cover and SR among littoral zones of three lakes #####
grp <- paste(macrophyte$Lake,macrophyte$Type, sep = '.')
# grp <- factor(grp,levels = c('L1.upper','L2.upper','L3.upper'
#                          ,'L1.middle','L2.middle','L3.middle')
#               )

grp <- factor(grp,levels = c('L1.upper','L1.middle'
                          ,'L2.upper','L2.middle'
                          ,'L3.upper','L3.middle')
              )

par(mfrow = c(2,1))
boxplot(cover~grp,lty=c(1,2),col='white')
boxplot(SR~grp,lty=c(1,2),col='white')
par(mfrow = c(1,1))

#   2.4 Barplot of macrophyte cover #####
# #barplot.macrophyte.cover <- function(){
#   par(mfrow=c(1,2))
#   barplot(average.cover, ylim = c(0,90)
#           ,names.arg = c('Lake 1', 'Lake 2', 'Lake 3')
#           ,ylab = 'Macrophyte cover (%)'
#           )
#   i=0
#   repeat{
#     i=i+1
#     arrows(x0=i,y0=average.cover[i],x1=i,y1=average.cover[i]+se[i],angle = 90,lwd=2)
#     if (i == 3) break
#   }
# 
#   plot(SR,ylim=c(0,12))
#   par(mfrow=c(1,1))
# #}
# #barplot.macrophyte.cover()
# 
#  # f.macrophyte.cover <- '221114 macrophyte cover.pptx'
#  # p.macrophyte.cover <- as.ggplot(~barplot.macrophyte.cover())
#  # topptx(p.macrophyte.cover,f.macrophyte.cover)


#   2.5 Comparison of macrophyte cover among lakes #####
title <- c('upper','middle')
i=0;result.F <- NULL
repeat{
  i=i+1
  data.grp <- macrophyte$Lake[macrophyte$Type == title[i]]
  data.cover <- cover[macrophyte$Type == title[i]]
  data.SR <- SR[macrophyte$Type == title[i]]
  
  t.cover <- Wilcox.test.grp(data.cover,data.grp)
  t.SR <- Wilcox.test.grp(data.SR,data.grp)
  
  # Cover , , ,  SR , , ,
  result <- cbind(t.cover,t.SR)
  
  # Upper
  # Middle
  result.F <- rbind(result.F,result)
  
  if(i==2)break
}
result.F
#write.csv(result.F,'221128 comparison of macrophyte cover and SR at site level among lake.csv')

#   2.6 Comparison of macrophyte weight among lakes
# wet.dry_weight <- macrophyte[,macrophyte.title=='weight']
# title.weight <- macrophyte$Lake[macrophyte$Wet_weight>0]
# wet.dry_weight <- wet.dry_weight[macrophyte$Wet_weight>0,]
Wilcox.test.grp(wet.dry_weight$Wet_weight,title.weight)
Wilcox.test.grp(wet.dry_weight$Dry_weight,title.weight)

# 3. BENTHOS SECTION#####
setwd('D:\\PhD\\20211204 Lake ecological study\\Statistical analysis')
benthos <- read.csv('Site-species data.csv')  # without the Qualitative data
ben.title <- c('Site','Lake','Site.number','Month','Type','Efforts',rep('Spe',66))
benthos.data <- benthos[,ben.title=='Spe']

#   3.1 SR and density of three lakes #####
# colsum <- colSums(benthos.data); effort <- sum(benthos$Efforts.m2.or.g.)
# Density <- sum(colsum)/effort                                            # cannot simply calculated anymore.. because I changed the effort of periphytic to dry weight
# colsum[colsum>0]=1;SR <- sum(colsum)
# 
# Density  # 1939.298 indi/m2
# SR       # 66 species

#   3.2 SR and species composition among lakes #####

# How many lakes
title <- Element.grp(benthos$Lake)
n <- length(title)

# transfer the first three letter of the species name into initials
spe.name <- colnames(benthos.data)
initial <- Left(spe.name,3)
order <- Element.grp(initial)# Which are higher orders

# Starts to calculate the SR, composition based on taxonomic and composition based on density
i=0;SR.result <- NULL;taxa.result.F.F <- NULL
repeat{
  i=i+1
  data <- benthos.data[benthos$Lake==title[i],]
  ben.col <- colSums(data)
  ben.col[ben.col>0]=1
  
  # # species accumulation curve
  # rarefaction <- specaccum(data, method = "rarefaction"
  #                          , permutations = 1000,conditioned =TRUE
  #                          , gamma = "chao1",  w = NULL)
  # plot(rarefaction
  #      , add = FALSE, random = FALSE
  #      ,col = "Blue", lty = 1       # the color of the accmulation curve
  #      ,ci = 2, ci.type = c("line") # setting of the ci (confidence interval)
  #      ,ci.col = "red", ci.lty = 2  # setting of the ci
  #      , ci.length = 10             # setting of the ci
  #      #,xlim=c(0,300), ylim=c(0,35)# limits of x- and y- axis
  #      , ylim=c(0,45)
  #      ,xvar = c("site")            # based on "individual" or 'site'
  #      ,main = title[i]
  # ) # ci - confidence interval)
  
# SR
  SR <- sum(ben.col)
  SR.result <- c(SR.result,SR)
  
# Species composition
  # Remove the non-existence species
  data <- initial[ben.col>0]
  
  # See how many taxa 
  t=0;taxa.result.F <- NULL
  repeat{
    t=t+1
    # Count how many initials as how many taxa under this order
    taxa.result <- length(data[data==order[t]])
    taxa.result.F <- c(taxa.result.F,taxa.result)
    if(t==length(order))break
  }
  taxa.result.F.F <- rbind(taxa.result.F.F,taxa.result.F)
  
  if (i==n) break
}
SR.result
#     Results of SR
# Lake 1: 40
# Lake 2: 42
# Lake 3: 39

colnames(taxa.result.F.F) <- order
rownames(taxa.result.F.F) <- title
taxa.result.F.F
#write.csv(taxa.result.F.F,'221129 Taxa composition among lakes.csv')

# barplot
color <- c('green','gray','orange','yellow','red','blue')

{
barplot(t(taxa.result.F.F)
        ,col=color
        ,ylim = c(0,50))
legend('topright',order
       ,col = color
       ,pch=15
       )

}
#   3.3 SR and density at site level #####

# Take out the upper data set
ben.lake.Upper <- benthos[benthos$Type=='Upper',]

# take out the Mid.P and Mid.F data set
ben.lake.Mid.F <- benthos[benthos$Type == 'Mid (F)',]
ben.lake.Mid.P <- benthos[benthos$Type == 'Mid (P)',]

ben.list <- list(ben.lake.Upper,ben.lake.Mid.F,ben.lake.Mid.P)

# Compare of SR and density among lakes
i=0;result.F <- NULL;mean.se.F <- NULL
repeat{
  i=i+1
  data.set <- ben.list[[i]]
  
  # remove the empty row
  data.set <- data.set[data.set$Efforts.m2.or.g.>0,]
  
  data <- data.set[,ben.title=='Spe']
 
  # Density
  density <- rowSums(data)/(data.set$Efforts.m2.or.g.)
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
  result.F <- rbind(result.F,result)    # Row names: Upper, Mid F, Mid P
  mean.se.F <- rbind(mean.se.F,mean.se) # Row names: Upper, Mid F, Mid P
  if(i==3) break
}
result.F  # col name: Lake Upper, Mid F, Mid P; row name: Density, SR
mean.se.F # col name: Lake (L1, L2, L3); row name: Upper, Mid F, Mid P
# write.csv(result.F,'221216 Macroinvertebrate wilcox.test of SR and density at site level.csv')
# write.csv(mean.se.F,'221216 Macroinvertebrate mean ± SE of upper midF midP among different lakes.csv')

# 4. NMDS and ANONSIM ####

#   4.1 NMDS calculation #####
ben.data <- benthos[,ben.title=="Spe"]
ben.rm.empty <- ben.data[rowSums(ben.data)>0,] 

ben.lg <- log(ben.rm.empty +1) # otherwise the NMDS cannot display properly
Jac.dist <- Jac.dis(ben.lg)# use 'BC distance'
nmds <- metaMDS(Jac.dist)
nmds$stress

# take out the coordinates of each sites
score <- scores(nmds)               
nmds.coord <- as.data.frame(score)  # convert coordinates to data frame

#   4.2 Plot NMDS #####

# Setting two groups
  # 4.2.1 Setting color factors based on the sampling zone
{
  grp1 <- benthos$Type[rowSums(ben.data)>0]
  main1 <- 'Sampling zone (Up, Mid P F)'

  color <- c('red','blue','green')
  Replace(grp1,color)
  col1 <- Replace(grp1,color)[[2]]
}
  # 4.2.2 Setting color factors based on the Lake
{
  grp2 <- benthos$Lake[rowSums(ben.data)>0]
  main2 <- 'Lake 1 2 3'

  color <- c('red','blue','green')
  Replace(grp2,color)
  col2 <- Replace(grp2,color)[[2]]
}
# 4.2.3 Setting color factors based on the Lake and sampling zone
{
  grp1 <- benthos$Type[rowSums(ben.data)>0]
  grp2 <- benthos$Lake[rowSums(ben.data)>0]

  main3 <- 'Lake 1 2 3 and Upper MidF MidP'
  color <- c('red','blue','green')
  Replace(grp2,color)
  col2 <- Replace(grp2,color)[[2]]
  
  pch <- c(15,16,17)
  Replace(grp1,pch)
  pch <- Replace(grp1,pch)[[2]]
}

#   4.3 ANOSIM analysis #####
main1;anosim(Jac.dist, grouping = grp1)
main2;anosim(Jac.dist, grouping = grp2)

#     4.3.1 Dissimilarity of different habitats #####
Lake.title <- (benthos[rowSums(ben.data)>0,])$Lake
Type.title <- (benthos[rowSums(ben.data)>0,])$Type

#       Dissimilarity of upper
ben.up.lg <- ben.lg[Type.title=='Upper',]
lake.title.up <- Lake.title[Type.title=='Upper']

#       Dissimilarity of MidF
ben.midF.lg <- ben.lg[Type.title=='Mid (F)',]
lake.title.midF <- Lake.title[Type.title=='Mid (F)']

#       Dissimilarity of MidP
ben.midP.lg <- ben.lg[Type.title=='Mid (P)',]
lake.title.midP <- Lake.title[Type.title=='Mid (P)']

#     4.3.2 ANOSIM of community dissimilarity among lakes #####
Janc.dis <- Jac.dis(ben.up.lg)  # BC and Jaccard distance are very similiar in terms of nMDS.
mod1 <- anosim(Janc.dis, grouping = lake.title.up)
mod1

Janc.dis <- Jac.dis(ben.midF.lg)
mod2 <- anosim(Janc.dis, grouping = lake.title.midF)
mod2

Janc.dis <- Jac.dis(ben.midP.lg)
mod3 <- anosim(Janc.dis, grouping = lake.title.midP)
mod3

par(mfrow = c(1,3))
plot(mod1);plot(mod2);plot(mod3)
par(mfrow = c(1,1))

#     4.4 Plot NMDS #####
#plot.nmds <- function()
{
  # par(mfrow=c(1,2))
  # 
  # #P1 - Littoral zone
  # plot(nmds.coord     # use the coordinate take out from nMDS
  #  ,col= col1      # separate the sites use different color based on grp
  #  ,pch=16
  #  ,main=main1
  #  )
  # legend('topleft',c('Upper','Lower F', 'Lower P')
  #        ,col = c('red','blue','green')
  #        ,pch=16)
  # text(x=-0.2,y=-0.3,'Stress = 0.17')
  # text(x=0.25,y=-0.3,'ANOSIM, R = 0.27, p = 0.001')
  # 
  # #P2 - Lake 
  # plot(nmds.coord     # use the coordinate take out from nMDS 
  #      ,col= col2      # separate the sites use different color based on grp
  #      ,pch=16
  #      ,main=main2
  # ) 
  # legend('topleft',c('Lake 1','Lake 2', 'Lake 3')
  #        ,col = c('red','blue','green')
  #        ,pch=16)
  # text(x=-0.2,y=-0.3,'Stress = 0.17')
  # text(x=0.25,y=-0.3,'ANOSIM, R = 0.30, p = 0.001')
  # 
  # par(mfrow=c(1,1)
}
#plot.nmds()
#
#f.nmds.plot <- '221129 nmds.pptx'
#p.nmds.plot <- as.ggplot(~plot.nmds())
#topptx(p.nmds.plot,f.nmds.plot)

#P3 - Lake and Littoral zone
  plot(nmds.coord     # use the coordinate take out from nMDS 
       ,col= col2      # separate the sites use different color based on grp
       ,pch=pch
       ,main=main3
  )
  legend('bottomleft',c('Lake 1','Lake 2','Lake 3'
                     ,'Upper','Lower F', 'Lower P')
         ,col = c('red','blue','green'
                  ,'black','black','black')
         ,pch = c(1,1,1,
                  15,16,17))

# 5. REGRESSION MODEL #####
# Data input
{
  setwd('D:\\PhD\\20211204 Lake ecological study\\Statistical analysis')
  env <- read.csv('Environmental data.csv')
  env.title <- c('Site', 'Lake', 'Site.number', 'Month', 'type', 'water transparency', 'macrophyte breath'
                 ,rep('substrate',4),rep('water qaulity',7)
                 ,rep('canopy',2))
  
  macrophyte <- read.csv('Macrophyte cover & weight.csv')
  macrophyte.title <- c('site','lake','site.number','month','type'
                        ,'wet_weight','dry_weight'
                        ,rep('species',15))
  
  benthos <- read.csv('Site-species data.csv')  # without the Qualitative data
  ben.title <- c('Site','Lake','Site.number','Month','Type','Efforts',rep('Spe',66))
}
# Separate data to Upper, midF and midP
{
  # Water quality data
  water.quality <- (cbind(env[,1:5],env[,env.title=='water qaulity']))[env$type=='upper',]
  # Macrophyte67t                
  macrophyte.up <- (cbind(macrophyte[,1:5],macrophyte[,macrophyte.title=='species']))[macrophyte$Type=='upper',]
  macrophyte.mid <- (cbind(macrophyte[,1:5],macrophyte[,macrophyte.title=='species']))[macrophyte$Type=='middle',]
  # Substrate
  substrate.up <- (cbind(env[,1:5],env[,env.title=='substrate']))[env$type=='upper',]
  substrate.mid <- (cbind(env[,1:5],env[,env.title=='substrate']))[env$type=='middle',]
  # Benthos
  ben.up <- benthos[benthos$Type=='Upper',]
  ben.midF <- benthos[benthos$Type=='Mid (F)',]
  ben.midP <- benthos[benthos$Type=='Mid (P)',]
}
#   5.1 Calculate distance #####
#     5.1.1 ben.up = water.quality + substrate.up + macrophyte.up #####
      # remove empty and standardize
      { # Remove empty row and take out only data
        judge <- ben.up[,ben.title == 'Spe'];judge <- rowSums(judge)
        ben.up <- ben.up[judge>0,]
        ben.up.data <- ben.up[,ben.title == 'Spe']
        #ben.up.data <- ben.up.data[,colSums(ben.up.data)>0]
        
        # Remove empty row, take out only data and standardize
        water.quality.up <- water.quality[judge>0,]
        water.quality.up.data <- water.quality.up[,-c(1:5)]
        water.quality.up.data.scale <- scale(water.quality.up.data,center = TRUE,scale = TRUE)
        
        # Apply correlation methods on water qualities
        #Correlation.plot(water.quality.up.data.scale)
        
        # rm pH.
        water.quality.up.data.scale <- water.quality.up.data.scale[,-4]
        
        # Remove empty row, take out only data and standardize
        substrate.up <- substrate.up[judge>0,]
        substrate.up.data <- substrate.up[,-c(1:5)]
        substrate.up.data.scale <- scale(substrate.up.data,center = TRUE,scale = TRUE)

        # Remove empty row, take out only data
        macrophyte.up <- macrophyte.up[judge>0,]
        macrophyte.up.data <- macrophyte.up[,-c(1:5)]
        macrophyte.up.data <- macrophyte.up.data[,colSums(macrophyte.up.data)>0]
        #macrophyte.up.data <- macrophyte.up.data[,colSums(macrophyte.up.data)>0]
      }    
      #calculate distance
      {
        #com.dist <- Jac.dis(ben.up.data)
        com.dist <- vegdist(ben.up.data,method = 'chao')
        wq.dist <- vegdist(water.quality.up.data.scale, method = 'euclidean')
        subs.dist <- vegdist(substrate.up.data.scale, method = 'euclidean')
        #plant.dist <- Jac.dis(macrophyte.up.data)
        plant.dist <- vegdist(macrophyte.up.data,method = 'bray')
        dis.up <- list(com.dist,wq.dist,subs.dist,plant.dist)
      } 
  #   5.1.2 ben.midF = water.quality + substrate.mid #####
      # remove empty and standardize
      {
        judge <- ben.midF[,ben.title == 'Spe'];judge <- rowSums(judge)
        ben.midF <- ben.midF[judge>0,]
        ben.midF.data <- ben.midF[,ben.title == 'Spe']
        ben.midF.data <- ben.midF.data[,colSums(ben.midF.data)>0]
        
        water.quality.midF <- water.quality[judge>0,]
        water.quality.midF.data <- water.quality.midF[,-c(1:5)]
        water.quality.midF.data.scale <- scale(water.quality.midF.data,center = TRUE,scale = TRUE)
        
        #Correlation.plot(water.quality.midF.data.scale)
        
        # remove pH
        water.quality.midF.data.scale <- water.quality.midF.data.scale[,-4]
        
        substrate.midF <- substrate.mid[judge>0,]
        substrate.midF.data <- substrate.midF[,-c(1:5)]
        substrate.midF.data.scale <- scale(substrate.midF.data,center = TRUE,scale = TRUE)
      }
      #calculate distance
      {
        #com.dist <- Jac.dis(ben.midF.data)
        com.dist <- vegdist(ben.midF.data,method = 'chao')
        wq.dist <- vegdist(water.quality.midF.data.scale, method = 'euclidean')
        subs.dist <- vegdist(substrate.midF.data.scale, method = 'euclidean')
        dis.midF <- list(com.dist,wq.dist,subs.dist)
      }
      
#     5.1.3 ben.midP = water.quality + macrophyte.mid #####
      # remove empty and standardize
      {
        judge <- ben.midP[,ben.title == 'Spe'];judge <- rowSums(judge)
        ben.midP <- ben.midP[judge>0,]
        ben.midP.data <- ben.midP[,ben.title == 'Spe']
        #ben.midP.data <- ben.midP.data[,colSums(ben.midP.data)>0]
        
        water.quality.midP <- water.quality[judge>0,]
        water.quality.midP.data <- water.quality.midP[,-c(1:5)]
        water.quality.midP.data.scale <- scale(water.quality.midP.data,center = TRUE,scale = TRUE)
        
        #Correlation.plot(water.quality.midP.data.scale)
        
        # remove pH
        water.quality.midP.data.scale <- water.quality.midP.data.scale[,-4]
        
        macrophyte.mid <- macrophyte.mid[judge>0,]
        macrophyte.mid.data <- macrophyte.mid[,-c(1:5)]
        macrophyte.mid.data <- macrophyte.mid.data[,colSums(macrophyte.mid.data)>0]
        
      }
      #calculate distance
      {
        #com.dist <- Jac.dis(ben.midP.data)
        com.dist <- vegdist(ben.midP.data,method = 'chao')
        wq.dist <- vegdist(water.quality.midP.data.scale, method = 'euclidean')
        #plant.dist <- Jac.dis(macrophyte.mid.data)
        plant.dist <- vegdist(macrophyte.mid.data,method = 'bray')
        dis.midP <- list(com.dist,wq.dist,plant.dist)
      }
      
#     5.2 Fit Multiple Regression on distance Matrices (MRM) #### 
# https://rdrr.io/cran/ecodist/man/MRM.html
dist <- list(dis.up,dis.midF,dis.midP)
dist.title <- c('ben.up = water.quality + substrate.up + macrophyte.up'
                ,' ben.midF = water.quality + substrate.mid'
                ,'ben.midP = water.quality + macrophyte.mid')
#lm.validation <-function(){
i=0; result.F <- NULL
n <- length(dist.title)
repeat{
  i=i+1
  data.list <- dist[[i]]

  length <- length(data.list)
  response <- data.list[[1]]
  
  par(mfrow=c(1,3))
  
  if(length==3){
    explan1 <- data.list[[2]];explan2 <- data.list[[3]]
    plot(response~explan1,main=dist.title[i]);plot(response~explan2)
    mod <- MRM(response ~ explan1 + explan2, nperm=999)
    lm <- gamlss(response ~ explan1 + explan2,family = 'BEINF')
    }
  if(length==4){
    explan1 <- data.list[[2]];explan2 <- data.list[[3]];explan3 <- data.list[[4]]
    plot(response~explan1,main=dist.title[i]);plot(response~explan2);plot(response~explan3)
    mod <- MRM(response ~ explan1 + explan2 + explan3, nperm=999)
    lm <- gamlss(response ~ explan1 + explan2 + explan3,family = 'BEINF')
  }
  
  par(mfrow=c(1,1))
  
  #plot(lm,main=dist.title[i])
  summary(lm)
  #wp(lm)
  Rsq(lm)
  result <- list(dist.title[i]
                 #,mod
                 ,summary(lm),Rsq(lm))

  if(i==1) result1 <- result
  if(i==2) result2 <- result
  if(i==3) result3 <- result
  
  if(i==n)break
}
result.F <- list(result1,result2,result3)
result.F
#}
# f.lm.validation <- '221204 lm.validation.pptx'
# p.lm.validation <- as.ggplot(~lm.validation())
# topptx(p.lm.validation,f.lm.validation)

#   5.3 Turnover and nested components of sørensen index in each habitats #####
# Sorensen dissimilarity = (b+c)/(2a+b+c), which would result lower value
# compare to Jaccard index, which are (b+c)/(a+b+c).
data.list <- list(ben.up.data,ben.midF.data,ben.midP.data)
title.list <- c('upper','midF','midP')
i=0; result.F <- NULL
par(mfrow=c(1,3))
repeat{
  i=i+1
  data <- data.list[[i]]
  data.binary <- data;data.binary[data.binary>0]=1
  
  jaccard.ab <- vegdist(data,method = 'jaccard')
  jaccard <- beta.pair(data.binary,index.family = "jaccard")
  total.dis <- jaccard$beta.jac
  turnover <- jaccard$beta.jtu
  nested <- jaccard$beta.jne
  
  t.ab.inci <-t.test(as.vector(jaccard.ab),as.vector(total.dis))
  t.ab.inci
  t.turn.nest <- t.test(as.vector(turnover),as.vector(nested))
  t.turn.nest
  # boxplot
  list <- list(jaccard.ab,total.dis,turnover,nested)
  grp <- c('Adun-Jaccard','Jaccard','Turnover','Nested')
  Boxplot.grp(list,grp,main = title.list[i]
              ,col = 'white',lty = 1)
  text(x=1,y=1,t.ab.inci$statistic)
  text(x=2,y=0.9,t.ab.inci$p.value)
  
  text(x=3,y=0.8,t.turn.nest$statistic)
  text(x=4,y=0.7,t.turn.nest$p.value)
  dissimilarity <- cbind(jaccard.ab,total.dis,turnover,nested)
  result <- Mean.SE(dissimilarity)
    
  # sorensen <- beta.multi(data.binary,index.family = "sorensen")
  # total.dis <- round(sorensen$beta.SOR,2)
  # turnover <- round(sorensen$beta.SIM,2)
  # nested <- round(sorensen$beta.SNE,2)
  # result <- c(title.list[i],total.dis,turnover,nested)
  # result.F <- rbind(result.F,result)
  if(i==3)break
}
par(mfrow=c(1,1))
#colnames(result.F) <- c('habitat','Total','turnover','nested')
#result.F
result

# 6. CCA #####
# Data input and remove empty
{
  setwd('D:\\PhD\\20211204 Lake ecological study\\Statistical analysis')
  
  benthos <- read.csv('Site-species data.csv')  # without the Qualitative data
  ben.title <- c('Site','Lake','Site.number','Month','Type','Efforts',rep('Spe',66))
  
  env <- read.csv('Environmental data.csv')
  env.title <- c('Site', 'Lake', 'Site.number', 'Month', 'type', 'water transparency', 'macrophyte breath'
                 ,rep('substrate',4),rep('water qaulity',7)
                 ,rep('canopy',2))
  
  macrophyte <- read.csv('Macrophyte cover & weight.csv')
  macrophyte.title <- c('site','lake','site.number','month','type'
                        ,'wet_weight','dry_weight'
                        ,rep('species',15))
}

{
  # Macrophyte                
  macrophyte.up <- (cbind(macrophyte[,macrophyte.title=='species']))[macrophyte$Type=='upper',]
  macrophyte.binary <- macrophyte.up; macrophyte.binary[macrophyte.binary>0]=1
  macrophyte.up <- macrophyte.up[,colSums(macrophyte.binary)>1]

  macrophyte.mid <- (cbind(macrophyte[,macrophyte.title=='species']))[macrophyte$Type=='middle',]
  macrophyte.binary <- macrophyte.mid; macrophyte.binary[macrophyte.binary>0]=1
  macrophyte.mid <- macrophyte.mid[,colSums(macrophyte.binary)>1]
  
  
  
  # Water quality data
  w.trans <- (env[,env.title== 'water transparency'])[env$type=='upper']
  m.breadth <- (env[,env.title== 'macrophyte breath'])[env$type=='upper']
  canopy.cover <- env$Modified.Caco
  water.quality <- (env[,env.title== 'water qaulity'])[env$type=='upper',]
  
  macrophyte.cover.up <- rowSums(macrophyte[,-c(1:7)])[env$type=='upper']
  macrophyte.cover.mid<- rowSums(macrophyte[,-c(1:7)])[env$type=='middle']
  
  water.quality.up <- cbind(w.trans, m.breadth,macrophyte.cover.up,  water.quality);water.quality.up.scale <- scale(water.quality.up,center=TRUE,scale=TRUE)
  water.quality.mid <- cbind(w.trans          ,macrophyte.cover.mid, water.quality);water.quality.mid.scale <- scale(water.quality.mid,center=TRUE,scale=TRUE)

  # Substrate
  substrate.up <- (env[,env.title=='substrate'])[env$type=='upper',]; substrate.up.scale <- scale(substrate.up,center=TRUE,scale=TRUE)
  substrate.mid <- (env[,env.title=='substrate'])[env$type=='middle',]; substrate.mid.scale <- scale(substrate.mid,center=TRUE,scale=TRUE)
  
  # Benthos
  ben.up <- benthos[benthos$Type=='Upper',];ben.up.data <- ben.up[,-c(1:6)]
  ben.midF <- benthos[benthos$Type=='Mid (F)',];ben.midF.data <- ben.midF[,-c(1:6)]
  ben.midP <- benthos[benthos$Type=='Mid (P)',];ben.midP.data <- ben.midP[,-c(1:6)]
}

ben.list <- list(ben.up,ben.midF,ben.midP)
ben.data.list <- list(ben.up.data,ben.midF.data,ben.midP.data)
explanatory.list <- list(
  list(water.quality.up.scale,substrate.up.scale,macrophyte.up)
  ,list(water.quality.mid.scale,substrate.mid)
  ,list(water.quality.mid.scale,macrophyte.mid)
)
title <- c('Upper','Mid F', 'Mid P')
# Columns that should be removed
  rm1 <- c(5,9,11,21,22,23,25) 
  # alias:'gravel.stone'[11], 'sp15'[25]
  # 'NH3'[5], 'sp11'[22], 'NaNO2'[7],'sp10'[21]
  # 'wooden debris'[9],'sp12'[23]
  
  rm2 <- c(8)  # 'silt and mud'
  rm3 <- c(5,13) # 'sp11','NH3-N'
  rm <- list(rm1,rm2,rm3)

#   6.1 CCA calculation #####

i=0;result.F <- NULL
par(mfrow=c(2,3))
repeat{
  i=i+1
  ben.data <- ben.data.list[[i]]
  exp <- explanatory.list[[i]]
  n <- length(exp)
  if(n==2) explanatory.original <- cbind(exp[[1]],exp[[2]])
  if(n==3) explanatory.original <- cbind(exp[[1]],exp[[2]],exp[[3]])
  
  ben.lg <- log(ben.data+1)
  # macrophyte.up.log <- log(macrophyte.up.data+1)
  decorana(ben.lg) # DCA analysis (<3 use RDA; 4 use CCA)
  
#   6.2 forward removal of alias values step by step #####
  # This step was conducted before the final plot
  
  colnames(explanatory.original)
  # cca analysis
  explanatory <- explanatory.original
  explanatory <- as.data.frame(explanatory)
  spe.cca <- cca(ben.lg~.,data=explanatory) # ben.lg with all the column in explanatory

  alias(spe.cca) # alias variables. i.e., empty values, values with less observations.
  vif.cca(spe.cca) # co-linear coeffecience (the small the better, need < 10, other wise remove some variables)
  explanatory <- explanatory.original
  
  # remove variables wit high collinearity
  explanatory <- explanatory.original
  explanatory <- explanatory[,-(rm[[i]])]
  
#   6.3 CCA analysis ######
  explanatory <- as.data.frame(explanatory)
  spe.cca <- cca(ben.lg~.,data=explanatory) # ben.lg with all the column in explanatory
  #summary(spe.cca) # summary of RDA analysis
  vif.cca(spe.cca) # co-linear coeffecience (the small the better, need < 10, other wise remove some varaibles)
  
#   6.4 CCA plot #####
    #       6.4.1 default plot ####
    pl <- plot(spe.cca,scaling=3,main=title[i]) 
    
    #       6.4.2 plot from a blank ####
    pl <- plot(spe.cca, type="none", scaling=3, correlation=TRUE)
    # 'arg'should be one of “text”, “points”, “none”其中的一个
    
    #           6.4.2.1 plot sites####
    lake.type <- (ben.list[[i]])$Lake
    pch <- Replace(lake.type,c(1,2,3))[[2]]
    points(pl, "site", pch=pch
           ,lwd=1                       # line width
           ,col=1,cex=1.2)
    
    #           6.4.2.2 plot species ####
    # Intrude the family grp to saperate the species by grp (ben.family.grp)
    species.name <-rownames(pl$species)
    ben.family.grp <- substring(species.name,1,3)# taking the first four character of speceis name
    ben.family.grp <- as.factor(ben.family.grp)  # translate string into factors
    
    points(pl, "species"      # species in pl  
           ,pch=21,cex=0.8
           , col=NULL,bg=ben.family.grp 
           # background based on 'ben.family.grp'
           # and no outline color because 'col=NULL'
    )
    
    #           6.4.2.3 plot variables ####
    text(spe.cca, display ="bp", scaling=3,col='lightblue') # plot arrows and variables
    # 'arg' should be one of “sites”, “species”, “wa”, “lc”, “bp”, “reg”, “cn”
    # bp=variables and arrow
    
    #           6.4.2.4 plot legends ####
    legend("bottomright"
           , levels(ben.family.grp)
           , pch=21, pt.bg=1:10
           , border = NULL
           , bty="n")             # legend has no border
    legend("topright"
           , Element.grp(lake.type)
           , pch=1:3,lwd=1
           , border = NULL
           , bty="n"              # legend has no border
           ) 


#   6.5 CCA statistical analysis #####
  sig.var <- anova.cca(spe.cca, by="term", step=1000)   # The sig. of every each variable
  sig.mod <- anova.cca(spe.cca, step=1000)            # Permutation test of result of CCA/RDA
  
  library(rdacca.hp)
  
  result <- list(title[i],sig.mod,sig.var)
  result.F <- list(result.F,result)
  result
  if(i==length(explanatory.list))break
}
par(mfrow=c(1,1))
result.F


# 7. Mid P: MACROPHYTE and BENTHOS #####
{
  setwd('D:\\PhD\\20211204 Lake ecological study\\Statistical analysis')
  env <- read.csv('Environmental data.csv')
  env.title <- c('Site', 'Lake', 'Site.number', 'Month', 'type', 'water transparency', 'macrophyte breath'
                 ,rep('substrate',4),rep('water qaulity',7)
                 ,rep('canopy',2))
  
  macrophyte <- read.csv('Macrophyte cover & weight.csv')
  macrophyte.title <- c('site','lake','site.number','month','type'
                        ,'wet_weight','dry_weight'
                        ,rep('species',15))
  
  benthos <- read.csv('Site-species data.csv')  # without the Qualitative data
  ben.title <- c('Site','Lake','Site.number','Month','Type','Efforts',rep('Spe',66))
}
#   7.1 indi/g and indi/cover% ##### 
# ben
ben.midP <- benthos[benthos$Type == 'Mid (P)',]
judge <- ben.midP$Efforts.m2

ben.up <- benthos[benthos$Type == 'Upper',]
# rm empty
ben.rm <- ben.midP[judge>0,]
ben.data <- ben.rm[,ben.title == 'Spe']

# 7.1.1 benthos density and SR
abundance.midP <- rowSums(ben.data)

ben.data.binary <- ben.data
ben.data.binary[ben.data.binary>0]=1
SR.midP <- rowSums(ben.data.binary)

# 7.1.2 macrophyte cover, weight
cover <- rowSums(macrophyte[,macrophyte.title=='species'])

cover.midP <- cover[macrophyte$Type=='middle'];cover.midP <- cover.midP[judge>0]
wet.weight <- (macrophyte$Wet_weight)[macrophyte$Type=='middle']; wet.weight<- wet.weight[judge>0]
dry.weight <- (macrophyte$Dry_weight)[macrophyte$Type=='middle']; dry.weight <- dry.weight[judge>0]

wet.w.cover <- wet.weight/cover.midP # g/cover%
dry.w.cover <- dry.weight/cover.midP # g/cover%

# # 7.1.3 Macroinvertebrate abundence and macrophyte cover
# plot(abundance.midP[abundance.midP<7000]~cover.midP[abundance.midP<7000])
# lm <- lm(abundance.midP[abundance.midP<7000]~cover.midP[abundance.midP<7000])
# summary(lm)
# abline(lm,col='red')
# par(mfrow=c(2,2));plot(lm);par(mfrow=c(1,1))
# cor.test(abundance.midP[abundance.midP<7000], cover.midP[abundance.midP<7000], method='spearman')

# 7.1.4 Macroinvertebrate density: how many individuals per macrophyte weight: indi/g
den.wet.w <- abundance.midP/wet.weight
den.dry.w <- abundance.midP/dry.weight
den.cover <- abundance.midP/(cover.midP*0.24)  # indi/m2

grp <- ben.rm$Lake
par(mfrow=c(1,3))
boxplot(den.wet.w~grp);boxplot(den.dry.w~grp);boxplot(den.cover~grp)
boxplot(wet.w.cover~grp);boxplot(dry.w.cover~grp)
par(mfrow=c(1,1))

Wilcox.test.grp(den.wet.w,grp)# All did not show significance p > 0.06
Wilcox.test.grp(den.dry.w,grp)# All did not show significance p > 0.08
Wilcox.test.grp(den.cover,grp)# All did not show significance p > 0.08

Wilcox.test.grp(wet.w.cover,grp) # only sig L1-L3: w=0, p=0.044444. Rest p>0.26
Wilcox.test.grp(dry.w.cover,grp) # All did not show significance p > 0.17

# 7.1.5 macroinvertebrate density and macrophyte in both upper and midP

#   Upper
ben.up.data <- ben.up[,ben.title == 'Spe']
abundance.up <- rowSums(ben.up.data)
den.up <- abundance.up/ben.up$Efforts.m2.or.g.
cover.up <- cover[macrophyte$Type=='upper']
plot(den.up[den.up<5000] ~ cover.up[den.up<5000],xlim=c(0,100))
#   midP
plot(den.dry.w ~ cover.midP)

lm.up <- gamlss(den.up[den.up<5000] ~ cover.up[den.up<5000])
lm.midP <- gamlss(den.dry.w ~ cover.midP)
summary(lm.up);summary(lm.midP)
par(mfrow=c(2,2))
plot(lm.up)
plot(lm.midP)
par(mfrow=c(1,1))

plot(den.up[den.up<5000] ~ cover.up[den.up<5000],xlim=c(0,100))
abline(lm.up,col='red')

plot(den.dry.w ~ cover.midP)
abline(lm.midP,col='red')

cor.test(den.up[den.up<5000], cover.up[den.up<5000],method='spearman')
cor.test(den.dry.w, cover.midP,method='spearman')
# 8. Macrophyte #####
{
  setwd('D:\\PhD\\20211204 Lake ecological study\\Statistical analysis')
  env <- read.csv('Environmental data.csv')
  env.title <- c('Site', 'Lake', 'Site.number', 'Month', 'type', 'water transparency', 'macrophyte breath'
                 ,rep('substrate',4),rep('water qaulity',7)
                 ,rep('canopy',2))
  
  macrophyte <- read.csv('Macrophyte cover & weight.csv')
  macrophyte.title <- c('site','lake','site.number','month','type'
                        ,'wet_weight','dry_weight'
                        ,rep('species',15))
  
}
#   8.1 Canopy and macrophyte cover #####
env.up <- env[env$type == 'upper',]
canopy <- env.up$Modified.Caco
breadth <- env.up$macrophyte.breadth

macrophyte.up <- macrophyte[macrophyte$Type == 'upper',]
macrophyte.data <- macrophyte.up[,macrophyte.title == 'species']
cover <- rowSums(macrophyte.data)
macrophyte.data[macrophyte.data>0]=1
SR <- rowSums(macrophyte.data)

judge <- env.up$Lake

plot(breadth[judge != 'L1']~canopy[judge != 'L1']) # breadth are all too small
plot(SR[judge != 'L1']~canopy[judge != 'L1'])      # clearly canopy can't affect SR
plot(cover~canopy)   # Canopy affect macrophyte cover

# fit linear model
lm <- gamlss(cover/100~canopy) # to make the distribution of cover within (0,1)
lm1 <- gamlss(cover/100~canopy,family = 'BE')
lm2 <- gamlss(cover/100~canopy,family = 'GB1')
lm3 <- gamlss(cover/100~canopy,family = 'LOGITNO')
AIC(lm,lm1,lm2,lm3)

summary(lm1);summary(lm1);summary(lm2)

summary(lm3);Rsq(lm3)

hist(cover)

plot(lm3);wp(lm3)

par(mfrow=c(1,1))
plot(cover/100~canopy,xlim=c(0,100),ylim=c(0,1))
abline(lm,col='red')
abline(lm1,col='blue')
abline(lm2,col='yellow')
abline(lm3,col='green')

# Correlation 
cor.test(cover,canopy
         ,method='spearman') # Method 'pearson' need data to follow normal distribution
