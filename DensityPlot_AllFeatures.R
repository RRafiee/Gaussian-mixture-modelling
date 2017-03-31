##################################################################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##################################################################################################
# Start
#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
# Density plot script written by Dr Reza Rafiee
# Research Associate, Northern Institute for Cancer Research, Newcastle University
# This script gets a csv file including age at diagnosis and other features for plotting density plots 

linear.density.plot <- function(present.score,not.present.score, condition, max.age)
  
{
  #max(c(present.score,not.present.score)) -> max.age
  score <- c(present.score,not.present.score)
  cat <- rep('a',length(score))
  cat[score<5] <- 'less.five'
  #cat[score<5&score>=3] <- 'three.to.five'
  cat[score<=max.age&score>=5] <- 'more.five'
  cat <- as.factor(cat)
  cat <- factor(cat,levels(cat)[c(1,3,2)])
  pres.fact <- factor(c(rep("Present",length(present.score)),rep("Absent",length(not.present.score))))
  counts <- table(pres.fact,cat)
  plot(density(present.score)$x,density(present.score)$y*14, col = "blue", type = "l",   #density(present.score)$y*13 is original, 15
       ylim=c(-0.2,max(density(present.score)$y)*14),  # original 
       #ylim=c(-0.2,4.00),   # maxium Y axis (3.59) is for the Chr17q.gain
       xlim=c(0,max.age), xlab = "Age (Years)", 
       bty="n",
       xaxt='n', ann=FALSE, yaxt='n'
  )
  ##rect(xleft=0, xright=5, ybottom=0, ytop=10, density=50, col="grey", border = NA, angle = 135)
  ##rect(xleft=5, xright=max.age, ybottom=0, ytop=10, density=10, col="grey", border = NA, angle = 135)
  lines(density(present.score)$x,density(present.score)$y*14, col = "blue")  #15
  axis(1,at=c(0,max.age),pos=0, tck=c(-0.1), labels = NA,lwd=1)
  axis(1,at=c(0,max.age),pos=0, tck=c(0.1), labels = NA,lwd=1)
  #lines(density(present.score), col = "blue")
  rug(present.score, col = "blue", side = 1, ticksize=0.2,lwd = 2)
  rug(not.present.score, col = "grey", side = 1, ticksize=0.2,lwd =0.5)
  abline(v=5, lty = 2)
  mtext(paste("",condition), side = 2, las = 1, cex = 0.75)
  
}

setwd("~/My Projects at NICR/2014/Infant/Our New Material")

#### density plots

# 31 March 2017
data_Infant16 <- read.csv("New under 5 vs over 5 NMB DH210716 for Reza density plots DH310317.csv", header=T)  # 
data <- data_Infant16

#### data stored in a list for each plot on the combined linear density plot

density.plots <- list()

###########

score <- data$Ageatdiagnosis
cat <- rep('a',length(score))

cat[score<5] <- 'less.five'
cat[score>=5] <- 'more.five'
cat <- as.factor(cat)
cat <- factor(cat,levels(cat)[c(1,3,2)])

present.score <- data_Infant16$Ageatdiagnosis[which(data_Infant16$Male=="1")]
not.present.score <- data_Infant16$Ageatdiagnosis[which(data_Infant16$Male=="2")]
density.plots[[1]]<-list(present.score,not.present.score)

present.score <- data_Infant16$Ageatdiagnosis[which(data_Infant16$STR=="1")]
not.present.score <- data_Infant16$Ageatdiagnosis[which(data_Infant16$STR=="2")]
density.plots[[2]]<-list(present.score,not.present.score)

present.score <- data_Infant16$Ageatdiagnosis[which(data_Infant16$M2.=="1")]  
not.present.score <- data_Infant16$Ageatdiagnosis[which(data_Infant16$M2.=="2")]
density.plots[[3]]<-list(present.score,not.present.score)

present.score <- data$Ageatdiagnosis[which(data_Infant16$CLA=="1")]
not.present.score <- data$Ageatdiagnosis[which(data_Infant16$CLA=="2")]
density.plots[[4]] <- list(present.score,not.present.score)

present.score <- data_Infant16$Ageatdiagnosis[which(data_Infant16$DN.MBEN=="1")]
not.present.score <- data_Infant16$Ageatdiagnosis[which(data_Infant16$DN.MBEN=="2")]
density.plots[[5]] <- list(present.score,not.present.score)

present.score <- data_Infant16$Ageatdiagnosis[which(data_Infant16$LCA=="1")]
not.present.score <- data_Infant16$Ageatdiagnosis[which(data_Infant16$LCA=="2")] 
density.plots[[6]] <- list(present.score,not.present.score)

present.score <- data_Infant16$Ageatdiagnosis[which(data_Infant16$WNT=="1")]
not.present.score <- data_Infant16$Ageatdiagnosis[which(data_Infant16$WNT=="2")]
density.plots[[7]] <- list(present.score,not.present.score)

present.score <- data_Infant16$Ageatdiagnosis[which(data_Infant16$SHH=="1")]
not.present.score <- data_Infant16$Ageatdiagnosis[which(data_Infant16$SHH=="2")]
density.plots[[8]] <- list(present.score,not.present.score)

present.score <- data_Infant16$Ageatdiagnosis[which(data_Infant16$G3=="1")]
not.present.score <- data_Infant16$Ageatdiagnosis[which(data_Infant16$G3=="2")]
density.plots[[9]] <- list(present.score,not.present.score)

present.score <- data_Infant16$Ageatdiagnosis[which(data_Infant16$G4=="1")]
not.present.score <- data_Infant16$Ageatdiagnosis[which(data_Infant16$G4=="2")]
density.plots[[10]] <- list(present.score,not.present.score)


present.score <- data_Infant16$Ageatdiagnosis[which(data_Infant16$MYC.amplification=="1")]
not.present.score <- data_Infant16$Ageatdiagnosis[which(data_Infant16$MYC.amplification=="2")]
density.plots[[11]] <- list(present.score,not.present.score)

present.score <- data_Infant16$Ageatdiagnosis[which(data_Infant16$MYCN.amplification=="1")]
not.present.score <- data_Infant16$Ageatdiagnosis[which(data_Infant16$MYCN.amplification=="2")]
density.plots[[12]] <- list(present.score,not.present.score)

present.score <- data_Infant16$Ageatdiagnosis[which(data_Infant16$TP53.mutation=="1")]
not.present.score <- data_Infant16$Ageatdiagnosis[which(data_Infant16$TP53.mutation=="2")]
density.plots[[13]]<-list(present.score,not.present.score)


present.score <- data_Infant16$Ageatdiagnosis[which(data_Infant16$i17q=="1")]
not.present.score <- data_Infant16$Ageatdiagnosis[which(data_Infant16$i17q=="2")]
density.plots[[14]]<-list(present.score,not.present.score)

# present.score <- data_Infant16$Ageatdiagnosis[which(data_Infant16$Under.5=="1")]  # Under 5
# not.present.score <- data_Infant16$Ageatdiagnosis[which(data_Infant16$Under.5=="2")]  # more than 5 years
# density.plots[[2]] <- list(present.score,not.present.score)



names(density.plots)<-c("Male", #1
                        "STR", #2
                        "M2+", #3 max(density(present.score)$y)*13 ===> 0.9822429
                        #"DN/MBEN", #4 max(density(present.score)$y)*13 ===> 0.9822429
                        "CLA", #4
                        "DN/MBEN", #5
                        "LCA", #6
                        "WNT", #7
                        "SHH", #8
                        "Grp3", #9
                        "Grp4", #10
                        "MYC amplification", #11
                        "MYCN amplification", #12
                        #"Male", #13 - max(density(present.score)$y)*13 ===> 1.750508
                        "TP53 Mutation", #13
                        "iso17q")    #14

names(density.plots)<- strrep(" ",1:14)    #14
##save(density.plots, file="~/My Projects at NICR/2014/Infant/Our New Material/density.plots_List20Var_14Jan16.rda")
##save(density.plots, file="~/My Projects at NICR/2014/Infant/Our New Material/density.plots_List16Var_08Jan16.rda")

###pdf(file="linear.density.plot.pdf")
#load("~/My Projects at NICR/2014/Infant/Our New Material/density.plots_List16Var_08Jan16.rda")
###data <- read.csv(file = "Infant data current 0-16 DH050116_Ch17.csv")
#load("~/My Projects at NICR/2014/Infant/Our New Material/density.plots_List20Var_14Jan16.rda")
#data <- read.csv(file = "FOR DENSITY PLOT FIG 1 ONLY_DATA INFIDELITY Copy of Infant data current 0-16 DH140116.csv")


max.age <- max(data_Infant16$Ageatdiagnosis, na.rm=TRUE)
#max.age <- max(data$Ageatdiagnosis)
number.of.graphs <- 14
#number.of.graphs <- 1

#pdf(file="linear.density.plot14Jan2016_Final.pdf")

par(mfrow=c(number.of.graphs,1),mar=c(0.5,10,0.5,5)) #original: mar=c(0.5,8,0.5,5)


for(i in 1:14)
{ 
  present.score <- as.numeric(density.plots[[i]][[1]])
  present.score <- present.score[!is.na(present.score)]
  
  not.present.score <- as.numeric(density.plots[[i]][[2]]) 
  not.present.score <- not.present.score[!is.na(not.present.score)]
  
  linear.density.plot(present.score ,not.present.score ,names(density.plots)[i],max.age)
  
}


#dev.off()

#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
# End
##################################################################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##################################################################################################

