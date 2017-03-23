##################################################################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##################################################################################################
# Start
#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
# Clustering 450K methylation CpGs probes written by Dr Reza Rafiee
# Research Associate, Northern Institute for Cancer Research, Newcastle University
# This script loads 20,000 methylation probes (from 450K methylation profiling) and doing clustering analysis
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 

library(mclust) # Gaussian Mixture Modelling package for Model-Based Clustering, Classification, and Density Estimation
library(scatterplot3d)
library(pheatmap)
library(apcluster) # Affinity Propagation Clustering

load("~/20KBetaValues_51InfantSHH.RData")    # 20,000 probes
length(colnames(BetaValues_51Samples_20K))

# Performs a principal components analysis on the given data matrix and returns the results as an object of class prcomp
PCA_Comp_Scaled_Centered <- prcomp(t(BetaValues_51Samples_20K), center = TRUE, scale=T)  # scale =T is better
summary(PCA_Comp_Scaled_Centered)

# Importance of components:
#                         PC1     PC2      PC3      PC4      PC5      PC6      PC7     PC8      PC9     PC10     PC11     PC12     PC13     PC14     PC15     PC16
# Standard deviation     56.400 45.3408 35.12111 32.24443 28.50901 27.18614 24.98240 24.3309 22.87971 21.83193 21.25324 20.75756 20.46585 18.90885 18.24397 17.89305
# Proportion of Variance  0.159  0.1028  0.06167  0.05199  0.04064  0.03695  0.03121  0.0296  0.02617  0.02383  0.02259  0.02154  0.02094  0.01788  0.01664  0.01601
# Cumulative Proportion   0.159  0.2618  0.32351  0.37550  0.41614  0.45309  0.48430  0.5139  0.54007  0.56390  0.58649  0.60803  0.62897  0.64685  0.66349  0.67950
#                         PC17     PC18     PC19     PC20     PC21     PC22     PC23     PC24     PC25    PC26     PC27     PC28     PC29     PC30     PC31
# Standard deviation     17.63880 17.53269 17.04100 16.66391 16.16575 15.92072 15.55023 15.35452 15.11828 14.7654 14.57411 14.51679 14.25707 14.21835 14.03399
# Proportion of Variance  0.01556  0.01537  0.01452  0.01388  0.01307  0.01267  0.01209  0.01179  0.01143  0.0109  0.01062  0.01054  0.01016  0.01011  0.00985
# Cumulative Proportion   0.69506  0.71043  0.72495  0.73883  0.75190  0.76457  0.77666  0.78845  0.79988  0.8108  0.82140  0.83194  0.84210  0.85221  0.86205
#                         PC32     PC33     PC34    PC35     PC36    PC37     PC38     PC39     PC40     PC41     PC42    PC43     PC44     PC45     PC46    PC47
# Standard deviation     13.89703 13.62558 13.54381 13.2693 13.12656 12.8847 12.78211 12.46655 12.22380 11.93376 11.72942 11.5741 11.41000 11.26610 10.89133 10.6728
# Proportion of Variance  0.00966  0.00928  0.00917  0.0088  0.00862  0.0083  0.00817  0.00777  0.00747  0.00712  0.00688  0.0067  0.00651  0.00635  0.00593  0.0057
# Cumulative Proportion   0.87171  0.88099  0.89017  0.8990  0.90758  0.9159  0.92405  0.93182  0.93930  0.94642  0.95330  0.9600  0.96650  0.97285  0.97878  0.9845
#                         PC48     PC49    PC50      PC51
# Standard deviation     10.44080 10.25791 9.81097 6.577e-14
# Proportion of Variance  0.00545  0.00526 0.00481 0.000e+00
# Cumulative Proportion   0.98993  0.99519 1.00000 1.000e+00
# 

# How many eigenvectors/components we need to use in this analysis?
# Taking the first k eigenvectors that capture at least 99% of the total variance

# creating a graphic visualization of the relationship between eigenvalues and number of factors

par(mfrow=c(1,1))
par(mar=c(5,4,4,5) + 0.1)
par(cex.axis=0.8)

plot(PCA_Comp_Scaled_Centered,main = "Variances vs. number of components", type = "l")

# biplot(prcomp(PCA_Comp_Scaled_Centered$x, scale = TRUE))
pairs(PCA_Comp_Scaled_Centered$x[,1:5], pch=19, col="black",log = "")

GMM_object_PCA <- Mclust(as.matrix(PCA_Comp_Scaled_Centered$x[,1:5]), G=1:4)
Best_Num_of_Clusters <- dim(GMM_object_PCA$z)[2]
cat("Model-based optimal number of clusters:", Best_Num_of_Clusters, "\n")
# model-based optimal number of clusters: 2 clusters


# Resampling-based Inference for Gaussian finite mixture models
# Bootstrap or jackknife estimation of standard errors and percentile bootstrap confidence intervals
# for the parameters of a Gaussian mixture model.
bootClust <- MclustBootstrap(GMM_object_PCA)
bootClust$modelName #[1] "VII"
summary(bootClust, what = "se")
summary(bootClust, what = "ci")

# plot(GMM_object_PCA)

plot(GMM_object_PCA, what = c("BIC", "classification", "uncertainty", "density"),
     dimens = NULL, xlab = NULL, ylab = NULL, ylim = NULL, addEllipses = TRUE, main = TRUE)

#sort(GMM_object_PCA$z[,1])
#sort(GMM_object_PCA$z[,2])

Probability_Assignment <- GMM_object_PCA$z
#
#         [,1]         [,2]
# NMB113 1.161098e-18 1.0000000000
# NMB138 1.087764e-13 1.0000000000
# NMB143 3.095577e-15 1.0000000000
# NMB200 9.673815e-01 0.0326185108
# ...

bestmodelBIC <- mclustBIC(GMM_object_PCA$data)
# Top 3 models based on the BIC criterion:
#   VII,2     VII,5     EII,5 
# -2562.430 -2565.113 -2569.015 


length(which(GMM_object_PCA$classification == 1)) # red colour, n=18
length(which(GMM_object_PCA$classification == 2)) # blue colour, n=33

which(GMM_object_PCA$classification == 1)
which(GMM_object_PCA$classification == 2)

# Group 1
# which(GMM_object_PCA$classification == 1)
# NMB200 NMB254 NMB272  NMB32 NMB324 NMB328 NMB363 NMB465 NMB471 NMB483 NMB497 NMB553 NMB594 NMB608 NMB621  NMB64 NMB676 NMB712 
# 4      5      6      7      8      9     10     15     17     22     27     29     32     33     35     36     43     45 

# Group 2
# which(GMM_object_PCA$classification == 2)
# NMB113 NMB138 NMB143 NMB364 NMB371 NMB379 NMB439 NMB466 NMB474 NMB477 NMB479 NMB482 NMB485 NMB486 NMB495 NMB496 NMB498 NMB554 NMB580 NMB612 NMB651 NMB667 NMB670 
# 1      2      3     11     12     13     14     16     18     19     20     21     23     24     25     26     28     30     31     34     37     38     39 
# NMB673 NMB674 NMB675 NMB690 NMB720 NMB726  NMB79 NMB798 NMB803 NMB873 
# 40     41     42     44     46     47     48     49     50     51 

#GMM_object_PCA$data[,1:3]

grpCols <- ifelse(GMM_object_PCA$classification == 1, "blue", "red")

# Add regression plane
shapes = c(16, 17, 18) 
my.lm <- lm(GMM_object_PCA$data[,1] ~ GMM_object_PCA$data[,2] + GMM_object_PCA$data[,3])
s3d <- scatterplot3d(GMM_object_PCA$data[,1:3], pch = 17, type = "h", color = grpCols,angle=-225, scale.y=0.9, col.grid="black", grid=TRUE) # angle=-225, 280, 55, 75
s3d$plane3d(my.lm)
text(s3d$xyz.convert(GMM_object_PCA$data[,1:3]), labels = rownames(PCA_Comp_Scaled_Centered$x), cex= 0.8, col = "black")

# Visualising heatmap of PCA components
pheatmap(GMM_object_PCA$data[,1:5],color = colorRampPalette(c("navy", "white", "firebrick3"))(100),clustering_method = "ward.D2")

write.csv(Probability_Assignment, file="~/Probabilities_Assignment.csv")


################################################## AP clustering ###############################
# Affinity Propagation clustering is used to confirm an aggreement/consensus with GMM+EM results

# Affinity propagation (AP) is a relatively new clustering algorithm that has been introduced by
# Brendan J. Frey and Delbert Dueck. The authors themselves describe affinity propagation as
# follows:
# “An algorithm that identifies exemplars among data points and forms clusters of data
# points around these exemplars. It operates by simultaneously considering all data
# point as potential exemplars and exchanging messages between data points until a
# good set of exemplars and clusters emerges.”

AP_object_PCA <- apcluster(negDistMat(r=2), PCA_Comp_Scaled_Centered$x[,1:5], q=0.0)
cat("affinity propogation optimal number of clusters:", length(AP_object_PCA@clusters), "\n")
plot(AP_object_PCA,GMM_object_PCA$data[,1:5])
heatmap(AP_object_PCA)
show(AP_object_PCA)
AP_object_PCA

# APResult object
# 
# Number of samples     =  51 
# Number of iterations  =  129 
# Input preference      =  -73914.21 
# Sum of similarities   =  -234507 
# Sum of preferences    =  -221742.6 
# Net similarity        =  -456249.6 
# Number of clusters    =  3 

# Exemplars:
#   NMB371 NMB498 NMB553
# Clusters:
#   Cluster 1, exemplar NMB371:
#   NMB143 NMB371 NMB439 NMB485 NMB486 NMB554 NMB612 NMB670 NMB673 NMB674 NMB675
# Cluster 2, exemplar NMB498:
#   NMB113 NMB138 NMB364 NMB379 NMB466 NMB474 NMB477 NMB479 NMB482 NMB495 NMB496 NMB498 NMB580 NMB651 NMB667 NMB690 NMB720 NMB726 NMB79 NMB798 NMB803 NMB873
# Cluster 3, exemplar NMB553:
#   NMB200 NMB254 NMB272 NMB32 NMB324 NMB328 NMB363 NMB465 NMB471 NMB483 NMB497 NMB553 NMB594 NMB608 NMB621 NMB64 NMB676 NMB712

length(AP_object_PCA@clusters[[1]]) # n1=11 (red colours)
length(AP_object_PCA@clusters[[2]]) # n2=22 (green colours)
# n1+n2 = 33 (Group 2)
length(AP_object_PCA@clusters[[3]]) # n3=18 (Group 1; blue colours)

# We could use all 20k probes in this clustering 
# AP_object_20k <- apcluster(negDistMat(r=2), t(BetaValues_51Samples_20K), q=0.01)
# cat("affinity propogation optimal number of clusters:", length(AP_object_20k@clusters), "\n")
# affinity propogation optimal number of clusters: 2

#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
# End
##################################################################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##################################################################################################
