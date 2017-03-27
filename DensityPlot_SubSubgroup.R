##################################################################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##################################################################################################
# Start
#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
# Density plot script written by Dr Reza Rafiee
# Research Associate, Northern Institute for Cancer Research, Newcastle University
# This script gets a csv file including age at diagnosis and subsubgroup for density plot 

trapezoidal_integration = function(x, f)
{
  ### 3 checks to ensure that the arguments are numeric and of equal lengths
  # check if the variable of integration is numeric
  if (!is.numeric(x))
  {
    stop('The variable of integration "x" is not numeric.')
  }
  # check if the integrand is numeric
  if (!is.numeric(f))
  {
    stop('The integrand "f" is not numeric.')
  }
  # check if the variable of integration and the integrand have equal lengths
  if (length(x) != length(f))
  {
    stop('The lengths of the variable of integration and the integrand do not match.')
  }
  ### finish checks
  # obtain length of variable of integration and integrand
  n = length(x)
  # integrate using the trapezoidal rule
  integral = 0.5*sum((x[2:n] - x[1:(n-1)]) * (f[2:n] + f[1:(n-1)]))
  # print the definite integral
  return(integral)
}

#=============================================================================================
#=============================================================================================
setwd("~/My Projects at NICR/2014/Infant/Our New Material")

## loading data
data0_5 <- read.csv("Infant_2SubgroupInSHH.csv", header=T) # Paper, sub-subgroup density plot, March 2017

Cohort_size <- nrow(data0_5) # n=51
Cohort_name <- paste("Age Distribution, Under 5, n=", Cohort_size)
lb <- 0   # lower bound of the age 
hb <- 5   # higher bound of the age 
 
data0_5_SHH_Grp1 <- data0_5 

Temp_Col <- data0_5_SHH_Grp1$SHH_Subsubgroup_March17
Temp_Col <- as.numeric(data0_5_SHH_Grp1$SHH_Subsubgroup_March17)
Temp_Col[Temp_Col[] == "1"] <- "Group 1"
Temp_Col[Temp_Col[] == "2"] <- "Group 2"
data0_5_SHH_Grp1$SHH_Subsubgroup_March17 <- Temp_Col

infant_cohort <- data0_5_SHH_Grp1
Cohort_size <- nrow(infant_cohort)

# Density plot Group 1 vs. Group 2 (SHH sub-subgroup)
ageGrp1vsGrp2 <- infant_cohort$Ageatdiagnosis
ageGrp1vsGrp2 <- ageGrp1vsGrp2[infant_cohort$SHH_Subsubgroup_March17 == "Group 1"| infant_cohort$SHH_Subsubgroup_March17 == "Group 2"]
tmp_group1 <- factor(infant_cohort$SHH_Subsubgroup_March17[infant_cohort$SHH_Subsubgroup_March17 == "Group 1" | infant_cohort$SHH_Subsubgroup_March17 == "Group 2"])

# obtain values of the kernel density estimates
density_var1 = density(ageGrp1vsGrp2[tmp_group1 == "Group 1"], na.rm = T, from = lb, to = hb)
density_var2 = density(ageGrp1vsGrp2[tmp_group1 == "Group 2"], na.rm = T, from = lb, to = hb)

###############################

kde = density_var1$y
support_var1 = density_var1$x
integral_kde = trapezoidal_integration(support_var1, kde)
kde_scaled = kde/integral_kde

kde2 = density_var2$y
support_var2 = density_var2$x
integral_kde2 = trapezoidal_integration(support_var2, kde2)
kde_scaled2 = kde2/integral_kde2

# number of points used in density plot
n_density1 = density_var1$n
n_density2 = density_var2$n

# bandwidth in density plot
bw_density1 = density_var1$bw
bw_density2 = density_var2$bw

# For SHH - Group 1 vs. Group 2
ifelse(max(kde_scaled2) > max(kde_scaled),yub <- max(kde_scaled2),yub <- max(kde_scaled))
plot(1, type="n", xlab="", ylab="", xlim=c(0, 5), ylim=c(0,1))
lines(support_var1, kde_scaled, col = 'black', lwd=3)  # 
lines(support_var2, kde_scaled2, col = 'black', lwd=3, lty="longdash") # 
axis(1, at=lb:hb)
rug(ageGrp1vsGrp2[tmp_group1 == "Group 1"], col="black", side = 1, ticksize=0.2,lwd = 2)
rug(ageGrp1vsGrp2[tmp_group1 == "Group 2"], col="black", side = 1, ticksize=0.2,lwd = 2, lty="longdash")
legend('topright', c("Group 1","Group 2"), 
       lty=1:2, col=c('black', 'black'), bty='n', cex=1.0)

#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
# End
##################################################################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##################################################################################################
