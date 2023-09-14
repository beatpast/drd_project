# step 0.Setting the environment
# To check the objects that are present in your workspace:
ls()
rm(list=ls())

# To know the current working directory

getwd()

# To set the new working directory

setwd("C:/Users/beatr/OneDrive/Documenti/dna_rna_dynamics/report")

## step 1.Load raw data with minfi and create an object called RGset storing the RGChannelSet object 
# checking if I have minfi package
library(minfi)
vignette("minfi")
# Loading input data
list.files("C:/Users/beatr/OneDrive/Documenti/dna_rna_dynamics/report/input")
SampleSheet <- read.table("C:/Users/beatr/OneDrive/Documenti/dna_rna_dynamics/report/input/Samplesheet_report_2023.csv",sep=",",header=T)
SampleSheet
# Setting the directory in which the input data are stored and load the samplesheet using the function read.metharray.sheet
baseDir <- ("C:/Users/beatr/OneDrive/Documenti/dna_rna_dynamics/report/input")
targets <- read.metharray.sheet(baseDir)
targets
# Creating an object of class RGChannelSet using the function read.metharray.exp
?read.metharray.exp
# create an object called RGset storing the RGChannelSet object
RGset <- read.metharray.exp(targets = targets)
# saving RGset
save(RGset,file="RGset.RData")
# step 2.Create the dataframes Red and Green to store the red and green fluorescences respectively 
# Let's explore the RGset object:
str(RGset)
?RGChannelSet
# extracting the Green and Red Channels using the functions getGreen and getRed
Red <- data.frame(getRed(RGset))
dim(Red)
head(Red)
Green <- data.frame(getGreen(RGset))
dim(Green)
head(Green)

## step 3.	Fill the following table: what are the Red and Green fluorescences for the address assigned to me (45652402) ? 
# The Illumina 450k manifest has been downloaded from this link: http://support.illumina.com/array/array_kits/infinium_humanmethylation450_beadchip_kit/downloads.html
# loading the cleaned version of this manifest file
load('C:/Users/beatr/OneDrive/Documenti/dna_rna_dynamics/Illumina450Manifest_clean.RData')


head(Illumina450Manifest_clean)
# I want to check the probe having the address assigned to me (45652402) (The first rowname in Red and Green objects)
Illumina450Manifest_clean[Illumina450Manifest_clean$AddressA_ID=="45652402",]
Illumina450Manifest_clean[Illumina450Manifest_clean$AddressB_ID=="45652402",]
# the probe is type I and it has code cg01707559, the address assigned to me is addressA (45652402) of this probe because if I check addressB is empty with this id code
# checking the probe 
Illumina450Manifest_clean[Illumina450Manifest_clean$IlmnID=="cg01707559",]
# the addressB of this probe Id is 64689504 and the color channel for this probe is red 
Red[rownames(Red)=="45652402",]
Red[rownames(Red)=="64689504",]
# checking green color
Green[rownames(Green)=="45652402",]
Green[rownames(Green)=="64689504",]
# these are out of band signals
# every probe contains data in both colors fluorescence intensity even if it was not emitted green. the output of this data is important for normalization(separate background signal from real signal)
# in type I probes the color doesn't give me info about the methylation status. AddressA (identify beads with binding unmethyl) and B(=identify beads with binding to methylated version of the genome)
# creating a subset with red fluoroscence values for my addressA (45652402) of the probe cg01707559
# Subsetting the Red and Green data frames for my addressA id
red_subset <- Red[rownames(Red) == "45652402", ]
green_subset <- Green[rownames(Green) == "45652402", ]
dim(red_subset)
red_subset
rownames(red_subset)
# Creating a new data frame with the desired columns
table_red_green <- data.frame(
  sample = colnames(red_subset),
  Red.fluor = c(red_subset[1,1], red_subset[1,2], red_subset[1,3], red_subset[1,4], red_subset[1,5], red_subset[1,6], red_subset[1,7], red_subset[1,8]),
  Green.fluor = c(green_subset[1,1], green_subset[1,2], green_subset[1,3], green_subset[1,4], green_subset[1,5], green_subset[1,6], green_subset[1,7], green_subset[1,8]),
  Type = "Type I",
  Color = "Red" #because color channel is red
)
table_red_green
## step 4.Create the object MSet.raw 
# Extracting methylated and unmethylated signals
# using the function MSet.raw
MSet.raw <- preprocessRaw(RGset)
MSet.raw
?MethylSet
save(MSet.raw,file="MSet_raw.RData")
Meth <- getMeth(MSet.raw)
str(Meth)
head(Meth)
Unmeth <- getUnmeth(MSet.raw)
str(Unmeth)
head(Unmeth)
Unmeth[rownames(Unmeth)=="cg01707559",]
Meth[rownames(Meth)=="cg01707559",]
# I get results for both unmeth and meth
## step 5.	Perform the following quality checks and provide a brief comment to each step:
# 5.1 QCplot
# First of all, considering the median of methylation and unmethylation channels for each sample, using the function getQC 
qc <- getQC(MSet.raw)
qc
plotQC(qc)
# All the samples have high median methylation and unmethylation signals, they have good quality

# Then,checking the control probes. The Illumina guide (GenomeStudio_control_probes.pdf) tells us the expected intensities for each type of control probes.

# To know the names and the numbers of the different types of control probes, apply the getProbeInfo function to the RGset object.
getProbeInfo(RGset, type = "Control")
df_TypeControl <- data.frame(getProbeInfo(RGset, type = "Control"))
str(df_TypeControl)
table(df_TypeControl$Type)
## 5.2 checking the intensity of negative controls using minfi
controlStripPlot(RGset, controls="NEGATIVE")
# Negative controls are all fine, as all below 1000 (log2(1000)=10) they should be below 10 in this graph
# Green and red labels are swapped because of the version of R
## 5.3 calculate detection p-values; for each sample, how many probes have a detection p-value higher than the threshold assigned to me(0.01)?
getwd()
?detectionP
# using The detectionP() functions.
detP <- detectionP(RGset) 
save(detP,file="detP.RData")
load("detP.RData") 
str(detP)
dim(detP)
head(detP)
# In the detP object rownames are the CpG probes names
# Now considering a detectionP threshold of 0.01
failed <- detP>0.01
head(failed)
dim(failed)
table(failed)
# summary(failed) is particularly useful, as it returns the number of failed (TRUE) and not failed (FALSE) positions for each sample
summary(failed) 
# creating the table
table_failed <- data.frame(
  failed_positions = colSums(failed == TRUE)
)
table_failed
## 6.	Calculate raw beta and M values and plot the densities of mean methylation values, dividing the samples in WT and MUT (suggestion: subset the beta and M values matrixes in order to retain WT or MUT subjects and apply the function mean to the 2 subsets). Do you see any difference between the two groups?
# Beta and M values
# We load the MSet object, which contains the methylation and unmethylation signals
load("C:/Users/beatr/OneDrive/Documenti/dna_rna_dynamics/report/MSet_raw.RData")
ls()
MSet.raw
# The functions getBeta and getM allow to retrieve beta and M values matrices:
beta <- getBeta(MSet.raw)
class(beta)
dim(beta)
head(beta)
summary(beta)
# In this dataset, minimum for beta values is 0 maximum is 1 We have also some NAs: these are the positions for which both the Methylation and Unmethylation values in the MSet.raw are equal to 0 (by default, in minfi the offset is set to 0)
# Computing M values,M values are more statistical reliable than beta values (std of m values is not effected by mean of m values)
M <- getM(MSet.raw)
dim(M)
head(M)
summary(M)
# In this dataset, minimum for M values is -Inf (Methylation value=0, Unmethylation value>0), maximum is +Inf (Methylation value>0, Unmethylation value=0). I have also some NAs: these are the positions for which both the Methylation and Unmethylation values in the MSet.raw are equal to 0 (by default, in minfi the offset is set to 0)

# Subsetting the RGset based on the wt and mut groups:
wt <- RGset@colData[RGset@colData$Group == "WT",]
mut <- RGset@colData[RGset@colData$Group == "MUT",]
wt
mut
# Storing the names of wt samples
samples_wt <- rownames(wt)
samples_wt
# Storing the names of mut samples
samples_mut <- rownames(mut)
samples_mut
# Subsetting the beta values into beta values for wt and beta values for mut
beta_wt <- beta[,samples_wt]
beta_mut <- beta[,samples_mut]  
beta_wt
beta_mut
# Calculating mean of beta values for wt samples and mean of beta values for mut samples
# To calculate the mean, we will use the mean() function. As we know that in our input there are some missing values, we will specify that NA values should be stripped (see this page https://www.statmethods.net/input/missingdata.html) using the na.rm=T argument.
# To calculate the mean of each row of the matrix, we will use the apply() function, which takes as input a dataframe or a matrix and returns a vector (see here for example: http://petewerner.blogspot.it/2012/12/using-apply-sapply-lapply-in-r.html). Note that the "MARGIN" argument in the function allows to specify that the function mean() should be applied to each row (set "MARGIN" to 1) or to each column (set "MARGIN" to 2) of the matrix/dataframe.

mean_of_beta_wt <- apply(beta_wt,1,mean,na.rm=T)
head(mean_of_beta_wt)
length(mean_of_beta_wt)
mean_of_beta_mut <-apply(beta_mut,1,mean,na.rm=T)
head(mean_of_beta_mut)
length(mean_of_beta_mut)
# Now I can calculate the density distribution
# removing missing values from the vector because density function doesn't accept those
mean_of_beta_wt <- mean_of_beta_wt[!is.na(mean_of_beta_wt)]

d_mean_of_beta_wt <- density(mean_of_beta_wt)
d_mean_of_beta_wt
mean_of_beta_mut <- mean_of_beta_mut[!is.na(mean_of_beta_mut)]
d_mean_of_beta_mut <- density(mean_of_beta_mut)
mean_of_beta_mut
# plotting the beta values
plot(d_mean_of_beta_wt, main="Density of Beta Values wt",col="orange")
lines(d_mean_of_beta_mut,col="red")
# the Beta values appear to be highly bimodal, because are compressed in the low (between 0 and 0.2) and high (between 0.8 and 1) methylation ranges. In other words, we have 2 peaks, one for low methylation values, near 0, and another peak for high methylation values,
# subsetting the m values into m values for wt and m values for mut
M_wt <- M[,samples_wt]
M_mut <- M[,samples_mut] 
M_wt
M_mut
# I apply the same steps on M values to plot M values for MT samples and WT samples:
mean_of_M_wt <- apply(M_wt,1,mean,na.rm=T)
mean_of_M_wt <- mean_of_M_wt[!is.na(mean_of_M_wt)]
mean_of_M_mut <- apply(M_mut,1,mean,na.rm=T)
mean_of_M_mut <- mean_of_M_mut[!is.na(mean_of_M_mut)]
d_mean_of_M_wt <- density(mean_of_M_wt)
d_mean_of_M_mut <- density(mean_of_M_mut)
# plotting the M values in the same plot
plot(d_mean_of_M_wt,main="Density of M Values wt",col="purple")
lines(d_mean_of_M_mut,col="red")
# checking boxplot for beta values
par(mfrow=c(1,2))
boxplot(beta_wt,main="beta values wt", col="purple")
boxplot(beta_mut,main="beta values wt", col="red")
boxplot(beta,main="beta values mut")
## Step 7.Normalize the data using the function assigned to each student and compare raw data and normalized data. Produce a plot with 6 panels in which, for both raw and normalized data, you show the density plots of beta mean values according to the chemistry of the probes, the density plot of beta standard deviation values according to the chemistry of the probes and the boxplot of beta values. Provide a short comment about the changes you observe. 
# subsetting the Illumina450Manifest_clean in two dataframes, containing only type I (dfI) or type II (dfII) probes
dfI <- Illumina450Manifest_clean[Illumina450Manifest_clean$Infinium_Design_Type=="I",]
dfI <- droplevels(dfI)
dim(dfI)
str(dfI)
dfII <- Illumina450Manifest_clean[Illumina450Manifest_clean$Infinium_Design_Type=="II",]
dfII <- droplevels(dfII)
# Now subsetting the beta matrix in order to retain only the rows whose name is in the first column of dfI...
beta_I <- beta[rownames(beta) %in% dfI$IlmnID,]
dim(beta_I)
# ... or in the first column of dfII
beta_II <- beta[rownames(beta) %in% dfII$IlmnID,]
dim(beta_II)

# For each probe in the mean_of_beta_I and mean_of_beta_II matrices, calculating the mean of beta values across the 8 samples...
mean_of_beta_I <- apply(beta_I,1,mean)
mean_of_beta_II <- apply(beta_II,1,mean)

# ... and then we calculate the density distribution of the 2 vectors of mean values:
d_mean_of_beta_I <- density(mean_of_beta_I,na.rm=T)
d_mean_of_beta_II <- density(mean_of_beta_II,na.rm=T)
# For Raw data, I already have the matrix of beta values, the d_mean_of_beta_I and d_mean_of_beta_II objects; I need to calculate the densities of the standard deviations, which can be calculated using the function sd():
sd_of_beta_I <- apply(beta_I,1,sd,na.rm=T)
sd_of_beta_II <- apply(beta_II,1,sd,na.rm=T)
d_sd_of_beta_I <- density(sd_of_beta_I,)
sd_of_beta_II <- sd_of_beta_II[!is.na(sd_of_beta_II)]
d_sd_of_beta_II <- density(sd_of_beta_II)
# normalisation with preprocessfunnnorm
# This function implements functional normalization preprocessing for Illumina methylation microarrays. Functional normalization extends the idea of quantile normalization by adjusting for known covariates measuring unwanted variation.
?preprocessFunnorm

RGset
preprocessFunnorm_results <- preprocessFunnorm(RGset)
save(preprocessFunnorm_results, file="preprocessFunnorm_results.R")
str(preprocessFunnorm_results)
class(preprocessFunnorm_results)
preprocessFunnorm_results
?GenomicRatioSet
# GetBeta is among the accessor functions of a GenomicRatioSet class object
beta_preprocessFunnorm <- getBeta(preprocessFunnorm_results)
head(beta_preprocessFunnorm)
# saving the file 
save(beta_preprocessFunnorm, file="beta_preprocessFunnorm.RData")
# Now dividing the beta_preprocessFunnorm matrix according to type I and type II probes, calculate the mean and the standartd deviation for each probe across the 8 samples and calculate the density distributions
beta_preprocessFunnorm_I <- beta_preprocessFunnorm[rownames(beta_preprocessFunnorm) %in% dfI$IlmnID,]
beta_preprocessFunnorm_II <- beta_preprocessFunnorm[rownames(beta_preprocessFunnorm) %in% dfII$IlmnID,]
mean_of_beta_preprocessFunnorm_I <- apply(beta_preprocessFunnorm_I,1,mean)
mean_of_beta_preprocessFunnorm_II <- apply(beta_preprocessFunnorm_II,1,mean)
d_mean_of_beta_preprocessFunnorm_I <- density(mean_of_beta_preprocessFunnorm_I,na.rm=T)
d_mean_of_beta_preprocessFunnorm_II <- density(mean_of_beta_preprocessFunnorm_II,na.rm=T)
sd_of_beta_preprocessFunnorm_I <- apply(beta_preprocessFunnorm_I,1,sd)
sd_of_beta_preprocessFunnorm_II <- apply(beta_preprocessFunnorm_II,1,sd)
d_sd_of_beta_preprocessFunnorm_I <- density(sd_of_beta_preprocessFunnorm_I,na.rm=T)
d_sd_of_beta_preprocessFunnorm_II <- density(sd_of_beta_preprocessFunnorm_II,na.rm=T)
# plotting all of them in the same page with par function
par(mfrow=c(2,3))
plot(d_mean_of_beta_I,col="blue",main="raw beta")
lines(d_mean_of_beta_II,col="red")
plot(d_sd_of_beta_I,col="blue",main="raw sd")
lines(d_sd_of_beta_II,col="red")
boxplot(beta,ylim=c(0,1))
plot(d_mean_of_beta_preprocessFunnorm_I,col="blue",main="preprocessFunnorm beta")
lines(d_mean_of_beta_preprocessFunnorm_II,col="red")
plot(d_sd_of_beta_preprocessFunnorm_I,col="blue",main="preprocessFunnorm sd")
lines(d_sd_of_beta_preprocessFunnorm_II,col="red")
boxplot(beta_preprocessFunnorm,ylim=c(0,1))
# It is easier to compare the plots if I have the same scales on x and y axes; to this aim, I can specify the xlim and ylim arguments:
par(mfrow=c(2,3))
plot(d_mean_of_beta_I,col="blue",main="raw beta",xlim=c(0,1),ylim=c(0,5))
lines(d_mean_of_beta_II,col="red")
plot(d_sd_of_beta_I,col="blue",main="raw sd",xlim=c(0,0.6),ylim=c(0,60))
lines(d_sd_of_beta_II,col="red")
boxplot(beta,ylim=c(0,1))
plot(d_mean_of_beta_preprocessFunnorm_I,col="blue",main="preprocessFunnorm beta",xlim=c(0,1),ylim=c(0,5))
lines(d_mean_of_beta_preprocessFunnorm_II,col="red")
plot(d_sd_of_beta_preprocessFunnorm_I,col="blue",main="preprocessFunnorm sd",xlim=c(0,0.6),ylim=c(0,60))
lines(d_sd_of_beta_preprocessFunnorm_II,col="red")
boxplot(beta_preprocessFunnorm,ylim=c(0,1))
# saving the plots in a pdf file
pdf("Plot_comparison_raw_preprocessFunnorm.pdf",height=7,width=15)
par(mfrow=c(2,3))
plot(d_mean_of_beta_I,col="blue",main="raw beta",xlim=c(0,1),ylim=c(0,5))
lines(d_mean_of_beta_II,col="red")
plot(d_sd_of_beta_I,col="blue",main="raw sd",xlim=c(0,0.6),ylim=c(0,60))
lines(d_sd_of_beta_II,col="red")
boxplot(beta,ylim=c(0,1))
plot(d_mean_of_beta_preprocessFunnorm_I,col="blue",main="preprocessFunnorm beta",xlim=c(0,1),ylim=c(0,5))
lines(d_mean_of_beta_preprocessFunnorm_II,col="red")
plot(d_sd_of_beta_preprocessFunnorm_I,col="blue",main="preprocessFunnorm sd",xlim=c(0,0.6),ylim=c(0,60))
lines(d_sd_of_beta_preprocessFunnorm_II,col="red")
boxplot(beta_preprocessFunnorm,ylim=c(0,1))
dev.off()
## Step 8.	Perform a PCA on the matrix of normalized beta values generated in step 7, after normalization. Comment the plot (do the samples divide according to the group? Do they divide according to the sex of the samples? Do they divide according to the batch, that is the column Sentrix_ID?). 
# using the function prcomp() to calculate the PCA on my matrix of beta values.
load("C:/Users/beatr/OneDrive/Documenti/dna_rna_dynamics/report/beta_preprocessFunnorm.RData")
# The matrix has samples in columns and CpG probes in rows: the prcomp() function should be applied to the transposed matrix, which is achieved using the t() function
?prcomp
pca_results <- prcomp(t(beta_preprocessFunnorm),scale=T)

# I can print and plot the variance accounted for each component
print(summary(pca_results))
plot(pca_results)
# Asking the str() of pca_object
str(pca_results)

# The principal components of interest are stored in the element named "x" of the list.
pca_results$x

# Plotting PC1 and PC2 and check if samples cluster according to some variable; the argument "cex" defines the size of the dots, the argument "pch" the dot type:
plot(pca_results$x[,1], pca_results$x[,2],cex=2,pch=2)

# labeling the dots using the text() function and setting the "labels" argument to the rownames of pca_results$x; the argument "pos" defines the position of the text:
?text
text(pca_results$x[,1], pca_results$x[,2],labels=rownames(pca_results$x),pos=1)

pheno <- read.csv("C:/Users/beatr/OneDrive/Documenti/dna_rna_dynamics/report/input/Samplesheet_report_2023.csv",header=T, stringsAsFactors=T)
str(pheno)
# In this dataset we have few samples, so it is ok to plot their names. When you work with larger datasets, on the contrary, it is easier to colour the dots according to some column of the pheno file of interest, for example Group:
pheno$Group
pheno$Sentrix_ID
# defining a palette of colours; each colour will be assigned to each levels of the factor according to the order of levels. For example, in this case Group A will be in orange, Group B in purple.adjusting the margins and add a legend to the plot.
# coloring according to the group wt and mut
levels(pheno$Group)
palette(c("yellow","orange"))
plot(pca_results$x[,1], pca_results$x[,2],cex=2,pch=2,col=pheno$Group,xlab="PC1",ylab="PC2",xlim=c(-1000,1000),ylim=c(-1000,1000))
text(pca_results$x[,1], pca_results$x[,2],labels=rownames(pca_results$x),cex=0.5,pos=1)
legend("bottomright",legend=levels(pheno$Group),col=c(1:nlevels(pheno$Group)),pch=2)

# coloring according to the sex of the subjects
levels(pheno$Sex)
palette(c("green","blue"))
plot(pca_results$x[,1], pca_results$x[,2],cex=2,pch=2,col=pheno$Sex,xlab="PC1",ylab="PC2",xlim=c(-1000,1000),ylim=c(-1000,1000))
text(pca_results$x[,1], pca_results$x[,2],labels=rownames(pca_results$x),cex=0.5,pos=1)
legend("bottomright",legend=levels(pheno$Sex),col=c(1:nlevels(pheno$Sex)),pch=2)

# Do they divide according to the batch, that is the column Sentrix_ID? No, because the levels of pheno$Sentrix_id is null, all the values of this column in the pheno dataframe are the same
levels(pheno$Sentrix_ID)
## step 9.Using the matrix of normalized beta values generated in step 7, identify differentially methylated probes between group WT and group MUT using the function assigned to me Mann-Whitney test (I tried also t-test)
# mannwhitney test
My_mannwhitney_function <- function(x) {
  mannwhitney <- wilcox.test(x~ pheno$Group )
  return(mannwhitney$p.value)
} 

beta_preprocessFunnorm

pValues_wilcox <- apply(beta_preprocessFunnorm,1, My_mannwhitney_function)
pValues_wilcox 
save(pValues_wilcox, file="pValues_wilcox.R")
plot(density(pValues_wilcox ))
hist(pValues_wilcox )
pValues_wilcox_df <- data.frame(pValues_wilcox)
pValues_wilcox_df
length(pValues_wilcox)
# creating a data.frame with all the beta values and the pValue column
final_wilcox <- data.frame(beta_preprocessFunnorm,pValues_wilcox)
head(final_wilcox)
dim(final_wilcox)

# ttest
My_ttest_function <- function(x) {
t_test <- t.test(x~ pheno$Group)
return(t_test$p.value)
} 
pValues_ttest <- apply(beta_preprocessFunnorm,1, My_ttest_function)
pValues_ttest
final_ttest <- data.frame(beta_preprocessFunnorm,pValues_ttest)
## step 10.	Apply multiple test correction and set a significant threshold of 0.05. How many probes do you identify as differentially methylated considering nominal pValues? How many after Bonferroni correction? How many after BH correction?
# applying the Benjamini & Hochberg and the Bonferroni corrections (I used both t test and Mann Withney test
# first with p-values from t test:
corrected_pValues_BH_t <- p.adjust(final_ttest$pValues_ttest,"BH")
corrected_pValues_Bonf_t <- p.adjust(final_ttest$pValues_ttest,"bonferroni")
final_ttest_corrected <- data.frame(final_ttest, corrected_pValues_BH_t, corrected_pValues_Bonf_t)
dim(final_ttest_corrected)
# We can visualize the distributions of the p-values and of the corrected p-values by boxplots:
colnames(final_wilcox_corrected)
boxplot(final_ttest_corrected[,9:11])

# How many probes survive the multiple test correction?
# How many probes do you identify as differentially methylated considering nominal pValues? 
dim(final_ttest_corrected[pValues_ttest<=0.05,])
#63964
# How many after Bonferroni correction? 
dim(final_ttest_corrected[corrected_pValues_BH_t<=0.05,])
# 194 probes survive after bonferroni correction
# How many after BH correction?
dim(final_ttest_corrected[corrected_pValues_Bonf_t<=0.05,])
#3 probes survive after BH correction

# now with Mann withney test p values 
corrected_pValues_BH_w <- p.adjust(final_wilcox$pValues_wilcox,"BH")
corrected_pValues_Bonf_w <- p.adjust(final_wilcox$pValues_wilcox,"bonferroni")
final_wilcox_corrected <- data.frame(final_wilcox, corrected_pValues_BH_w, corrected_pValues_Bonf_w)
dim(final_wilcox_corrected)
# I can visualize the distributions of the p-values and of the corrected p-values by boxplots:
colnames(final_wilcox_corrected)
boxplot(final_wilcox_corrected[,9:11])
# bonferroni correction is really strict, none of the p-values in final_mann were significant after correction.
# How many probes survive the multiple test correction?
# How many probes do you identify as differentially methylated considering nominal pValues? 
dim(final_wilcox_corrected[pValues_wilcox<=0.05,])
# 57261
# How many after Bonferroni correction? 
dim(final_wilcox_corrected[corrected_pValues_BH_w<=0.05,])
# 0 probes survive after bonferroni correction
# How many after BH correction?
dim(final_wilcox_corrected[corrected_pValues_Bonf_w<=0.05,])
# 0 probes survive after BH correction
# step 11.	Produce a volcano plot and a Manhattan plot of the results of differential methylation analysis 
## 11.1 Volcano plots
# First of all, calculating the difference between the average of group WT values and the average of group MUT values. To this aim, first creating two matrixes containing the beta-values of group WT and group MUT samples, and then calculating the mean within each group for each row.
pheno$Group
dim(final_wilcox_corrected)
beta_values <- beta_preprocessFunnorm
beta_groupWT <- beta_values[,pheno$Group=="WT"]
beta_groupMUT
dim(beta_groupWT)
mean_beta_groupWT <- apply(beta_groupWT,1,mean)
beta_groupMUT <- beta_values[,pheno$Group=="MUT"]
mean_beta_groupMUT <- apply(beta_groupMUT,1,mean)
mean_beta_groupWT
# calculating the difference between average values:
delta <- mean_beta_groupMUT - mean_beta_groupWT
head(delta)

# creating dataframe with two columns, one containing the delta values and the other with the -log10 of p-values 
# using p-values from t-test
toVolcPlot <- data.frame(delta, -log10(final_ttest$pValues_ttest))
head(toVolcPlot)
dim(toVolcPlot)
plot(toVolcPlot[,1], toVolcPlot[,2])
beta_preprocessFunnorm
# Let's improve a bit the plot: http://www.statmethods.net/advgraphs/parameters.html
# adding a threshold for pvalue significance (for example, 0.01) by using the abline function

plot(toVolcPlot[,1], toVolcPlot[,2],pch=16,cex=0.5)
-log10(0.01)
?abline
abline(h=-log10(0.01),col="red")

#  highlight the probes (that is, the points), that have a nominal pValue<0.01 and a delta > 0.1)
plot(toVolcPlot[,1], toVolcPlot[,2],pch=16,cex=0.5)
toHighlight <- toVolcPlot[toVolcPlot[,1]>0.1 & toVolcPlot[,2]>(-log10(0.01)),]
head(toHighlight)
points(toHighlight[,1], toHighlight[,2],pch=16,cex=0.7,col="yellow")

# If I want to highlight the points with an absolute delta > 0.01
plot(toVolcPlot[,1], toVolcPlot[,2],pch=16,cex=0.5)
toHighlight <- toVolcPlot[abs(toVolcPlot[,1])>0.1 & toVolcPlot[,2]>(-log10(0.01)),]
head(toHighlight)
points(toHighlight[,1], toHighlight[,2],pch=16,cex=0.7,col="red")

# Volcano plot using p-values from mann-whitney test
toVolcPlot_w <- data.frame(delta, -log10(final_wilcox$pValues_wilcox))
head(toVolcPlot_w)
dim(toVolcPlot_w)
plot(toVolcPlot_w[,1], toVolcPlot_w[,2])
## 11.2 Manhattan plot
install.packages("qqman", version = "0.1.4")
library(qqman)

# First I have to annotate the dataframe, that is add genome annotation information for each cpg probe. I will use the Illumina450Manifest_clean object 
load('C:/Users/beatr/OneDrive/Documenti/dna_rna_dynamics/report/Illumina450Manifest_clean.RData')
# I want to merge on the basis of the CpG probes, but unfortunately in the final_ttest_corrected object the CpG probes are stored in the rownames, not in a column. I can overcome this issue as follows:
final_ttest_1 <- data.frame(rownames(final_ttest),final_ttest)
head(final_ttest_1)
colnames(final_ttest_1)
colnames(final_ttest_1)[1] <- "IlmnID"
colnames(final_ttest_1)

final_ttest_annotated <- merge(final_ttest_1, Illumina450Manifest_clean,by="IlmnID")
dim(final_ttest_1)
dim(final_ttest_annotated)
str(final_ttest_annotated)
# dataframe is automathically reordered on the basis of alphabetical order of the column used for merging
head(final_ttest_annotated)
# Now I can create the input for the Manhattan plot analysis. The input object should contain 4 info: probe, chromosome, position on the chromosome and p-value. I will select these columns in the final_ttest_corrected_annotated object
input_Manhattan <- final_ttest_annotated[colnames(final_ttest_annotated) %in% c("IlmnID","CHR","MAPINFO","pValues_ttest")]
dim(input_Manhattan)
head(input_Manhattan)
str(input_Manhattan)
levels(input_Manhattan$CHR)
# It is better to reorder the levels of the CHR
order_chr <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y")
input_Manhattan$CHR <- factor(input_Manhattan$CHR,levels=order_chr )
levels(input_Manhattan$CHR)
# The function that I will use is "manhattan"
?manhattan
# As you see, the column "CHR" should be numeric --> we will convert factors to numbers
input_Manhattan$CHR <- as.numeric(input_Manhattan$CHR)
table(input_Manhattan$CHR)
# finally I can produce our Manhattan plot
manhattan(input_Manhattan, snp="IlmnID",chr="CHR", bp="MAPINFO", p="pValues_ttest" )
#let's try for mann whitney test p values
final_wilcox_1 <- data.frame(rownames(final_wilcox),final_wilcox)
head(final_wilcox_1)
colnames(final_wilcox_1)
colnames(final_wilcox_1)[1] <- "IlmnID"
colnames(final_wilcox_1)

final_wilcox_annotated <- merge(final_wilcox_1, Illumina450Manifest_clean,by="IlmnID")
dim(final_wilcox_1)
dim(final_wilcox_annotated)
str(final_wilcox_annotated)
# dataframe is automathically reordered on the basis of alphabetical order of the column used for merging
head(final_wilcox_annotated)
# Now I can create the input for the Manhattan plot analysis. The input object should contain 4 info: probe, chromosome, position on the chromosome and p-value. I will select these columns in the final_ttest_corrected_annotated object
input_Manhattan <- final_wilcox_annotated[colnames(final_wilcox_annotated) %in% c("IlmnID","CHR","MAPINFO","pValues_wilcox")]
dim(input_Manhattan)
head(input_Manhattan)
str(input_Manhattan)
levels(input_Manhattan$CHR)
# It is better to reorder the levels of the CHR
order_chr <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y")
input_Manhattan$CHR <- factor(input_Manhattan$CHR,levels=order_chr )
levels(input_Manhattan$CHR)
# The function that I will use is "manhattan"
?manhattan
# As you see, the column "CHR" should be numeric --> we will convert factors to numbers
input_Manhattan$CHR <- as.numeric(input_Manhattan$CHR)
table(input_Manhattan$CHR)
# finally I can produce our Manhattan plot
manhattan(input_Manhattan, snp="IlmnID",chr="CHR", bp="MAPINFO", p="pValues_wilcox" )
save(final_wilcox, file="final_wilcox.R")
