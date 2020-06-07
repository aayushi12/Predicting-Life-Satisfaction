#Group 6
#Aayushi Pandey (2038981/u848144)
#Shreshtha Sharma (2037590/u434770)
#Nasharetty Wiesken (2011820/u876834)
#Carolyn Landbrug (2045788/u901160)

rm(list = ls(all = TRUE)) # Clear workspace

set.seed(235711)

dataDir <- "../data/"
fileName <- "wvs_data.rds"

##Data cleaning

data <- readRDS(paste0(dataDir, fileName))
subset_data <- data[,c("V23","V4","V5","V6","V7","V8","V9","V10","V11","V45","V46","V47","V48","V49","V50","V51","V52","V53","V54","V55","V57","V58","V59","V75","V102","V104","V143","V153","V181","V182","V183","V184","V185","V186","V188","V189","V190","V191","V194","V195","V232","V233","V238","V239","V242","V240","V248","V250","V170","V203","V204","V205","V206","V207","V208","V209","V210")]

library('naniar')
library('dplyr')

subset_data$V181[subset_data$V181 == -3] <- 5
subset_data$V182[subset_data$V182 == -3] <- 5

na_value <- c(-3)
subset_data_na <- subset_data %>% replace_with_na_all(condition = ~.x %in% na_value)

subset_data_na <- na.omit(subset_data_na)

na_values <- c(-1,-2,-4,-5)
subset_data_na <- subset_data_na %>% replace_with_na_all(condition = ~.x %in% na_values)

#Finally we have 12243 observations after dealing with -3, and setting all other missing values to NA

#Univariate Outlier Analysis

#On Age V242, Rest of all is a fixed scale 1-10

#Median Absolute Deviation Method 

madOutliers <- function(x, cut = 2.5, na.rm = TRUE) {
  ## Compute the median and MAD of x:
  mX   <- median(x, na.rm = na.rm)
  madX <- mad(x, na.rm = na.rm)
  
  ## Return row indices of observations for which |T_MAD| > cut:
  which(abs(x - mX) / madX > cut)
} 

madOutliers(subset_data_na$V242, cut = 2.5) # cut = 2.5: 2402 2439 2733 3157 7129 7445

indices = list(madOutliers(subset_data_na$V242, cut = 2.5))

#Setting the univariate outliers to NA
for (i in indices){
  subset_data_na$V242[i] = NA
}

#Missing data descriptives
library(mice) 
library(MASS)

pm <- colMeans(is.na(subset_data_na))*100
pm #percentage missing for each column 

range(pm) #0.0245038 14.2122029
mean(pm) #4.545813
median(pm) #4.141142

#Covariance coverage descriptives
cc <- md.pairs(subset_data_na)$rr / nrow(subset_data_na) 
cc
range(cc) #0.7797109 0.9997550

pat <- cc <= 0.8
apply(pat, 1, function(x) names(x)[x])

#(V153,V186), (V186,V203) are the problematic pairs

##Multiple imputation
library(miceadds)
source("studentFunctions.R")
source("miPredictionRoutines.R")

meth        <- rep("pmm", ncol(subset_data_na)) 
names(meth) <- colnames(subset_data_na)

subset_data_na$V240 <- as.factor(subset_data_na$V240)
subset_data_na$V57 <- as.factor(subset_data_na$V57)
subset_data_na$V250 <- as.factor(subset_data_na$V250)

meth["V240"]    <- "logreg" #Gender
meth["V57"] <- "polyreg"  #Marital status (1,2,3,4,5,6)
meth["V250"] <- "logreg" #Yes/No Answer: Do you live with your parents

#Generate a predictor matrix
predMat <- quickpred(subset_data_na, mincor = 0.2, include = "V240")

# Impute missing using the predictor matrix from above:
miceOut <- mice(data            = subset_data_na,
                m               = 20,
                maxit           = 10,
                method          = meth,
                predictorMatrix = predMat,
                seed            = 235711)
#Remove the factor variables from each imputed dataset (because we don't want
#categorical variables in our data when checking for multivariate outliers):

miceOut2 <-
  subset_datlist(datlist = miceOut,
                 select  = setdiff(colnames(subset_data_na), c("V240", "V250","V57")),
                 toclass = "mids")

impList <- complete(miceOut2, "all") #This list is created for checking multivariate outliers

#Checking the quality of imputation: Performing convergence checks

plot(miceOut,"V232", layout = c(2,1))
plot(miceOut,"V242", layout = c(2,1))
plot(miceOut, "V57",layout = c(2,1))

#Sanity checks: imputed vs observed values

densityplot(miceOut,~V232)
densityplot(miceOut,~V242)
densityplot(miceOut,~V57)

#Carrying out Multivariate Outlier Analysis on the 20 imputed datasets

mcdMahalanobis <- function(data, prob, ratio = 0.75, seed = NULL) {
  ## Set a seed, if one is provided:
  if(!is.null(seed)) set.seed(seed)
  
  ## Compute the MCD estimates of the mean and covariance matrix:
  stats <- cov.mcd(data, quantile.used = floor(ratio * nrow(data)))
  
  ## Compute robust squared Mahalanobis distances
  md <- mahalanobis(x = data, center = stats$center, cov = stats$cov)
  
  ## Find the cutoff value:
  crit <- qchisq(prob, df = ncol(data))
  
  ## Return row indices of flagged observations:
  which(md > crit)
}

outlier_list <- vector(mode = "list", length = 20)

dataset_number = c('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20')

for (i in dataset_number ){
  print(i)
  index = as.integer(i)
  data = impList[index][[i]]
  outliers = mcdMahalanobis(data = data[ ,!(colnames(data) %in% c("V102","V4"))] , prob = 0.999, seed = 235711) #Some columns are excluded because they have IQR 0
  outlier_list[[index]] = outliers
}

mv_outliers = unlist(outlier_list, recursive = TRUE, use.names = TRUE)

unique_mv_outliers = unique(mv_outliers) #1642

count_outlier <- vector(mode = "integer", length = length(unique_mv_outliers))
#Count outlier keeps count in how many datasets was observation flagged as outlier
tab = table(mv_outliers)

for(i in 1:length(unique_mv_outliers)){
  row = unique_mv_outliers[i]
  count = tab[names(tab) == row]
  count_outlier[i] = count
}

hist(count_outlier, main = "", ylab = "Number of observations", xlab = "Number of Datasets", breaks = 5, col = "green", border = "blue")
length(unique_mv_outliers[count_outlier == 20]) #99.9 - 893 
length(unique_mv_outliers[count_outlier >= 15]) #99.9 - 964

#We can remove the rows which are labelled as multivariate outliers(in atleast 15 out of 20 datasets)
#Then we can report our analysis twice with and without deletion

indices_mv_outliers = unique_mv_outliers[count_outlier >= 15]
#This leaves us with 12243 - 964 = 11279 observations

#Doing multiple imputation after removing multivariate outliers

miceOut_mv <- subset_datlist(datlist = miceOut, # We're using the original imputations
                             subset  = setdiff(1 : nrow(subset_data_na), indices_mv_outliers),
                             toclass = "mids")

impList_mv <- complete(miceOut_mv, "all")
impList <- complete(miceOut, "all") #Resetting implist according to original imputation


###The data cleaning part is finished##





##Inferential Modeling##

# Variable information:
# social conservative predictors: V45 (When jobs are scarce men should have more right to do job than women)
# V47(when woman earns more money than her husband it is almost certain to cause problems)
# V203-210 
# demographics variables/confounders: Age V242 (continuous), sex V240 (1male), education V248 (no education/low to high level 9 levels), Martial status V57 

library("ggplot2")
library('dplyr')
library("miceadds")
select <- dplyr::select

#Exploratory Data Analysis

# Significant positive correlation between age and level of happiness, although small
cor.test(subset_data_na$V10,subset_data_na$V242)

#From the graph we can see that feeling of happiness is evenly distributed over sex
#p <- ggplot(subset_data_na, aes(V10, fill=V240)) +
  #geom_bar()
#p <- p + labs(x = "Feeling of Happiness")
#p <- p + labs(y = "Count")
#p <- p + labs(title =  "Distribution of feeling of happiness over sex")
#p <- p + scale_fill_discrete(name = "Sex", labels = c("Male", "Female", "NA"))
#p

# Significant negative Correlation between feeling of happiness and education - significant though small
cor.test(subset_data_na$V10,subset_data_na$V248)

# significant positive correlation between feeling of happiness and health status 
cor.test(subset_data_na$V10,subset_data_na$V11)

# Some conservative attitudes apppear to be related to sex 
p<- ggplot(subset_data_na, aes(V45, fill=V240)) +geom_bar(position = "dodge")
p <- p + labs(x = "1 - Agree, 2- Neither, 3- Disagree")
p <- p + labs(y = "Count")
p <- p + scale_fill_discrete(name = "Sex", labels = c("Male", "Female","NA"))
p <- p + labs(title =  "Men have more right to jobs when scarce")
p

p<- ggplot(subset_data_na, aes(V47, fill=V240)) +geom_bar(position = "dodge")
p <- p + labs(x = "1 - Agree, 2- Neither, 3- Disagree")
p <- p + labs(y = "Count")
p <- p + scale_fill_discrete(name = "Sex", labels = c("Male", "Female","NA"))
p <- p + labs(title =  "Women earning more cause household issues")
p


#We will compare 2 models
#model 1: all social conservative predictors
#model 2: social conservative predictors + demographics confounding variables and their interactions


#Modelling with data in which multivariate outliers are removed

#model 1: restricted model
fit_intercept <- with(data=miceOut_mv, exp=lm(V10 ~ 1))
fit1 <- with(data=miceOut_mv, exp=lm(V10 ~ V45 + V47+ V203 +V204 +
                                       V205 + V206 + V207+ V208+ V209 + V210))

pool1 <- pool(fit1)
summary(pool1)

#how much variation in level of happiness is explained by the model
#R^2 values
fit1_r <- pool.r.squared(fit1) #0.02732127


#test if the R^2 > 0
f_r_test1 <- pool.compare(fit1, fit_intercept) 

f_r_test1$Dm     # Test statistic 19.98763
f_r_test1$pvalue # P-Value=0

#model 2: Full model

fit2 <- with(data=miceOut_mv, exp=lm(V10 ~ V45+V47+ V203 +V204 
                                     +V205+ V206 + V207+ V208+ V209 + V210+ V242+ V248 + V240 + V240*V47 + V240*V45 + V11))
pool2 <- pool(fit2)
summary(pool2)
fit2_r <- pool.r.squared((fit2)) #0.1912979

#test if the R^2 > 0
f_r_test2 <- pool.compare(fit2, fit_intercept) 
f_r_test2$Dm     # Test statistic  125.1444
f_r_test2$pvalue # P-Value=0

#Both models R^2 is different from zero with the inclusion of predictors. 

#compare R^2 to quantify the additional variation explained by the demographic variables
fit2_r - fit1_r # 0.1639766

### Model comparison ###
# We compare pooled estimates from model 1 and 2
# we test if change in R^2= 0.1639766 represents a significantly greater degree of explained variation

f_test <- pool.compare(fit2,fit1) 
f_test$Dm # test statistics 325.221
f_test$pvalue # pvalue= 0 

#The increase is significantly greater than zero.


#Modelling with data in which multivariate outliers aren't removed

#model 1: restricted model
fit_intercept <- with(data=miceOut, exp=lm(V10 ~ 1))
fit1 <- with(data=miceOut, exp=lm(V10 ~ V45 + V47+ V203 +V204 +
                                    V205 + V206 + V207+ V208+ V209 + V210))

pool1 <- pool(fit1)
summary(pool1)

#how much variation in level of happiness is explained by the model
#R^2 values
fit1_r <- pool.r.squared(fit1) #0.02184623

#test if the R^2 > 0
f_r_test1 <- pool.compare(fit1, fit_intercept) 

f_r_test1$Dm     # Test statistic 19.48286
f_r_test1$pvalue # P-Value=0

#model 2: Full model

fit2 <- with(data=miceOut, exp=lm(V10 ~ V45+V47+ V203 +V204 
                                  +V205+ V206 + V207+ V208+ V209 + V210+ V242+ V248 + V240 + V240*V47 + V240*V45+ V11))
pool2 <- pool(fit2)
summary(pool2)
fit2_r <- pool.r.squared((fit2)) #0.1902551

#test if the R^2 > 0
f_r_test2 <- pool.compare(fit2, fit_intercept) 
f_r_test2$Dm     # Test statistic 140.3137
f_r_test2$pvalue # P-Value=0

#Both models R^2 is different from zero with the inclusion of predictors. 

#compare R^2 to quantify the additional variation explained by the demographic variables
fit2_r - fit1_r #0.1684089

### Model comparison ###
# We compare pooled estimates from model 1 and 2
# we test if change in R^2= 0.1684089 represents a significantly greater degree of explained variation

f_test <- pool.compare(fit2,fit1) 
f_test$Dm # test statistics 367.6707
f_test$pvalue # pvalue= 0 significantly 

#The increase is significantly greater than zero.







##Predictive Modelling##

#EDA#

library("ggplot2")

# A great part of the sample is married: makes sense to include family related satisfaction indicators
#p <- ggplot(subset_data_na, aes(V57, fill=V240)) +geom_bar(position = "dodge")
#p <- p + labs(x = "Marital Status")
#p <- p + labs(y = "Count")
#p <- p + scale_fill_discrete(name = "Sex", labels = c("Male", "Female", "NA"))
#p <- p + labs(title =  "Distribution of marital status in sample")
#p

# Satisfaction relation with gender: not a clear pattern
#p<- ggplot(subset_data_na, aes(V23, fill=V240)) +geom_bar(position = "dodge")
#p <- p + labs(x = "Satisfaction with life")
#p <- p + labs(y = "Count")
#p <- p + scale_fill_discrete(name = "Sex", labels = c("Male", "Female", "NA"))
#p <- p + labs(title =  "Distribution of sex over satisfaction")
#p


# Satisfaction relation with marital status : not a clear pattern
#p<- ggplot(subset_data_na, aes(V23, fill=V57)) +geom_bar(position = "dodge")
#p <- p + labs(x = "Satisfaction with life")
#p <- p + labs(y = "Count")
#p <- p + scale_fill_discrete(name = "Marital Status", labels = c("Married", "Living Together", "Divorced","Separated","Widowed","Single"))
#p <- p + labs(title =  "Distribution of Marital Status over satisfaction")
#p


#Correlations

#Financial satisfaction of household: 0.49
cor.test(x = subset_data_na$V23, y = subset_data_na$V59, use = "pairwise.complete.obs") #0.49, Significant

#How much do you trust your family: -0.096
cor.test(x = subset_data_na$V23, y = subset_data_na$V102, use = "pairwise.complete.obs")

#How many children do you have: -0.021, Rejected  at 1% significance
cor.test(x = subset_data_na$V23, y = subset_data_na$V58, use = "pairwise.complete.obs")

#Family felt unsafe from crime in own house: 0.149
cor.test(x = subset_data_na$V23, y = subset_data_na$V189, use = "pairwise.complete.obs")

#Feeling of happiness: -0.489
cor.test(x = subset_data_na$V23, y = subset_data_na$V10, use = "pairwise.complete.obs")

#State of health: -0.351
cor.test(x = subset_data_na$V23, y = subset_data_na$V11, use = "pairwise.complete.obs")

#Freedom of choice and control over life: 0.443
cor.test(x = subset_data_na$V23, y = subset_data_na$V55, use = "pairwise.complete.obs")

#Being sucessful is important to this person: -0.044
cor.test(x = subset_data_na$V23, y = subset_data_na$V75, use = "pairwise.complete.obs")

#How much do you trust people you know personally: -0.126
cor.test(x = subset_data_na$V23, y = subset_data_na$V104, use = "pairwise.complete.obs")

#Thinking about meaning and purpose of life: -0.031
cor.test(x = subset_data_na$V23, y = subset_data_na$V143, use = "pairwise.complete.obs")

#Nature of tasks: routine vs creative: 0.153
cor.test(x = subset_data_na$V23, y = subset_data_na$V232, use = "pairwise.complete.obs")

#Nature of tasks: independence : 0.217
cor.test(x = subset_data_na$V23, y = subset_data_na$V233, use = "pairwise.complete.obs")

#Age: -0.037
cor.test(x = subset_data_na$V23, y = subset_data_na$V242, use = "pairwise.complete.obs")

#Education: 0.156
cor.test(x = subset_data_na$V23, y = subset_data_na$V248, use = "pairwise.complete.obs")

#Worry: Losing my job: 0.107
cor.test(x = subset_data_na$V23, y = subset_data_na$V181, use = "pairwise.complete.obs")

#Social class: -0.260
cor.test(x = subset_data_na$V23, y = subset_data_na$V238, use = "pairwise.complete.obs")

#Secure in neighbourhood: -0.203
cor.test(x = subset_data_na$V23, y = subset_data_na$V170, use = "pairwise.complete.obs")

#Family gone without enough food: 0.201
cor.test(x = subset_data_na$V23, y = subset_data_na$V188, use = "pairwise.complete.obs")

#Family gone without needed medicine: 0.211
cor.test(x = subset_data_na$V23, y = subset_data_na$V190, use = "pairwise.complete.obs")

#Family gone without cash income: 0.2715
cor.test(x = subset_data_na$V23, y = subset_data_na$V191, use = "pairwise.complete.obs")

#Interactions

#State of health and age: 0.325301, Sigificant
cor.test(x = subset_data_na$V11, y = subset_data_na$V242, use = "pairwise.complete.obs")

#State of health and gender: Men healthier
p<- ggplot(subset_data_na, aes(V11, fill=V240)) +geom_bar(position = "dodge")
p <- p + labs(x = "State of health")
p <- p + labs(y = "Count")
p <- p + scale_fill_discrete(name = "Sex", labels = c("Male", "Female", "NA"))
p <- p + labs(title =  "Distribution of sex over state of health")
p

#Financial situation of household with education: 0.091, Significant
cor.test(x = subset_data_na$V248, y = subset_data_na$V59, use = "pairwise.complete.obs")

#Financial situation of household with social class: -0.3128789, Significant
cor.test(x = subset_data_na$V238, y = subset_data_na$V59, use = "pairwise.complete.obs")

#Financial situation of household with V188,V190, V191 >0.2 and significant
cor.test(x = subset_data_na$V188, y = subset_data_na$V59, use = "pairwise.complete.obs")
cor.test(x = subset_data_na$V190, y = subset_data_na$V59, use = "pairwise.complete.obs")
cor.test(x = subset_data_na$V191, y = subset_data_na$V59, use = "pairwise.complete.obs")


#Full Model, Restricted Model and both with interactions 
#We will first work on the data which has multivariate outliers removed 

#Model 1: Full Model
pred_fit <- lm.mids(V23 ~ V59 + V102 + V182 + V188 + V189 + V190 + V191 + V10+ V11+ V55+ V75+ V104+ V143+ V232 + V233 + V238 + V170 + V240 +V242+ V248+ V57 + V181, data = miceOut_mv) 
pred_poolFit <- pool(pred_fit)
summary(pred_poolFit)


# Model 2: Restricted Model
pred_fit_1 <- lm.mids(V23 ~ V59 + V188 + V190 + V191 + V10+ V11+ V55+ V233 + V238 + V170 + V240 + V57, data = miceOut_mv) 
pred_poolFit_1 <- pool(pred_fit_1)
summary(pred_poolFit_1)


# Model 3: Full Model + Interactions
pred_fit_2 <- lm.mids(V23 ~ V59 + V102 + V182 + V188 + V189 + V190 + V191 + V10+ V11+ V55+ V75+ V104+ V143+ V232 + V233 + V238 + V170 + V240 +V242+ V248+ V57 + V181 + V11*V242 + V11*V240 + V59*V248 + V59*V238 + V59*V188 + V59*V190 + V59*V191, data = miceOut_mv) 
pred_poolFit_2 <- pool(pred_fit_2)
summary(pred_poolFit_2)


# Model 4: Restricted Model + Interactions
pred_fit_3 <- lm.mids(V23 ~ V59 + V188 + V190 + V191 + V10+ V11+ V55+ V233 + V238 + V170 + V240 + V57 + V11*V242 + V11*V240 + V59*V248 + V59*V238 + V59*V188 + V59*V190 + V59*V191, data = miceOut_mv) 
pred_poolFit_3 <- pool(pred_fit_3)
summary(pred_poolFit_3)


#CV
library(MLmetrics)

#MI-Based Prediction Cross-Validation 

#Split the multiply imputed data into training, validation, and testing sets:
set.seed(235711)

n = nrow(impList_mv$`1`)
index <- sample(
  c(rep("train", 6000), rep("valid", 3279), rep("test", n - 9279))
)

impList_2 <- splitImps(imps = impList_mv, index = index)

# Define some models to compare:
mods <- c("V23 ~ V59 + V102 + V182 + V188 + V189 + V190 + V191 + V10+ V11+ V55+ V75+ V104+ V143+ V232 + V233 + V238 + V170 + V240 +V242+ V248+ V57 + V181",
          "V23 ~ V59 + V188 + V190 + V191 + V10+ V11+ V55+ V233 + V238 + V170 + V240 + V57",
          "V23 ~ V59 + V102 + V182 + V188 + V189 + V190 + V191 + V10+ V11+ V55+ V75+ V104+ V143+ V232 + V233 + V238 + V170 + V240 +V242+ V248+ V57 + V181 + V11*V242 + V11*V240 + V59*V248 + V59*V238 + V59*V188 + V59*V190 + V59*V191",
          "V23 ~ V59 + V188 + V190 + V191 + V10+ V11+ V55+ V233 + V238 + V170 + V240 + V57 + V11*V242 + V11*V240 + V59*V248 + V59*V238 + V59*V188 + V59*V190 + V59*V191"
)


#Merge the MI training and validations sets:
index2   <- gsub(pattern = "valid", replacement = "train", x = index)
impList_3 <- splitImps(impList_mv, index2)

# K-Fold Cross-Validation:

# Conduct 10-fold cross-validation in each multiply imputed dataset:
tmp <- sapply(impList_3$train, cv.lm, K = 10, models = mods, seed = 235711)

## Aggregate the MI-based CVEs:
cve <- rowMeans(tmp)
cve

## Refit the winning model and compute test-set MSEs:
fits <- lapply(X   = impList_3$train,
               FUN = function(x, mod) lm(mod, data = x),
               mod = mods[which.min(cve)])
mse <- mseMi(fits = fits, newData = impList_3$test)

mse


#Predictive Modelling without removing Multivariate Outliers

#Model 1: Full Model
pred_fit <- lm.mids(V23 ~ V59 + V102 + V182 + V188 + V189 + V190 + V191 + V10+ V11+ V55+ V75+ V104+ V143+ V232 + V233 + V238 + V170 + V240 +V242+ V248+ V57 + V181, data = miceOut) 
pred_poolFit <- pool(pred_fit)
summary(pred_poolFit)


# Model 2: Restricted Model
pred_fit_1 <- lm.mids(V23 ~ V59 + V188 + V190 + V191 + V10+ V11+ V55+ V233 + V238 + V170 + V240 + V57, data = miceOut) 
pred_poolFit_1 <- pool(pred_fit_1)
summary(pred_poolFit_1)


# Model 3: Full Model + Interactions
pred_fit_2 <- lm.mids(V23 ~ V59 + V102 + V182 + V188 + V189 + V190 + V191 + V10+ V11+ V55+ V75+ V104+ V143+ V232 + V233 + V238 + V170 + V240 +V242+ V248+ V57 + V181 + V11*V242 + V11*V240 + V59*V248 + V59*V238 + V59*V188 + V59*V190 + V59*V191, data = miceOut) 
pred_poolFit_2 <- pool(pred_fit_2)
summary(pred_poolFit_2)


# Model 4: Restricted Model + Interactions
pred_fit_3 <- lm.mids(V23 ~ V59 + V188 + V190 + V191 + V10+ V11+ V55+ V233 + V238 + V170 + V240 + V57 + V11*V242 + V11*V240 + V59*V248 + V59*V238 + V59*V188 + V59*V190 + V59*V191, data = miceOut) 
pred_poolFit_3 <- pool(pred_fit_3)
summary(pred_poolFit_3)


#CV

#MI-Based Prediction Cross-Validation 

#Split the multiply imputed data into training, validation, and testing sets:
set.seed(235711)

n = nrow(impList$`1`)
index <- sample(
  c(rep("train", 7000), rep("valid", 3243), rep("test", n - 10243))
)


impList_2 <- splitImps(imps = impList, index = index)

# Define some models to compare:
mods <- c("V23 ~ V59 + V102 + V182 + V188 + V189 + V190 + V191 + V10+ V11+ V55+ V75+ V104+ V143+ V232 + V233 + V238 + V170 + V240 +V242+ V248+ V57 + V181",
          "V23 ~ V59 + V188 + V190 + V191 + V10+ V11+ V55+ V233 + V238 + V170 + V240 + V57",
          "V23 ~ V59 + V102 + V182 + V188 + V189 + V190 + V191 + V10+ V11+ V55+ V75+ V104+ V143+ V232 + V233 + V238 + V170 + V240 +V242+ V248+ V57 + V181 + V11*V242 + V11*V240 + V59*V248 + V59*V238 + V59*V188 + V59*V190 + V59*V191",
          "V23 ~ V59 + V188 + V190 + V191 + V10+ V11+ V55+ V233 + V238 + V170 + V240 + V57 + V11*V242 + V11*V240 + V59*V248 + V59*V238 + V59*V188 + V59*V190 + V59*V191"
)


#Merge the MI training and validations sets:
index2   <- gsub(pattern = "valid", replacement = "train", x = index)
impList_3 <- splitImps(impList, index2)

# K-Fold Cross-Validation:

# Conduct 10-fold cross-validation in each multiply imputed dataset:
tmp <- sapply(impList_3$train, cv.lm, K = 10, models = mods, seed = 235711)

# Aggregate the MI-based CVEs:
cve <- rowMeans(tmp)
cve

## Refit the winning model and compute test-set MSEs:
fits <- lapply(X   = impList_3$train,
               FUN = function(x, mod) lm(mod, data = x),
               mod = mods[which.min(cve)])
mse <- mseMi(fits = fits, newData = impList_3$test)

mse




