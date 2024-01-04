#*********** FINAL PROJECT *************#
  
#set the working directory
setwd('')

#install and load necessary packages
library(haven)

#load the data
data <- read_sav("00_raw_data\\CY07_MSU_STU_QQQ.sav")

#Choose Italy with 11785 observations, drop other countries
new_data <- subset(data, CNT == 'ITA')

#dependent variable = science literacy score: calculate
#mean value of Plausible values
# Calculate the mean of numeric columns for each row
columns_to_include <- colnames(new_data)[grep("^PV\\d+SCIE$", colnames(new_data))]
new_data$scie_score <- rowMeans(new_data[columns_to_include])

#save the tidy data
write.csv(new_data, file = '01_tidy_data\\italy_data', row.names = FALSE)

rm(list = ls())

library(readr)
df <- read.csv("01_tidy_data\\italy_data")

#I. DEVELOP THE LINEAR REGRESSION MODEL:
##1. Clean the data: remove all missing values
data <- df[c('scie_score', 'ST004D01T', 'ST011Q04TA', 
               'ST034Q01TA', 'ST059Q03TA', 'USESCH', 'IC150Q03HA', 
               'ESCS', 'ST160Q02IA')]

cleaned_data <- na.omit(data)
#transform dependent variable to its log
cleaned_data$log_scie_score <- log(cleaned_data$scie_score)

##2. summary statistics
library(dplyr)
cleaned_data %>%
  select(log_scie_score, ST059Q03TA, USESCH, ESCS) %>%
  summary()

##calculate sd
selected_columns <- c('log_scie_score', 'ST059Q03TA', 'USESCH', 'ESCS')
standard_deviations <- sapply(cleaned_data[selected_columns], sd)
print(standard_deviations)

#Factor variables: ST004D01T, ST011Q04TA, ST034Q01TA, IC150Q03HA, ST160Q02IA
ST004D01T_counts <- table(cleaned_data$ST004D01T)
print(ST004D01T_counts)

ST011Q04TA_counts <- table(cleaned_data$ST011Q04TA)
print(ST011Q04TA_counts)

ST034Q01TA_counts <- table(cleaned_data$ST034Q01TA)
print(ST034Q01TA_counts)

IC150Q03HA_counts <- table(cleaned_data$IC150Q03HA)
print(IC150Q03HA_counts)

ST160Q02IA_counts <- table(cleaned_data$ST160Q02IA)
print(ST160Q02IA_counts)

attach(cleaned_data)

#converts quantitative variables into qualitative variables

ST004D01T <- as.factor(ST004D01T)
ST011Q04TA <- as.factor(ST011Q04TA)
ST034Q01TA <- as.factor(ST034Q01TA)
IC150Q03HA <- as.factor(IC150Q03HA)
ST160Q02IA <- as.factor(ST160Q02IA)

##simple linear regression

lm.fit <- lm(log_scie_score ~ ST004D01T + ST011Q04TA + ST034Q01TA +
               ST059Q03TA + USESCH + IC150Q03HA + ESCS + ST160Q02IA)
summary(lm.fit)

# Model diagnostics

#1. Check for mean 0, homoschedasticity (constant variance) and independence

## Mean == 0?
###Option 1:
plot(lm.fit$fitted.values, lm.fit$residuals, xlab = "Fitted values",
     ylab = "Residuals", main = "Plot of the residuals")

abline(h=0, lwd = 2)
round(mean(lm.fit$residuals), 3)

###Option 2:
residuals <- lm.fit$residuals # Calculate the residuals
t.test(residuals, mu = 0) # Test if the mean of residuals is significantly different from zero
#=> both p-value = 1 and the confidence interval, which spans from -1.570725 to 1.570725
# suggest that the mean of the residuals is not significantly different from zero
#=> aligns with the assumption of mean zero in linear regression.

#2. To check for the independence, you have to see if residuals are randomly
#distributed in the plot or if they exhibit any pattern

sigma.hat =summary(lm.fit)$sigma
sigma.hat #standard deviation of the residuals

abline(h=sigma.hat, col="red", lty=2, lwd = 2)
abline(h=-sigma.hat, col="red", lty=2, lwd = 2)
abline(h=2*sigma.hat, col="green", lty=2, lwd = 2)
abline(h=-2*sigma.hat, col="green", lty=2, lwd = 2)
legend(c(-20,20), legend=c("+/- residual standard deviation", "+/- 2 residual 
                           #standard deviation"),
       col = c("red", "green"), pch=20, bty="n",
       x.intersp = 0.5, y.intersp = 0.5)

##Option 2:
# Install and load the 'lmtest' package for the 'dwtest' function
library(lmtest)

# Perform the Durbin-Watson test
dw_test <- dwtest(lm.fit)

# Print the test results
print(dw_test)

#=> Conclusion: p-value < 2.2e-16 suggests strong evidence against 
#the null hypothesis of no autocorrelation

#3. Check for homoskedasticity: How many points lie outside the two green lines.
#Check for homoskedasticity: can use the test (instead):
library(lmtest)
bptest(lm.fit)

#Ho: the variance of the residuals is constant (homoscedasticity)
#H1: the variance of the residuals is not constant (heteroscedasticity)
#=> The p-value is very small (2.2e-16 < 0.05) => reject the null hypothesis 
#which indicates that there is enough evidence to suggest that the variance 
#of the residuals is not constant

# Using sandwich package for robust standard errors
library(sandwich)

robust_model <- coeftest(lm.fit, vcov = vcovHC(lm.fit, type = "HC1"))
robust_model

#4. Check for normality assumption

hist(lm.fit$residuals, breaks=20, main = "Histogram of the residuals",
     xlab="Residuals", freq=F,ylab="Density", col = "lightblue")

# Observed (empirical) quantiles of the residuals
qqnorm(lm.fit$residuals,cex=0.6, main= "", xlab="Quantiles", ylab="Residuals")
qqline(lm.fit$residuals, lwd = 2, col = "red") 

emp = quantile(lm.fit$residuals, probs = seq(0,1, by=0.1))
norm = qnorm(seq(0,1, by=0.1))

points(norm, emp, col = "green", pch = 20)

#5. Check for multicollinearity
library(car)
# Calculate VIF for each predictor variable
vif_values <- vif(lm.fit)

# Print the VIF values
print(vif_values)

###VIF values around 1 indicate low multicollinearity.
##VIF values between 1 and 5 suggest moderate multicollinearity.
##VIF values above 5 suggest high multicollinearity.


# Install and load the required package
library(lmtest)

# Perform Ramsey RESET test
ramsey_test <- resettest(lm.fit, power = 3)  # You can adjust the 'power' argument

# Print the test results
print(ramsey_test)


#II. Export results
library(broom)
library(openxlsx)

# Tidy the model results
tidy_results <- tidy(robust_model)

# Define the path where you want to save the Excel file
excel_file_path <- "04_report\\linear_regression_results.xlsx"

# Write the tidy results to an Excel file
write.xlsx(tidy_results, excel_file_path, sheetName = "LinearRegressionResults", rowNames = FALSE)
