#
# Regression analysis for comorbidities study
#

library('stargazer')
library('MASS')
library('arm')
library('foreign')
library('plm')

##################################################
# LOAD AND CLEAN DATA
##################################################
data <- read.csv('data/regression_data.csv')
cost_data <- read.csv('data/cost_estimation.csv')
data <- merge(data, cost_data, by.x=c("member_no", "year"), by.y=c("membno", "year"))

data$zipcode[data$zipcode==''] <- NA
data$zipcode[data$zipcode=='<NA>'] <- NA
data$zipcode[data$zipcode=='NA'] <- NA
data$zipcode <- factor(data$zipcode)

data$sex[data$sex==''] <- NA
data$sex[data$sex=='<NA>'] <- NA
data$sex[data$sex=='NA'] <- NA
data$sex <- factor(data$sex)

# make the number of comorbidities a factor
no_comorbidities.f <- factor(data$no_comorbidities)
data <- cbind(no_comorbidities.f, data)

##  make dummy variables factors
for (column in colnames(data)[c(c(14:83), c(87:131))]) {
    data[column] <- factor(data[[column]])
}

# make panel data
panel_data <- data[data$count>1,]

####################################################
# CREATE FORMULAS
####################################################

formula.base <- as.formula('log_costs ~ no_comorbidities')
formula.additive <- as.formula('log_costs ~ no_comorbidities + age + sex')
formula.base.factor <- as.formula('log_costs ~ no_comorbidities.f')

# adding in cost center and comorbidity type variables 
ccc <- paste('~ . + ', paste(colnames(data)[c(14:83)], collapse = "+"))
ccc <- paste(paste(ccc, '+'),
                 paste(colnames(data)[c(87:131)], collapse = "+"))

formula.ccc <- update.formula(formula.base, as.formula(ccc))


# adding interaction terms for age and sex dummies
# age
ccc_age <- paste('~ . + age:', paste(colnames(data)[c(14:83)], collapse = "+age:"))
ccc_age <- paste(paste(ccc_age, '+'),
                 paste(colnames(data)[c(87:131)], collapse = "+"))

formula.ccc_age <- update.formula(formula.base, as.formula(ccc_age))


###################################################
# MODELLING
###################################################

# base model
model.base <- lm(formula.base, data=data)
model.base_pd <- lm(formula.base, data=panel_data)

# factor model :: there's a nice plot out of this one
model.base.factor <- lm(formula.base.factor, data=data)

# within estimation
model.base.within <- plm(formula.base, 
                         data=data, 
                         index=c("member_no", "year"), 
                         model="within")
model.ccc.within <- plm(formula.ccc, 
                        data=data, 
                        index=c("member_no", "year"), 
                        model="within")
model.ccc_age.within <- plm(formula.ccc_age, 
                            data=data, 
                            index=c("member_no", "year"), 
                            model="within")

###################################################
# OUTPUT
###################################################

png('../output/coefplot.png', width=10, height=15, units='cm', res=300)
coefplot(model.base.factor, vertical = FALSE, xlab = "Number of comorbidities",
         ylab = "Coefficient size",
         c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12",
           "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23",
           "24", "25", "26", "27", "28", "29", "30", "31", "32", "33", "34",
           "36"))

stargazer(model.base, model.base_pd, model.base.within,
          title = "Comparing bias amongst models", keep = "no_comorbidities",
          covariate.labels = "Number of comorbidities",
          dep.var.labels = c("OLS (full sample)", "OLS (panel sample)", "Within estimator"), 
          type = "html",
          #notes = "Comparing results across models and datasets. The full sample uses all patients over the three year period, the panel sample uses only patients that visited the hospital in at least two years. The within estimator uses the panel sample"
          # notes.append=TRUE,
          out = '../output/regression_bias.htm')

stargazer(model.base.within, model.ccc.within, 
          model.ccc_age.within,
          title = "Within estimator results", keep = "no_comorbidities",
          covariate.labels = "Number of comorbidities",
          dep.var.labels = c("Base", "Dummies"), type = "html",
          out = '../output/regression_within.htm')

# EOF
