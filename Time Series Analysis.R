rm(list=ls())     # Clean memory
graphics.off()    # Close graphs
cat("\014")       # Clear Console

library("vars")
library("urca")
library("aTSA")
library("beepr")
library("plotrix")
library("tseries")

invest <- read.table("Real Gross Private Domestic Investment.csv",sep=";", header=TRUE)
cons <- read.table("Real Personal Consumption Expenditures.csv",sep=";", header=TRUE)
gnp <- read.table("Real gross national product per capita.csv",sep=";", header=TRUE)

dates <- as.Date(invest[,1], "%Y-%m-%d")

tsInvest <- as.ts(invest[,2], start = invest[1,1], end= invest[121,1], frequency = 4)
tsCons <- as.ts(cons[,2], start = invest[1,1], end= invest[121,1], frequency = 4)
tsGnp <- as.ts(gnp[,2], start = invest[1,1], end= invest[121,1], frequency = 4)

###### Order of integration ########

#### Plots of time series & first differences #####
plot(dates, tsInvest, type = "l", xlab = "Date", ylab = "Investment")
plot(dates[-1], diff(tsInvest), type = "l", xlab = "Date", ylab = "First difference of Investment")

plot(dates, tsCons, type = "l", xlab = "Date", ylab = "Consumption")
plot(dates[-1], diff(tsCons), type = "l", xlab = "Date", ylab = "First difference of consumption")

plot(dates, tsGnp, type = "l", xlab = "Date", ylab = "GNP")
plot(dates[-1], diff(tsGnp), type = "l", xlab = "Date", ylab = "First difference of GNP")

#### Plot of Log values #####
plot(dates, log(tsInvest), type = "l", xlab = "Date", ylab = "Investment")
plot(dates[-1], diff(log(tsInvest)), type = "l", xlab = "Date", ylab = "First difference of Investment") 

plot(dates, log(tsCons), type = "l", xlab = "Date", ylab = "Consumption")
plot(dates[-1], diff(log(tsCons)), type = "l", xlab = "Date", ylab = "First difference of consumption")

plot(dates, log(tsGnp), type = "l", xlab = "Date", ylab = "GNP")
plot(dates[-1], diff(log(tsGnp)), type = "l", xlab = "Date", ylab = "First difference of GNP")

##### ACFs #######

##### Investment #####
acf(tsInvest, main = "ACF Investment") # Slowly decreasing
acf(diff(tsInvest), main = "ACF Difference of Investment") #Two significant lags
#potentially log
acf(log(tsInvest), main = "ACF Logged Investment")
acf(diff(log(tsInvest)), main = "ACF Difference of Logged Investment")

##### Consumption ####
acf(tsCons, main = "ACF Consumption") # Slowly decreasing
acf(diff(tsCons), main = "ACF Difference of Consumption") #Five significant lags
#potentially log
acf(log(tsCons), main = "ACF Logged Consumption") # slowly decreasing
acf(diff(log(tsCons)), main = "ACF Difference of Logged Consumption") # 6 significant lags

#### GNP ####
acf(tsGnp, main = "ACF GNP") # Slowly decreasing
acf(diff(tsGnp), main = "ACF Difference of GNP") #Two significant lags
#potentially log
acf(log(tsGnp), main = "ACF Logged GNP") # slowly decreasing 
acf(diff(log(tsGnp)), main = "ACF Difference of Logged GNP") # two significant lags

# Dickey Fuller Test
###### Investment ######
aTSA::adf.test(log(tsInvest), nlag = 1)
aTSA::adf.test(diff(log(tsInvest)), nlag = 1)

###### Consumption #######
aTSA::adf.test(log(tsCons), nlag = 1)
aTSA::adf.test(diff(log(tsCons)), nlag = 1)

###### GNP ######
aTSA::adf.test(log(tsGnp), nlag = 1)
aTSA::adf.test(diff(log(tsGnp)), nlag = 1)

### all are difference stationary ###

# Augmented Dickey Fuller Test
###### Investment ######
aTSA::adf.test(log(tsInvest))
aTSA::adf.test(diff(log(tsInvest)))

###### Consumption #######
aTSA::adf.test(log(tsCons))
aTSA::adf.test(diff(log(tsCons)))

###### GNP ######
aTSA::adf.test(log(tsGnp))
aTSA::adf.test(diff(log(tsGnp)))

### all are difference stationary ###

# Phillips-Perron Test
PP.test(log(tsInvest))
PP.test(log(tsCons))
PP.test(log(tsGnp))

## Do not reject null of unit root, so we take first differences p-value >= 0.66 for all series

PP.test(diff(log(tsInvest)))
PP.test(diff(log(tsCons)))
PP.test(diff(log(tsGnp))) 

## Null of unit root rejected, p-value = 0.01 for all logged series in differences

inv <- diff(log(tsInvest))
con <- diff(log(tsCons))
gn <- diff(log(tsGnp))

# VAR
mat <- cbind(inv, con, gn)
series <- ts(mat)
#VAR(y = series , type = "const")

# IC for which lag length
VARselect(y = series, lag.max=10, type = "const")
# Autocorrelation for which lag length
#residuals <- residuals(VAR(y = series , type = "const", lag.max = 1))
# Sequential test on the nullity of coefficient matrix
p <- 9
optimal.lag.lr <- 0
optimal.lag.wald <- 0
for(i in p:1){
  ## LR
  trash.unrest <- summary(VAR(y = series , type = "const", p = i))
  det.un <- det(trash.unrest$covres)
  trash.re <- summary(VAR(y = series , type = "const", p = i-1))
  det.re <- det(trash.re$covres)
  t <- length(series)-i
  lr <- t*(log(det.re)-log(det.un))
  critical.value <- qchisq(0.95, df = 9) # 3 I(1) variables, so n^2 is 9
  if(lr > critical.value){
    optimal.lag.lr = i
    break
  } 
  
  ## Wald
  x <- as.matrix(VAR(y = series , type = "const", p = i)$datamat[,-1:-3])
  q <- (t(x)%*%x)/t
  c1 <- matrix(0,nrow =9, ncol = 3+9*(i-1))
  c2 <- diag(9)
  c <- cbind(c1,c2)
  c.small = rep(0, 9)
  omega <- trash.unrest$covres
  pi <- c(trash.unrest$varresult$inv$coefficients[,1],trash.unrest$varresult$con$coefficients[,1],trash.unrest$varresult$gn$coefficients[,1])
  
  #wald <- t*t(c%*%pi-c.small)%*%solve((c%*%(omega%x%solve(q))%*%t(c)))%*%(c%*%pi-c.small)
  #if(wald > critical.value){
   # optimal.lag.wald = i
    #break
  #}
}

#library(VAR.etp)
#VAR.Wald(series, restrict = )

#### optimal length: 3
###ACF plots in residuals
var3 <- summary(VAR(y = series , type = "const", p = 3))
varmodel3 <- VAR(y = series , type = "const", p = 3)
acf(var3$varresult$inv$residuals, main = "ACF plot Residuals of Investment Equation")
acf(var3$varresult$con$residuals,main = "ACF plot Residuals of Consumption Equation")
acf(var3$varresult$gn$residuals,main = "ACF plot Residuals of GNP Equation")

## stability of F (first write in companion form)
f <- rbind(var3$varresult$inv$coefficients[1:9,1],var3$varresult$con$coefficients[1:9,1],var3$varresult$gn$coefficients[1:9,1])
ident <- diag(1,3)
ident0 <- diag(0,3)
Fmatrix <- rbind(f,cbind(ident, ident0,ident0), cbind(ident0,ident,ident0))

plot(1, type="n", xlab="", ylab="", xlim=c(-1, 1), ylim=c(-1, 1), asp=1)
draw.circle(0,0,1)
lines(x = eigen(Fmatrix)$values, type = "p")

#check for stability with R function
vars::roots(VAR(y = series , type = "const", p = 3), modulus = FALSE)

### Test for weak stationarity of the residuals 
# -> together with stability, implies stationarity of y_t

var1 <- summary(VAR(y = series , type = "const", p = 1))
acf(var1$varresult$inv$residuals, main = "ACF plot Residuals of Investment Equation")
acf(var1$varresult$con$residuals,main = "ACF plot Residuals of Consumption Equation")
acf(var1$varresult$gn$residuals,main = "ACF plot Residuals of GNP Equation")
#beep(0)

F2matrix <- rbind(var1$varresult$inv$coefficients[1:3,1],var1$varresult$con$coefficients[1:3,1],var1$varresult$gn$coefficients[1:3,1])
plot(1, type="n", xlab="", ylab="", xlim=c(-1, 1), ylim=c(-1, 1), asp=1)
draw.circle(0,0,1)
points(x = eigen(F2matrix)$values, y = c(0,0,0))
vars::roots(VAR(y = series , type = "const", p = 1), modulus = FALSE)

##### Plot true vs model

plot(x = ts(mat[-1:-3,1], start = c(1988,1), end = c(2018,1), frequency = 4), type = "l", main = "Estimated vs Actual for Investment", xlab = "Date", ylab = "Value")
lines(ts(varmodel3$varresult$inv$fitted.values, start = c(1988,1), end = c(2018,1), frequency = 4), col = 2)

plot(x = ts(mat[-1:-3,2], start = c(1988,1), end = c(2018,1), frequency = 4), type = "l", main = "Estimated vs Actual for Consumption", xlab = "Date", ylab = "Value")
lines(ts(varmodel3$varresult$con$fitted.values, start = c(1988,1), end = c(2018,1), frequency = 4), col = 2)

plot(x = ts(mat[-1:-3,3], start = c(1988,1), end = c(2018,1), frequency = 4), type = "l", main = "Estimated vs Actual for GNP", xlab = "Date", ylab = "Value")
lines(ts(varmodel3$varresult$gn$fitted.values, start = c(1988,1), end = c(2018,1), frequency = 4), col = 2)


###########################################
####### Granger causality VAR #############
###########################################
restr.var.3nognp <- summary(VAR(y = ts(cbind(inv,con)), type = "const", p = 3))
restr.var.3noinv <- summary(VAR(y = ts(cbind(con, gn)), type = "const", p = 3))
restr.var.3nocon <- summary(VAR(y = ts(cbind(inv,gn)), type = "const", p = 3))

sse.full.inv <- sum((var3$varresult$inv$residuals)^2)
sse.full.con <- sum((var3$varresult$con$residuals)^2)
sse.full.gnp <- sum((var3$varresult$gn$residuals)^2)

sse.rest.nognp.inv <- sum((restr.var.3nognp$varresult$inv$residuals)^2)
sse.rest.nognp.con <- sum((restr.var.3nognp$varresult$con$residuals)^2)
sse.rest.noinv.con <- sum((restr.var.3noinv$varresult$con$residuals)^2)
sse.rest.noinv.gnp <- sum((restr.var.3noinv$varresult$gn$residuals)^2)
sse.rest.nocon.inv <- sum((restr.var.3nocon$varresult$inv$residuals)^2)
sse.rest.nocon.gnp <- sum((restr.var.3nocon$varresult$gn$residuals)^2)

critic.val.f <- qf(.95, df1 = 3, df2 = 110)
print(critic.val.f)

partial.f.gnp.improves.inv <- ((sse.rest.nognp.inv - sse.full.inv)/(3))/((sse.full.inv)/110)
print(partial.f.gnp.improves.inv)
df(partial.f.gnp.improves.inv, df1 = 3, df2 = 110)
partial.f.gnp.improves.con <- ((sse.rest.nognp.con - sse.full.con)/(3))/((sse.full.con)/110)
print(partial.f.gnp.improves.con)
df(partial.f.gnp.improves.con, df1 = 3, df2 = 110)
partial.f.con.improves.inv <- ((sse.rest.nocon.inv - sse.full.inv)/(3))/((sse.full.inv)/110)
print(partial.f.con.improves.inv)
df(partial.f.con.improves.inv, df1 = 3, df2 = 110)
partial.f.con.improves.gnp <- ((sse.rest.nocon.gnp - sse.full.gnp)/(3))/((sse.full.gnp)/110)
print(partial.f.con.improves.gnp)
df(partial.f.con.improves.gnp, df1 = 3, df2 = 110)
partial.f.inv.improves.con <- ((sse.rest.noinv.con - sse.full.con)/(3))/((sse.full.con)/110)
print(partial.f.inv.improves.con)
df(partial.f.inv.improves.con, df1 = 3, df2 = 110)
partial.f.inv.improves.gnp <- ((sse.rest.noinv.gnp - sse.full.gnp)/(3))/((sse.full.gnp)/110)
print(partial.f.inv.improves.gnp)
df(partial.f.inv.improves.gnp, df1 = 3, df2 = 110)

###########################################
###########################################
###########################################

## Impulse Response Analysis
# Look at how a shock in income (GNP) affects savings (investment) and spending (consumption)
irf_gn <- irf(varmodel3, impulse = "gn", cumulative = TRUE)
# If we assume a closed economy, all changes in savings + spending = changes in income (CLOSE ENOUGH!) "cool finding" - (Verkade, 2018)
plot(irf_gn, main = "Impulse Response shocking GNP")
irf(varmodel3, impulse = "gn", n.ahead = 8)

# We look at shock in investment
irf_in <- irf(varmodel3, impulse = "inv", cumulative = TRUE, n.ahead = 8)
plot(irf(varmodel3, impulse = "inv", cumulative = TRUE, n.ahead = 8),main = "Impulse Response shocking Investment")

# We look at shock in consumption
irf(varmodel3, impulse = "con", cumulative = TRUE, n.ahead = 8)
plot(irf(varmodel3, impulse = "con", cumulative = TRUE, n.ahead = 8), main = "Impulse Response shocking Consumption")

# Forecasting the var model using predict
predict.varmodel3 <- predict(varmodel3, n.ahead = 8, ci = 0.95)
fanchart(predict.varmodel3, colors = c("orange", "blue", "red", "green", "yellow"), plot.type = "multiple")
plot(predict.varmodel3)

# Forecasting how it's asked in the assignment
h <- 8
# Deleting the last 8 observations
con.in = con[1:(120-h)]
inv.in = inv[1:(120-h)]
gn.in = gn[1:(120-h)]
# Re-estimating the model
# VAR
mat.in <- cbind(inv.in, con.in, gn.in)
series.in <- ts(mat.in)
varmodel3.in <- VAR(y = series.in , type = "const", p = 3)
# Performing one-step ahead forecasts
# Using the predict function
predict.varmodel3.in <- predict(varmodel3.in, n.ahead = h, ci = 0.95)
plot(x = inv[113:120], type = "l", main = "Forecast vs Actual for Investment", xlab = "Forecast Horizon", ylab = "Value")
lines(predict.varmodel3.in$fcst$inv.in[,1], col = 2)

plot(x = con[113:120], type = "l", main = "Forecast vs Actual for Consumption", xlab = "Forecast Horizon", ylab = "Value")
lines(predict.varmodel3.in$fcst$con.in[,1], col = 2)

plot(x = gn[113:120], type = "l", main = "Forecast vs Actual for GNP", xlab = "Forecast Horizon", ylab = "Value")
lines(predict.varmodel3.in$fcst$gn.in[,1], col = 2)

# Using our own stuff
# Using for loops
inv.in.coeffs <- varmodel3.in$varresult$inv.in$coefficients
con.in.coeffs <- varmodel3.in$varresult$con.in$coefficients
gn.in.coeffs <- varmodel3.in$varresult$gn.in$coefficients

in.coeffs <- rbind(inv.in.coeffs, con.in.coeffs, gn.in.coeffs)

in.coeffs.first.lag <- in.coeffs[,1:3]
in.coeffs.second.lag <- in.coeffs[,4:6]
in.coeffs.third.lag <- in.coeffs[,7:9]
constant <- in.coeffs[,10]

for (i in 1:h) {
  dummy <- in.coeffs.first.lag %*% mat.in[length(mat.in[,1]),] + in.coeffs.second.lag %*% mat.in[length(mat.in[,1])-1,] + in.coeffs.third.lag %*% mat.in[length(mat.in[,1])-2,] + constant
  mat.in <- rbind(mat.in,t(dummy))
}

plot(x = inv[113:120], type = "l", main = "Forecast vs Actual for Investment", xlab = "Forecast Horizon", ylab = "Value")
lines(mat.in[113:120,1], col = 2)

plot(x = con[113:120], type = "l", main = "Forecast vs Actual for Consumption", xlab = "Forecast Horizon", ylab = "Value")
lines(mat.in[113:120,2], col = 2)

plot(x = gn[113:120], type = "l", main = "Forecast vs Actual for GNP", xlab = "Forecast Horizon", ylab = "Value")
lines(mat.in[113:120,3], col = 2)

# In sample prediction vs outcome
plot(x = inv[4:120], type = "l")
lines(c(varmodel3.in$varresult$inv.in$fitted.values, mat.in[113:120,1]), col = 2)

plot(x = con[4:120], type = "l")
lines(c(varmodel3.in$varresult$con.in$fitted.values, mat.in[113:120,2]), col = 2)

plot(x = gn[4:120], type = "l")
lines(c(varmodel3.in$varresult$gn.in$fitted.values, mat.in[113:120,3]), col = 2)

#### Now onto Cointegration - Engle-Granger approach

loggedSeries <- matrix(c(log(tsInvest), log(tsCons), log(tsGnp)), nrow =121, ncol = 3)
colnames(loggedSeries) <- c("investment", "consumption", "gnp")

####################################################
#### Change every variable to be dependent variable 
####################################################
logged.regr <- lm(loggedSeries[,2] ~ loggedSeries[,1]+loggedSeries[,3])

#normalSeries <- matrix(c(tsInvest, tsCons, tsGnp), nrow =121, ncol = 3)
#colnames(normalSeries) <- c("investment", "consumption", "gnp")
#normal.regr <- lm(normalSeries[,3] ~ normalSeries[,1:2])
## Get residuals first
logged.resid <- residuals(logged.regr)
#normal.resid <- residuals(normal.regr)

## Phillips-Ouliaris test (po test package in r needs model, not residuals so this is fine)
# H0: No cointegration
po.test(loggedSeries)
# p-value >= 0.15 so do not reject null of no cointegration

## Phillips-Perron test
# H0: Unit root
PP.test(diff(diff(logged.resid)))
PP.test(diff(logged.resid))
PP.test(logged.resid)
# p-value = 0.1669 so do not reject null of unit root

## DF test
aTSA::adf.test(diff(diff(logged.resid)), nlag = 1)
aTSA::adf.test(diff(logged.resid), nlag = 1)
aTSA::adf.test(logged.resid, nlag = 1)
# p-value's differ, but do not reject if drift and trend included

## ADF test
aTSA::adf.test(diff(diff(logged.resid)))
aTSA::adf.test(diff(logged.resid))
aTSA::adf.test(logged.resid)
# same story as above

### Johansen's Approach  (https://stats.stackexchange.com/questions/186208/johansen-test-for-cointegration)

# Use VARselect to determine lag order p of VAR in levels
VARselect(y = loggedSeries, lag.max=10, type = "const")
# VARselect selects 2 lags for all IC, hence p=2 

## Trace test
trace_test.diffed <- ca.jo(loggedSeries, type="trace", K=2, ecdet = "const", spec="transitory")
summary(trace_test.diffed)
# 2 cointegrating relations by the Trace test (p < 0.01)!!! 

# Max Eigen test
eigen_test.diffed <- ca.jo(loggedSeries, type="eigen", K=2, ecdet = "const", spec="transitory")
summary(eigen_test.diffed)
# 2 cointegrating relations by the Maximum Eigenvalue Test (p < 0.01)!!! 

### So we have two cointegrating relations, now we need to represent our model using VECM
fit.vecm2 <- cajorls(trace_test.diffed, r = 2)
fit.vecm2

fit.vecm2.eigen <- cajorls(eigen_test.diffed, r= 2)
fit.vecm2.eigen
# Need to test for Granger-causality again 
## VERY MUCH UNSURE WTF IS GOING ON FUCK IT WE CAN DO IT IN EVIEWS
unres.vecm2 <- cajools(trace_test.diffed, r = 2)
vecm2.sse.full <- sum(unres.vecm2$residuals^2)
vecm2.sse.res <- sum((fit.vecm2$rlm$residuals)^2)


vecm2.crit.val <- qf(.95, df1 = 2, df2 = 111)
print(vecm2.crit.val)

vecm2.partial.f <- ((vecm2.sse.res - vecm2.sse.full)/(2))/((vecm2.sse.full)/111)
print(vecm2.partial.f)
df(vecm2.partial.f, df1 = 2, df2 = 111)

