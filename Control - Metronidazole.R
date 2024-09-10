# Control Antibiotic - METRONIDAZOLE

# Load Packages
library(fpp2)
library(ggplot2)
library(tidyverse)

# Load Main Sheets (Data for both and Precovid also)
precovidmet <- read.csv("precovidmet.csv", as.is = T)
metboth <- read.csv("metboth.csv", as.is = T)
view(precovidmet)
view(metboth)

# Convert needed data columns to time series 
met <- ts(precovidmet, start = c(2016, 5), frequency = 12)
premet <- ts(precovidmet$Metronidazole, start = c(2016, 5), frequency = 12)
patientsmet <- ts(precovidmet$Patientsmillion, start = c(2016, 5), frequency = 12)
bothmet <- ts(metboth$Metronidazole, start = c(2016, 5), frequency = 12)

# Plot control prescription
autoplot(premet, main = "Pre-Covid Metronidazole prescription (UK) time series")
autoplot(bothmet, main = "Metronidazole prescription (UK) time series")

#lag plot and autocorrelation
premetlag <- window(premet, start = c(2016, 5))
gglagplot(premetlag)
ggAcf(premetlag, lag = 44)

premetdiff <- diff(premet)
autoplot(premetdiff, main = "Pre-Covid Metronidazole prescription (UK) monthly differences")
ggAcf(premetdiff, lag = 44)

bothmetlag <- window(bothmet, start = c(2016, 5))
gglagplot(bothmetlag)
ggAcf(bothmetlag, lag = 61)

bothmetdiff <- diff(bothmet)
autoplot(bothmetdiff, main = "Metronidazole prescription (UK) monthly differences")
ggAcf(bothmetdiff, lag = 61)

# Box-Cox Transformation
(lambdaprecovidmet <- BoxCox.lambda(premetlag))
autoplot(BoxCox(premetlag, lambdaprecovidmet))

# Residual checks
preresmet <- residuals(naive(premet))
autoplot(preresmet) + xlab(" ") + ylab(" ") +
  ggtitle("Residuals from naive method")
gghistogram(preresmet) + ggtitle("Histogram of residuals")
ggAcf(preresmet) + ggtitle("ACF of residuals")

# Portmanteau Tests
Box.test(preresmet, lag = 10, fitdf = 0)
Box.test(preresmet, lag = 10, fitdf = 0, type = "Lj" )
checkresiduals(naive(premet))

#Basic principle forecast
autoplot(premet) +
  autolayer(meanf(premet, h = 17),
            series = "Mean", PI = FALSE) +
  autolayer(naive(premet, h = 17),
            series = "Naive", PI = FALSE) +
  autolayer(snaive(premet, h = 17),
            series = "Seasonal Naive", PI = FALSE) +
  autolayer(snaive(premet)) +
  ggtitle("Forecasts for Metronidazole Prescription") +
  xlab("Year") + ylab("Metronidazole Prescription") +
  guides(color=guide_legend(title="Forecast"))

# Variable Plot Comparison
autoplot(met[, c("Metronidazole", "Patientsmillion")])

#Correlation Plot
GGally::ggpairs(as.data.frame(met[, c("Metronidazole", "Patientsmillion")]))

met %>% 
  as.data.frame() %>% 
  ggplot(aes(x = Metronidazole, y = Patientsmillion)) +
  ylab("Patients (millions)") +
  xlab("Metronidazole Prescription") +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE)
tslm(Patientsmillion ~ Metronidazole, data = met)

# Fit Linear regression plot for both variables
fit.met <- tslm(
  Metronidazole ~ patientsmet,
  data=met
)

autoplot(met[, "Metronidazole"], series = "Data") +
  autolayer(fitted(fit.met), series = "Fitted")

cbind(Data = met[, "Metronidazole"],
      Fitted = fitted(fit.met)) %>% 
  as.data.frame() %>% 
  ggplot(aes(x=Data, y=Fitted)) +
  geom_point() +
  geom_abline(intercept=0, slope = 1)
checkresiduals(fit.met)

df <- as.data.frame(met)
df[, "Residuals"] <- as.numeric(residuals((fit.met)))
ggplot(df, aes(premet, Residuals)) +
  geom_point()
cbind(Fitted= fitted(fit.met),
      Residuals= residuals(fit.met)) %>% 
  as.data.frame() %>% 
  ggplot(aes(x=Fitted, y=Residuals)) + geom_point()
summary(fit.met)
CV(fit.met)

# Time series decomposition plot
metfc <- window(premet, start = c(2016, 5))
fit.metfc <- tslm(premet ~ trend + season)
metfcast <- forecast(fit.metfc, h = 17)
autoplot(metfcast, PI = FALSE) +
  autolayer(bothmet, colour = FALSE) +
  ggtitle("Forecast of Metronidazole prescription using Time series decomposition") +
  xlab("Year") + ylab(" ")

premet %>% decompose(type="multiplicative") %>% 
  autoplot() + xlab("Year")

#Export New Prediction
write.csv(metfcast, file = "metfc.csv")

# Load change sheet for differences
metfc <- read.csv("metfc2.csv", as.is = T)

percentchange <- ts(metfc$PercentageChange, start = 2020, frequency = 12)
percentdecrease <- ts(metfc$PercentageDecrease, start = 2020, frequency = 12)
percentincrease <- ts(metfc$PercentageIncrease, start = 2020, frequency = 12)
metchange <- ts(metfc, start = 2020, frequency = 12)

autoplot(percentchange, series = "Point Forecast") +
  autolayer(percentdecrease, series = "95% Interval (Low)") +
  autolayer(percentincrease, series = "95% Interval (High)") +
  geom_hline(aes(yintercept = 0))
