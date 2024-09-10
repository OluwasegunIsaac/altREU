library(fpp2)
library(ggplot2)
library(tidyverse)

aziboth <- read.csv("Azidata/Aziboth.csv", as.is = T)
view(aziboth)

  
azibothdata <- ts(aziboth, start = c(2016, 5), frequency = 12)
Per1000Patients <- ts(aziboth$Per1000Patients, start = c(2016, 5), frequency = 12)
Patients <- ts(aziboth$Patients, start = c(2016, 5), frequency = 12)
Patientsmillion <- ts(aziboth$Patientsmillion, start = c(2016, 5), frequency = 12)
Azithromycin2 <- ts(aziboth$Azithromycin, start = c(2016, 5), frequency = 12)
STARPU <- ts(aziboth$STARPU, start = c(2016, 5), frequency = 12)
Per1000STARPU <- ts(aziboth$Per1000STARPU, start = c(2016, 5), frequency = 12)
Raw <- ts(aziboth$Rawazipatients, start = c(2016, 5), frequency = 12)

#Insert abline v to aviboth
plot(Azithromycin2, main = "Azithromycin dispensation across all CCGs in the UK")
abline(v = 2020, col = "red", lty = 2, lwd = 1.5)
abline(h = 71381, col = "red", lty = 2, lwd = 1.5)

#Plot Populationpermillion with time
Azipatientsm <- ts(aziboth$Patientsmillion, start = c(2016, 5), frequency = 12)

autoplot(Patientsmillion, main = "Patient Population (millions) across all CCGs in the UK") +
  ylab("Patient population(millions)")

autoplot(STARPU, main = "STARPU across all CCGs in the UK") +
  ylab(" ")

autoplot(Raw, main = "STARPU across all CCGs in the UK") +
  ylab(" ")

autoplot(Per1000STARPU, main = "STARPU across all CCGs in the UK") +
  ylab(" ")

autoplot(azibothdata[, c("Azithromycin", "Per1000Patients", "Patientsmillion")], facet = TRUE) +
  ylab("Patient population(millions), Azi per 100 Patients, Azi Prescription")

plot(Azipatientsm, main = "Patients across all CCGs in the UK")
abline(v = 2020, col = "red", lty = 2, lwd = 1.5)
abline(h = 60.32292, col = "red", lty = 2, lwd = 1.5)

#Seasonal Plot
autoplot(Azithromycin) +
  ggtitle("Azithromycin prescription across UK") +
  ylab("Prescription") +
  xlab("Year")

autoplot(Per1000Patients) +
  ggtitle("Azithromycin prescription per 1000 patients across UK") +
  ylab("Prescription per 1000 Patients") +
  xlab("Year")

ggseasonplot(Azithromycin, year.labels = TRUE, year.labels.left = TRUE) +
  ylab("Prescription") +
  ggtitle("Seasonal plot: Azithromycin Prescription")
ggseasonplot(Azithromycin, polar = TRUE) +
  ylab("Prescription") +
  ggtitle("Seasonal plot: Azithromycin Prescription")

ggseasonplot(Per1000Patients, year.labels = TRUE, year.labels.left = TRUE) +
  ylab("Prescription per 1000 Patients") +
  ggtitle("Seasonal plot: Azithromycin Prescription per 1000 Patients")
ggseasonplot(Per1000Patients, polar = TRUE) +
  ylab("Prescription per 1000 Patients") +
  ggtitle("Seasonal plot: Azithromycin Prescription per 1000 Patients")

ggseasonplot(Azipatientsm, year.labels = TRUE, year.labels.left = TRUE) +
  ylab("Patients Population in millions") +
  ggtitle("Seasonal plot: Patient Population, UK")
ggseasonplot(Azipatientsm, polar = TRUE) +
  ylab("Patients Population in millions") +
  ggtitle("Seasonal plot: Patient Population, UK")

#Correlation Plot
GGally::ggpairs(as.data.frame(aziboth[, c("Azithromycin", "Patientsmillion", "STARPUmillion")]))

#Lag plots
Azilag <- window(Azithromycin, start = c(2016, 5))
gglagplot(Azilag)
ggAcf(Azilag, lag = 61)
frequency(Azithromycin)
Patientslag <- window(Patients, start = c(2016, 5))
gglagplot(Patientslag)
ggAcf(Patientslag, lag = 61)

patientsdiff <- diff(Patients)
autoplot(Patients)
ggAcf(patientsdiff)

#Applying the Mean, Naive, and Seasonal Naive method
#Set training date from 1992 to 2007
Azi2 <- window(Azithromycin, start=c(2016, 5), end=c(2020))
#Plot FC
autoplot(Azi2)

autoplot(Azithromycin) +
  autolayer(meanf(Azi2, h = 16),
            series = "Mean", PI = FALSE) +
  autolayer(naive(Azi2, h = 16),
            series = "Naive", PI = FALSE) +
  autolayer(snaive(Azi2, h = 16),
            series = "Seasonal Naive", PI = FALSE) +
  autolayer(snaive(Azi2)) +
  ggtitle("Forecasts for Azithromycin Prescription") +
  xlab("Year") + ylab("Azithromycin Prescription") +
  guides(color=guide_legend(title="Forecast"))

bothres <- residuals(naive(Azithromycin))
autoplot(bothres) + xlab(" ") + ylab(" ") +
  ggtitle("Residuals from naive method")
gghistogram(bothres) + ggtitle("Histogram of residuals")
ggAcf(bothres, lag = 61) + ggtitle("ACF of residuals")

Box.test(bothres, lag = 10, fitdf = 0)
Box.test(bothres, lag = 10, fitdf = 0, type = "Lj" )
checkresiduals(naive(Azithromycin))

snaive(Azi2)
autoplot(snaive(Azi2))

azibothdata %>% 
  as.data.frame() %>% 
  ggplot(aes(x = Azithromycin, y = Patientsmillion)) +
  ylab("Patients (millions)") +
  xlab("Azithromycin Prescription") +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE)
tslm(Patientsmillion ~ Azithromycin, data = precovid)

fc <- window(azibothdata, start = c(2016, 5))
view(fc)
fit.fc <- tslm(fc ~ Patients)
fcast <- forecast(fit.fc)
autoplot(fcast)
