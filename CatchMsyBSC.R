setwd("C:/Users/Ren/Desktop/BSC")
library(fishmethods)
mydata<-read.csv("BSCdata.csv")

plot(mydata$Year,mydata$Catch,type="o",pch=16,xlab="Year",ylab="Harvest (MT)")


catchmsy <- dget("catchmsyFUNCTION.R")
outpt<-catchmsy(year=mydata$Year,
catch=mydata$Catch,catchCV=NULL,
catargs=list(dist="none",low=0,up=Inf,unit="MT"),
l0=list(low=0.75,up=0.99,step=0.01),
lt=list(low=0.25,up=0.75,refyr=2015),sigv=0,
k=list(dist="unif",low=52488,up=200000,mean=0,sd=0),
#r=list(dist="unif",low=0.01,up=0.3,mean=0,sd=0),
r=list(dist="unif",low=0.8,up=1.8,mean=0,sd=0),
M=list(dist="unif",low=1.12,up=1.12,mean=1.12,sd=0.00),
nsims=3000)

outpt$Initial
outpt$Estimates
outpt$Parameters

values<-outpt$Values
valuesFin<-values[values$likeli==1,]
head(valuesFin)
median(valuesFin$l0)
mean(valuesFin$l0)

##########################################################
setwd("E:/Research UCSB/Indonesia Project/CatchMSY")
library(fishmethods)
library(xlsx)
fishdata <- read.xlsx("SampleDataFin.xlsx", 1)

#catchmsy <- dget("catchmsy.R")

outpt<-catchmsy(year=fishdata$year,
catch=fishdata$catch,catchCV=NULL,
catargs=list(dist="none",low=0,up=Inf,unit="MT"),
l0=list(low=0.8,up=0.8,step=0),
lt=list(low=0.1,up=0.6,refyr=2002),sigv=0,
k=list(dist="unif",low=4333,up=433300,mean=0,sd=0),
r=list(dist="unif",low=0.015,up=0.1,mean=0,sd=0),
M=list(dist="unif",low=0.18,up=0.18,mean=0.00,sd=0.00),nsims=30000)

outpt$Estimates


###alternative method when everything is working fine
setwd("E:/Research UCSB/Indonesia Project/CatchMSY")
library(fishmethods)
fishdata <- read.csv("skipjackdata.csv")
#plot(fishdata)
catchmsy <- dget("catchmsy.R")
head(fishdata)
#data(lingcod)
outpt<-catchmsy(year=fishdata$year,
catch=fishdata$catch,catchCV=NULL,
catargs=list(dist="none",low=0,up=Inf,unit="MT"),
l0=list(low=0.75,up=0.99,step=0),
lt=list(low=0.3,up=0.75,refyr=2014),sigv=0,
k=list(dist="unif",low=496682*2,up=496682*20,mean=0,sd=0),
r=list(dist="unif",low=0.2,up=1,mean=0,sd=0),
M=list(dist="unif",low=0,up=0.3,mean=0.00,sd=0.00),
nsims=10000)#30000)

#medium resilience based on fishbase

outpt$Initial
outpt$Estimate
outpt$Parameters

values<-outpt$Values
valuesFin<-values[values$likeli==1,]
median(valuesFin$l0)
mean(valuesFin$l0)



###applied to real data [BSC data]
library(fishmethods)
setwd("E:/Research UCSB/Indonesia Project/CatchMSY")
library(xlsx)
BSCdata <- read.xlsx("BSCdata.xlsx", 1)
catchmsyPT <- dget("catchmsyPT.R")
outpt<-catchmsyPT(year=BSCdata$Year,
catch=BSCdata$Catch,catchCV=NULL,
catargs=list(dist="none",low=0,up=Inf,unit="MT"),
l0=list(low=0.50,up=0.99,step=0.01),
lt=list(low=0.01,up=0.5,refyr=2015),sigv=0,
k=list(dist="unif",low=52488,up=5248800,mean=0,sd=0),
r=list(dist="unif",low=0.01,up=0.5,mean=0,sd=0),
M=list(dist="unif",low=0,up=0.3,mean=0.00,sd=0.00),
nsims=30000)

outpt$Initial
outpt$Estimates
outpt$Parameters

values<-outpt$Values
valuesFin<-values[values$likeli==1,]
head(valuesFin)
median(valuesFin$l0)
mean(valuesFin$l0)