# BSC Indonesia
# 28 Apr 2017
# Author: Ren Cabral
# 29 August 2017 - revised






BSCFUNCTION<-function(BerriedPolicyCompliance,PricePremium,TrawlBanCompliance, SizeLCompliance, Efforts,PLOTRESULTS,SIZELIMIT,SIZELIMITpolicy,TRAWLBANpolicy,OPENACCESSpolicy){
  
  Nmat<-180 #number of months in age-structured model
  AdjustTime<- 12*200 #months (12 months times x years)
  ProjectionTime<- 12*20 #months (12 months times x years)
  leadT<-0
  ModelTimeRange<- -AdjustTime:ProjectionTime #stabilizing time and projection time
  
  age<-c(1:Nmat)
  data <- data.frame(age=age)
  data$CW<-187*(1-exp(-1.13*((data$age/12)+0.0038)))
  maturity<- 1/(12*(1+exp(12.5*(1-(data$CW/118.98)))))
  data$m<-0.345
  data$m[6:length(age)]<-maturity[6:length(age)]
  natmort<-exp(-0.0934)
  eggs<-((19.128*data$CW)-1517.1)*1000
  data$eggs<-eggs*(eggs>=0)
  data$meggs<-data$m*data$eggs
  data$weight<-0.00007*(data$CW^2.997)
  K<-160960#1058543 #can be derived using catch-MSY method
  
  #Price<- 24967000/12500 #USD per ton, using data provided by the kkp and conversion of 12,500IDR/USD
  if (PricePremium==1){
    Price<- ((156.25*data$CW)+6875)*1000/12500
  }else{
    Price<- (50000000/12500)*(data$CW>0)
  }
  
  #source of Price: MMAF-National Value
  
  #Source of cost: SFP 2014
  #CostTripTrap<- 5000#165139
  #TrapCostPerF<- CostTripTrap*123910*15*0.7/12500 #497267*123910*15*0.7/12500 #trap per trip #revised 1/25/17 using KKP data (165,139 per trip)
  #NetCostPerF<- CostTripTrap*123910*15*0.2/12500 #497267*123910*15*0.2/12500
  #TrawlCostPerF<- CostTripTrap*123910*15*0.1/12500 #497267*123910*15*0.1/12500
  
#  CostMonthTrap<-(222000*15*123910)/(10*12500) # IDR per boat per day * no. of days per month * number of traps/ (30 traps * conversion to $)
#  CostMonthNet<-(146000*15*123910)/(10*12500) 
#  CostMonthTrawl<-(342000*15*123910)/(10*12500)   
  
#  TrapCostPerF<-  CostMonthTrap*0.7
#  NetCostPerF<- CostMonthNet*0.2
#  TrawlCostPerF<-  CostMonthTrawl*0.1
  
  TrapCostPerF<-(222000*52488000*.7)/(12*12500*10*Efforts[1]) #10 kg of harvest per day, 52,488MTx1000kg/MT, 12 is harvest per month
  NetCostPerF<- (146000*52488000*.2)/(12*12500*10*Efforts[2])
  TrawlCostPerF<- (343000*52488000*.1)/(12*12500*10*Efforts[3])  
  
  
  #POLICY PARAMETER HERE-------------------------
  BerriedPolicyCompliance<-BerriedPolicyCompliance #max is 1, 0 means no Berried Policy
  epsilon<-0.1#0.05 #for the open access scenario
  #epsilon<-0.5
  
  #SizeLimit<- 100 #unit is mm
  SizeLimit<-SIZELIMIT
  SizeLimitCompliance<-SizeLCompliance  #max is 1 (full compliance). 0 means no compliance
  
  PopINI<-matrix(1000000, ncol = 1, nrow = Nmat)
  
  Leslie <- matrix(0, ncol = Nmat, nrow = Nmat)
  x<-rep(natmort,Nmat)
  y <- diag(x)
  Leslie[c(2:Nmat),]<-y[c(1:Nmat-1),]
  Leslie[1,]<-data$meggs
  
  ##What is Rsurv?
  ###########----------------
  Biomass<-vector()
  SSB<-vector()
  Pop<-PopINI
  #GK<-7.85e8 ###ADJUST THIS TO GET BIOMASS == K


# ----GET GK
# for (GK in seq(1.197e8,1.2e8,by=0.0001e8)){
#   count<-0
#   RsurvVEC<-vector()
#   for (month in ModelTimeRange){
#     count<-count+1
#     Pop2<-floor(Leslie%*%Pop)
#     RsurvVEC[count]<-GK/Pop2[1,1] #output is 2.022027e-06
#     
#     #if (month%%6==0){
#     Pop2[1,1]<-GK
#     #} else
#     #{Pop2[1,1]<-0}
#     
#     Pop<-Pop2
#     Biomass[count]<-sum(Pop*data$weight)/(1000*1000) #from grams to mt
#     SSB[count]<-sum(floor(data$m*Pop)*data$weight)/(1000*1000)
#   }
#   indicator<-tail(Biomass,n=1)
#   if (indicator >K){break}
# }
# indicator
# GK

  GK<-1.1988e8
  count<-0
  RsurvVEC<-vector()
  for (month in ModelTimeRange){
    count<-count+1
    Pop2<-floor(Leslie%*%Pop)
    RsurvVEC[count]<-GK/Pop2[1,1] #output is 2.022027e-06
    
    #if (month%%6==0){
    Pop2[1,1]<-GK
    #} else
    #{Pop2[1,1]<-0}
    
    Pop<-Pop2
    Biomass[count]<-sum(Pop*data$weight)/(1000*1000) #from grams to mt
    SSB[count]<-sum(floor(data$m*Pop)*data$weight)/(1000*1000)
  }
  tail(Biomass,n=1)
  #plot(Biomass)
  SSBk<-tail(SSB,n=1)
  #########------------------
  Rsurv<-tail(RsurvVEC,n=1)   #Rsurv<-2.022027e-06
  
  ######--- MAIN MODEL HERE
  Strap<-1/(1+exp(25*(1-(data$CW/105))))
  Sgnet<-1/(1+exp(10*(1-(data$CW/90))))
  Strawl<-1/(1+exp(10*(1-(data$CW/80))))
  
  alpha<-1 ##ADJUSTABLE
  Biomass<-vector()
  SSB<-vector()
  TotHarvest<-vector()
  TrapHarvest<-vector()
  NetHarvest<-vector()
  TrawlHarvest<-vector()
  TotProfit<-vector()
  TrapProfit<-vector()
  NetProfit<-vector()
  TrawlProfit<-vector()
  
  Pop<-PopINI
  PopBerriedPolicy<-PopINI*0
  
  
  Etrap<-Efforts[1]#0.01#0.105
  Egnet<-Efforts[2]#0.002#0.025
  Etrawl<-Efforts[3]#0.00065#0.0114
  
  #Etrap0<-Etrap
  #Egnet0<-Egnet
  #Etrawl0<-Etrawl
  
  EtrapVec<-vector()
  EgnetVec<-vector()
  EtrawlVec<-vector()
  
  #catch limit
  data$CatchLimit<- (1-SizeLimitCompliance)
  data$CatchLimit[data$CW>=SizeLimit] = 1
  
  
  #data$CatchLimit[9:length(age)]<-1
  
  BerriedPolicyEfficiency=0
  
  count<-0
  for (month in ModelTimeRange){
    count<-count+1
    
    if (BerriedPolicyCompliance>0){
      if (month>0){
        BerriedPolicyEfficiency=BerriedPolicyCompliance
      }else {}
    }
    
    Ftrap=exp(-Etrap*Strap)
    Fgnet=exp(-Egnet*Sgnet)
    Ftrawl=exp(-Strawl*Etrawl)
    F=Ftrap*Fgnet*Ftrawl #F=exp(-Etrap*Strap-Egnet*Sgnet-Strawl*Etrawl)
    M=rep(natmort,Nmat)
    
    LeslieModel <- matrix(0, ncol = Nmat, nrow = Nmat)
    
    Mort<-M*F
    DiagMort <- diag(Mort)
    LeslieModel[c(2:Nmat),]<-DiagMort[c(1:Nmat-1),]
    LeslieModelBerried<-(LeslieModel>0)*1 
    LeslieModelBerried[1,]<-data$eggs
    LeslieModel[1,]<-data$meggs
    
    #harvest
    Harvest<-Pop*((1-Mort)*(1-F)/ ((1-M)+(1-F)))
    
    #BERRIED POLICY
    PopBerriedPolicy<-floor(Harvest*data$m*BerriedPolicyEfficiency)
    Harvest<-Harvest-PopBerriedPolicy
    
    TotHarvest[count]<-sum(data$weight*Harvest)/(1000*1000)
    
    #population transition
    Pop2<-floor((LeslieModel%*%Pop) + (LeslieModelBerried%*%PopBerriedPolicy))
    
    ##RECRUITMENT
    if ((Pop2[1,1]*Rsurv) < (GK)){
      Pop2[1,1]<-((GK^alpha+(GK-(Rsurv*Pop2[1,1]))^alpha)^(1/alpha))
    } else {
      Pop2[1,1]<-GK
    }
    
    #RECRUITMENT
    #if (month%%6==0){
    #Pop2[1,1]<-GK
    #} else
    #{Pop2[1,1]<-0}
    
    ##RECRUITMENT based on BH stock-recruitment relationship (assumption)
    #SSB<-sum(floor(data$m*Pop2)*data$weight)/(1000*1000)
    #varB<-0.2 #SSB level relative to SSBk at 1/2 Rmax 
    #Pop2[1,1]<-((GK/(varB*SSBk))*SSB)/(1+(SSB/(varB*SSBk)))
    ##Reference: https://notendur.hi.is/gunnar/kennsla/fii/lorna/fish480.pdf
    
    #PLOT the Stock-recruitment function
    if (FALSE){
      
      GK<-100 #max recruitment
      varB<-0.7 #fraction of SSBk where R is 1/2
      SSBk<-100 #SSB at K
      count<-0
      R<-vector()
      SSBLen<-1:150
      for (SSB in SSBLen){
        count<-count+1
        #R[count]<-((GK/(0.66*varB*SSBk))*SSB)/(1+(SSB/(varB*SSBk)))
        R[count]<-((GK/(varB*SSBk))*SSB)/(1+(SSB/(varB*SSBk)))
      }
      plot(SSBLen,R)
      
      
    }
    
    FracTrap<-(1-Ftrap)/((1-Ftrap)+(1-Fgnet)+(1-Ftrawl))
    FracNet<-(1-Fgnet)/((1-Ftrap)+(1-Fgnet)+(1-Ftrawl))
    FracTrawl<-(1-Ftrawl)/((1-Ftrap)+(1-Fgnet)+(1-Ftrawl))
    
    is.nan.data.frame <- function(x)
      do.call(cbind, lapply(x, is.nan))
    FracTrap[is.nan.data.frame(FracTrap)] <- 0
    FracNet[is.nan.data.frame(FracNet)] <- 0
    FracTrawl[is.nan.data.frame(FracTrawl)] <- 0
    
    TrapHarvest[count]<-sum(data$weight*Harvest*FracTrap)/(1000*1000)
    NetHarvest[count]<-sum(data$weight*Harvest*FracNet)/(1000*1000)
    TrawlHarvest[count]<-sum(data$weight*Harvest*FracTrawl)/(1000*1000)
    
    TrapProfit[count]<-(sum(data$weight*Harvest*FracTrap*Price)/(1000*1000)) - (TrapCostPerF*Etrap)
    NetProfit[count]<-(sum(data$weight*Harvest*FracNet*Price)/(1000*1000)) - (NetCostPerF*Egnet)
    TrawlProfit[count]<-(sum(data$weight*Harvest*FracTrawl*Price)/(1000*1000)) - (TrawlCostPerF*Etrawl)
    
    #TrapProfit[count]<-(TrapHarvest[count]*Price) - (TrapCostPerF*Etrap/Etrap0)
    #NetProfit[count]<-(NetHarvest[count]*Price) - (NetCostPerF*Egnet/Egnet0)
    #TrawlProfit[count]<-(TrawlHarvest[count]*Price) - (TrawlCostPerF*Etrawl/Etrawl0)
    
    Pop<-Pop2
    Biomass[count]<-sum(Pop*data$weight)/(1000*1000) #from grams to mt
    
    EtrapVec[count]<-Etrap
    EgnetVec[count]<-Egnet
    EtrawlVec[count]<-Etrawl
    
    #---POLICIES
    
    if (month>=0){
      #---OPEN ACCESS
      if (OPENACCESSpolicy==1){
        Etrap<-Etrap+(epsilon*(((sum(data$weight*Harvest*FracTrap*Price)/(1000*1000))/(TrapCostPerF*Etrap))-1)*Etrap)
        Egnet<-Egnet+(epsilon*(((sum(data$weight*Harvest*FracNet*Price)/(1000*1000))/(NetCostPerF*Egnet))-1)*Egnet)
        Etrawl<-Etrawl+(epsilon*(((sum(data$weight*Harvest*FracTrawl*Price)/(1000*1000))/(TrawlCostPerF*Etrawl))-1)*Etrawl)
      }
      
      #---TRAWL BAN
      if (TRAWLBANpolicy==1){
        
        Etrawl<-Efforts[3]*(1-TrawlBanCompliance)
        #print(Etrawl)
      } #no need to worry about OA. It should not be OA when there is a trawl ban
      
      #Etrawl<-0}
      
    }
    
    if (month == 0){
      #---SIZE LIMIT
      if (SIZELIMITpolicy==1){
        Strap<-Strap*data$CatchLimit
        Sgnet<-Sgnet*data$CatchLimit
        Strawl<-Strawl*data$CatchLimit}
    }
  }
  
  if (PLOTRESULTS==1){
    
    par(mfrow=c(3,4))
    
    
    
    
    plot(-leadT:ProjectionTime,Biomass[(length(Biomass)-leadT-ProjectionTime):length(Biomass)], xlab="Time (month)",ylab="Biomass (mt)",main="Total Biomass",type="o")
    plot(Pop,xlab="Age (month)",ylab="Individuals")
    plot(Pop*data$weight/(1000*1000), xlab="Age (month)",ylab="Biomass (mt)")
    plot(data$CW,Pop*data$weight/(1000*1000), xlab="CW (mm)",ylab="Biomass (mt)")
    
    plot(-leadT:ProjectionTime,TotHarvest[(length(Biomass)-leadT-ProjectionTime):length(Biomass)], xlab="Time (month)", ylab="Harvest (mt)",main="Total Harvest",type="o")
    plot(-leadT:ProjectionTime,TrapHarvest[(length(Biomass)-leadT-ProjectionTime):length(Biomass)], xlab="Time (month)", col="blue",ylab="Harvest (mt)",main="Trap Harvest",type="o")
    plot(-leadT:ProjectionTime,NetHarvest[(length(Biomass)-leadT-ProjectionTime):length(Biomass)], xlab="Time (month)", col="red",ylab="Harvest (mt)",main="Gillnet Harvest",type="o")
    plot(-leadT:ProjectionTime,TrawlHarvest[(length(Biomass)-leadT-ProjectionTime):length(Biomass)], xlab="Time (month)", col="darkgreen", ylab="Harvest (mt)",main="Trawl Harvest",type="o")
    
    plot(-leadT:ProjectionTime,TrapProfit[(length(Biomass)-leadT-ProjectionTime):length(Biomass)], xlab="Time (month)", ylab="Profit (USD)",col="blue",main="Trap Profit",type="o")
    plot(-leadT:ProjectionTime,NetProfit[(length(Biomass)-leadT-ProjectionTime):length(Biomass)], xlab="Time (month)", ylab="Profit (USD)",col="red",main="Gillnet Profit",type="o")
    plot(-leadT:ProjectionTime,TrawlProfit[(length(Biomass)-leadT-ProjectionTime):length(Biomass)], xlab="Time (month)", ylab="Profit (USD)",col="darkgreen",main="Trawl Profit",type="o")
    plot(-leadT:ProjectionTime,EtrapVec[(length(Biomass)-leadT-ProjectionTime):length(Biomass)], xlab="Time (month)", ylab="Effort",col="blue",main="Effort",ylim=c(0,max(EtrapVec,EgnetVec,EtrawlVec)),type="o")
    points(-leadT:ProjectionTime,EgnetVec[(length(Biomass)-leadT-ProjectionTime):length(Biomass)], xlab="Time (month)", col="red",type="o")
    points(-leadT:ProjectionTime,EtrawlVec[(length(Biomass)-leadT-ProjectionTime):length(Biomass)], xlab="Time (month)", col="darkgreen",type="o")
  }
  
  newList <- list("CW"=data$CW,"BIOM"=Pop*data$weight/(1000*1000),"Biomass"=Biomass,"TrapHarvestENDYEAR"=data$weight*Harvest*FracTrap/(1000*1000),"NetHarvestENDYEAR"=data$weight*Harvest*FracNet/(1000*1000),"TrawlHarvestENDYEAR"=data$weight*Harvest*FracTrawl/(1000*1000),"TrapProfit"=TrapProfit,"NetProfit"=NetProfit,"TrawlProfit"=TrawlProfit,"TrapHarvest"=TrapHarvest,"NetHarvest"=NetHarvest,"TrawlHarvest"=TrawlHarvest,"ModelTimeRange"=ModelTimeRange)
  return(newList)
  
}

#Efforts<-c(0.015,0.0041,0.002)
#Efforts<-c(0.004,0.001,0.0005)

#------------------------------------------------------------------------------------------------------------------------------------


#CHECKER
Efforts<-c(0.0301,0.00793,0.0038)
ProfitListCheck<-BSCFUNCTION(BerriedPolicyCompliance=0,PricePremium=0,TrawlBanCompliance=1, SizeLCompliance=1, Efforts,PLOTRESULTS=1,SIZELIMIT=100,SIZELIMITpolicy=1,TRAWLBANpolicy=1,OPENACCESSpolicy=0)
sum(ProfitListCheck$TrapHarvestENDYEAR)
sum(ProfitListCheck$NetHarvestENDYEAR)
sum(ProfitListCheck$TrawlHarvestENDYEAR)

#------------------plot for paper OPEN ACCESS
library(ggplot2)
library(gridExtra)

Efforts<-c(0.0301,0.00793,0.0038)
OASL<-BSCFUNCTION(BerriedPolicyCompliance=0,PricePremium=0,TrawlBanCompliance=1, SizeLCompliance=1, Efforts,PLOTRESULTS=0,SIZELIMIT=100,SIZELIMITpolicy=1,TRAWLBANpolicy=0,OPENACCESSpolicy=1)
OATB<-BSCFUNCTION(BerriedPolicyCompliance=0,PricePremium=0,TrawlBanCompliance=1, SizeLCompliance=1, Efforts,PLOTRESULTS=0,SIZELIMIT=100,SIZELIMITpolicy=0,TRAWLBANpolicy=1,OPENACCESSpolicy=1)
OABOTH<-BSCFUNCTION(BerriedPolicyCompliance=0,PricePremium=0,TrawlBanCompliance=1, SizeLCompliance=1, Efforts,PLOTRESULTS=0,SIZELIMIT=100,SIZELIMITpolicy=1,TRAWLBANpolicy=1,OPENACCESSpolicy=1)

time<-OASL$ModelTimeRange #we can use this in all inquiry below

#this is Open Access - profit
dPTrapOASL<-(OASL$TrapProfit[which(time==max(time),)] - OASL$TrapProfit[which(time==0,)])*100 / OASL$TrapProfit[which(time==0,)]
dPNetOASL<-(OASL$NetProfit[which(time==max(time),)] - OASL$NetProfit[which(time==0,)])*100 / OASL$NetProfit[which(time==0,)]
dPTrawlOASL<-(OASL$TrawlProfit[which(time==max(time),)] - OASL$TrawlProfit[which(time==0,)])*100 / OASL$TrawlProfit[which(time==0,)]

dPTrapOATB<-(OATB$TrapProfit[which(time==max(time),)] - OATB$TrapProfit[which(time==0,)])*100 / OATB$TrapProfit[which(time==0,)]
dPNetOATB<-(OATB$NetProfit[which(time==max(time),)] - OATB$NetProfit[which(time==0,)])*100 / OATB$NetProfit[which(time==0,)]
dPTrawlOATB<-(OATB$TrawlProfit[which(time==max(time),)] - OATB$TrawlProfit[which(time==0,)])*100 / OATB$TrawlProfit[which(time==0,)]

dPTrapOABOTH<-(OABOTH$TrapProfit[which(time==max(time),)] - OABOTH$TrapProfit[which(time==0,)])*100 / OABOTH$TrapProfit[which(time==0,)]
dPNetOABOTH<-(OABOTH$NetProfit[which(time==max(time),)] - OABOTH$NetProfit[which(time==0,)])*100 / OABOTH$NetProfit[which(time==0,)]
dPTrawlOABOTH<-(OABOTH$TrawlProfit[which(time==max(time),)] - OABOTH$TrawlProfit[which(time==0,)])*100 / OABOTH$TrawlProfit[which(time==0,)]

dPOA <- data.frame(Gear=c("Trap","Gillnet","Trawl","Trap","Gillnet","Trawl","Trap","Gillnet","Trawl"),
                 Policy=c("Size Limit","Size Limit","Size Limit","Trawl Ban","Trawl Ban","Trawl Ban","Size Limit and Trawl Ban","Size Limit and Trawl Ban","Size Limit and Trawl Ban"),
                 dP=c(dPTrapOASL,dPNetOASL,dPTrawlOASL,dPTrapOATB,dPNetOATB,dPTrawlOATB,dPTrapOABOTH,dPNetOABOTH,dPTrawlOABOTH))

#this is open access - catch
dCTrapOASL<-(OASL$TrapHarvest[which(time==max(time),)] - OASL$TrapHarvest[which(time==0,)])*100 / OASL$TrapHarvest[which(time==0,)]
dCNetOASL<-(OASL$NetHarvest[which(time==max(time),)] - OASL$NetHarvest[which(time==0,)])*100 / OASL$NetHarvest[which(time==0,)]
dCTrawlOASL<-(OASL$TrawlHarvest[which(time==max(time),)] - OASL$TrawlHarvest[which(time==0,)])*100 / OASL$TrawlHarvest[which(time==0,)]

dCTrapOATB<-(OATB$TrapHarvest[which(time==max(time),)] - OATB$TrapHarvest[which(time==0,)])*100 / OATB$TrapHarvest[which(time==0,)]
dCNetOATB<-(OATB$NetHarvest[which(time==max(time),)] - OATB$NetHarvest[which(time==0,)])*100 / OATB$NetHarvest[which(time==0,)]
dCTrawlOATB<-(OATB$TrawlHarvest[which(time==max(time),)] - OATB$TrawlHarvest[which(time==0,)])*100 / OATB$TrawlHarvest[which(time==0,)]

dCTrapOABOTH<-(OABOTH$TrapHarvest[which(time==max(time),)] - OABOTH$TrapHarvest[which(time==0,)])*100 / OABOTH$TrapHarvest[which(time==0,)]
dCNetOABOTH<-(OABOTH$NetHarvest[which(time==max(time),)] - OABOTH$NetHarvest[which(time==0,)])*100 / OABOTH$NetHarvest[which(time==0,)]
dCTrawlOABOTH<-(OABOTH$TrawlHarvest[which(time==max(time),)] - OABOTH$TrawlHarvest[which(time==0,)])*100 / OABOTH$TrawlHarvest[which(time==0,)]


dCOA <- data.frame(Gear=c("Trap","Gillnet","Trawl","Trap","Gillnet","Trawl","Trap","Gillnet","Trawl"),
                   Policy=c("Size Limit","Size Limit","Size Limit","Trawl Ban","Trawl Ban","Trawl Ban","Size Limit and Trawl Ban","Size Limit and Trawl Ban","Size Limit and Trawl Ban"),
                   dC=c(dCTrapOASL,dCNetOASL,dCTrawlOASL,dCTrapOATB,dCNetOATB,dCTrawlOATB,dCTrapOABOTH,dCNetOABOTH,dCTrawlOABOTH))


#this is open access - biomass
dBOASL<-(OASL$Biomass[which(time==max(time),)] - OASL$Biomass[which(time==0,)])*100 / OASL$Biomass[which(time==0,)]
dBOATB<-(OATB$Biomass[which(time==max(time),)] - OATB$Biomass[which(time==0,)])*100 / OATB$Biomass[which(time==0,)]
dBOABOTH<-(OABOTH$Biomass[which(time==max(time),)] - OABOTH$Biomass[which(time==0,)])*100 / OABOTH$Biomass[which(time==0,)]

dBOA<-data.frame(Policy=c("Size Limit","Trawl Ban","Size Limit and Trawl Ban"),
                 dB=c(dBOASL,dBOATB,dBOABOTH))



dPOA$Policy <- factor(dPOA$Policy, levels = dPOA$Policy)
dPOA$Gear <- factor(dPOA$Gear, levels = dPOA$Gear)
dCOA$Policy <- factor(dCOA$Policy, levels = dCOA$Policy)
dCOA$Gear <- factor(dCOA$Gear, levels = dCOA$Gear)
dBOA$Policy <- factor(dBOA$Policy, levels = dBOA$Policy)

P1_Legend<-ggplot(data=dPOA, aes(x=Policy, y=dP, fill=Gear)) +
  geom_bar(stat="identity", position=position_dodge())+
  scale_fill_brewer(palette="Paired")+
  theme_minimal()

library(gridExtra)
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
legend <- get_legend(P1_Legend)



P1<-ggplot(data=dPOA, aes(x=Policy, y=dP, fill=Gear)) +
  geom_bar(stat="identity", position=position_dodge())+
  scale_fill_brewer(palette="Paired")+
  xlab("") +
  ylab("% Change in Profit") +
  guides(fill=FALSE)+
  theme_minimal()+
  theme(axis.text.x = element_blank())

P2<-ggplot(data=dCOA, aes(x=Policy, y=dC, fill=Gear)) +
  geom_bar(stat="identity", position=position_dodge())+
  scale_fill_brewer(palette="Paired")+
  xlab("") +
  guides(fill=FALSE)+
  ylab("% Change in Catch") +
  theme_minimal()+
  theme(axis.text.x = element_blank())

P3<-ggplot(data=dBOA, aes(x=Policy, y=dB)) +
  geom_bar(stat="identity", position=position_dodge())+
  scale_fill_brewer(palette="Paired")+
  xlab("") +
  ylab("% Change in Biomass") +
  theme_minimal()

grid.arrange(arrangeGrob(arrangeGrob(P1,P2,P3)),legend,ncol=2,widths=c(5/6,1/6))#,ncol=2,widths=c(5/6,1/6)))
dPOA
dCOA
dBOA


#------------------plot for paper Constant effort
library(ggplot2)
library(gridExtra)

Efforts<-c(0.0301,0.00793,0.0038)
CESL<-BSCFUNCTION(BerriedPolicyCompliance=0,PricePremium=0,TrawlBanCompliance=1, SizeLCompliance=1, Efforts,PLOTRESULTS=0,SIZELIMIT=100,SIZELIMITpolicy=1,TRAWLBANpolicy=0,OPENACCESSpolicy=0)
CETB<-BSCFUNCTION(BerriedPolicyCompliance=0,PricePremium=0,TrawlBanCompliance=1, SizeLCompliance=1, Efforts,PLOTRESULTS=0,SIZELIMIT=100,SIZELIMITpolicy=0,TRAWLBANpolicy=1,OPENACCESSpolicy=0)
CEBOTH<-BSCFUNCTION(BerriedPolicyCompliance=0,PricePremium=0,TrawlBanCompliance=1, SizeLCompliance=1, Efforts,PLOTRESULTS=0,SIZELIMIT=100,SIZELIMITpolicy=1,TRAWLBANpolicy=1,OPENACCESSpolicy=0)

time<-CESL$ModelTimeRange #we can use this in all inquiry below

#this is CE - profit
dPTrapCESL<-(CESL$TrapProfit[which(time==max(time),)] - CESL$TrapProfit[which(time==0,)])*100 / CESL$TrapProfit[which(time==0,)]
dPNetCESL<-(CESL$NetProfit[which(time==max(time),)] - CESL$NetProfit[which(time==0,)])*100 / CESL$NetProfit[which(time==0,)]
dPTrawlCESL<-(CESL$TrawlProfit[which(time==max(time),)] - CESL$TrawlProfit[which(time==0,)])*100 / CESL$TrawlProfit[which(time==0,)]

dPTrapCETB<-(CETB$TrapProfit[which(time==max(time),)] - CETB$TrapProfit[which(time==0,)])*100 / CETB$TrapProfit[which(time==0,)]
dPNetCETB<-(CETB$NetProfit[which(time==max(time),)] - CETB$NetProfit[which(time==0,)])*100 / CETB$NetProfit[which(time==0,)]
dPTrawlCETB<-(CETB$TrawlProfit[which(time==max(time),)] - CETB$TrawlProfit[which(time==0,)])*100 / CETB$TrawlProfit[which(time==0,)]

dPTrapCEBOTH<-(CEBOTH$TrapProfit[which(time==max(time),)] - CEBOTH$TrapProfit[which(time==0,)])*100 / CEBOTH$TrapProfit[which(time==0,)]
dPNetCEBOTH<-(CEBOTH$NetProfit[which(time==max(time),)] - CEBOTH$NetProfit[which(time==0,)])*100 / CEBOTH$NetProfit[which(time==0,)]
dPTrawlCEBOTH<-(CEBOTH$TrawlProfit[which(time==max(time),)] - CEBOTH$TrawlProfit[which(time==0,)])*100 / CEBOTH$TrawlProfit[which(time==0,)]


dPCE <- data.frame(Gear=c("Trap","Gillnet","Trawl","Trap","Gillnet","Trawl","Trap","Gillnet","Trawl"),
                   Policy=c("Size Limit","Size Limit","Size Limit","Trawl Ban","Trawl Ban","Trawl Ban","Size Limit and Trawl Ban","Size Limit and Trawl Ban","Size Limit and Trawl Ban"),
                   dP=c(dPTrapCESL,dPNetCESL,dPTrawlCESL,dPTrapCETB,dPNetCETB,dPTrawlCETB,dPTrapCEBOTH,dPNetCEBOTH,dPTrawlCEBOTH))

#this is CE - catch
dCTrapCESL<-(CESL$TrapHarvest[which(time==max(time),)] - CESL$TrapHarvest[which(time==0,)])*100 / CESL$TrapHarvest[which(time==0,)]
dCNetCESL<-(CESL$NetHarvest[which(time==max(time),)] - CESL$NetHarvest[which(time==0,)])*100 / CESL$NetHarvest[which(time==0,)]
dCTrawlCESL<-(CESL$TrawlHarvest[which(time==max(time),)] - CESL$TrawlHarvest[which(time==0,)])*100 / CESL$TrawlHarvest[which(time==0,)]

dCTrapCETB<-(CETB$TrapHarvest[which(time==max(time),)] - CETB$TrapHarvest[which(time==0,)])*100 / CETB$TrapHarvest[which(time==0,)]
dCNetCETB<-(CETB$NetHarvest[which(time==max(time),)] - CETB$NetHarvest[which(time==0,)])*100 / CETB$NetHarvest[which(time==0,)]
dCTrawlCETB<-(CETB$TrawlHarvest[which(time==max(time),)] - CETB$TrawlHarvest[which(time==0,)])*100 / CETB$TrawlHarvest[which(time==0,)]

dCTrapCEBOTH<-(CEBOTH$TrapHarvest[which(time==max(time),)] - CEBOTH$TrapHarvest[which(time==0,)])*100 / CEBOTH$TrapHarvest[which(time==0,)]
dCNetCEBOTH<-(CEBOTH$NetHarvest[which(time==max(time),)] - CEBOTH$NetHarvest[which(time==0,)])*100 / CEBOTH$NetHarvest[which(time==0,)]
dCTrawlCEBOTH<-(CEBOTH$TrawlHarvest[which(time==max(time),)] - CEBOTH$TrawlHarvest[which(time==0,)])*100 / CEBOTH$TrawlHarvest[which(time==0,)]

dCCE <- data.frame(Gear=c("Trap","Gillnet","Trawl","Trap","Gillnet","Trawl","Trap","Gillnet","Trawl"),
                   Policy=c("Size Limit","Size Limit","Size Limit","Trawl Ban","Trawl Ban","Trawl Ban","Size Limit and Trawl Ban","Size Limit and Trawl Ban","Size Limit and Trawl Ban"),
                   dC=c(dCTrapCESL,dCNetCESL,dCTrawlCESL,dCTrapCETB,dCNetCETB,dCTrawlCETB,dCTrapCEBOTH,dCNetCEBOTH,dCTrawlCEBOTH))


#this is CE - biomass
dBCESL<-(CESL$Biomass[which(time==max(time),)] - CESL$Biomass[which(time==0,)])*100 / CESL$Biomass[which(time==0,)]
dBCETB<-(CETB$Biomass[which(time==max(time),)] - CETB$Biomass[which(time==0,)])*100 / CETB$Biomass[which(time==0,)]
dBCEBOTH<-(CEBOTH$Biomass[which(time==max(time),)] - CEBOTH$Biomass[which(time==0,)])*100 / CEBOTH$Biomass[which(time==0,)]

dBCE<-data.frame(Policy=c("Size Limit","Trawl Ban","Size Limit and Trawl Ban"),
                 dB=c(dBCESL,dBCETB,dBCEBOTH))



dPCE$Policy <- factor(dPCE$Policy, levels = dPCE$Policy)
dPCE$Gear <- factor(dPCE$Gear, levels = dPCE$Gear)
dCCE$Policy <- factor(dCCE$Policy, levels = dCCE$Policy)
dCCE$Gear <- factor(dCCE$Gear, levels = dCCE$Gear)
dBCE$Policy <- factor(dBCE$Policy, levels = dBCE$Policy)

P1_Legend<-ggplot(data=dPCE, aes(x=Policy, y=dP, fill=Gear)) +
  geom_bar(stat="identity", position=position_dodge())+
  scale_fill_brewer(palette="Paired")+
  theme_minimal()

library(gridExtra)
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
legend <- get_legend(P1_Legend)



P1<-ggplot(data=dPCE, aes(x=Policy, y=dP, fill=Gear)) +
  geom_bar(stat="identity", position=position_dodge())+
  scale_fill_brewer(palette="Paired")+
  xlab("") +
  ylab("% Change in Profit") +
  guides(fill=FALSE)+
  theme_minimal()+
  theme(axis.text.x = element_blank())

P2<-ggplot(data=dCCE, aes(x=Policy, y=dC, fill=Gear)) +
  geom_bar(stat="identity", position=position_dodge())+
  scale_fill_brewer(palette="Paired")+
  xlab("") +
  guides(fill=FALSE)+
  ylab("% Change in Catch") +
  theme_minimal()+
  theme(axis.text.x = element_blank())

P3<-ggplot(data=dBCE, aes(x=Policy, y=dB)) +
  geom_bar(stat="identity", position=position_dodge())+
  scale_fill_brewer(palette="Paired")+
  xlab("") +
  ylab("% Change in Biomass") +
  theme_minimal()

grid.arrange(arrangeGrob(arrangeGrob(P1,P2,P3)),legend,ncol=2,widths=c(5/6,1/6))#,ncol=2,widths=c(5/6,1/6)))
dPCE
dCCE
dBCE



###--------------The difference HERE

library(ggplot2)
library(gridExtra)

#this is diff - profit
dPTrapDIFFSL<-(CESL$TrapProfit[which(time==max(time),)] - OASL$TrapProfit[which(time==max(time),)]) / CESL$TrapProfit[which(time==0,)]
dPNetDIFFSL<-(CESL$NetProfit[which(time==max(time),)] - OASL$NetProfit[which(time==max(time),)]) / CESL$NetProfit[which(time==0,)]
dPTrawlDIFFSL<-(CESL$TrawlProfit[which(time==max(time),)] - OASL$TrawlProfit[which(time==max(time),)]) / CESL$TrawlProfit[which(time==0,)]

dPTrapDIFFTB<-(CETB$TrapProfit[which(time==max(time),)] - OATB$TrapProfit[which(time==max(time),)]) / CETB$TrapProfit[which(time==0,)]
dPNetDIFFTB<-(CETB$NetProfit[which(time==max(time),)] - OATB$NetProfit[which(time==max(time),)]) / CETB$NetProfit[which(time==0,)]
dPTrawlDIFFTB<-(CETB$TrawlProfit[which(time==max(time),)] - OATB$TrawlProfit[which(time==max(time),)]) / CETB$TrawlProfit[which(time==0,)]

dPTrapDIFFBOTH<-(CEBOTH$TrapProfit[which(time==max(time),)] - OABOTH$TrapProfit[which(time==max(time),)]) / CEBOTH$TrapProfit[which(time==0,)]
dPNetDIFFBOTH<-(CEBOTH$NetProfit[which(time==max(time),)] - OABOTH$NetProfit[which(time==max(time),)]) / CEBOTH$NetProfit[which(time==0,)]
dPTrawlDIFFBOTH<-(CEBOTH$TrawlProfit[which(time==max(time),)] - OABOTH$TrawlProfit[which(time==max(time),)]) / CEBOTH$TrawlProfit[which(time==0,)]


dPDIFF <- data.frame(Gear=c("Trap","Gillnet","Trawl","Trap","Gillnet","Trawl","Trap","Gillnet","Trawl"),
                   Policy=c("Size Limit","Size Limit","Size Limit","Trawl Ban","Trawl Ban","Trawl Ban","Size Limit and Trawl Ban","Size Limit and Trawl Ban","Size Limit and Trawl Ban"),
                   dP=c(dPTrapDIFFSL,dPNetDIFFSL,dPTrawlDIFFSL,dPTrapDIFFTB,dPNetDIFFTB,dPTrawlDIFFTB,dPTrapDIFFBOTH,dPNetDIFFBOTH,dPTrawlDIFFBOTH))

#this is CE - catch
dCTrapDIFFSL<-(CESL$TrapHarvest[which(time==max(time),)] - OASL$TrapHarvest[which(time==max(time),)]) / OASL$TrapHarvest[which(time==0,)]
dCNetDIFFSL<-(CESL$NetHarvest[which(time==max(time),)] - OASL$NetHarvest[which(time==max(time),)]) / OASL$NetHarvest[which(time==0,)]
dCTrawlDIFFSL<-(CESL$TrawlHarvest[which(time==max(time),)] - OASL$TrawlHarvest[which(time==max(time),)]) / OASL$TrawlHarvest[which(time==0,)]

dCTrapDIFFTB<-(CETB$TrapHarvest[which(time==max(time),)] - OATB$TrapHarvest[which(time==max(time),)]) / OATB$TrapHarvest[which(time==0,)]
dCNetDIFFTB<-(CETB$NetHarvest[which(time==max(time),)] - OATB$NetHarvest[which(time==max(time),)]) / OATB$NetHarvest[which(time==0,)]
dCTrawlDIFFTB<-(CETB$TrawlHarvest[which(time==max(time),)] - OATB$TrawlHarvest[which(time==max(time),)]) / OATB$TrawlHarvest[which(time==0,)]

dCTrapDIFFBOTH<-(CEBOTH$TrapHarvest[which(time==max(time),)] - OABOTH$TrapHarvest[which(time==max(time),)]) / OABOTH$TrapHarvest[which(time==0,)]
dCNetDIFFBOTH<-(CEBOTH$NetHarvest[which(time==max(time),)] - OABOTH$NetHarvest[which(time==max(time),)]) / OABOTH$NetHarvest[which(time==0,)]
dCTrawlDIFFBOTH<-(CEBOTH$TrawlHarvest[which(time==max(time),)] - OABOTH$TrawlHarvest[which(time==max(time),)]) / OABOTH$TrawlHarvest[which(time==0,)]

dCDIFF <- data.frame(Gear=c("Trap","Gillnet","Trawl","Trap","Gillnet","Trawl","Trap","Gillnet","Trawl"),
                   Policy=c("Size Limit","Size Limit","Size Limit","Trawl Ban","Trawl Ban","Trawl Ban","Size Limit and Trawl Ban","Size Limit and Trawl Ban","Size Limit and Trawl Ban"),
                   dC=c(dCTrapDIFFSL,dCNetDIFFSL,dCTrawlDIFFSL,dCTrapDIFFTB,dCNetDIFFTB,dCTrawlDIFFTB,dCTrapDIFFBOTH,dCNetDIFFBOTH,dCTrawlDIFFBOTH))

#this is diff - biomass
dBDIFFSL<-(CESL$Biomass[which(time==max(time),)] - OASL$Biomass[which(time==max(time),)]) / CESL$Biomass[which(time==0,)]
dBDIFFTB<-(CETB$Biomass[which(time==max(time),)] - OATB$Biomass[which(time==max(time),)]) / CETB$Biomass[which(time==0,)]
dBDIFFBOTH<-(CEBOTH$Biomass[which(time==max(time),)] - OABOTH$Biomass[which(time==max(time),)]) / CEBOTH$Biomass[which(time==0,)]

dBDIFF<-data.frame(Policy=c("Size Limit","Trawl Ban","Size Limit and Trawl Ban"),
                 dB=c(dBDIFFSL,dBDIFFTB,dBDIFFBOTH))

dPDIFF$Policy <- factor(dPDIFF$Policy, levels = dPDIFF$Policy)
dPDIFF$Gear <- factor(dPDIFF$Gear, levels = dPDIFF$Gear)
dCDIFF$Policy <- factor(dCDIFF$Policy, levels = dCDIFF$Policy)
dCDIFF$Gear <- factor(dCDIFF$Gear, levels = dCDIFF$Gear)
dBDIFF$Policy <- factor(dBDIFF$Policy, levels = dBDIFF$Policy)

P1_Legend<-ggplot(data=dPDIFF, aes(x=Policy, y=dP, fill=Gear)) +
  geom_bar(stat="identity", position=position_dodge())+
  scale_fill_brewer(palette="Paired")+
  theme_minimal()

library(gridExtra)
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
legend <- get_legend(P1_Legend)



P1<-ggplot(data=dPDIFF, aes(x=Policy, y=dP, fill=Gear)) +
  geom_bar(stat="identity", position=position_dodge())+
  scale_fill_brewer(palette="Paired")+
  xlab("") +
  ylab(expression(paste(Delta, "Profit relative to present")))+
#    "Profit difference relative to present") +
  guides(fill=FALSE)+
  theme_minimal()+
  theme(axis.text.x = element_blank())

P2<-ggplot(data=dCDIFF, aes(x=Policy, y=dC, fill=Gear)) +
  geom_bar(stat="identity", position=position_dodge())+
  scale_fill_brewer(palette="Paired")+
  xlab("") +
  guides(fill=FALSE)+
  ylab(expression(paste(Delta, "Catch relative to present")))+
  #  ylab("Catch difference relative to present") +
  theme_minimal()+
  theme(axis.text.x = element_blank())

P3<-ggplot(data=dBDIFF, aes(x=Policy, y=dB)) +
  geom_bar(stat="identity", position=position_dodge())+
  scale_fill_brewer(palette="Paired")+
  xlab("") +
#  ylab("Biomass difference relative to present") +
  ylab(expression(paste(Delta, "Biomass relative to present")))+
  theme_minimal()

grid.arrange(arrangeGrob(arrangeGrob(P1,P2,P3)),legend,ncol=2,widths=c(5/6,1/6))#,ncol=2,widths=c(5/6,1/6)))
dPDIFF
dCDIFF
dBDIFF

######----------------------------COMPLIANCE

#Profit - varying compliance - 
TrapProfit<-vector()
NetProfit<-vector()
TrawlProfit<-vector()
#Efforts<-c(0.102,0.024,0.011)
count<-0
for (i in seq(0,1,0.1)){
  count<-count+1
  ProfitListCheck<-BSCFUNCTION(BerriedPolicyCompliance=0,PricePremium=0,TrawlBanCompliance=i,SizeLCompliance=i, Efforts,PLOTRESULTS=0,SIZELIMIT=100,SIZELIMITpolicy=1,TRAWLBANpolicy=1,OPENACCESSpolicy=0)
  TrapProfit[count]<-tail(ProfitListCheck$TrapProfit,1)
  NetProfit[count]<-tail(ProfitListCheck$NetProfit,1)
  TrawlProfit[count]<-tail(ProfitListCheck$TrawlProfit,1)
}

TotalCatch<-(TrapProfit+NetProfit+TrawlProfit)-(TrapProfit[1]+NetProfit[1]+TrawlProfit[1])
plot(seq(0,1,0.1),TotalCatch, xlim=c(0,1), ylim=c(min(TotalCatch,na.rm=TRUE),max(TotalCatch,na.rm=TRUE)),xlab="Compliance",ylab="Additional Profit (US$)",pch=1,col="black")

TrapOnly<-(TrapProfit)-(TrapProfit[1])
plot(seq(0,1,0.1),TrapOnly, xlim=c(0,1), ylim=c(min(TrapOnly,na.rm=TRUE),max(TrapOnly,na.rm=TRUE)),xlab="Compliance",ylab="Additional Profit for trap fishers (US$)",pch=1,col="black")



plot(seq(0,1,0.1),TrapProfit, xlim=c(0,1), ylim=c(min(TrapProfit,NetProfit,TrawlProfit,na.rm=TRUE),max(TrapProfit,NetProfit,TrawlProfit,na.rm=TRUE)),xlab="Compliance",ylab="Profit (US$)",pch=1,col="blue")
points(seq(0,1,0.1),NetProfit,pch=0,col="red")
points(seq(0,1,0.1),TrawlProfit,pch=2,col="darkgreen")
legend("topright", c("Trap", "Gillnet","Trawl"),pch=c(1,0,2),col = c("blue","red","darkgreen"), bty = "n")
#ProfitListCheck<-BSCFUNCTION(SizeLCompliance=1, Efforts,PLOTRESULTS=1,SIZELIMIT=100,SIZELIMITpolicy=1,TRAWLBANpolicy=0,OPENACCESSpolicy=1)


##-----------------------Plotting the growth curves

a<-seq(1,180)

CWfemale<-187*(1-exp(-1.13*((a/12)+0.0038)))
CWmale<-185*(1-exp(-1.26*((a/12)+0.0034)))

plot(a,CWfemale,pch=16, col="blue",xlab="month",ylab="Carapace width (mm)")
points(a,CWmale,pch=15, col="red")



#tradeoff
TRAPHARV<-vector()
NETHARV<-vector()
count<-0
for (a in seq(0,0.1,0.001)){
  for (b in seq(0,0.1,0.001)){
    count<-count+1
    Efforts<-c(a,b,0)
    ProfitListCheck<-BSCFUNCTION(BerriedPolicyCompliance=0,PricePremium=0,TrawlBanCompliance=1, SizeLCompliance=1, Efforts,PLOTRESULTS=0,SIZELIMIT=100,SIZELIMITpolicy=1,TRAWLBANpolicy=1,OPENACCESSpolicy=0)
    TRAPHARV[count]<-tail(ProfitListCheck$TrapHarvest, n=1)
    NETHARV[count]<-tail(ProfitListCheck$NetHarvest, n=1)
  }
}

write.csv(TRAPHARV, file = "C:/Users/Ren/Desktop/BSC/TRAPHARV.csv")
write.csv(NETHARV, file = "C:/Users/Ren/Desktop/BSC/NETHARV.csv")

TRAPHARV <- read.csv(file = "C:/Users/Ren/Desktop/BSC/TRAPHARV.csv", header=TRUE)
NETHARV <- read.csv(file = "C:/Users/Ren/Desktop/BSC/NETHARV.csv", header=TRUE)


tiff('C:/Users/Ren/Desktop/BSC/Fig7.tiff')
plot(TRAPHARV$x*12,NETHARV$x*12,xlab="Trap harvest (MT)",ylab="Gillnet harvest (MT)")
#lines(c(36741.6,100000),c(10497.6,10497.6),col="red",lwd=2)
#lines(c(36741.6,36741.6),c(10497.6,100000),col="red",lwd=2)

lines(c(0,100000),c(10497.6,10497.6),col="red",lwd=2)
lines(c(36741.6,36741.6),c(0,100000),col="red",lwd=2)

dev.off()




#------------------------

#ProfitListCheck<-BSCFUNCTION(Efforts<-c(0.012,0.0032,0.0016),PLOTRESULTS=1,SIZELIMIT=100,SIZELIMITpolicy=0,TRAWLBANpolicy=0,OPENACCESSpolicy=0)


#PROFIT - Varying SIZE LIMIT, no trawl ban
TrapProfit<-vector()
NetProfit<-vector()
TrawlProfit<-vector()
Efforts<-c(0.102,0.024,0.011)
for (i in seq(50,150,10)){
  ProfitListCheck<-BSCFUNCTION(BerriedPolicyCompliance=0,PricePremium=0,SizeLCompliance=1, Efforts,PLOTRESULTS=0,SIZELIMIT=i,SIZELIMITpolicy=1,TRAWLBANpolicy=0,OPENACCESSpolicy=0)
  TrapProfit[i]<-tail(ProfitListCheck$TrapProfit,1)
  NetProfit[i]<-tail(ProfitListCheck$NetProfit,1)
  TrawlProfit[i]<-tail(ProfitListCheck$TrawlProfit,1)
}
plot(TrapProfit, xlim=c(50,150), ylim=c(0,max(TrapProfit,NetProfit,TrawlProfit,na.rm=TRUE)),xlab="Size Limit (mm)",ylab="Profit (US$)",pch=1,col="blue")
points(NetProfit,pch=0,col="red")
points(TrawlProfit,pch=2,col="darkgreen")
legend("topright", c("Trap", "Gillnet","Trawl"),pch=c(1,0,2),col = c("blue","red","darkgreen"), bty = "n")


#HARVEST - varying SIZE LIMIT, open access 
TrapHarvest<-vector()
NetHarvest<-vector()
TrawlHarvest<-vector()
Efforts<-c(0.102,0.024,0.011)
for (i in seq(50,150,10)){
  ProfitListCheck<-BSCFUNCTION(BerriedPolicyCompliance=0,PricePremium=0,SizeLCompliance=1, Efforts,PLOTRESULTS=0,SIZELIMIT=i,SIZELIMITpolicy=1,TRAWLBANpolicy=0,OPENACCESSpolicy=1)
  TrapHarvest[i]<-sum(ProfitListCheck$TrapHarvest)
  NetHarvest[i]<-sum(ProfitListCheck$NetHarvest)
  TrawlHarvest[i]<-sum(ProfitListCheck$TrawlHarvest)
}
plot(TrapHarvest, xlim=c(50,150), ylim=c(0,max(NetHarvest,TrapHarvest,na.rm=TRUE)),xlab="Size Limit (mm)",ylab="Harvest (MT)",pch=1,col="blue")
points(NetHarvest,pch=0,col="red")
points(TrawlHarvest,pch=2,col="darkgreen")
legend("topright", c("Trap", "Gillnet","Trawl"),pch=c(1,0,2),col = c("blue","red","darkgreen"), bty = "n")




#varying Compliance SIZE LIMIT Harvest
TrapHarvest<-vector()
NetHarvest<-vector()
TrawlHarvest<-vector()
Efforts<-c(0.102,0.024,0.011)
count<-0
for (i in seq(0,1,0.1)){
  count<-count+1
  ProfitListCheck<-BSCFUNCTION(BerriedPolicyCompliance=1,PricePremium=0, SizeLCompliance=i, Efforts,PLOTRESULTS=0,SIZELIMIT=120,SIZELIMITpolicy=1,TRAWLBANpolicy=0,OPENACCESSpolicy=1)
  TrapHarvest[count]<-sum(ProfitListCheck$TrapHarvest)
  NetHarvest[count]<-sum(ProfitListCheck$NetHarvest)
  TrawlHarvest[count]<-sum(ProfitListCheck$TrawlHarvest)
}
plot(seq(0,1,0.1),TrapHarvest, xlim=c(0,1), ylim=c(0,max(TrawlHarvest,NetHarvest,TrapHarvest,na.rm=TRUE)),xlab="Size Limit Compliance",ylab="Harvest (MT)",pch=1,col="blue")
points(seq(0,1,0.1),NetHarvest,pch=0,col="red")
points(seq(0,1,0.1),TrawlHarvest,pch=2,col="darkgreen")
legend("topright", c("Trap", "Gillnet","Trawl"),pch=c(1,0,2),col = c("blue","red","darkgreen"), bty = "n")






if (FALSE){
  #---- open access
  SIZELIMIT<-100 #in mm
  OPENACCESSpolicy<- 1	#1 for yes, 0 for no
  SIZELIMITpolicy<- 0	#1 for yes, 0 for no
  TRAWLBANpolicy<- 0    	#1 for yes, 0 for no
  PLOTRESULTS<- 1		#0 for yes, 0 for no
  
  ProfitList<-BSCFUNCTION(Efforts,PLOTRESULTS,SIZELIMIT,SIZELIMITpolicy,TRAWLBANpolicy,OPENACCESSpolicy)
  
  SUMMARY<-data.frame(ProfitList$ModelTimeRange)
  colnames(SUMMARY)<-"Month"
  SUMMARY$OATRAPprofit<-ProfitList$TrapProfit
  SUMMARY$OANetProfit<-ProfitList$NetProfit
  SUMMARY$OATrawlProfit<-ProfitList$TrawlProfit
  
  plot(ProfitList$CW,ProfitList$BIOM)
  
  BIOMSUMMARY<-data.frame(ProfitList$CW)
  CATCHSUMMARY<-data.frame(ProfitList$CW)
  
  BIOMSUMMARY$OAbiom<-ProfitList$BIOM
  CATCHSUMMARY$OATrapHarvest<-ProfitList$TrapHarvest
  CATCHSUMMARY$OANetHarvest<-ProfitList$NetHarvest
  CATCHSUMMARY$OATrawlHarvest<-ProfitList$TrawlHarvest
  
  #---size limit
  ProfitList_SizeL<-BSCFUNCTION(Efforts,PLOTRESULTS=1,SIZELIMIT=100,SIZELIMITpolicy=1,TRAWLBANpolicy=0,OPENACCESSpolicy=0)
  BIOMSUMMARY$SizeLbiom<-ProfitList_SizeL$BIOM - BIOMSUMMARY$OAbiom
  CATCHSUMMARY$SLTrapHarvest<-ProfitList_SizeL$TrapHarvest - CATCHSUMMARY$OATrapHarvest
  CATCHSUMMARY$SLNetHarvest<-ProfitList_SizeL$NetHarvest - CATCHSUMMARY$OANetHarvest
  CATCHSUMMARY$SLTrawlHarvest<-ProfitList_SizeL$TrawlHarvest - CATCHSUMMARY$OATrawlHarvest
  
  SUMMARY$SLTRAPprofit<-ProfitList_SizeL$TrapProfit
  SUMMARY$SLNetProfit<-ProfitList_SizeL$NetProfit
  SUMMARY$SLTrawlProfit<-ProfitList_SizeL$TrawlProfit
  
  #---trawl ban
  ProfitList_TrawlB<-BSCFUNCTION(Efforts,PLOTRESULTS=1,SIZELIMIT=100,SIZELIMITpolicy=0,TRAWLBANpolicy=1,OPENACCESSpolicy=0)
  BIOMSUMMARY$TrawlBbiom<-ProfitList_TrawlB$BIOM - BIOMSUMMARY$OAbiom
  CATCHSUMMARY$TBTrapHarvest<-ProfitList_TrawlB$TrapHarvest - CATCHSUMMARY$OATrapHarvest
  CATCHSUMMARY$TBNetHarvest<-ProfitList_TrawlB$NetHarvest - CATCHSUMMARY$OANetHarvest
  CATCHSUMMARY$TBTrawlHarvest<-ProfitList_TrawlB$TrawlHarvest - CATCHSUMMARY$OATrawlHarvest
  
  SUMMARY$TBTRAPprofit<-ProfitList_TrawlB$TrapProfit
  SUMMARY$TBNetProfit<-ProfitList_TrawlB$NetProfit
  SUMMARY$TBTrawlProfit<-ProfitList_TrawlB$TrawlProfit
  
  #--- trawl ban and size limit
  ProfitList_SizeLTrawlB<-BSCFUNCTION(Efforts,PLOTRESULTS=1,SIZELIMIT=100,SIZELIMITpolicy=1,TRAWLBANpolicy=1,OPENACCESSpolicy=0)
  BIOMSUMMARY$SizeLTrawlBbiom<-ProfitList_SizeLTrawlB$BIOM - BIOMSUMMARY$OAbiom
  CATCHSUMMARY$SLTBTrapHarvest<-ProfitList_SizeLTrawlB$TrapHarvest - CATCHSUMMARY$OATrapHarvest
  CATCHSUMMARY$SLTBNetHarvest<-ProfitList_SizeLTrawlB$NetHarvest - CATCHSUMMARY$OANetHarvest
  CATCHSUMMARY$SLTBTrawlHarvest<-ProfitList_SizeLTrawlB$TrawlHarvest - CATCHSUMMARY$OATrawlHarvest
  
  SUMMARY$SLTBTRAPprofit<-ProfitList_SizeLTrawlB$TrapProfit
  SUMMARY$SLTBNetProfit<-ProfitList_SizeLTrawlB$NetProfit
  SUMMARY$SLTBTrawlProfit<-ProfitList_SizeLTrawlB$TrawlProfit
  
  #write.csv(BIOMSUMMARY, file = "C:/Users/Ren/Desktop/BSC/BSCBiomSummary23Jan17.csv",row.names=FALSE)
  #write.csv(CATCHSUMMARY, file = "C:/Users/Ren/Desktop/BSC/BSCCatchSummary24Jan17.csv",row.names=FALSE)
  
  SUMMARYsub<-SUMMARY[SUMMARY$Month==0 | SUMMARY$Month==max(SUMMARY$Month),]*12
  write.csv(SUMMARYsub, file = "C:/Users/Ren/Desktop/BSC/BSCProfitSummary24Jan17.csv",row.names=FALSE)
  
  SUMMARYsub<-SUMMARYsub[,-1]
  #barplot(as.matrix(SUMMARYsub), beside=T,legend=c("Current", "Future (in year 20)"), names.arg=c("Trap","Gillnet","Trawl","Trap","Gillnet","Trawl","Trap","Gillnet","Trawl","Trap","Gillnet","Trawl"))
  barplot(as.matrix(SUMMARYsub), beside=T, names.arg=c("Trap","Gillnet","Trawl","Trap","Gillnet","Trawl","Trap","Gillnet","Trawl","Trap","Gillnet","Trawl"),col=c("aquamarine3","brown1"))
  title(ylab="Profit (US$)")
  text(x=5,y=6e7,"Open Access")
  text(x=15,y=6e7,"Size Limit (100mm CW)")
  text(x=23,y=6e7,"Trawl Ban")
  text(x=33,y=6e7,"Size Limit & Trawl Ban")
  legend("topright", c("Current", "Future (in year 20)"),fill = c("aquamarine3","brown1"), bty = "n")
  
}