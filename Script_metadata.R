#### Script of the complete datatable ####
#script is based on work of Olivier Beauchard and Elisabeth Debusschere
#Jerico-Next LifeWatch cruise by RV Simon Stevin, led by VLIZ (Jonas Mortelmans)
#May 8-11 2017
t=read.csv("data/Cruise_data_2017_V11_complete_Flowcam_NGS.csv",h=T,row.names=1)

#### 1.  Complete df with additional parameters ####

#make DayorNight numeric
t$DayorNight<-ifelse(t$DayorNight=="Night","0","1")

#add extra ratios of nutrients to the dataframe
t$N.Si<-(t$Nutr_NO2+t$Nutr_NO3+t$Nutr_NOX)/t$Nutr_SiO2
t$N.P<-(t$Nutr_NO2+t$Nutr_NO3+t$Nutr_NOX)/t$Nutr_PO4
t$N.P.Si<-(t$Nutr_NO2+t$Nutr_NO3+t$Nutr_NOX)/t$Nutr_PO4/t$Nutr_SiO2
# Sensors present during the trip

#### 2. add units to each parameter ####
#1.Flowcam: phytoplankton groups (100 - 300µm), processed by Luz Amadei Martinéz
#2.Zooscan: zooplankton groups (300 - 2000µm), processed by Jonas Mortelmans
#3.FCM (vliz): size groups of Phytoplankton (<800µm), processed by Machteld Rijkeboer
#4.Pigment data transformed via chemtax to phytoplankton groups, processed by Reinhoud de Blok
#5.Abiotic parameters derived from the underway system of the RV Simon Stevin, processed by...

#dataframe:
##rows are stations
##columns are parameters measured at that station

#units of all the parameters:
names(t)
#StartDate and EndDate in UTC
#latitude and longitude in decimal degrees
#nutriens in µmol per liter
#pimgents in µg per liter
#zooplankton: ind per m³
#UW_navDepth50kHz in meters
#UW_Salinity in PSU
#UW_Temperatuur in °C
#BO2_curvelmean_ss: mean surface current velocity (m/s)
#flowcam sammples: cells per liter
#Chemtax groups: µg chlorophyll a per liter
#Number of day: number of day within the sampling cruise
#Dayor NIght: time of sampling at that station
#fcm data of rws, CNRS-LOG and VLIZ: cells per ml and total fluorescence red per ml
#extra ratios of nutrients were calculated:
##N.Si; N.P; N.P.Si no units here


#ctd data was skipped, due to the low number of sampled station
#pigment data was used in the format of chemtax
#Turb_Secchi; UW_id_nav"; "UW_OctansHeading";"UW_OdomDepth200khz";"UW_CourseOverGround" ; "UW_SpeedOverGround" ; "UW_OdomDepth33khz";"UW_Speedlog"; "UW_FLRTchla"  ;" "UW_RelativeHumidity;"UW_WindDirection"  ;"UW_WindSpeed" ; "UW_AirPressure";"UW_SoundVelocity"; "UW_WaterFlow"  are not relevant for our research
#only fcm of vliz was used since this had the highest number of station sampled.
#wisp, frrf data was skipped

#### 3. rename parameters ####
#create final dataframe with all relevant parameters per station
df=t[,c(3:4, 5:10,168:170,72,75:76, 127,166:167,43:64,107:112, 130:157)]
names(df)
#rename all parameters
names(df)[3] <- "NH4"
names(df)[4] <- "NO2"
names(df)[5] <- "NO3"
names(df)[6] <- "NOX"
names(df)[7] <- "PO4"
names(df)[8] <- "SiO2"
names(df)[12] <- "Depth"
names(df)[13] <- "Salinity"
names(df)[14] <- "Temperature"
names(df)[15] <- "Current Velocity"

names(df)[18] <- "Amphipoda"
names(df)[19] <- "Annelida"
names(df)[20] <- "Anomura"
names(df)[21] <- "Appendicularia"
names(df)[22] <- "Brachyura.zoe"
names(df)[23] <- "Branchiopoda"
names(df)[24] <- "Calanoidea"
names(df)[25] <- "Caridae.zoe"
names(df)[26] <- "Chaetognatha"
names(df)[27] <- "Cirripeda.cypris"
names(df)[28] <- "Cirripeda.nauplius"
names(df)[29] <- "Cnidaria"
names(df)[30] <- "Ctenophora"
names(df)[31] <- "Cumacea"
names(df)[32] <- "Echinodermata"
names(df)[33] <- "Harpacticoida"
names(df)[34] <- "Mollusca"
names(df)[35] <- "Mysidae"
names(df)[36] <- "Noctiluca"
names(df)[37] <- "Pisces.egg"
names(df)[38] <- "Pisces.larvae"
names(df)[39] <- "Porcellanidae.zoe"
#zooplankton unit needs to be rescaled to cells per L (now it's cells per 1000L)
df[,18:39]<-df[,18:39]/1000
#overview of the newly formed dataframe
names(df)
#1:2 lat & lon (decimal degrees)
#3:11 nutrients (µmol per liter)
#12:15 environmental parameters 
#16:17 cruise info
#18:39 Zooplankton (cells per liter)
#40:45 chemtax based phytoplankton groups (µg chlorophyll a per liter)
#46:73 flowcam based phytoplankton groups (cells per liter)

####4. clean dataframe ####
#remove NA rows
df<-na.omit(df) #from 44 stations to 38 stations

write.csv(df,"stationdf.csv")
#first column are the station numbers, I added manually the column name "Station" after export
