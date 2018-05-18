rm(list=ls())
require(ade4)
require(maps)
require(mapdata)
require(dplyr)
setwd("~/Jerico next 2017 cruise/data")
#Data
t=read.csv("Cruise_data_2017_V9_complete_Flowcam_NGS.csv",h=T,row.names=1)

##chek correlation with flowcam noctiluca and zooscan noctiluca
colnames(t)
#both sensors give a concnetration, therefore pearson correlation is used
#first deleta all NA rows
test<-t[,c(147,61)]
test<-na.omit(test)
cor(test$Noctiluca,test$Zoopl_Noctiluca)

