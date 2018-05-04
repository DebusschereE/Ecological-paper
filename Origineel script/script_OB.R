require(ade4)
require(maps)
require(mapdata)
require(packfor)
install.packages("packfor", repos="http://R-Forge.R-project.org") 
#Data
t=read.csv("t.csv",h=T,row.names=1)
colnames(t)

#Geographic coordinates
xy=t[,2:3]

#Estuaries
est=t$est

#Distances
dis=t[,4:7]

#Abiotic
abio=t[,8:14]

#Zoo
zoo=t[,15:ncol(t)]

#Mapping the stations
#When I see such a pattern,
#I really wonder how stations
#could have been ascribed to
#a respective "estuary"!...
par(mfrow=c(2,2))
for(i in 1:4){
  plot(xy,col="grey",main=levels(est)[i],cex.main=1.5)
  points(xy[est==levels(est)[i],],pch=20,cex=1.5)
  map("worldHires",add=T)
}

#Mapping abiotic data
#Large white squares, low values
#Small squares, intermediate values
#Large squares, high values
#Remark: a priori, is "relative humidity"
#a pertinent variable for aquatic organisms?
par(mfrow=c(3,3),mar=c(2,2,3,1))
for(i in 1:ncol(abio)){
  plot(xy,type="n",main=colnames(abio)[i],cex.main=1.5)
  map("worldHires",add=T)
  s.value(xy,scalewt(abio)[,i],cleg=0,add.p=T)
}

#PCA of zoo
#Keep 3 axes
#"scale=F", no reduction, better for
#matrices with a unique measurement unit
#"scale=T", systematically in the opposite case
#Here, keep 3 axes
pcaz=dudi.pca(log(zoo+1),scale=F)

#Eigenvalues expressed in %
100*pcaz$eig/sum(pcaz$eig)

#Zoo variables
s.arrow(pcaz$co,clab=1)

#Zoo stations
s.label(pcaz$li)
s.class(pcaz$li,est)
s.label(pcaz$li,add.p=T)

#Together
par(mfrow=c(2,1))
s.arrow(pcaz$co,clab=1)
s.class(pcaz$li,est)

#Correlation of abiotic variables?
cor(abio,pcaz$li)
par(mfrow=c(3,3))
for(i in 1:ncol(abio)){
  s.value(pcaz$li,scalewt(abio[,i]),
          sub=colnames(abio)[i],csub=3,possub="topleft")
}

#The same with axes 1 and 3
#NO2 seems to be the only one correlated to Axis 3
par(mfrow=c(3,3))
for(i in 1:ncol(abio)){
  s.value(pcaz$li[,c(1,3)],scalewt(abio[,i]),
          sub=colnames(abio)[i],csub=3,possub="topleft")
}

#PCA on Instrumental Variables (PCAIV)
#In this numerical context, you can say that PCAIV <=> RDA
#Keep 3 axes
iv=pcaiv(pcaz,abio)
iv

#Proportion of explained variance (R squared)
sum(iv$eig)/sum(pcaz$eig)

#Global summary (axes 1 and 2)
plot(iv)

#Global summary (axes 1 and 3)
plot(iv,yax=3)

#Predictions
iv$li

#Observations
iv$ls

#Discripancies
#ls = "lignes suppl?mentaires
#They results from the passive projection
#of the raw zoo data onto the axes (arrow tips)
#Arrow length = lack of fitting = residual
s.match(iv$li,iv$ls)

#Abio variables
#Salinity is the most influential variable
s.corcircle(iv$cor)

#Permutation test
test=randtest(iv,999)
test
plot(test)

#Mapping the three axes
#This enables to check a possible spatialization
#of the three gradients
par(mfrow=c(3,1),mar=c(2,2,3,1))
for(i in 1:3){
  plot(xy,type="n",main=paste ("Axis",i,sep=" "),cex.main=1.5)
  map("worldHires",add=T)
  s.value(xy,iv$li[,i],cleg=0,add.p=T)
  box()
}

#Do distances to the estuaries explain
#this correlation between abio and zoo?
iv.dis=pcaiv(iv,dis)
plot(iv.dis)

#Distance act significantly
randtest(iv.dis,999)

par(mfrow=c(2,1))
s.corcircle(iv.dis$cor,box=T)
#Shledt and Meuse are far from Channel and Thames (left)
#and reciprocally (right)
s.class(iv.dis$li,est)

#Abiotic variables
iv.abio=pcaiv(dudi.pca(abio,scan=F),dis,scan=F)
sup=t(iv.abio$tab)%*%as.matrix(iv.dis$l1)
s.arrow(sup)

#Mapping
par(mfrow=c(2,1),mar=c(2,2,3,1))
for(i in 1:2){
  plot(xy,type="n",main=paste ("Axis",i,sep=" "),cex.main=1.5)
  map("worldHires",add=T)
  s.value(xy,iv.dis$li[,i],cleg=0,add.p=T)
}

#I do not feel confortable with
#this distance approach because
#(1) I do not agree with the "estuarine" grouping
#(2) I do not know how these distances are computed
#(3) Marine connectivity can be very complex (currents)
#(4) Every ecological pattern can results from
#    processes expressed at different scales
#I suggest to simply forget these distances
#and to find other distances that could,
#a posteriori, explain the observed ecological patterns

require(adegraphics)
require(spdep)
require(adespatial)
require(vegan)
require(packfor)


#A relatively recent development generating
#spatial predictors hierarchically
#organized from large to small scale
#They are build simply from the xy geographic coordinates
#A specifically weighted distance matrix among stations is computed
#and spatial predictors are derived
#They are called Moran's Eigenvector Maps (MEM)
#They are perfectly independent (orthogonal)

#Choosing a connection network among stations
cn=chooseCN(xy,res="listw",type=1,plot.nb=F)
cn<-connection.network
#Triangulation is the simplest one
plot(cn,xy)
map("worldHires",add=T)

#Note that some distance crossing lands could be removed
#You can remove them
nb=tri2nb(xy)
#The map becomes interactive
#Remove the links crossing the lands
#by clicking on the connected stations
nb=edit(nb,xy) ##R loopt hierop vast

#Plot again
cn=nb2listw(nb,style="W")
plot(cn,xy)
map("worldHires",add=T)

#MEM computation
#There are 38 stations, so 38-1=37 spatial predictors
#They represent all the possible scales
#that can be investigated in this data set
umem=scores.listw(cn)
#Test the significance of Moran's I for MEMs
moran.randtest <- function(x, listw, nrepet = 999, ... ){
  
  if(missing(listw) & inherits(x, "orthobasisSp"))
    
    if(!is.null(attr(x, "listw")))
      
      listw <- attr(x, "listw")
    
    
    
    x <- as.matrix(x)
    
    if(!is.numeric(x))
      
      stop("x should contain only numeric values")
    
    
    
    if(NCOL(x) == 1){
      
      res <- moran.mc(x, listw = listw, nsim = nrepet, zero.policy = TRUE)
      
      res <- as.randtest(obs = res$statistic, sim = res$res[-(nrepet + 1)], call = match.call(), ...)
      
    } else {
      
      res <- apply(x, 2, moran.mc, listw = listw, nsim = nrepet, zero.policy = TRUE)
      
      res <- as.krandtest(obs = sapply(res, function(x) x$statistic), sim = sapply(res, function(x) x$res[-(nrepet + 1)]), call = match.call(), ...)
      
    }
    
    
    
    return(res)
    
    
    
}
memtest=moran.randtest(umem,cn,999,"two-sided")
#Keeping only MEMS with significant values
Umem=umem[,memtest[[7]]<0.05]

#Test of zoo vs xy correlation
randtest(pcaiv(pcaz,xy,scan=F),999)
#Test of abio vs xy correlation
randtest(pcaiv(dudi.pca(abio,scan=F),xy,scan=F),999)
#They are strongly correlated to xy
#suggesting that processes expressed at larger
#scale than the scale of the studied area
#affect the data
#Hence, this influence must be removed
#by detrending data from xy (removing xy effect)
#through Orthogonal PCAIV
zoo.o=pcaivortho(pcaz,xy,scan=F)
abio.o=pcaivortho(dudi.pca(abio,scan=F),xy,scan=F)

#The detrended PCAIV
iv=pcaiv(zoo.o,abio.o$tab)
randtest(iv,999)
#https://cran.r-project.org/web/packages/adespatial/vignettes/tutorial.html
#Finding spatial predictors of iv
#28 % of iv variance is explained by space (last line, AdjR2Cum)
##is anders bij mij???????????????????
rda.mem=rda(iv$tab,Umem)
R2a.mem=RsquareAdj(rda.mem)$adj.r.squared
fw.mem=forward.sel(iv$tab,Umem,adjR2thresh=R2a.mem,99999)
fw.mem
fw.mem=fw.mem[order(fw.mem$order),]
#F values optellen? voor 28% te krijgen 10+6+12


#Selected predictors ##for each station
pre.mem=data.frame(Umem[,colnames(Umem)%in%fw.mem$variables==T])
colnames(pre.mem)=fw.mem$variables
pre.mem=data.frame(pre.mem[,order(fw.mem$R2,decreasing=T)])

#pre.mem=Umem[,colnames(Umem)%in%fw.mem$variables==T]
#pre.mem=pre.mem[,order(fw.mem$R2,decreasing=T)]
#Represent them
#These predictors evidence processes
#(1) dominantly along waves of circa 50 km (MEM5, best predictor)
#    (alternation of white and black squares)
#(2) slightly along waves of similar length (MEM4, MEM7 and MEM8) 
##wat bedoelt hij hiermee????
detach('package:adegraphics')
detach('package:spdep')
detach('package:adespatial')
detach('package:vegan')
detach('package:packfor')

par(mfrow=c(2,2),mar=c(3,3,3,3))
for(i in 1:ncol(pre.mem)){
  plot(xy,type="n",main=colnames(pre.mem)[i],
       cex.main=1.5,bty="n",xaxt="n",yaxt="n")
  map("worldHires",add=T)
  s.value(xy,pre.mem[,i],cleg=0,csize=1.5,add.p=T)
}

#PCAIV on MEMs
iv.mem=pcaiv(iv,pre.mem)
plot(iv.mem)

#Projection of abiotic variables
w=pcaiv(dudi.pca(abio.o$tab,scan=F),pre.mem,scan=F)
#The data (iv.abio$tab) are projected onto the iv.mem system
sup=t(w$tab)%*%as.matrix(iv.mem$l1)
s.arrow(sup)

#The spatialization of the ecological pattern is mostly
#explained along a first axis by MEM5 and MEM7
#They explain the increasing concentration of
#Zoopl_Appendicularia and Zoopl_Noctiluca densities
#with decreasing salinity and increasing N concentration
#in coastal areas, from the center of the area to the coasts
par(mfrow=c(2,3),mar=c(3,3,3,3))
for(i in 1:2){
  plot(xy,type="n",main=colnames(pre.mem)[i],
       cex.main=1.5,bty="n",xaxt="n",yaxt="n")
  map("worldHires",add=T)
  s.value(xy,pre.mem[,i],cleg=0,csize=1.5,add.p=T)
}
s.corcircle(iv.mem$cor,clab=2)
s.arrow(sup,clab=2,xlim=c(-50,40))
s.arrow(iv.mem$c1,clab=2,xlim=c(-1.2,0.4))

#You should obtain the following graph

