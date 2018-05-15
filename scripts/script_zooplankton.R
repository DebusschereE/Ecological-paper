rm(list=ls())
require(ade4)
require(maps)
require(mapdata)
require(dplyr)

#Data
t=read.csv("data/Cruise_data_2017_V10_complete_Flowcam_NGS.csv",h=T,row.names=1)

t$N.Si<-(t$Nutr_NO2+t$Nutr_NO3+t$Nutr_NOX)/t$Nutr_SiO2
t$N.P<-(t$Nutr_NO2+t$Nutr_NO3+t$Nutr_NOX)/t$Nutr_PO4
t$N.P.Si<-(t$Nutr_NO2+t$Nutr_NO3+t$Nutr_NOX)/t$Nutr_PO4/t$Nutr_SiO2

colnames(t)
#select dataframe geographic coordinates, zooplankton, abiotic parameters
#
df=t[,c(3:4, 5:10, 166:168,72,75:76 ,127, 43:64)]
#remove NA rows for zooplankton
df<-na.omit(df)
#Zoo
colnames(df)
zoo=df[,16:ncol(df)]
#check for high number of zeros per column
#colSums(zoo == 0)
#colnames(zoo)
#zoo<-zoo[,c(2,4,5,7, 13:16,19:20)]
#Geographic coordinates
xy=df[,c(2,1)]

#Abiotic
abio=df[,3:15]


#Mapping abiotic data
#Large white squares, low values
#Small squares, intermediate values
#Large squares, high values
#Remark: a priori, is "relative humidity"
#a pertinent variable for aquatic organisms?
par(mfrow=c(5,3),mar=c(2,2,3,1))
for(i in 1:ncol(abio)){
  plot(xy,type="n",main=colnames(abio)[i],cex.main=1.5)
  map("worldHires",add=T)
  s.value(xy,scalewt(abio)[,i],cleg=0,add.p=T)
}



#zoo
par(mfrow=c(5,3),mar=c(2,2,3,1))
for(i in 1:ncol(zoo)){
  plot(xy,type="n",main=colnames(zoo)[i],cex.main=1.5)
  map("worldHires",add=T)
  s.value(xy,scalewt(zoo)[,i],cleg=0,add.p=T)
}

#PCA of zoo
#Keep 3 axes
#"scale=F", no reduction, better for
#matrices with a unique measurement unit
#"scale=T", systematically in the opposite case
#Here, keep 3 axes
pcaz=dudi.pca(log(zoo+1),scale=F, scannf=FALSE, nf=3)

#Eigenvalues expressed in %
100*pcaz$eig/sum(pcaz$eig)
dev.off()
#Zoo variables
s.arrow(pcaz$co,clab=0.5)

#Zoo stations
s.label(pcaz$li)
s.label(pcaz$li,add.p=T)

#Correlation of abiotic variables?
cor(abio,pcaz$li)
##Salinity has the highest correlation with axis 1
##NO2 seems to be correlated with Axis 2
##
#visualisations:
#axes 1 and 2:
par(mfrow=c(3,3))
for(i in 1:ncol(abio)){
  s.value(pcaz$li,scalewt(abio[,i]),
          sub=colnames(abio)[i],csub=3,possub="topleft")
}

#The same with axes 1 and 3
par(mfrow=c(3,3))

for(i in 1:ncol(abio)){
  s.value(pcaz$li[,c(1,3)],scalewt(abio[,i]),
          sub=colnames(abio)[i],csub=3,possub="topleft")
}

#PCA on Instrumental Variables (PCAIV)
#In this numerical context, you can say that PCAIV <=> RDA
#Keep 3 axes
iv=pcaiv(pcaz,abio,scannf=FALSE, nf=3)
iv

#Proportion of explained variance (R squared)
sum(iv$eig)/sum(pcaz$eig)
##60% explained
#Global summary (axes 1 and 2)
plot(iv)

#Global summary (axes 1 and 3)
plot(iv,yax=3)
dev.off()
scatter(iv,clab.col =0.5)
#Predictions
iv$li

#Observations
iv$ls

#Discripancies
#ls = "lignes suppl?mentaires
#They results from the passive projection
#of the raw zoo data onto the axes (arrow tips)
#Arrow length = lack of fitting = residual
dev.off()
s.match(iv$li,iv$ls, clab=0.5)
s.corcircle(iv$cor,clab=1)
s.arrow(iv$c1,clab=0.55,xlim=c(-0.1,0.2))
s.arrow(iv$c1,clab=0.55,xlim=c(-0.1,0.2),ylim=c(-0.1,0.2))



#Abio variables
#Salinity is the most influential variable
s.corcircle(iv$cor)

dev.off()
#Permutation test
test=randtest(iv,999)
test
#significant permutation test: siginficant relationship between zooplankton and abiotic parameters
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
##axis1: determined by salinity;depth, current velocity (white boxes) and T, NH4 and SiO2 black boxes)
##axis 2 : nutrients at the english coast
##axis 3: also by nutriens, (temp and current velocity)

###spatial analysis ####

require(adegraphics)
require(spdep)
require(adespatial)
require(vegan)
install.packages("packfor", repos="http://R-Forge.R-project.org") 
require(packfor)
require(ade4)

#A relatively recent development generating
#spatial predictors hierarchically
#organized from large to small scale
#They are build simply from the xy geographic coordinates
#A specifically weighted distance matrix among stations is computed
#and spatial predictors are derived
#They are called Moran's Eigenvector Maps (MEM)
#They are perfectly independent (orthogonal)
library(spdep)

#Choosing a connection network among stations
load("data/chooseCN.Rdata")
cn=chooseCN(xy,res="listw",type=1,plot.nb=F)
cn
#cn<-connection.network
#Triangulation is the simplest one
dev.off()
detach('package:purrr')
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
cnzoo<-cn
save(cnzoo, file="cnzoo.RData")
load("data/cnzoo.RData")
cn<-cnzoo
plot(cn,xy)
map("worldHires",add=T)
library(adespatial)
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
names(memtest)
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
randtest(iv,999)##significant
#https://cran.r-project.org/web/packages/adespatial/vignettes/tutorial.html
#Finding spatial predictors of iv
#43.5 % of iv variance is explained by space (last line, AdjR2Cum)
#7 spatial predictors
rda.mem=rda(iv$tab,Umem)
summary(rda.mem)
R2a.mem=RsquareAdj(rda.mem)$adj.r.squared
fw.mem=forward.sel(iv$tab,Umem,adjR2thresh=R2a.mem,99999)
fw.mem
fw.mem=fw.mem[order(fw.mem$order),]##selection of MEM 


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
require(purrr)
par(mfrow=c(4,2),mar=c(3,3,3,3))
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
dev.off()

#The spatialization of the ecological pattern is mostly
#explained along a first axis by MEM5 and MEM7
#They explain the increasing concentration of
#Zoopl_Appendicularia and Zoopl_Noctiluca densities
#with decreasing salinity, depth and increasing N concentration and T
#in coastal areas, from the center of the area to the coasts
##a second pattern expressed along the second axis (MEM4 and MEM3), explain the decreasing current velocity, and align with the echinodermata whilst harpacticoida like high current velocity.

par(mfrow=c(2,3),mar=c(3,3,3,3))
for(i in 1:7){
  plot(xy,type="n",main=colnames(pre.mem)[i],
       cex.main=1.5,bty="n",xaxt="n",yaxt="n")
  map("worldHires",add=T)
  s.value(xy,pre.mem[,i],cleg=0,csize=1.5,add.p=T)
}
dev.off()
s.corcircle(iv.mem$cor,clab=1)
s.arrow(sup,clab=1,xlim=c(-50,40))
s.arrow(iv.mem$c1,clab=0.8,xlim=c(-0.7,1.3))









# Function chooseCN
#####################


#' Function to choose a connection network
#'
#' The function \code{chooseCN} is a simple interface to build a connection
#' network (CN) from xy coordinates. The user chooses from 6 types of graph and
#' one additional weighting scheme.  \code{chooseCN} calls functions from
#' appropriate packages, handles non-unique coordinates and returns a
#' connection network either with classe \code{nb} or \code{listw}. For graph
#' types 1-4, duplicated locations are not accepted and will issue an error.
#'
#' There are 7 kinds of graphs proposed: \cr Delaunay triangulation (type 1)\cr
#' Gabriel graph (type 2)\cr Relative neighbours (type 3)\cr Minimum spanning
#' tree (type 4)\cr Neighbourhood by distance (type 5)\cr K nearests neighbours
#' (type 6)\cr Inverse distances (type 7)\cr
#'
#' The last option (type=7) is not a true neighbouring graph: all sites are
#' neighbours, but the spatial weights are directly proportional to the
#' inversed spatial distances.\cr Also not that in this case, the output of the
#' function is always a \code{listw} object, even if \code{nb} was
#' requested.\cr
#'
#' The choice of the connection network has been discuted on the adegenet
#' forum. Please search the archives from adegenet website (section 'contact')
#' using 'graph' as keyword.
#'
#' @param xy an matrix or data.frame with two columns for x and y coordinates.
#' @param ask a logical stating whether graph should be chosen interactively
#' (TRUE,default) or not (FALSE). Set to FALSE if \code{type} is provided.
#' @param type an integer giving the type of graph (see details).
#' @param result.type a character giving the class of the returned object.
#' Either "nb" (default) or "listw", both from \code{spdep} package. See
#' details.
#' @param d1 the minimum distance between any two neighbours. Used if
#' \code{type=5.}
#' @param d2 the maximum distance between any two neighbours. Used if
#' \code{type=5}. Can also be a character: "dmin" for the minimum distance so
#' that each site has at least one connection, or "dmax" to have all sites
#' connected (despite the later has no sense).
#' @param k the number of neighbours per point. Used if \code{type=6}.
#' @param a the exponent of the inverse distance matrix. Used if \code{type=7}.
#' @param dmin the minimum distance between any two distinct points. Used to
#' avoid infinite spatial proximities (defined as the inversed spatial
#' distances). Used if \code{type=7}.
#' @param plot.nb a logical stating whether the resulting graph should be
#' plotted (TRUE, default) or not (FALSE).
#' @param edit.nb a logical stating whether the resulting graph should be
#' edited manually for corrections (TRUE) or not (FALSE, default).
#' @param check.duplicates a logical indicating if duplicate coordinates should be detected; this can be an issue for some graphs; TRUE by default.
#'
#' @return Returns a connection network having the class \code{nb} or
#' \code{listw}. The xy coordinates are passed as attribute to the created
#' object.
#'
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#'
#' @seealso \code{\link{spca}}
#' @keywords spatial utilities
#' @examples
#'
#' \dontrun{
#' data(nancycats)
#'
#' par(mfrow=c(2,2))
#' cn1 <- chooseCN(nancycats@@other$xy,ask=FALSE,type=1)
#' cn2 <- chooseCN(nancycats@@other$xy,ask=FALSE,type=2)
#' cn3 <- chooseCN(nancycats@@other$xy,ask=FALSE,type=3)
#' cn4 <- chooseCN(nancycats@@other$xy,ask=FALSE,type=4)
#' par(mfrow=c(1,1))
#' }
#'
#' @export chooseCN
#' @importFrom spdep "tri2nb" "gabrielneigh" "graph2nb" "relativeneigh" "dnearneigh" "knearneigh" "knn2nb" "nb2listw" "mat2listw" "listw2mat" "lag.listw" "card"
#' @import ade4
#'
chooseCN <- function(xy, ask = TRUE, type = NULL, result.type = "nb",
                     d1 = NULL, d2 = NULL, k = NULL, a = NULL,
                     dmin = NULL, plot.nb = TRUE, edit.nb = FALSE,
                     check.duplicates = TRUE){
  
  if(is.data.frame(xy)) xy <- as.matrix(xy)
  if(ncol(xy) != 2) stop("xy does not have two columns.")
  if(any(is.na(xy))) stop("NA entries in xy.")
  result.type <- tolower(result.type)
  if(is.null(type) & !ask) stop("Non-interactive mode but no graph chosen; please provide a value for 'type' argument.")
  
  ## if(!require(spdep, quietly=TRUE)) stop("spdep library is required.")
  
  res <- list()
  
  if(!is.null(d2)){
    if(d2=="dmin"){
      tempmat <- as.matrix(dist(xy))
      d2min <- max(apply(tempmat, 1, function(r) min(r[r>1e-12])))
      d2min <- d2min * 1.0001 # to avoid exact number problem
      d2 <- d2min
    } else if(d2=="dmax"){
      d2max <- max(dist(xy))
      d2max <- d2max * 1.0001 # to avoid exact number problem
      d2 <- d2max
    }
  } # end handle d2
  
  d1.first <- d1
  d2.first <- d2
  k.first <- k
  
  ## handle type argument
  if(!is.null(type)){
    type <- as.integer(type)
    if(type < 1 |type > 7) stop("type must be between 1 and 7")
    ask <- FALSE
  }
  
  ## check for uniqueness of coordinates
  if(check.duplicates && any(xyTable(xy)$number>1)){ # if duplicate coords
    DUPLICATE.XY <- TRUE
  } else {
    DUPLICATE.XY <- FALSE
  }
  
  
  ## if(is.null(type) & !ask) { type <- 1 }
  
  ### begin large while ###
  chooseAgain <- TRUE
  while(chooseAgain){
    # re-initialisation of some variables
    d1 <- d1.first
    d2 <- d2.first
    k <- k.first
    
    ## read type from console
    if(ask){
      temp <- TRUE
      while(temp){
        cat("\nChoose a connection network:\n")
        cat("\t Delaunay triangulation (type 1)\n")
        cat("\t Gabriel graph (type 2)\n")
        cat("\t Relative neighbours (type 3)\n")
        cat("\t Minimum spanning tree (type 4)\n")
        cat("\t Neighbourhood by distance (type 5)\n")
        cat("\t K nearest neighbours (type 6)\n")
        cat("\t Inverse distances (type 7)\n")
        cat("Answer: ")
        
        type <- as.integer(readLines(con = getOption('adegenet.testcon'), n = 1))
        temp <- type < 1 |type > 7
        if(temp) cat("\nWrong answer\n")
        
        if(type %in% 1:4 & DUPLICATE.XY){
          cat("\n\n== PROBLEM DETECTED ==")
          cat("\nDuplicate locations detected\nPlease choose another graph (5-7) or add random noise to locations (see ?jitter).\n")
          temp <- TRUE
        }
        
      } # end while
    }
    ##
    
    ## warning about duplicate xy coords
    if(type %in% 1:4 & DUPLICATE.XY){
      stop("Duplicate locations detected and incompatible with graph type 1-4.\nPlease choose another graph (5-7) or add random noise to locations (see ?jitter).")
    }
    
    ## graph types
    ## type 1: Delaunay
    if(type==1){
      ## if(!require(tripack, quietly=TRUE)) stop("tripack library is required.")
      cn <- tri2nb(xy)
    }
    
    # type 2: Gabriel
    if(type==2){
      cn <- gabrielneigh(xy)
      cn <- graph2nb(cn, sym=TRUE)
    }
    
    ## type 3: Relative neighbours
    if(type==3){
      cn <- relativeneigh(xy)
      cn <- graph2nb(cn, sym=TRUE)
    }
    
    ## type 4: Minimum spanning tree
    if(type==4){
      cn <- ade4::mstree(dist(xy)) # there is also a spdep::mstree
      cn <- neig2nb(cn)
    }
    
    ## type 5: Neighbourhood by distance
    if(type==5){
      if(is.null(d1) |is.null(d2)){
        tempmat <- as.matrix(dist(xy))
        d2min <- max(apply(tempmat, 1, function(r) min(r[r>1e-12])))
        d2min <- d2min * 1.0001 # to avoid exact number problem
        d2max <- max(dist(xy))
        d2max <- d2max * 1.0001 # to avoid exact number problem
        dig <- options("digits")
        options("digits=5")
        cat("\n Enter minimum distance: ")
        d1 <- as.numeric(readLines(con = getOption('adegenet.testcon'), n = 1))
        cat("\n Enter maximum distance \n(dmin=", d2min, ", dmax=", d2max, "): ")
        d2 <- readLines(con = getOption('adegenet.testcon'), n = 1)
        ## handle character
        if(d2=="dmin") {
          d2 <- d2min
        } else if(d2=="dmax") {
          d2 <- d2max
        } else {
          d2 <- as.numeric(d2)
        }
        ## restore initial digit option
        options(dig)
      }
      # avoid that a point is its neighbour
      dmin <- mean(dist(xy))/100000
      if(d1<dmin) d1 <- dmin
      if(d2<d1) stop("d2 < d1")
      cn <- dnearneigh(x=xy, d1=d1, d2=d2)
    }
    
    ## type 6: K nearests
    if(type==6){
      if(is.null(k)) {
        cat("\n Enter the number of neighbours: ")
        k <- as.numeric(readLines(con = getOption('adegenet.testcon'), n = 1))
      }
      cn <- knearneigh(x=xy, k=k)
      cn <- knn2nb(cn, sym=TRUE)
    }
    
    ## type 7: inverse distances
    if(type==7){
      if(is.null(a)) {
        cat("\n Enter the exponent: ")
        a <- as.numeric(readLines(con = getOption('adegenet.testcon'), n = 1))
      }
      cn <- as.matrix(dist(xy))
      if(is.null(dmin)) {
        cat("\n Enter the minimum distance \n(range = 0 -", max(cn),"): ")
        dmin <- as.numeric(readLines(con = getOption('adegenet.testcon'), n = 1))
      }
      if(a<1) { a <- 1 }
      thres <- mean(cn)/1e8
      if(dmin < thres) dmin <- thres
      cn[cn < dmin] <- dmin
      cn <- 1/(cn^a)
      diag(cn) <- 0
      cn <- prop.table(cn,1)
      plot.nb <- FALSE
      edit.nb <- FALSE
      result.type <- "listw"
    } # end type 7
    
    ## end graph types
    
    if(ask & plot.nb) {
      plot(cn,xy)
      cat("\nKeep this graph (y/n)? ")
      ans <- tolower(readLines(con = getOption('adegenet.testcon'), n=1))
      if(ans=="n") {chooseAgain <- TRUE} else {chooseAgain <- FALSE}
    }
    else if(plot.nb){
      plot(cn,xy)
      chooseAgain <- FALSE
    }
    else {chooseAgain <- FALSE}
    
  }
  ### end large while
  
  if(edit.nb) {cn <- edit(cn,xy)}
  
  if(result.type == "listw") {
    if(type!=7) {
      cn <- nb2listw(cn, style="W", zero.policy=TRUE)
    } else {
      cn <- mat2listw(cn)
      cn$style <- "W"
    }
  }
  
  res <- cn
  
  attr(res,"xy") <- xy
  
  return(res)
  
} # end chooseCN
save("chooseCN", file="chooseCN.Rdata")