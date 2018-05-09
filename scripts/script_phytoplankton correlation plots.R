### Phytoplankton ####
t=read.csv("Cruise_data_2017_V9_complete_Flowcam_NGS.csv",h=T,row.names=1)
install.packages("easypackages")
library(easypackages)

# Required packages
libraries("Hmisc", "corrplot", "PerformanceAnalytics", "ggplot2", "vegan", "tidyr", "data.table", "vegan","dplyr", "knitr", "magrittr", "checkmate","Hmisc")

#select dataframe: chemtax, ngs, flowcam, fcm (rws, cnrs, vliz) parameters
colnames(t)
df=t[,c(83:112,130:165)]

colnames(df)
### compare NGS and Chemtax ####

dff<-df[,c(25:30, 60:66)]
#remove NA rows for plankton

res <- rcorr(as.matrix(dff,type="spearman"))
res$r
res$P

corrplot(res$r, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45)
corrplot(res$r, type="upper", order="hclust", 
         p.mat = res$P, sig.level = 0.05, insig = "blank")
chart.Correlation(dff, histogram=TRUE, pch=19)
col<- colorRampPalette(c("blue", "white", "red"))(20)
heatmap(x = res$r, col = col, symm = TRUE, main = "Correlation heatmap")


