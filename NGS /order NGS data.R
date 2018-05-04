rm(list=ls())

# Load packages
p = c("ggplot2", "RColorBrewer", "dplyr", "plotly", "vegan", "plyr", "sp", "ggmap",  "broom", "rgdal", "xlsx")
# lapply(p, install.packages)
lapply(p, require, character = T)
rm(p)

# Load NGS data ####
setwd("~/Jerico next 2017 cruise/NGS ")
otu_18s = read.table("otu_table_18S_max1.txt", header = T, sep = "")
otu_18s[, c(2:53)] = as.data.frame(lapply(otu_18s[, c(2:53)], as.numeric))
rownames(otu_18s) = otu_18s[ ,1]
otu_18s = otu_18s[ ,-1]
otu_18s.sub = otu_18s[,colSums(otu_18s) > 1000]
otu_18s.sub$OTUId <- as.factor(rownames(otu_18s.sub))
otu_18s.sub = otu_18s.sub [ , c(48, 2:47)]
otu_16s = read.table("otu_table_16S_max1.txt", header = T, sep = "")
otu_16s[ , c(2:53)] = lapply(otu_16s[ , c(2:53)], as.numeric)
tax_18s = read.table("tax18.csv", header = T, sep = ";")
tax_16s = read.table("tax16.csv", header = T, sep = ";")
# Data manipulation
merged_18s = merge(otu_18s.sub, tax_18s, by = "OTUId")
rm(otu_18s)
rm(tax_18s)
rm(otu_18s.sub)
head(merged_18s[0,])
#names(merged_18s)[names(merged_18s) == "R1"] = "st.43"
#names(merged_18s)[names(merged_18s) == "R2"] = "st.14"
names(merged_18s)[names(merged_18s) == "R3"] = "st.35"
#names(merged_18s)[names(merged_18s) == "R4"] = "st.24"
names(merged_18s)[names(merged_18s) == "R5"] = "st.19" 
names(merged_18s)[names(merged_18s) == "R6"] = "st.7"
#names(merged_18s)[names(merged_18s) == "R7"] = "st.16"
#names(merged_18s)[names(merged_18s) == "R8"] = "st.10"
names(merged_18s)[names(merged_18s) == "R9"] = "st.22.26.1"
#names(merged_18s)[names(merged_18s) == "R10"] = "st.21"
names(merged_18s)[names(merged_18s) == "R11"] = "st.40"
names(merged_18s)[names(merged_18s) == "R12"] = "st.41"
names(merged_18s)[names(merged_18s) == "R13"] = "st.25"
names(merged_18s)[names(merged_18s) == "R14"] = "st.38"
names(merged_18s)[names(merged_18s) == "R15"] = "st.6"
names(merged_18s)[names(merged_18s) == "R16"] = "st.33"
names(merged_18s)[names(merged_18s) == "R17"] = "st.44"
names(merged_18s)[names(merged_18s) == "R18"] = "st.32"
names(merged_18s)[names(merged_18s) == "R19"] = "st.11"
names(merged_18s)[names(merged_18s) == "R20"] = "st.5"
names(merged_18s)[names(merged_18s) == "R21"] = "st.28"
names(merged_18s)[names(merged_18s) == "R22"] = "st.37"
names(merged_18s)[names(merged_18s) == "R23"] = "st.34"
names(merged_18s)[names(merged_18s) == "R24"] = "st.29"
names(merged_18s)[names(merged_18s) == "R25"] = "st.13"
names(merged_18s)[names(merged_18s) == "R26"] = "st.31"
names(merged_18s)[names(merged_18s) == "R27"] = "st.36"
names(merged_18s)[names(merged_18s) == "R28"] = "st.9"
names(merged_18s)[names(merged_18s) == "R29"] = "st.3"
names(merged_18s)[names(merged_18s) == "R30"] = "st.27"
names(merged_18s)[names(merged_18s) == "R31"] = "st.42"
names(merged_18s)[names(merged_18s) == "R32"] = "st.8"
names(merged_18s)[names(merged_18s) == "R33"] = "st.17"
names(merged_18s)[names(merged_18s) == "R34"] = "st.2"
names(merged_18s)[names(merged_18s) == "R35"] = "st.12"
names(merged_18s)[names(merged_18s) == "R36"] = "st.1"
names(merged_18s)[names(merged_18s) == "R37"] = "st.22.26.3"
names(merged_18s)[names(merged_18s) == "R38"] = "st.30"
names(merged_18s)[names(merged_18s) == "R39"] = "st.15"
names(merged_18s)[names(merged_18s) == "R40"] = "st.39"
names(merged_18s)[names(merged_18s) == "R41"] = "st.23"
names(merged_18s)[names(merged_18s) == "R42"] = "st.4"
names(merged_18s)[names(merged_18s) == "R43"] = "st.20"
# name instead of number 
merged_18s = merged_18s[ , c("OTUId", "Empire", "Kingdom", "Phylum", "Classes", "Order",  "Class", "Family", "Genus", "st.1", "st.2", "st.3", "st.4", "st.5", "st.6", "st.7", "st.8", "st.9", "st.11", "st.12", "st.13", "st.15", "st.17", "st.19", "st.20", "st.22.26.1", "st.23", "st.25", "st.29", "st.27", "st.28", "st.22.26.3", "st.30", "st.31", "st.32", "st.33", "st.34", "st.35", "st.36", "st.37", "st.38", "st.39", "st.40", "st.41", "st.42", "st.44")]
#divide by filtered volume (stickertje and excel sheet give different volumes)
merged_18s$st.1 = merged_18s$st.1/200
merged_18s$st.2 = merged_18s$st.2/200
merged_18s$st.3 = merged_18s$st.3/100
merged_18s$st.4 = merged_18s$st.4/100
merged_18s$st.5 = merged_18s$st.5/75
merged_18s$st.6 = merged_18s$st.6/150
merged_18s$st.7 = merged_18s$st.7/100
merged_18s$st.8 = merged_18s$st.8/200
merged_18s$st.9 = merged_18s$st.9/200
merged_18s$st.10 = merged_18s$st.10/300
merged_18s$st.11 = merged_18s$st.11/100
merged_18s$st.12 = merged_18s$st.12/125
merged_18s$st.13 = merged_18s$st.13/125
merged_18s$st.14 = merged_18s$st.14/300
merged_18s$st.15 = merged_18s$st.15/300
merged_18s$st.16 = merged_18s$st.16/200
merged_18s$st.17 = merged_18s$st.17/400
merged_18s$st.19 = merged_18s$st.19/75
merged_18s$st.20 = merged_18s$st.20/300
merged_18s$st.21 = merged_18s$st.21/250
merged_18s$st.22.26.1 = merged_18s$st.22.26.1/200
merged_18s$st.23 = merged_18s$st.23/300
merged_18s$st.24 = merged_18s$st.24/200
merged_18s$st.25 = merged_18s$st.25/200
merged_18s$st.29 = merged_18s$st.29/150
merged_18s$st.27 = merged_18s$st.27/200
merged_18s$st.28 = merged_18s$st.28/200
merged_18s$st.22.26.3 = merged_18s$st.22.26.3/200
merged_18s$st.30 = merged_18s$st.30/200
merged_18s$st.31 = merged_18s$st.31/400
merged_18s$st.32 = merged_18s$st.32/100
merged_18s$st.33 = merged_18s$st.33/200
merged_18s$st.34 = merged_18s$st.34/200
merged_18s$st.35 = merged_18s$st.35/100
merged_18s$st.36 = merged_18s$st.36/200
merged_18s$st.37 = merged_18s$st.37/100
merged_18s$st.38 = merged_18s$st.38/50
merged_18s$st.39 = merged_18s$st.39/200
merged_18s$st.40 = merged_18s$st.40/100
merged_18s$st.41 = merged_18s$st.41/100
merged_18s$st.42 = merged_18s$st.42/100
merged_18s$st.43 = merged_18s$st.43/200
merged_18s$st.44 = merged_18s$st.44/200
# Subset 18s NGS on CHEMTAX groups ####
# Phaeocystis
Phaeocystis = subset(merged_18s, (Phylum == "Haptophyta"))
#Phaeocystis = Phaeocystis[ , c(10:52)]
Phaeocystis = Phaeocystis[ , c(10:46)]
Phaeocystis = as.data.frame(colSums(Phaeocystis))
names(Phaeocystis)[names(Phaeocystis) == "colSums(Phaeocystis)"] = "Phaeocystis_NGS"
# Diatoms                                                          
Diatoms = subset(merged_18s, (Classes == "Bacillariophyta"))
Diatoms = Diatoms[ , c(10:46)]
Diatoms = as.data.frame(colSums(Diatoms))
names(Diatoms)[names(Diatoms) == "colSums(Diatoms)"] = "Diatoms_NGS"
# Chlorophytes
Chlorophytes = subset(merged_18s, (Phylum == "Chlorophyta"))
Chlorophytes = Chlorophytes[ , c(10:46)]
Chlorophytes = as.data.frame(colSums(Chlorophytes))
names(Chlorophytes)[names(Chlorophytes) == "colSums(Chlorophytes)"] = "Chlorophytes_NGS"
# Dinoflagellates
Dinoflagellates = subset(merged_18s, (Phylum == "Dinophyta"))
Dinoflagellates = Dinoflagellates[ , c(10:46)]
Dinoflagellates = as.data.frame(colSums(Dinoflagellates))
names(Dinoflagellates)[names(Dinoflagellates) == "colSums(Dinoflagellates)"] = "Dinoflagellates_NGS"
# Euglenophytes
Euglenophytes = subset(merged_18s, (Order == "Euglyphida"))
Euglenophytes = Euglenophytes[ , c(10:46)]
Euglenophytes = as.data.frame(colSums(Euglenophytes))
names(Euglenophytes)[names(Euglenophytes) == "colSums(Euglenophytes)"] = "Euglenophytes_NGS"
# Cryptophytes
Cryptophytes = subset(merged_18s, (Phylum == "Cryptophyta"))
Cryptophytes = Cryptophytes[ , c(10:46)]
Cryptophytes = as.data.frame(colSums(Cryptophytes))
names(Cryptophytes)[names(Cryptophytes) == "colSums(Cryptophytes)"] = "Cryptophytes_NGS"
# Subset 16s NGS on CHEMTAX groups ####
merged_16s = merge(otu_16s, tax_16s, by = "OTUId")
rm(otu_16s)
rm(tax_16s)
head(merged_16s[0,])
names(merged_16s)[names(merged_16s) == "R1"] = "st.43"
names(merged_16s)[names(merged_16s) == "R2"] = "st.14"
names(merged_16s)[names(merged_16s) == "R3"] = "st.35"
names(merged_16s)[names(merged_16s) == "R4"] = "st.24"
names(merged_16s)[names(merged_16s) == "R5"] = "st.19" 
names(merged_16s)[names(merged_16s) == "R6"] = "st.7"
names(merged_16s)[names(merged_16s) == "R7"] = "st.16"
names(merged_16s)[names(merged_16s) == "R8"] = "st.10"
names(merged_16s)[names(merged_16s) == "R9"] = "st.22.26.1"
names(merged_16s)[names(merged_16s) == "R10"] = "st.21"
names(merged_16s)[names(merged_16s) == "R11"] = "st.40"
names(merged_16s)[names(merged_16s) == "R12"] = "st.41"
names(merged_16s)[names(merged_16s) == "R13"] = "st.25"
names(merged_16s)[names(merged_16s) == "R14"] = "st.38"
names(merged_16s)[names(merged_16s) == "R15"] = "st.6"
names(merged_16s)[names(merged_16s) == "R16"] = "st.33"
names(merged_16s)[names(merged_16s) == "R17"] = "st.44"
names(merged_16s)[names(merged_16s) == "R18"] = "st.32"
names(merged_16s)[names(merged_16s) == "R19"] = "st.11"
names(merged_16s)[names(merged_16s) == "R20"] = "st.5"
names(merged_16s)[names(merged_16s) == "R21"] = "st.28"
names(merged_16s)[names(merged_16s) == "R22"] = "st.37"
names(merged_16s)[names(merged_16s) == "R23"] = "st.34"
names(merged_16s)[names(merged_16s) == "R24"] = "st.29"
names(merged_16s)[names(merged_16s) == "R25"] = "st.13"
names(merged_16s)[names(merged_16s) == "R26"] = "st.31"
names(merged_16s)[names(merged_16s) == "R27"] = "st.36"
names(merged_16s)[names(merged_16s) == "R28"] = "st.9"
names(merged_16s)[names(merged_16s) == "R29"] = "st.3"
names(merged_16s)[names(merged_16s) == "R30"] = "st.27"
names(merged_16s)[names(merged_16s) == "R31"] = "st.42"
names(merged_16s)[names(merged_16s) == "R32"] = "st.8"
names(merged_16s)[names(merged_16s) == "R33"] = "st.17"
names(merged_16s)[names(merged_16s) == "R34"] = "st.2"
names(merged_16s)[names(merged_16s) == "R35"] = "st.12"
names(merged_16s)[names(merged_16s) == "R36"] = "st.1"
names(merged_16s)[names(merged_16s) == "R37"] = "st.22.26.3"
names(merged_16s)[names(merged_16s) == "R38"] = "st.30"
names(merged_16s)[names(merged_16s) == "R39"] = "st.15"
names(merged_16s)[names(merged_16s) == "R40"] = "st.39"
names(merged_16s)[names(merged_16s) == "R41"] = "st.23"
names(merged_16s)[names(merged_16s) == "R42"] = "st.4"
names(merged_16s)[names(merged_16s) == "R43"] = "st.20"
# name instead of number 
merged_16s = merged_16s[ , c("OTUId", "Domain", "Phylum", "Class", "Order", "Family", "Genus", "st.1", "st.2", "st.3", "st.4", "st.5", "st.6", "st.7", "st.8", "st.9", "st.10", "st.11", "st.12", "st.13", "st.14", "st.15", "st.16", "st.17", "st.19", "st.20", "st.21", "st.22.26.1", "st.23", "st.24", "st.25", "st.29", "st.27", "st.28", "st.22.26.3", "st.30", "st.31", "st.32", "st.33", "st.34", "st.35", "st.36", "st.37", "st.38", "st.39", "st.40", "st.41", "st.42", "st.43", "st.44")]
merged_16s$st.1 = merged_16s$st.1/200
merged_16s$st.2 = merged_16s$st.2/200
merged_16s$st.3 = merged_16s$st.3/100
merged_16s$st.4 = merged_16s$st.4/100
merged_16s$st.5 = merged_16s$st.5/75
merged_16s$st.6 = merged_16s$st.6/150
merged_16s$st.7 = merged_16s$st.7/100
merged_16s$st.8 = merged_16s$st.8/200
merged_16s$st.9 = merged_16s$st.9/200
merged_16s$st.10 = merged_16s$st.10/300
merged_16s$st.11 = merged_16s$st.11/100
merged_16s$st.12 = merged_16s$st.12/125
merged_16s$st.13 = merged_16s$st.13/125
merged_16s$st.14 = merged_16s$st.14/300
merged_16s$st.15 = merged_16s$st.15/300
merged_16s$st.16 = merged_16s$st.16/200
merged_16s$st.17 = merged_16s$st.17/400
merged_16s$st.19 = merged_16s$st.19/75
merged_16s$st.20 = merged_16s$st.20/300
merged_16s$st.21 = merged_16s$st.21/250
merged_16s$st.22.26.1 = merged_16s$st.22.26.1/200
merged_16s$st.23 = merged_16s$st.23/300
merged_16s$st.24 = merged_16s$st.24/200
merged_16s$st.25 = merged_16s$st.25/200
merged_16s$st.29 = merged_16s$st.29/150
merged_16s$st.27 = merged_16s$st.27/200
merged_16s$st.28 = merged_16s$st.28/200
merged_16s$st.22.26.3 = merged_16s$st.22.26.3/200
merged_16s$st.30 = merged_16s$st.30/200
merged_16s$st.31 = merged_16s$st.31/400
merged_16s$st.32 = merged_16s$st.32/100
merged_16s$st.33 = merged_16s$st.33/200
merged_16s$st.34 = merged_16s$st.34/200
merged_16s$st.35 = merged_16s$st.35/100
merged_16s$st.36 = merged_16s$st.36/200
merged_16s$st.37 = merged_16s$st.37/100
merged_16s$st.38 = merged_16s$st.38/50
merged_16s$st.39 = merged_16s$st.39/200
merged_16s$st.40 = merged_16s$st.40/100
merged_16s$st.41 = merged_16s$st.41/100
merged_16s$st.42 = merged_16s$st.42/100
merged_16s$st.43 = merged_16s$st.43/200
merged_16s$st.44 = merged_16s$st.44/200
Cyanobacteria = subset(merged_16s, (Class == "Cyanobacteria"))
Cyanobacteria = Cyanobacteria[ , c(8:50)]
Cyanobacteria = as.data.frame(colSums(Cyanobacteria))
names(Cyanobacteria)[names(Cyanobacteria) == "colSums(Cyanobacteria)"] = "Cyanobacteria_NGS"
# Join data into one dataframe ####
Phaeocystis$rn <- rownames(Phaeocystis)
Diatoms$rn <- rownames(Diatoms)
Chlorophytes$rn <- rownames(Chlorophytes)
Cryptophytes$rn <- rownames(Cryptophytes)
Dinoflagellates$rn <- rownames(Dinoflagellates)
Euglenophytes$rn <- rownames(Euglenophytes)
Cyanobacteria$rn <- rownames(Cyanobacteria)
NGS_groups <- join_all(list(Phaeocystis,Diatoms,Chlorophytes,Cryptophytes,Dinoflagellates,Euglenophytes,Cyanobacteria), by = 'rn', type = 'full')



