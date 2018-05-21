# Prepare and load datasets
# You can prepare your datasets in this script
# You will need a dataframe with your responses variable and your drivers
# If you work from multivariate data you need to transform them into univariate data using PCA or DCA for instance

# You can use our script with our Example data, but we also include below the code we used to clean it
# This code is in the 'prepare_datasets.R' script

## Example data ##
response <- read.delim(paste(getwd(),"/Input/Data/Bruel_etal_2018_response_example.txt",sep=""))
driver <- read.delim(paste(getwd(),"/Input/Data/Bruel_etal_2018_driver_example.txt",sep=""))

## Below is an example of code to tranform your data from multi- to uni-variate ##

#### Response(s) ####
# Load responses
clado <- read.delim(paste(getwd(),"/Input/Data/JOU16-02.clado.txt",sep=""))
clado$Year <- round(clado$Year,0)
colnames(clado) <- c("Year", "LK","SC","DL","EC","EL","BL","Eurycercus.spp.",
                    "C..rectirostris","A..harpae","A..elongata","A..affinis",
                    "A..quadrangularis","A..guttata","A..rectangula","A..costata",
                    "L..leydigi","G..testudinaria","R..falcata","M..dispar" ,"A..excisa",
                    "A..exigua","A..nana","P..trigonellus","P..uncinatus",
                    "C..sphaericus","C..gibbus","P..pigra","BYL")

# Ordination method
# If your data are already uni-variate you can skip this part but make sure to name your data file as requested below
str(clado);dim(clado)
clado$Year<-as.numeric(clado$Year)
datacp<-clado[,-c(1)]
row.names(datacp)<-as.numeric(clado$Year)
# ACP
datacp2<-decostand(datacp,"hellinger",na.rm=TRUE)
acp<-dudi.pca(datacp2,scannf=F,nf=10)
plot(acp$li[,1:2], type="l")
points(acp$li[,1:2], pch=20)
text(acp$li[,1:2], labels = rownames(acp$li), cex=0.6, pos = 3)
s.arrow(acp$c1, lab = colnames(datacp), grid = F)
summary(acp)
scores.pca<-as.data.frame(acp$li)
scores.species.pca<-as.data.frame(acp$co)

# DCA
datadca <- decorana(datacp)
plot(datadca)
plot(scores(datadca), asp=1, type="b", main="DCA")
scores.dca<-as.data.frame(datadca$rproj)
scores.species.dca<-as.data.frame(datadca$cproj)

# Two ordinations methods are proposed here
# You need to make sure you want to keep the same transformation as us ("hellinger")
# You can make your decision on which ordination method to use based on the percentage of variability explained on each axis
# Note that Principal Curves can be an interesting alternative when the species are distributed along a single gradient (see Bennion et al. 2015)

# Select the output of one of the ordination method
# Your responses data file must be named 'response' to run the code as smoothly as possible
response <- cbind(clado$Year,scores.pca)
colnames(response) <- c("Year",colnames(scores.pca))

#### Driver(s) ####
# Always check that the import worked. I use a specific format from mac + comma as separator
# Adapt the read.delim fonction
tp <- read.delim(paste(getwd(),"/Input/Data/JOU13-02_TP_sediment.txt",sep=""))
sata <- read.delim(paste(getwd(),"/Input/Data/Buntgen_2006_SAT_anomaly.txt",sep=""))

# If your data are multivariate, you can use ordination methods to transform them
# If the correspondance between x-axis of responses and drivers is not exactly the same, sub-sample the datasets
# To run the code, you need a value of driver(s) for each response

driver <- as.data.frame(matrix(rep(NA, nrow(clado)*3), ncol=3))
colnames(driver) <- c("Year", "TP", "SAT")
driver$Year <- response$Year
for (i in 1:nrow(driver)) {
  driver$TP[i] <- tp$TP[which(abs(driver$Year[i]-tp$Year)==min(abs(driver$Year[i]-tp$Year), na.rm=T))]
  driver$SAT[i] <- sata$SATa[which(abs(driver$Year[i]-sata$Year)==min(abs(driver$Year[i]-sata$Year), na.rm=T))]
}

#### Clean environment ####
rm(acp,datacp,datacp2,datadca,scores.pca,scores.species.pca,scores.dca,scores.species.dca,
   clado,tp,sata,i)

