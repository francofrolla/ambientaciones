
## 5. MULTIVARIATE SITE CLASSIFICATION ---------------------------------------------------

# After all the variables have been processed from step 1 up to the interpolation with the same prediction grid (step 5),
# the different data sets obtained should be concatenated using the cbind function.
# Below the R code was deactivated because a new database that has been concatenated is used.

# Pred <- cbind(PredECa30[,1:3], PredECa90[,3], PredElev[,3],PredSd[,3])
# names(Pred)[3]<-paste("ECa30")
# names(Pred)[4]<-paste("ECa90")
# names(Pred)[5]<-paste("Elev")
# names(Pred)[6]<-paste("Sd")

Pred <-read.table("C:\\.....\\Pred.txt", header = TRUE)

# Spatial principal component analysis (MULTISPATI-PCA)
pca <- dudi.pca(Pred[,3:6], center=T,scannf = FALSE,  nf = 5)

cord_1 <- coordinates(Pred[,1:2])
gri_1 <- dnearneigh(cord_1,0,25)
lw_1 <- nb2listw(gri_1, style = "W")

ms <- multispati(pca, lw_1, scannf = F, nfposi = 5)
s.arrow(ms$c1,xax = 1, yax = 2, clabel = 1)

# Extraction of spatial principal components
sPC <- ms$li[,1:4]
PredMA <- cbind(Pred,sPC) ;PredMA

#  Fuzzy k-means cluster analysis
MC_2<-cmeans(PredMA[,7:8],2,100,method="cmeans",m=1.3)
MC_3<-cmeans(PredMA[,7:8],3,100,method="cmeans",m=1.3)
MC_4<-cmeans(PredMA[,7:8],4,100,method="cmeans",m=1.3)

# Indices for selecting the number of classes: two (I2MC), three (I3MC) and four (I4MC)
I2MC <- fclustIndex(MC_2,PredMA[,7:8], index=c("xie.beni", "fukuyama.sugeno",
"partition.coefficient", "partition.entropy"))

I3MC <- fclustIndex(MC_3,PredMA[,7:8], index=c("xie.beni", "fukuyama.sugeno",
"partition.coefficient", "partition.entropy"))

I4MC <- fclustIndex(MC_4,PredMA[,7:8], index=c("xie.beni", "fukuyama.sugeno",
"partition.coefficient", "partition.entropy"))

Indices0 <- cbind(I2MC,I3MC,I4MC)

XieBeni <-Indices0[1,]
FukSug <-Indices0[2,]
PartCoef_1 <-Indices0[3,]
PartCoef <- 1/PartCoef_1
PartEntr <-Indices0[4,]

Indices <- as.data.frame(rbind(XieBeni,FukSug,PartCoef,PartEntr))
Indices

# Summary indices
XieBeniMax<-max(Indices[1,])
FukSugMax<-max(Indices[2,])
PartCoefMax<-max(Indices[3,])
PartEntrMax<-max(Indices[4,])

XieBeniN<- XieBeni/XieBeniMax
FukSugN<- FukSug/FukSugMax
PartCoefN<- PartCoef/PartCoefMax
PartEntrN<-PartEntr/PartEntrMax

IndicesN <- as.data.frame(rbind(XieBeniN,FukSugN,PartCoefN,PartEntrN))
IndicesN2 <- (IndicesN)^2

Indice2MC <- sqrt(sum(IndicesN2[,1]))
Indice3MC <- sqrt(sum(IndicesN2[,2]))
Indice4MC<- sqrt(sum(IndicesN2[,3]))

# Summary indices for selection of two, three or four management zones
Indice2MC; Indice3MC; Indice4MC

# Maps with management classes delimited
MC_22 <-as.data.frame(MC_2$cluster)
MC_33 <-as.data.frame(MC_3$cluster)
MC_44 <-as.data.frame(MC_4$cluster)

baseMC <- cbind(PredMA[,1:2],MC_22,MC_33,MC_44)

coordinates(baseMC) <- ~x+y
gridded(baseMC) <- T
spplot(baseMC["MC_2$cluster"],col.regions=gray.colors(2),colorkey = F)

spplot(baseMC["MC_3$cluster"],col.regions=gray.colors(10),colorkey = F)

spplot(baseMC["MC_4$cluster"],col.regions=gray.colors(4),colorkey = F)



## 6. SMOOTHING OF CLASSIFICATION RESULTS ------------------------------------------------

# Median filter function
smooth <-function(mytable,mywindow){
  newtable<-matrix(1:(dim(mytable)[1]*dim(mytable)[2]),dim(mytable)[1],dim(mytable)[2])
      vecinity<-function(pos) {
        col=as.integer((pos-1)/nrow(newtable))+1
        row=pos-((nrow(newtable)*col)-nrow(newtable))
        if (is.na(mytable[row,col])) NA else{
        myrow1<-ifelse(row-mywindow<1,1,row-mywindow)
        mycol1<-ifelse(col-mywindow<1,1,col-mywindow)
        myrow2<-ifelse(row+mywindow>dim(newtable)[1],row,row+mywindow)
        mycol2<-ifelse(col+mywindow>dim(newtable)[2],col,col+mywindow)

        neighbor<-na.omit(as.vector(mytable[myrow1:myrow2,mycol1:mycol2]))
        round(median(neighbor),digits=0)
      }}
      as.matrix(apply(newtable,c(1,2),vecinity))}

# Function to obtain a matrix
obtainM <- function(mytable){
  x<-as.numeric(names(table(mytable$x)))
  y<-as.numeric(names(table(mytable$y)))
  myframe <- matrix(1:(length(x)*length(y)), length(x), length(y))
  position<-function(pos) {
    col=as.integer((pos-1)/nrow(myframe))+1
    row=pos-((nrow(myframe)*col)-nrow(myframe))
    myindex=which(mytable$x==x[row] & mytable$y==y[col],arr.ind=T)
    if(length(myindex)==0) return(NA) else mytable[myindex,3]
  }
 thematrix<-as.matrix(apply(myframe,c(1,2),position))
 rownames(thematrix)<-x
 colnames(thematrix)<-y
 thematrix}

base0 <- cbind(PredMA[,1:2],MC_22)
datafilter <- obtainM(base0)

# Windows 5 x 5
smoot5x5 <- smooth(datafilter,5)

# Windows 7 x 7
smoot7x7 <- smooth(datafilter,7)

# Windows 9 x 9
smoot9x9 <- smooth(datafilter,9)

par(mfrow=c(2,2))
image(datafilter, main= "Original Zonification", axes = FALSE, xlab="",ylab="",col=palette(c("grey94","grey34")))
image(smoot5x5, main= "Median Filter 5 x 5",axes = FALSE, xlab="",ylab="",col=palette(c("grey94","grey34")))
image(smoot7x7, main= "Median Filter 7 x 7",axes = FALSE, xlab="",ylab="",col=palette(c("grey94","grey34")))
image(smoot9x9, main= "Median Filter 9 x 9",axes = FALSE, xlab="",ylab="",col=palette(c("grey94","grey34")))

# Data set with smoothing classification
# Function to ransform a matrix to table
MtoT <- function(mymatrix){
position <- function(ij){
data.frame(x=rownames(mymatrix)[ij[1]],y=colnames(mymatrix)[ij[2]],z=mymatrix[ij[1],ij[2]])
}
myindex <- which(!is.na(mymatrix),arr.ind=T);rownames(myindex)=NULL
b <- apply(myindex,1,position)
b <- do.call("rbind",b);b}

# New data set
base1 <- as.data.frame(smoot9x9)
base2 <- MtoT(base1)
base2[order(base2[,1], base2[,2]),]
PredMA[order(PredMA[,1], PredMA[,2]),]
Finalbase <- cbind(PredMA[,1:6],base2[,3])
names(Finalbase)[7]<-paste("Zone")



## 7. VALIDATION OF MANAGEMENT ZONES------------------------------------------------------

# Load data set
Sample <-read.table("C:\\.....\\Sample.txt", header = TRUE)
Sample$Zone<-as.factor(Sample$Zone)

# SOM data ---------------------------------
# Model with exponencial spatial correlation
mod1_som <-gls(SOM~1+Zone
,correlation=corExp(form=~as.numeric(as.character(X))+as.numeric(as.character(Y))
,metric="euclidean"
,nugget=FALSE)
,method="REML"
,na.action=na.omit
,data=Sample)

# Model with exponential spatial correlation and nugget effect
mod2_som <-gls(SOM~1+Zone
,correlation=corExp(form=~as.numeric(as.character(X))+as.numeric(as.character(Y))
,metric="euclidean"
,nugget=TRUE)
,method="REML"
,na.action=na.omit
,data=Sample)

# Model with spherical spatial correlation
mod3_som <-gls(SOM~1+Zone
,correlation=corSpher(form=~as.numeric(as.character(X))+as.numeric(as.character(Y))
,metric="euclidean"
,nugget=FALSE)
,method="REML"
,na.action=na.omit
,data=Sample)

# Model with spherical spatial correlation and nugget effect
mod4_som <-gls(SOM~1+Zone
,correlation=corSpher(form=~as.numeric(as.character(X))+as.numeric(as.character(Y))
,metric="euclidean"
,nugget=TRUE)
,method="REML"
,na.action=na.omit
,data=Sample)

# Model of independent errors
mod5_som <-gls(SOM~1+Zone
,method="REML"
,na.action=na.omit
,data=Sample)

# Selecting spatial correlation model using the Akaike information criterion
AICmod1_som <- AIC(mod1_som)
AICmod2_som <- AIC(mod2_som)
AICmod3_som <- AIC(mod3_som)
AICmod4_som <- AIC(mod4_som)
AICmod5_som <- AIC(mod5_som)

AICmod1_som
AICmod2_som
AICmod3_som
AICmod4_som
AICmod5_som

# Summary of selected model (SOM)
summary(mod1_som)
SOMmeans <- summary(lsmeans(mod1_som,"Zone")); SOMmeans


# Clay data --------------------------------
# Model with exponencial spatial correlation
mod1_clay <-gls(Clay~1+Zone
,correlation=corExp(form=~as.numeric(as.character(X))+as.numeric(as.character(Y))
,metric="euclidean"
,nugget=FALSE)
,method="REML"
,na.action=na.omit
,data=Sample)

# Model with exponential spatial correlation and nugget effect
mod2_clay <-gls(Clay~1+Zone
,correlation=corExp(form=~as.numeric(as.character(X))+as.numeric(as.character(Y))
,metric="euclidean"
,nugget=TRUE)
,method="REML"
,na.action=na.omit
,data=Sample)

# Model with spherical spatial correlation
mod3_clay <-gls(Clay~1+Zone
,correlation=corSpher(form=~as.numeric(as.character(X))+as.numeric(as.character(Y))
,metric="euclidean"
,nugget=FALSE)
,method="REML"
,na.action=na.omit
,data=Sample)

# Model with spherical spatial correlation and nugget effect
mod4_clay <-gls(Clay~1+Zone
,correlation=corSpher(form=~as.numeric(as.character(X))+as.numeric(as.character(Y))
,metric="euclidean"
,nugget=TRUE)
,method="REML"
,na.action=na.omit
,data=Sample)

# Model of independent errors
mod5_clay <-gls(Clay~1+Zone
,method="REML"
,na.action=na.omit
,data=Sample)

# Selecting spatial correlation model using the Akaike information criterion
AICmod1_clay <- AIC(mod1_clay)
AICmod2_clay <- AIC(mod2_clay)
AICmod3_clay <- AIC(mod3_clay)
AICmod4_clay <- AIC(mod4_clay)
AICmod5_clay <- AIC(mod5_clay)

AICmod1_clay
AICmod2_clay
AICmod3_clay
AICmod4_clay
AICmod5_clay

# Summary of selected model (Clay)
summary(mod3_clay)
Claymeans <- summary(lsmeans(mod3_clay,"Zone")); Claymeans


# Wheat yield data --------------------------------

# Load data set
WheatYield <-read.table("C:\\.....\\WheatYield.txt", header = TRUE)

# Making random sample (n = 1000) on the original wheat yiel data (n = 5982)
set.seed(20)
SampleWY <- WheatYield[sample(1:nrow(WheatYield), 1000, replace=FALSE),1:4]
SampleWY$Zone<-as.factor(SampleWY$Zone)

# Model with exponencial spatial correlation
mod1_Wy <-gls(Wy~1+Zone
,correlation=corExp(form=~as.numeric(as.character(X))+as.numeric(as.character(Y))
,metric="euclidean"
,nugget=FALSE)
,method="REML"
,na.action=na.omit
,data=SampleWY)

# Model with exponential spatial correlation and nugget effect
mod2_Wy <-gls(Wy~1+Zone
,correlation=corExp(form=~as.numeric(as.character(X))+as.numeric(as.character(Y))
,metric="euclidean"
,nugget=TRUE)
,method="REML"
,na.action=na.omit
,data=SampleWY)

# AModel with spherical spatial correlation
mod3_Wy <-gls(Wy~1+Zone
,correlation=corSpher(form=~as.numeric(as.character(X))+as.numeric(as.character(Y))
,metric="euclidean"
,nugget=FALSE)
,method="REML"
,na.action=na.omit
,data=SampleWY)

# Model with spherical spatial correlation and nugget effect
mod4_Wy <-gls(Wy~1+Zone
,correlation=corSpher(form=~as.numeric(as.character(X))+as.numeric(as.character(Y))
,metric="euclidean"
,nugget=TRUE)
,method="REML"
,na.action=na.omit
,data=SampleWY)

# Model of independent errors
mod5_Wy <-gls(Wy~1+Zone
,method="REML"
,na.action=na.omit
,data=SampleWY)

# Selecting spatial correlation model using the Akaike information criterion
AICmod1_Wy <- AIC(mod1_Wy)
AICmod2_Wy <- AIC(mod2_Wy)
AICmod3_Wy <- AIC(mod3_Wy)
AICmod4_Wy <- AIC(mod4_Wy)
AICmod5_Wy <- AIC(mod5_Wy)

AICmod1_Wy
AICmod2_Wy
AICmod3_Wy
AICmod4_Wy
AICmod5_Wy

# Summary of selected model (Yield)
summary(mod1_Wy)
Wymeans <- summary(lsmeans(mod1_Wy,"Zone")); Wymeans

# Mean differences of soil and yield variables for the delineated management zones
par(mfrow=c(1,3))
attach(SOMmeans)
SOMmean <-by(lsmean,Zone,mean)
so <- barplot(SOMmean,xlab="Management Zone", ylab="SOM (%)",col=c("grey94","grey34")
,ylim=c(3,max(SOMmean+SOMmean*0.1)),xpd=F)
letters = c("a","b")
text(x=so,y=SOMmean+SOMmean*0.05,label=letters,cex = 1)

attach(Claymeans)
Claymean <-by(lsmean,Zone,mean)
cl <-barplot(Claymean,xlab="Management Zone", ylab="Clay (%)",col=c("grey94","grey34")
,ylim=c(15,max(Claymean+Claymean*0.1)),xpd=F)
letters = c("a","a")
text(x=cl,y=Claymean+Claymean*0.05,label=letters,cex = 1)

attach(Wymeans)
Wymean <-by(lsmean,Zone,mean)
Wyse<-by(SE,Zone,mean)
wy <- barplot(Wymean,xlab="Management Zone", ylab="Wheat Yield (t/ha)",col=c("grey94","grey34")
,ylim=c(3,max(Wymean+Wymean*0.1)),xpd=F)
letters = c("a","b")
text(x=wy,y=Wymean+Wymean*0.05,label=letters,cex = 1)