source("http://zzlab.net/GAPIT/gapit_functions.txt") 
library(dplyr)
library(tidyverse)
library(BGLR)
setwd("/home/peter.schmuker/LEEROY/")
#kernel code from A guide for kernel generalized regression methods for genomic-enabled prediction
Kernel_computation=function(X,name, degree, nL){
  
  p=ncol(X)
  
  x1=X
  
  x2=X
  
  d=degree
  
  ############Polynomial kernel##################
  
  K.Polynomial=function(x1, x2=x1, gamma=1, b=0, d=3)
    
  { (gamma*(as.matrix(x1)%*%t(x2))+b)^d}
  
  ############Sigmoid kernel####################
  
  K.Sigmoid=function(x1,x2=x1, gamma=1, b=0)
    
  { tanh(gamma*(as.matrix(x1)%*%t(x2))+b) }
  
  ############Gaussian kernel##################
  
  l2norm=function(x){sqrt(sum(x^2))}
  
  K.Gaussian=function(x1,x2=x1, gamma=1){
    
    exp(-gamma*outer(1:nrow(x1<- as.matrix(x1)), 1:ncol(x2<- t(x2)),
                     
                     Vectorize(function(i, j) l2norm(x1[i,]-x2[,j])^2)))}
  
  ##########Arc-cosine kernel with 1 hidden layer
  
  K.AK1_Final<-function(x1,x2){
    
    n1<-nrow(x1)
    
    n2<-nrow(x2)
    
    x1tx2<-x1%*%t(x2)
    
    norm1<-sqrt(apply(x1,1,function(x) crossprod(x)))
    
    norm2<-sqrt(apply(x2,1,function(x) crossprod(x)))
    
    costheta = diag(1/norm1)%*%x1tx2%*%diag(1/norm2)
    
    costheta[which(abs(costheta)>1,arr.ind = TRUE)] = 1
    
    theta<-acos(costheta)
    
    normx1x2<-norm1%*%t(norm2)
    
    J = (sin(theta)+(pi-theta)*cos(theta))
    
    AK1 = 1/pi*normx1x2x2*J
    
    AK1<-AK1/median(AK1)
    
    colnames(AK1)<-rownames(x2)
    
    rownames(AK1)<-rownames(x1)
    
    return(AK1)
    
  }
  
  ####Kernel Arc-Cosine with deep=L#########
  
  AK_L_Final<-function(AK1,nL){
    
    n1<-nrow(AK1)
    
    n2<-ncol(AK1)
    
    AKl1 = AK1
    
    for (l in 1:nL){
      
      AKAK<-tcrossprod(diag(AKl1),diag(AKl1))
      
      costheta<-AKl1*(AKAK^(-1/2))
      
      costheta[which(costheta>1,arr.ind = TRUE)] = 1
      
      theta<-acos(costheta)
      
      AKl<-(1/pi)*(AKAK^(1/2))*(sin(theta)+(pi-theta)*cos(theta))
      
      AKl1 = AKl
      
    }
    
    AKl<-AKl/median(AKl)
    
    rownames(AKl)<-rownames(AK1)
    
    colnames(AKl)<-colnames(AK1)
    
    return(AKl)
    
  }
  
  ########Exponencial Kernel############
  
  K.exponential=function(x1,x2=x1, gamma=1){
    
    exp(-gamma*outer(1:nrow(x1<- as.matrix(x1)), 1:ncol(x2<- t(x2)),
                     
                     Vectorize(function(i, j) l2norm(x1[i,]-x2[,j]))))}
  
  if (name=='Linear') {
    
    K=X%*%t(X)/p
    
  } else if (name=='Polynomial') {
    
    K=K.Polynomial(x1=x1, x2=x1, gamma=1/p, b=0, d=d)
    
  } else if (name=='Sigmoid') {
    
    K=K.Sigmoid(x1=x1, x2=x1, gamma=1/p, b=0)
    
  }else if (name=='Gaussian') {
    
    K=K.Gaussian(x1=x1, x2=x1, gamma=1/p)
    
  } else if (name=='AK1') {
    
    K= K.AK1_Final(x1=x1, x2=x1)
    
  } else if (name=='AKL') {
    
    AK1=K.AK1_Final(x1=x1, x2=x1)
    
    K=AK_L_Final(AK1=AK1,nL=nL)
    
  } else {
    
    K=K.exponential(x1=x1,x2=x1,gamma=1/p)
    
  }
  
}
trait22<-read.csv("trait22.csv",header=T)
map22<-read.csv("map22.csv",header=T)
gene22<-read.csv("gene22.csv",header=T)
kin22<-read.csv("kinship22.csv",header=F)
colnames(gene22)[1] = "Taxa"
print(gene22)
gene22 <- gene22[match(trait22$Taxa, gene22$Taxa), ]    
print(gene22)
gene22n<-gene22
rownames(gene22n)<-NULL

fixpls<-column_to_rownames(gene22n, var='Taxa')
sigkernel=Kernel_computation(X=fixpls,name='Sigmoid', degree=NULL, nL=NULL)
polykernel=Kernel_computation(X=fixpls,name='Polynomial', degree=2, nL=NULL)

fixpls<-fixpls/2
D<-as.matrix(dist(fixpls,method="euclidean"))^2
D<-D/mean(D)
   
K01<-exp(-.01*D)
K1<-exp(-1*D)
K3<-exp(-2*D)
df2 = data.frame(model = "starter", cor=0, rmse = 1)
trait22$Yield<-scale(trait22$Yield)

for(i in 1:50) {
n=nrow(gene22)
print(n)
testing=sample(n,round(n/5),replace=F)
training=-testing

y<-trait22[,6]

yNA<-y
yNA[testing]<-NA

training=-testing
myGAPIT=GAPIT(
  Y=trait22[training,c(1,6)],  SNP.MAF =.1,
  GD=gene22,
  GM=map22,
  PCA.total=3,
  model=c("BLINK"),file.output = F,
  Multiple_analysis=F)


GPCA<-myGAPIT$PCA
gwas<-myGAPIT$GWAS
low5<-gwas %>% slice_min(P.value, n = 5)
low5<-low5$SNP
new_df5 = subset(gene22, select = low5)
new_df5$Taxa<-gene22$Taxa

colnames(GPCA)[1] = "Taxa"
print(GPCA)

print(new_df5)
myQTN5<-inner_join(GPCA, new_df5, by="Taxa")
low3<-gwas %>% slice_min(P.value, n = 3)
low3<-low3$SNP
new_df3 = subset(gene22, select = low3)
new_df3$Taxa<-gene22$Taxa
#new_df3$Taxa<-gene22$Taxa
myQTN3<-inner_join(GPCA, new_df3, by="Taxa")
low1<-gwas %>% slice_min(P.value, n = 1)
low1<-low1$SNP
new_df1 = subset(gene22, select = low1)
new_df1$Taxa<-gene22$Taxa
myQTN1<-inner_join(GPCA, new_df1, by="Taxa")

myGAPIT5 <-GAPIT(
  Y=trait22[training,c(1,6)],
  CV=myQTN5, 
  model="GLM", 
  SNP.test=FALSE,file.output = F,
  memo="MAS")
myPred2=myGAPIT5$Pred

ff<-inner_join(myPred2,trait22, by="Taxa")
ff <- ff[match(trait22$Taxa, ff$Taxa), ] 

df2 <- df2 %>% 
  add_row(model = 'QTN5', cor=cor(ff$Prediction[testing],ff$Yield[testing]), rmse =  (mean(abs(ff$Prediction[testing]-ff$Yield[testing]))))
myGAPIT3 <-GAPIT(
  Y=trait22[training,c(1,6)],
  CV=myQTN3, 
  model="GLM", 
  SNP.test=FALSE,file.output = F,
  memo="MAS")
myPred2=myGAPIT3$Pred

ff<-inner_join(myPred2,trait22, by="Taxa")
ff <- ff[match(trait22$Taxa, ff$Taxa), ] 
df2 <- df2 %>% 
  add_row(model = 'QTN3', cor=cor(ff$Prediction[testing],ff$Yield[testing]), rmse =  (mean(abs(ff$Prediction[testing]-ff$Yield[testing]))))

myGAPIT1 <-GAPIT(
  Y=trait22[training,c(1,6)],
  CV=myQTN1, 
  model="GLM", 
  SNP.test=FALSE,file.output = F,
  memo="MAS")
myPred2=myGAPIT1$Pred

ff<-inner_join(myPred2,trait22, by="Taxa")
ff <- ff[match(trait22$Taxa, ff$Taxa), ] 
df2 <- df2 %>% 
  add_row(model = 'QTN1', cor=cor(ff$Prediction[testing],ff$Yield[testing]), rmse =  (mean(abs(ff$Prediction[testing]-ff$Yield[testing]))))
myGAPIT=GAPIT(
  Y=trait22[training,c(1,6)], 
  GD=gene22,
  GM=map22,
  PCA.total=3, 
  model=c("gBLUP"),file.output = F,
  Multiple_analysis=F)

myPred2=myGAPIT$Pred

ff<-inner_join(myPred2,trait22, by="Taxa")
ff <- ff[match(trait22$Taxa, ff$Taxa), ] 
df2 <- df2 %>% 
  add_row(model = 'gBLUP', cor=cor(ff$Prediction[testing],ff$Yield[testing]), rmse =  (mean(abs(ff$Prediction[testing]-ff$Yield[testing]))))

myGAPITgblup1=GAPIT(
  Y=trait22[training,c(1,6)], 
  GD=gene22,
  GM=map22,
  CV=myQTN1, 
  model=c("gBLUP"),file.output = F,
  Multiple_analysis=F)

myPred2=myGAPITgblup1$Pred

ff<-inner_join(myPred2,trait22, by="Taxa")
ff <- ff[match(trait22$Taxa, ff$Taxa), ] 
df2 <- df2 %>% 
  add_row(model = 'gBLUP qtn1', cor=cor(ff$Prediction[testing],ff$Yield[testing]), rmse =  (mean(abs(ff$Prediction[testing]-ff$Yield[testing]))))
myGAPITgblup3=GAPIT(
  Y=trait22[training,c(1,6)], 
  GD=gene22,
  GM=map22,
  CV=myQTN3, 
  model=c("gBLUP"),file.output = F,
  Multiple_analysis=F)
myPred2=myGAPITgblup3$Pred

ff<-inner_join(myPred2,trait22, by="Taxa")
ff <- ff[match(trait22$Taxa, ff$Taxa), ] 
df2 <- df2 %>% 
  add_row(model = 'gBLUP qtn3', cor=cor(ff$Prediction[testing],ff$Yield[testing]), rmse =  (mean(abs(ff$Prediction[testing]-ff$Yield[testing]))))

myGAPITgblup5=GAPIT(
  Y=trait22[training,c(1,6)], 
  GD=gene22,
  GM=map22,
  CV=myQTN5, 
  model=c("gBLUP"),file.output = F,
  Multiple_analysis=F)
myPred2=myGAPITgblup5$Pred

ff<-inner_join(myPred2,trait22, by="Taxa")
ff <- ff[match(trait22$Taxa, ff$Taxa), ] 
df2 <- df2 %>% 
  add_row(model = 'gBLUP qtn5', cor=cor(ff$Prediction[testing],ff$Yield[testing]), rmse =  (mean(abs(ff$Prediction[testing]-ff$Yield[testing]))))

myGAPIT=GAPIT(
  Y=trait22[training,c(1,6)], 
  GD=gene22,
  GM=map22,
  PCA.total=3, 
  model=c("cBLUP"),file.output = F,
  Multiple_analysis=F)
myPred2=myGAPIT$Pred

ff<-inner_join(myPred2,trait22, by="Taxa")
ff <- ff[match(trait22$Taxa, ff$Taxa), ] 
df2 <- df2 %>% 
  add_row(model = 'cblup', cor=cor(ff$Prediction[testing],ff$Yield[testing]), rmse =  (mean(abs(ff$Prediction[testing]-ff$Yield[testing]))))
myGAPITcBLUP1=GAPIT(
  Y=trait22[training,c(1,6)], 
  GD=gene22,
  GM=map22,
  CV=myQTN1,
  model=c("cBLUP"),file.output = F,
  Multiple_analysis=F)
myPred2=myGAPITcBLUP1$Pred

ff<-inner_join(myPred2,trait22, by="Taxa")
ff <- ff[match(trait22$Taxa, ff$Taxa), ] 
df2 <- df2 %>% 
  add_row(model = 'cblup qtn1', cor=cor(ff$Prediction[testing],ff$Yield[testing]), rmse =  (mean(abs(ff$Prediction[testing]-ff$Yield[testing]))))
myGAPITcBLUP3=GAPIT(
  Y=trait22[training,c(1,6)], 
  GD=gene22,
  GM=map22,
  CV=myQTN3,
  model=c("cBLUP"),file.output = F,
  Multiple_analysis=F)
myPred2=myGAPITcBLUP3$Pred
ff<-inner_join(myPred2,trait22, by="Taxa")
ff <- ff[match(trait22$Taxa, ff$Taxa), ] 
df2 <- df2 %>% 
  add_row(model = 'cblup qtn3', cor=cor(ff$Prediction[testing],ff$Yield[testing]), rmse =  (mean(abs(ff$Prediction[testing]-ff$Yield[testing]))))

myGAPITcBLUP5=GAPIT(
  Y=trait22[training,c(1,6)], 
  GD=gene22,
  GM=map22,
  CV=myQTN5, 
  model=c("cBLUP"),file.output = F,
  Multiple_analysis=F)
myPred2=myGAPITcBLUP5$Pred

ff<-inner_join(myPred2,trait22, by="Taxa")
ff <- ff[match(trait22$Taxa, ff$Taxa), ] 
df2 <- df2 %>% 
  add_row(model = 'cblup qtn5', cor=cor(ff$Prediction[testing],ff$Yield[testing]), rmse =  (mean(abs(ff$Prediction[testing]-ff$Yield[testing]))))

ETA1<-list(list(K=K01,model='RKHS'))
fm01<-BGLR(y=yNA,ETA=ETA1,
           nIter=10000,burnIn=5000,df0=5,S0=2)
df2 <- df2 %>% 
  add_row(model = 'bandwidth .01', cor=cor(fm01$yHat[testing],y[testing]), rmse =  (mean(abs(fm01$yHat[testing]-y[testing]))))
ETA2<-list(list(K=K1,model='RKHS'))
fm01<-BGLR(y=yNA,ETA=ETA2,
          nIter=10000,burnIn=5000,df0=5,S0=2)
df2 <- df2 %>% 
  add_row(model = 'bandwidth 1', cor=cor(fm01$yHat[testing],y[testing]), rmse =  (mean(abs(fm01$yHat[testing]-y[testing]))))
ETA2<-list(list(K=K3,model='RKHS'))
fm01<-BGLR(y=yNA,ETA=ETA2,
          nIter=10000,burnIn=5000,df0=5,S0=2)
df2 <- df2 %>% 
  add_row(model = 'bandwidth 1.5', cor=cor(fm01$yHat[testing],y[testing]), rmse =  (mean(abs(fm01$yHat[testing]-y[testing]))))

ETA3<-list(list(K=sigkernel,model='RKHS'))
fm01<-BGLR(y=yNA,ETA=ETA3,
          nIter=10000,burnIn=5000,df0=5,S0=2)
df2 <- df2 %>% 
  add_row(model = 'sigmoid', cor=cor(fm01$yHat[testing],y[testing]), rmse =  (mean(abs(fm01$yHat[testing]-y[testing]))))

ETA3<-list(list(K=polykernel,model='RKHS'))
fm01<-BGLR(y=yNA,ETA=ETA3,
          nIter=10000,burnIn=5000,df0=5,S0=2)
df2 <- df2 %>% 
  add_row(model = 'polynomial', cor=cor(fm01$yHat[testing],y[testing]), rmse =  (mean(abs(fm01$yHat[testing]-y[testing]))))

n1 <- names(myQTN1)
n1<-n1[-1]
f1 <- as.formula(paste("~ ", paste(n1[!n1 %in% new_df5], collapse = " + ")))
n3 <- names(myQTN3)
n3<-n3[-1]
f3 <- as.formula(paste("~ ", paste(n3[!n3 %in% new_df5], collapse = " + ")))
n5 <- names(myQTN5)
n5<-n5[-1]
f5 <- as.formula(paste("~ ", paste(n5[!n5 %in% new_df5], collapse = " + ")))

ETA1<-list(fixed=list(f1,data=myQTN1,model='FIXED'),list(K=K01,model='RKHS'))
ETA3<-list(fixed=list(f3,data=myQTN3,model='FIXED'),list(K=K01,model='RKHS'))
ETA5<-list(fixed=list(f5,data=myQTN5,model='FIXED'),list(K=K01,model='RKHS'))

fm01<-BGLR(y=yNA,ETA=ETA1,
           nIter=10000,burnIn=5000,df0=5,S0=2)
df2 <- df2 %>% 
  add_row(model = 'bandwidth 01 qtn1', cor=cor(fm01$yHat[testing],y[testing]), rmse =  (mean(abs(fm01$yHat[testing]-y[testing]))))

fm01<-BGLR(y=yNA,ETA=ETA3,
           nIter=10000,burnIn=5000,df0=5,S0=2)
df2 <- df2 %>% 
  add_row(model = 'bandwidth 01 qtn3', cor=cor(fm01$yHat[testing],y[testing]), rmse =  (mean(abs(fm01$yHat[testing]-y[testing]))))

fm01<-BGLR(y=yNA,ETA=ETA5,
           nIter=10000,burnIn=5000,df0=5,S0=2)
df2 <- df2 %>% 
  add_row(model = 'bandwidth 01 qtn5', cor=cor(fm01$yHat[testing],y[testing]), rmse =  (mean(abs(fm01$yHat[testing]-y[testing]))))


ETA1<-list(fixed=list(f1,data=myQTN1,model='FIXED'),list(K=K1,model='RKHS'))
ETA3<-list(fixed=list(f3,data=myQTN3,model='FIXED'),list(K=K1,model='RKHS'))
ETA5<-list(fixed=list(f5,data=myQTN5,model='FIXED'),list(K=K1,model='RKHS'))

fm01<-BGLR(y=yNA,ETA=ETA1,
           nIter=10000,burnIn=5000,df0=5,S0=2)
df2 <- df2 %>% 
  add_row(model = 'bandwidth 1 qtn1', cor=cor(fm01$yHat[testing],y[testing]), rmse =  (mean(abs(fm01$yHat[testing]-y[testing]))))

fm01<-BGLR(y=yNA,ETA=ETA3,
           nIter=10000,burnIn=5000,df0=5,S0=2)
df2 <- df2 %>% 
  add_row(model = 'bandwidth 1 qtn3', cor=cor(fm01$yHat[testing],y[testing]), rmse =  (mean(abs(fm01$yHat[testing]-y[testing]))))

fm01<-BGLR(y=yNA,ETA=ETA5,
           nIter=10000,burnIn=5000,df0=5,S0=2)
df2 <- df2 %>% 
  add_row(model = 'bandwidth 1 qtn5', cor=cor(fm01$yHat[testing],y[testing]), rmse =  (mean(abs(fm01$yHat[testing]-y[testing]))))


ETA1<-list(fixed=list(f1,data=myQTN1,model='FIXED'),list(K=K3,model='RKHS'))
ETA3<-list(fixed=list(f3,data=myQTN3,model='FIXED'),list(K=K3,model='RKHS'))
ETA5<-list(fixed=list(f5,data=myQTN5,model='FIXED'),list(K=K3,model='RKHS'))

fm01<-BGLR(y=yNA,ETA=ETA1,
           nIter=10000,burnIn=5000,df0=5,S0=2)
df2 <- df2 %>% 
  add_row(model = 'bandwidth 1.5 qtn1', cor=cor(fm01$yHat[testing],y[testing]), rmse =  (mean(abs(fm01$yHat[testing]-y[testing]))))

fm01<-BGLR(y=yNA,ETA=ETA3,
           nIter=10000,burnIn=5000,df0=5,S0=2)
df2 <- df2 %>% 
  add_row(model = 'bandwidth 1.5 qtn3', cor=cor(fm01$yHat[testing],y[testing]), rmse =  (mean(abs(fm01$yHat[testing]-y[testing]))))

fm01<-BGLR(y=yNA,ETA=ETA5,
           nIter=10000,burnIn=5000,df0=5,S0=2)
df2 <- df2 %>% 
  add_row(model = 'bandwidth 1.5 qtn5', cor=cor(fm01$yHat[testing],y[testing]), rmse =  (mean(abs(fm01$yHat[testing]-y[testing]))))


ETA1<-list(fixed=list(f1,data=myQTN1,model='FIXED'),list(K=sigkernel,model='RKHS'))
ETA3<-list(fixed=list(f3,data=myQTN3,model='FIXED'),list(K=sigkernel,model='RKHS'))
ETA5<-list(fixed=list(f5,data=myQTN5,model='FIXED'),list(K=sigkernel,model='RKHS'))

fm01<-BGLR(y=yNA,ETA=ETA1,
           nIter=10000,burnIn=5000,df0=5,S0=2)
df2 <- df2 %>% 
  add_row(model = 'sigmoid qtn1', cor=cor(fm01$yHat[testing],y[testing]), rmse =  (mean(abs(fm01$yHat[testing]-y[testing]))))

fm01<-BGLR(y=yNA,ETA=ETA3,
           nIter=10000,burnIn=5000,df0=5,S0=2)
df2 <- df2 %>% 
  add_row(model = 'sigmoid qtn3', cor=cor(fm01$yHat[testing],y[testing]), rmse =  (mean(abs(fm01$yHat[testing]-y[testing]))))

fm01<-BGLR(y=yNA,ETA=ETA5,
           nIter=10000,burnIn=5000,df0=5,S0=2)
df2 <- df2 %>% 
  add_row(model = 'sigmoid qtn5', cor=cor(fm01$yHat[testing],y[testing]), rmse =  (mean(abs(fm01$yHat[testing]-y[testing]))))


ETA1<-list(fixed=list(f1,data=myQTN1,model='FIXED'),list(K=polykernel,model='RKHS'))
ETA3<-list(fixed=list(f3,data=myQTN3,model='FIXED'),list(K=polykernel,model='RKHS'))
ETA5<-list(fixed=list(f5,data=myQTN5,model='FIXED'),list(K=polykernel,model='RKHS'))

fm01<-BGLR(y=yNA,ETA=ETA1,
           nIter=10000,burnIn=5000,df0=5,S0=2)
df2 <- df2 %>% 
  add_row(model = 'polynomial qtn1', cor=cor(fm01$yHat[testing],y[testing]), rmse =  (mean(abs(fm01$yHat[testing]-y[testing]))))

fm01<-BGLR(y=yNA,ETA=ETA3,
           nIter=10000,burnIn=5000,df0=5,S0=2)
df2 <- df2 %>% 
  add_row(model = 'polynomial qtn3', cor=cor(fm01$yHat[testing],y[testing]), rmse =  (mean(abs(fm01$yHat[testing]-y[testing]))))

fm01<-BGLR(y=yNA,ETA=ETA5,
           nIter=10000,burnIn=5000,df0=5,S0=2)
df2 <- df2 %>% 
  add_row(model = 'polynomial qtn5', cor=cor(fm01$yHat[testing],y[testing]), rmse =  (mean(abs(fm01$yHat[testing]-y[testing]))))}




write.csv(df2, 'Yieldloopsecond25.csv')
