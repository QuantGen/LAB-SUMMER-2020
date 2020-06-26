##########################################################################################
#Example 1

#Loading software
library(BGLR)

#Loading data
data(wheat)
y<-wheat.Y
X<-wheat.X 
X<-scale(X)/sqrt(ncol(X))

#Linear predictor 
ETA<-list(list(X=X,model="BRR"))

#Model fitting
set.seed(123) 
fm<-Multitrait(y=y,ETA=ETA,nIter=10000,burnIn=5000)


#Note equivalent code to specify structures
#ETA<-list(list(X=X,model="BRR",Cov=list(type="UN")))
#residual<-list(type="UN")

#Model fitting
#set.seed(123) 
#fm<-Multitrait(y=y,ETA=ETA,resCov=residual,nIter=10000,burnIn=5000)

#Residual covariance matrix
fm$resCov

#Genetic covariance matrix
fm$ETA[[1]]$Cov

#Marker effects
fm$ETA[[1]]$beta

unlink("*.dat")

#Some plots

par(mfrow=c(2,2))

for(j in 1:4)
{
	plot(y[,j],fm$yHat[,j],xlab="y",ylab=expression(hat(y)),main=paste("Environment ", j))
}

par(mfrow=c(2,2))

for(j in 1:4)
{
	hist(y[,j]-fm$yHat[,j],xlab="Residuals",main=paste("Environment ", j))
}

#Plots of posterior distributions

R<-read.table(file="R.dat",header=FALSE)
element<-1
par(mfrow=c(1,2))
plot(R[,element],type="b",ylab=expression(R[11]),xlab="Thinned iteration")
plot(density(R[,element]),xlab=expression(R[11]),
     ylab=expression(p(group("",R[11],"|")*data)),
     main="")

##########################################################################################
#Example 2

library(BGLR)
data(wheat)
K<-wheat.A
y<-wheat.Y

ETA<-list(list(K=K,model="RKHS"))

#Fit model
set.seed(123)
fm<-Multitrait(y=y,ETA=ETA,nIter=10000,burnIn=5000)

#Residual covariance matrix
fm$resCov

#Genetic covariance matrix
fm$ETA[[1]]$Cov

#Random effects
fm$ETA[[1]]$u


##########################################################################################
#Missing values
#Example 3

rm(list=ls())

library(BGLR)
data(wheat)

#Compute genomic relationship matrix
M<-scale(wheat.X,center=TRUE)
K1<-tcrossprod(M)/ncol(M)
#Relationship matrix derived from pedigree
K2<-wheat.A

y<-wheat.Y

#Define linear predictor
ETA1<-list(mar=list(K=K1,model="RKHS"),
           ped=list(K=K2,model="RKHS"))

#Sets for training and testing
set.seed(1)
n<-nrow(y)*ncol(y)
whichNa<-sample(1:n,size=as.integer(0.20*n),replace=FALSE)

#Copy of y
yNa<-y

#Assign missing values
yNa[whichNa]<-NA

fm1<-Multitrait(y=yNa,ETA=ETA1,nIter=10000,burnIn=5000)
unlink("*.dat")

#Missing values in environment 4
whichNa4<-fm1$missing_records[fm1$patterns[,4]]

#Predictions in training and testing
par(mar=c(5,5,2,2))
plot(yNa[,4],fm1$yHat[,4],xlab=expression(y[4]),
     ylab=expression(hat(y)[4]))
points(y[whichNa4,4],fm1$yHat[whichNa4,4],col="red",pch=19)
rtrn<-cor(yNa[,4],fm1$yHat[,4],use="complete.obs")
rtst<-cor(y[whichNa4,4],fm1$yHat[whichNa4,4])
txt1<-paste("Training, r=",round(rtrn,2))
txt2<-paste("Testing,  r=",round(rtst,2))
legend("bottomright",legend=c(txt1,txt2),pch=c(1,19),
       col=c("black","red"),bty="n")

##########################################################################################
#Example 4

library(BGLR)
data(simulated3t)

y<-as.matrix(simulated3t.pheno[,1:3])
g<-as.matrix(simulated3t.pheno[,4:6])
cov(g)
y<-scale(y,center=TRUE,scale=FALSE)
y.orig<-y
	
X<-simulated3t.X
X<-scale(X)/sqrt(ncol(X))

ETA1<-list(list(X=X,model="SpikeSlab",
		        inclusionProb=list(probIn=rep(1/100,ncol(y)),
		        counts=rep(1E6,ncol(y)))))

#Fit the model

fm1<-Multitrait(y=y,ETA=ETA1,nIter=1000,burnIn=500)

#Residual covariance, UN
fm1$resCov

#Compare against the TRUE residual covariance matrix
#6.0 6.0 1.0
#6.0 8.0 2.0  
#1.0 2.0 1.0

#Genetic co-variance
crossprod(fm1$ETA[[1]]$beta)/fm1$ETA[[1]]$p

#Compare against the TRUE genetic co-variance matrix
#1.00   0.34   0.07
#0.34   1.00   0.21 
#0.07   0.21   1.00 
	
#Covariance matrix for b, UN
fm1$ETA[[1]]$Cov
