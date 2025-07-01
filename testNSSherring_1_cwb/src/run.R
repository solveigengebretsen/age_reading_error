
##Install stockassessment library with age-reading error correction
##devtools::install_github("fishfollower/SAM/stockassessment", ref = "ageconfusion")

library(stockassessment)

## Function not included in current ageconfusion branch of stockassessment
caytable <- function(fit, fleet=which(fit$data$fleetTypes==0)){
   getfleet <- function(f){
     idx <- fit$conf$keyLogFsta[f,]+2    
     F <- cbind(NA,exp(t(fit$pl$logF)))[,idx]
     F[is.na(F)] <- 0
     M <- fit$data$natMor
     N <- exp(t(fit$pl$logN))
     F/(F+M)*N*(1-exp(-F-M))
   }
   ret <- Reduce("+",lapply(fleet,getfleet)) 
   colnames(ret) <- fit$conf$minAge:fit$conf$maxAge
   rownames(ret) <- fit$data$years
   return(ret)
}



cn<-read.ices("data/cn.dat")
cw<-read.ices("data/cw.dat")
dw<-read.ices("data/dw.dat")
lw<-read.ices("data/lw.dat")
mo<-read.ices("data/mo.dat")
nm<-read.ices("data/nm.dat")
pf<-read.ices("data/pf.dat")
pm<-read.ices("data/pm.dat")
sw<-read.ices("data/sw.dat")
lf<-read.ices("data/lf.dat")
surveys<-read.ices("data/survey.dat")


cnWithAgeError<-cn
tmp<-as.matrix(read.table("data/NSS_skjell_bestemodell_min2.csv",sep=",",header=T,row.names=1))
                                    
attr(cnWithAgeError,"ageConfusion")<-tmp


## Set up and estimation without age-readings errors
datWithoutAgeError<-setup.sam.data(surveys=surveys,
                    residual.fleet=cnWithoutAgeError, 
                    prop.mature=mo, 
                    stock.mean.weight=sw, 
                    catch.mean.weight=cw, 
                    dis.mean.weight=dw, 
                    land.mean.weight=lw,
                    prop.f=pf, 
                    prop.m=pm, 
                    natural.mortality=nm, 
                    land.frac=lf)

confWithoutAgeError<-loadConf(datWithoutAgeError,"conf/model.cfg", patch=TRUE)
confWithoutAgeError$fixVarToWeight<-0 ## Necesarry correction for current ageconfusion branch of stockassessment
confWithoutAgeError$fracMixN<-0       ## Necesarry correction for current ageconfusion branch of stockassessment

parWithoutAgeError<-defpar(datWithoutAgeError,confWithoutAgeError)
fitWithoutAgeError<-sam.fit(datWithoutAgeError,confWithoutAgeError,parWithoutAgeError)
if(fitWithoutAgeError$opt$convergence!=0) stop("Model did not converge.")



## Set up and estimation with age-reading errors on catch data
datWithAgeError<-setup.sam.data(surveys=surveys,
                    residual.fleet=cnWithAgeError, 
                    prop.mature=mo, 
                    stock.mean.weight=sw, 
                    catch.mean.weight=cw, 
                    dis.mean.weight=dw, 
                    land.mean.weight=lw,
                    prop.f=pf, 
                    prop.m=pm, 
                    natural.mortality=nm, 
                    land.frac=lf)

confWithAgeError<-loadConf(datWithAgeError,"conf/model.cfg", patch=TRUE)
confWithAgeError$fixVarToWeight<-0 ## Necesarry correction for current ageconfusion branch of stockassessment
confWithAgeError$fracMixN<-0       ## Necesarry correction for current ageconfusion branch of stockassessment

parWithAgeError<-defpar(datWithAgeError,confWithAgeError)
fitWithAgeError<-sam.fit(datWithAgeError,confWithAgeError,parWithAgeError)
if(fitWithAgeError$opt$convergence!=0) stop("Model did not converge.")




caaWithoutAgeError<-caytable(fitWithoutAgeError)
caaWithAgeError<-caytable(fitWithAgeError)

pdf("caa4years.pdf")
par(mfrow=c(2,2))
for (year in c("2016","2017","2018","2019")) {
  plot(2:12,cn[year,],type="n",xlab="age",ylab="catch (millions)",main=year,ylim=c(0,600))
  points(2:12,cn[year,],col="black")
  lines(2:12,caaWithoutAgeError[year,],lwd=2,col="black")
  lines(2:12,caaWithAgeError[year,],lwd=2,col="red")
  if (year=="2016") {
    legend("topleft",
           legend=c("data","SAM fit without age-reading error",
                    "SAM fit with age-reading error"),

           col=c("black","black","red"),lwd=c(1,2,2),lty=c(NA,1,1),pch=c(1,NA,NA))
  }
}
dev.off()

pdf("SAMwithAgeErrors.pdf",width=8,height=4)
par(mfrow=c(1,2),mar=c(4,4,3,1))
year<-"2018"
plot(2:12,cn[year,],type="n",xlab="age",ylab="",main=year,ylim=c(0,600))
title(ylab="catch (millions)",line=2)
  points(2:12,cn[year,],col="black")
  lines(2:12,caaWithoutAgeError[year,],lwd=2,col="black")
  lines(2:12,caaWithAgeError[year,],lwd=2,col="red")
    legend("topleft",
           legend=c("data","SAM fit without age-reading error",
                    "SAM fit with age-reading error"),
           col=c("black","black","red"),lwd=c(1,1.9,1.9),lty=c(NA,1,1),pch=c(1,NA,NA),
           cex=c(0.85,0.85,0.85))

ssbplot(fitWithoutAgeError, addCI=TRUE)   ### SSB in thousand tonnes
lines(1988:2024,ssbtable(fitWithAgeError)[,"Estimate"],col="red",lwd=2)
dev.off()
