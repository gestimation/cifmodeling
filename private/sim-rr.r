
#########################################################
library(mets)
library(prodlim)
data(bmt)

cif1 <- cbind(c(0,10,20,100),c(0,0.15,0.2,0.22))
#cif2 <- cbind(c(0,10,20,100),c(0,0.1,0.15,0.18))
#cif2 <- cbind(c(0,10,20,100),c(0,0.3,0.35,0.38))
#cif2 <- cbind(c(0,10,20,100),c(0,0.5,0.55,0.58))
cif2 <- cbind(c(0,10,20,100),c(0,0.7,0.75,0.78))

#plotl(cif1,ylim=c(0,1))
#lines(cif2,col=2)

fcif <- function(cif,rr,type=c("cif","cloglog","logistic"))
{# {{{
	mcif <- max(cif[,2])
	if (type[1]=="cif") mcif <- mcif*rr 
	if (type[1]=="cloglog") mcif <- 1- exp(-mcif*rr)
	if (type[1]=="logistic") mcif <- mcif* rr/(1 + mcif * rr)
	return(mcif)
}
# }}}

fcif(cif1,c(1,1.1,1.2))

#lrr1=c(-0.69314718056,-0.69314718056);lrr2=c(-0.69314718056,-0.69314718056);cens=NULL
#lrr1=c(-0.69314718056,0);lrr2=c(0,-0.69314718056);cens=NULL

#lrr1=c(-0.69314718056,-0.69314718056);lrr2=c(0,-0.69314718056);cens=NULL
lrr1=c(0,-0.69314718056);lrr2=c(0,-0.69314718056);cens=NULL
type1=c("cif","cloglog","logistic")
		  type2=c("cif","cloglog","logistic")
		  type1=type1[1]
		  type2=type2[1]
		  n <- 100

		  #simRR <- function(n,lrr1=c(-0.69314718056,-0.69314718056),lrr2=c(-0.69314718056,-0.69314718056),cens=NULL,
#simRR <- function(n,lrr1=c(0,-0.693147180560),lrr2=c(-0.69314718056,-0.69314718056),cens=NULL,
#simRR <- function(n,lrr1=c(-0.693147180560,-0.69314718056),lrr2=c(0,-0.69314718056),cens=NULL,
simRR <- function(n,lrr1=c(0,-0.69314718056),lrr2=c(0,-0.69314718056),cens=NULL,
		  type1=c("cif","cloglog","logistic"),type2=c("cif","cloglog","logistic")
		  ) {# {{{
### A binary, L binary
A <- rbinom(n,1,0.5)
L <- rbinom(n,1,0.5)
###
rr1 <- exp(cbind(A,L) %*% lrr1)
rr2 <- exp(cbind(A,L) %*% lrr2)
## model is fine
f1 <- fcif(cif1,max(rr1),type=type1[1])
f2 <- fcif(cif2,max(rr2),type=type2[1])
mmm<- f1+f2
mcif1 <- fcif(cif1,rr1,type=type1[1])
mcif2 <- fcif(cif2,rr2,type=type2[1])
if (mmm>1) warning(" models not satisfying sum <=1\n")
###
T1 <- simsubdist(cif1,rr1,type=type1[1])
T2 <- simsubdist(cif2,rr2,type=type2[1])
###
dies <- rbinom(n,1,mcif1+mcif2)
sel1 <- rbinom(n,1,mcif2/(mcif1+mcif2))+1
epsilon  <- dies*(sel1)
T1$epsilon <- epsilon
###
T1$A <- A
T1$L <- L
## times given 
T1$time <- T1$timecause
T1$time2 <- T2$timecause
T1$status <- epsilon
T1 <- dtransform(T1,time=100,epsilon==0)
T1 <- dtransform(T1,status=0,epsilon==0)
###
T1 <- dtransform(T1,time=time2,epsilon==2)
T1 <- dtransform(T1,status=2,epsilon==2)
data <- T1

if (!is.null(cens))  {
	cc <- rexp(n)/cens
	data$status <- ifelse(data$time<cc,data$status,0)
	data$time <- pmin(data$time,cc)
}

###
return(data)
}
# }}}

T1 <- simRR(2000000,cens=0.01)
#T1 <- simRR(2000000,cens=0.01,type1="logistic",type2="cif")
#T1 <- simRR(2000000,cens=0.01,type1="logistic",type2="logistic")
#T1 <- simRR(100000,cens=0.01,type1="cloglog",type2="logistic")

#####################################################
### check it out ####################################
#####################################################
dtable(T1,~status)

#lrr1=c(-0.69314718056,-0.69314718056);lrr2=c(-0.69314718056,-0.69314718056)
#lrr1=c(-0.69314718056,-0.69314718056);lrr2=c(0,-0.69314718056)
#lrr1=c(0,-0.69314718056);lrr2=c(-0.69314718056,-0.69314718056)
lrr1=c(0,-0.69314718056);lrr2=c(0,-0.69314718056)
cifs <- prodlim(Hist(time,status)~A+L,T1)
###
newd <- data.frame(expand.grid(A=0:1,L=0:1))
rr1 <- c(exp(as.matrix(newd) %*% lrr1))
rr2 <- c(exp(as.matrix(newd) %*% lrr2))
###
cifm1 <- cbind(cif1[,1],cif1[,2] %o% rr1)
cifm2 <- cbind(cif2[,1],cif2[,2] %o% rr2)
###
#par(mfrow=c(1,2))
#plot(cifs,cause=1,ylim=c(0,0.3)); 
#matlines(cifm1[,1],cifm1[,-1],col=1,lwd=2)
###
#plot(cifs,cause=2,ylim=c(0,0.7))
#matlines(cifm2[,1],cifm2[,-1],col=1,lwd=2)

#write.csv(T1, "test.csv", row.names = FALSE)
write.csv(T1, "RR_10051005_2075.csv", row.names = FALSE)
getwd()


