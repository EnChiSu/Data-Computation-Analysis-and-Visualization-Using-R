par(mfrow=c(3,2))
#正確版
fx <- function(x){
  return((2/pi)^0.5*exp(-(x^2)/2))
}

b <- function(){
  k <- rbinom(1,1,0.5)
  if(k==0) return(-1)
  else return(1)
}

set.seed(10)
n <- 10000
g <- rep(0,n)
count <- 0

for(i in 1:n){
  x <- rexp(1)
  y <- runif(1,0,(2*exp(1)/pi)^0.5*exp(-x))
  if (y<=fx(x)){
    g[i]<- x*b()
    count <- count+1
  }else{
    g[i] <- NA
  }
}

k<-g[!is.na(g)]

x=seq(-5,5,by=0.1)
hist(k,breaks=seq(-5,5,by=1),freq=F,ylim=c(0,1),main="Standard Normal Distribution\n Bin width=1")
lines(x,dnorm(x),col="red")
hist(k,breaks=seq(-5,5,by=0.5),freq=F,ylim=c(0,1),main="Standard Normal Distribution\n Bin width=0.5")
lines(x,dnorm(x),col="red")
hist(k,breaks="FD",freq=F,ylim=c(0,1),main="Standard Normal Distribution\n Bin width=h=2¡ÑIQR¡Ñn^(-1/3)=0.12")  #FD=(quantile(k,0.75)-quantil(k,0.25))*10000^(-1/3)*2
lines(x,dnorm(x),col="red")
hist(k,breaks=seq(-5,5,by=0.01),freq=F,ylim=c(0,1),main="Standard Normal Distribution\n Bin width=0.01")
lines(x,dnorm(x),col="red")
hist(k,breaks=seq(-5,5,by=0.001),freq=F,ylim=c(0,1),main="Standard Normal Distribution\n Bin width=0.001")
lines(x,dnorm(x),col="red")
hist(k,breaks=seq(-5,5,by=0.0001),freq=F,ylim=c(0,1),main="Standard Normal Distribution\n Bin width=0.0001")
lines(x,dnorm(x),col="red")
ap_rate= count/n

#Beta n=m版本
set.seed(10)
n <- 10000
g <- rep(0,n)
count <- 0

for(i in 1:n){
  y1<-runif(1,0,1)
  y2<-runif(1,0,1)
  if (y2<=4*(y1*(1-y1))){
    g[i]<-y1
    count <- count+1
  }else{
    g[i] <- NA
  }
}

k<-g[!is.na(g)]

x=seq(0,1,by=0.001)
hist(k,breaks="FD",freq=F,ylim=c(0,3),main="Beta Distribution n=m")
lines(x,5.93*x*(1-x),col="red")
ap_rate= count/n



#錯誤版
acceptance <- function(fx){
  while(TRUE){
    x <- rexp(1)
    y <- runif(1,0,(2*exp(1)/pi)^0.5*exp(-x))
    if (y<=fx(x)) return(x)
  }
}

fx <- function(x){
    return((2/pi)^0.5*exp(-x^2/2))
}

b <- function(){
  k <- rbinom(1,1,0.5)
  if(k==0) return(-1)
  else return(1)
}

set.seed(100)
n <- 10000
g <- rep(0,n)

for(i in 1:n){
  g[i] <- acceptance(fx)*b()
}

z=seq(-5,5,by=0.001)
hist(g,breaks=seq(-5,5,by=1),freq=F,ylim=c(0,1),main="¿ù»~ª©")
lines(z,dnorm(z),col="yellow")
lines(z,(2/pi)^0.5*exp(-z^2/2),col="red")
