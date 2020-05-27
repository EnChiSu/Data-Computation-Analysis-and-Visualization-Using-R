#步驟0: 參數設定
library(mvtnorm)
mu <- 4
r <- 0.5; r_x <- r; r_y <- r_x
alpha <- 16
lamda <- 4
s <- 1
size <- 80
range <- 2 

#Gamma轉Gaussian的function建立
Ax <- 1+(r_x/6)^4;      Ay <- 1+(r_y/6)^4
Bx <- r_x/6-(r_x/6)^3;   By <- r_y/6-(r_y/6)^3
Cx <- ((r_x/6)^2)/3;     Cy <- ((r_y/6)^2)/3
f <- function(rho_uv, rho_xy, Ax, Bx, Cx, Ay, By, Cy){    
  return((Ax*Ay-3*Ax*Cy-3*Cx*Ay+9*Cx*Cy)*rho_uv+2*Bx*By*rho_uv^2+6*Cx*Cy*rho_uv^3-rho_xy)
}

#一、range=1的80x80網格往下迭代1000張
rg1 <- list()
for (q in 1:1000){
  p <- rmvnorm(80, matrix(0,nrow = 80, ncol = 1), diag(80))
  
  #步驟3: 將模擬好的Gaussian random field轉回Gamma random field
  u <- matrix(NA, nrow = size, ncol = size)
  alpha <- 16
  lamda <- 4
  for (i in 1:80){
    for (j in 1:80){
      u[i,j] <- alpha/lamda*(1-1/9/alpha+p[i,j]*(1/9/alpha)^0.5)^3
    }
  }
  rg1[[q]] <- u
}
heatmap(rg1[[1]], Colv = NA, Rowv = NA)
#3D plot
# Library
library(plotly)
# Plot
m <- plot_ly(z = rg1[[1]], type = "surface")
m

#二、range=2的80x80網格往下迭代1000張

#A、產生range=2下的gaussian random field的sigma(如此就不用每迭代一張都重算一次)
#2.1 通過y(h)求解gamma matrix
sigma_gamma_11 <- matrix(0, nrow = size, ncol = size)
for (i in 1:80){
  for (j in 1:80){
    if(abs(i-j)<range){
      sigma_gamma_11[i,j] <- (s*exp(-3*abs(i-j)/range))/s^2 
    }else{
      sigma_gamma_11[i,j] <- 0
    }
  }
}
sigma_gamma_22<-sigma_gamma_11

sigma_gamma_12<-matrix(NA,nrow=size,ncol=size)
for (i in 1:80){
  for (j in 1:80){
    if((1^2+(i-j)^2)^0.5>=range){
      sigma_gamma_12[i,j] <- 0
    }else{
      sigma_gamma_12[i,j] <-  s*exp(-3*(1^2+(i-j)^2)^0.5/range)/s^2
    }
  }
}
sigma_gamma_21<-sigma_gamma_12

#2.2 通過rho_xy公式轉化為Gaussian matrix
sigma_gaussian_11 <- matrix(0,80,80)
for (i in 1:80){
  for (j in 1:80){
    if(sigma_gamma_11[i,j]==1){
      sigma_gaussian_11[i,j] <- 1
    }
    else if(sigma_gamma_11[i,j]==0){
      sigma_gaussian_11[i,j] <- 0
    }
    else{
      sigma_gaussian_11[i,j] <- uniroot(f, c(-1,1), Ax=Ax, Bx=Bx,Cx=Cx, Ay=Ay, By=By, Cy=Cy, rho_xy=sigma_gamma_11[i,j])$root
    }
  }
}
sigma_gaussian_22<-sigma_gaussian_11

sigma_gaussian_12 <- matrix(0,80,80)
for (i in 1:80){
  for (j in 1:80){
    if(sigma_gamma_12[i,j]==1){
      sigma_gaussian_12[i,j] <- 1
    }
    else if(sigma_gamma_12[i,j]==0){
      sigma_gaussian_12[i,j] <- 0
    }
    else{
      sigma_gaussian_12[i,j] <- uniroot(f, c(-1,1), Ax=Ax, Bx=Bx,Cx=Cx, Ay=Ay, By=By, Cy=Cy, rho_xy=sigma_gamma_12[i,j])$root
    }
  }
}
sigma_gaussian_21<-sigma_gaussian_12

#2.3 計算mu_star和sigma_star
#mu_star<- sigma_gaussian_12%*% solve(sigma_gaussian_22)%*% w[1,]
sigma_star <- sigma_gaussian_22 - sigma_gaussian_21%*% solve(sigma_gaussian_11) %*% t(sigma_gaussian_21) 


#B、將做好的Gaussian random field的sigma
rg2 <- list()
for (q in 1:1000){
  #步驟1: 第一列Gaussian Random Field模擬
  #library(mvtnorm)
  w <- matrix(NA, nrow = size, ncol = size)
  w[1,] <- rmvnorm(80,mean=c(0),sigma=diag(1))
  # 步驟2: 將模擬好的第一列Gaussian random field往下拓展79列
  for(i in 1:79){
    mu_star<- sigma_gaussian_21%*% solve(sigma_gaussian_11)%*% w[i,]
    w[i+1,]<-rmvnorm(1,mu_star,sigma_star)
  }
  #步驟3: 將模擬好的Gaussian random field轉回Gamma random field
  x <- matrix(NA, nrow = size, ncol = size)
  alpha <- 16
  lamda <- 4
  for (i in 1:80){
    for (j in 1:80){
      x[i,j] <- alpha/lamda*(1-1/9/alpha+w[i,j]*(1/9/alpha)^0.5)^3
    }
  }
  rg2[[q]] <- x
}
heatmap(rg2[[1]], Colv = NA, Rowv = NA)
#3D plot
# Library
library(plotly)
# Plot
n <- plot_ly(z = rg2[[1]], type = "surface")
n 

#步驟四：畫圖
#par(mfrow=c(2,2))

#5.1 驗證range
semi<-matrix(data=NA,nrow=310497600,ncol=2)
ecd<-rep(NA,112)
count=0
for(x1 in 1:80){
  for(y1 in 1:80){
    for(x2 in x1:80){
      for(y2 in y1:80){
        count=count+1
        semi[count,1]<-((x1-x2)^2+(y1-y2)^2)^(1/2)
        semi[count,2]<-(1/2)*(x[x2,y2]-x[x1,y1])^2
      }
    }
  }
}
semi[,1]<-ceiling(semi[,1])
for(h in 0:100){
  sum=0
  count=0
  for(i in 1:10497600){
    if(semi[i,1]==h){
      sum=sum+semi[i,2]
      count=count+1
    }
    ecd[h]<- sum/count
  }
}
plot(ecd)

#5.2驗證1000張Gamma的參數
z <- c()
sq_error <- c()
cube_error <- c()
N<-1000
for (i in 1:N){
  z[i] <- rg2[[i]][80,75]
}
mu_hat <- mean(z); mu_hat

for (i in 1:N){
  sq_error[i] <- (rg2[[i]][80,75]-mu_hat)^2
}
sigma_z <- sqrt(sum(sq_error)/(N-1)); sigma_z

for (i in 1:N){
  cube_error[i] <- (rg2[[i]][80,75]-mu_hat)^3
}
skew <- N/(N-1)/(N-2)/sigma_z^3*sum(cube_error); skew

