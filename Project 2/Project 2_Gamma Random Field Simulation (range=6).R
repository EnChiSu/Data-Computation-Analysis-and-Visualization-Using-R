## gamma random filed模擬1000張

#步驟0: 參數設定
library(mvtnorm)
mu <- 4
r <- 0.5; r_x <- r; r_y <- r_x
alpha <- 16
lamda <- 4
s <- 1
size <- 80
range <- 6

#Gamma轉Gaussian的function建立
Ax <- 1+(r_x/6)^4;      Ay <- 1+(r_y/6)^4
Bx <- r_x/6-(r_x/6)^3;   By <- r_y/6-(r_y/6)^3
Cx <- ((r_x/6)^2)/3;     Cy <- ((r_y/6)^2)/3
f <- function(rho_uv, rho_xy, Ax, Bx, Cx, Ay, By, Cy){    
  return((Ax*Ay-3*Ax*Cy-3*Cx*Ay+9*Cx*Cy)*rho_uv+2*Bx*By*rho_uv^2+6*Cx*Cy*rho_uv^3-rho_xy)
}

#三、range=6的80x80網格往下迭代1000張

#A、產生range=6下的gaussian random field的sigma(如此就不用每迭代一張都重算一次)
#2.1 通過y(h)求解gamma matrix
#2.11 計算sigma gamma矩陣當中間隔距離為0的sigma gamma
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
sigma_gamma_33<-sigma_gamma_11
sigma_gamma_44<-sigma_gamma_11
sigma_gamma_44<-sigma_gamma_11
sigma_gamma_66<-sigma_gamma_11

#2.12 計算sigma gamma矩陣當中間隔距離為1的sigma gamma
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
sigma_gamma_23<-sigma_gamma_12
sigma_gamma_32<-sigma_gamma_12
sigma_gamma_34<-sigma_gamma_12
sigma_gamma_43<-sigma_gamma_12
sigma_gamma_45<-sigma_gamma_12
sigma_gamma_54<-sigma_gamma_12
sigma_gamma_56<-sigma_gamma_12
sigma_gamma_65<-sigma_gamma_12

#2.13 計算sigma gamma矩陣當中間隔距離為2的sigma gamma
sigma_gamma_13<-matrix(NA,nrow=size,ncol=size)
for (i in 1:80){
  for (j in 1:80){
    if(sqrt(2^2+(i-j)^2) >= range){
      sigma_gamma_13[i,j] <- 0
    }else{
      sigma_gamma_13[i,j] <-  s*exp(-3*(sqrt(2^2+(i-j)^2))/range)/s^2
    }
  }
}
sigma_gamma_31<-sigma_gamma_13
sigma_gamma_24<-sigma_gamma_13
sigma_gamma_42<-sigma_gamma_13
sigma_gamma_35<-sigma_gamma_13
sigma_gamma_53<-sigma_gamma_13
sigma_gamma_46<-sigma_gamma_13
sigma_gamma_64<-sigma_gamma_13

#2.14 計算sigma gamma矩陣當中間隔距離為3的sigma gamma
sigma_gamma_14<-matrix(NA,nrow=size,ncol=size)
for (i in 1:80){
  for (j in 1:80){
    if(sqrt(3^2+(i-j)^2) >= range){
      sigma_gamma_14[i,j] <- 0
    }else{
      sigma_gamma_14[i,j] <-  s*exp(-3*(sqrt(3^2+(i-j)^2))/range)/s^2
    }
  }
}
sigma_gamma_41<-sigma_gamma_14
sigma_gamma_25<-sigma_gamma_14
sigma_gamma_52<-sigma_gamma_14
sigma_gamma_36<-sigma_gamma_14
sigma_gamma_63<-sigma_gamma_14

#2.15 計算sigma gamma矩陣當中間隔距離為4的sigma gamma
sigma_gamma_15<-matrix(NA,nrow=size,ncol=size)
for (i in 1:80){
  for (j in 1:80){
    if(sqrt(4^2+(i-j)^2) >= range){
      sigma_gamma_15[i,j] <- 0
    }else{
      sigma_gamma_15[i,j] <-  s*exp(-3*(sqrt(4^2+(i-j)^2))/range)/s^2
    }
  }
}
sigma_gamma_51<-sigma_gamma_15
sigma_gamma_26<-sigma_gamma_15
sigma_gamma_62<-sigma_gamma_15

#2.16 計算sigma gamma矩陣當中間隔距離為5的sigma gamma
sigma_gamma_16<-matrix(NA,nrow=size,ncol=size)
for (i in 1:80){
  for (j in 1:80){
    if(sqrt(5^2+(i-j)^2) >= range){
      sigma_gamma_16[i,j] <- 0
    }else{
      sigma_gamma_16[i,j] <-  s*exp(-3*(sqrt(5^2+(i-j)^2))/range)/s^2
    }
  }
}
sigma_gamma_61<-sigma_gamma_16

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
sigma_gaussian_33<-sigma_gaussian_11
sigma_gaussian_44<-sigma_gaussian_11
sigma_gaussian_55<-sigma_gaussian_11
sigma_gaussian_66<-sigma_gaussian_11

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
sigma_gaussian_23<-sigma_gaussian_12
sigma_gaussian_32<-sigma_gaussian_12
sigma_gaussian_34<-sigma_gaussian_12
sigma_gaussian_43<-sigma_gaussian_12
sigma_gaussian_45<-sigma_gaussian_12
sigma_gaussian_54<-sigma_gaussian_12
sigma_gaussian_56<-sigma_gaussian_12
sigma_gaussian_65<-sigma_gaussian_12

sigma_gaussian_13 <- matrix(0,80,80)
for (i in 1:80){
  for (j in 1:80){
    if(sigma_gamma_13[i,j]==1){
      sigma_gaussian_13[i,j] <- 1
    }
    else if(sigma_gamma_13[i,j]==0){
      sigma_gaussian_13[i,j] <- 0
    }
    else{
      sigma_gaussian_13[i,j] <- uniroot(f, c(-1,1), Ax=Ax, Bx=Bx,Cx=Cx, Ay=Ay, By=By, Cy=Cy, rho_xy=sigma_gamma_13[i,j])$root
    }
  }
}
sigma_gaussian_31<-sigma_gaussian_13
sigma_gaussian_24<-sigma_gaussian_13
sigma_gaussian_42<-sigma_gaussian_13
sigma_gaussian_35<-sigma_gaussian_13
sigma_gaussian_53<-sigma_gaussian_13
sigma_gaussian_46<-sigma_gaussian_13
sigma_gaussian_64<-sigma_gaussian_13

sigma_gaussian_14 <- matrix(0,80,80)
for (i in 1:80){
  for (j in 1:80){
    if(sigma_gamma_14[i,j]==1){
      sigma_gaussian_14[i,j] <- 1
    }
    else if(sigma_gamma_14[i,j]==0){
      sigma_gaussian_14[i,j] <- 0
    }
    else{
      sigma_gaussian_14[i,j] <- uniroot(f, c(-1,1), Ax=Ax, Bx=Bx,Cx=Cx, Ay=Ay, By=By, Cy=Cy, rho_xy=sigma_gamma_14[i,j])$root
    }
  }
}
sigma_gaussian_41<-sigma_gaussian_14
sigma_gaussian_25<-sigma_gaussian_14
sigma_gaussian_52<-sigma_gaussian_14
sigma_gaussian_36<-sigma_gaussian_14
sigma_gaussian_63<-sigma_gaussian_14

sigma_gaussian_15 <- matrix(0,80,80)
for (i in 1:80){
  for (j in 1:80){
    if(sigma_gamma_15[i,j]==1){
      sigma_gaussian_15[i,j] <- 1
    }
    else if(sigma_gamma_15[i,j]==0){
      sigma_gaussian_15[i,j] <- 0
    }
    else{
      sigma_gaussian_15[i,j] <- uniroot(f, c(-1,1), Ax=Ax, Bx=Bx,Cx=Cx, Ay=Ay, By=By, Cy=Cy, rho_xy=sigma_gamma_15[i,j])$root
    }
  }
}
sigma_gaussian_51<-sigma_gaussian_15
sigma_gaussian_26<-sigma_gaussian_15
sigma_gaussian_62<-sigma_gaussian_15

sigma_gaussian_16 <- matrix(0,80,80)
for (i in 1:80){
  for (j in 1:80){
    if(sigma_gamma_16[i,j]==1){
      sigma_gaussian_16[i,j] <- 1
    }
    else if(sigma_gamma_16[i,j]==0){
      sigma_gaussian_16[i,j] <- 0
    }
    else{
      sigma_gaussian_16[i,j] <- uniroot(f, c(-1,1), Ax=Ax, Bx=Bx,Cx=Cx, Ay=Ay, By=By, Cy=Cy, rho_xy=sigma_gamma_16[i,j])$root
    }
  }
}
sigma_gaussian_61<-sigma_gaussian_16

#2.3 計算mu_star和sigma_star
library(MASS)
sigma_star <- sigma_gaussian_66 - cbind(sigma_gaussian_65,sigma_gaussian_64,sigma_gaussian_63,sigma_gaussian_62,sigma_gaussian_61)%*% ginv(rbind(cbind(sigma_gaussian_55,sigma_gaussian_54,sigma_gaussian_53,sigma_gaussian_52,sigma_gaussian_51),cbind(sigma_gaussian_45,sigma_gaussian_44,sigma_gaussian_43,sigma_gaussian_42,sigma_gaussian_41),cbind(sigma_gaussian_35,sigma_gaussian_34,sigma_gaussian_33,sigma_gaussian_32,sigma_gaussian_31),cbind(sigma_gaussian_25,sigma_gaussian_24,sigma_gaussian_23,sigma_gaussian_22,sigma_gaussian_21),cbind(sigma_gaussian_15,sigma_gaussian_14,sigma_gaussian_13,sigma_gaussian_12,sigma_gaussian_11))) %*%  rbind(sigma_gaussian_56,sigma_gaussian_46,sigma_gaussian_36,sigma_gaussian_26,sigma_gaussian_16)

#B、根據前面做好的sigma gaussian去產生gaussain random field的模擬
rg6 <- list()
for (q in 1:100){
  #步驟1: 前五列Gaussian Random Field模擬
  #library(mvtnorm)
  w <- matrix(NA, nrow = size, ncol = size)
  w[1,] <- rmvnorm(80,mean=c(0),sigma=diag(1))
  mu_star<- sigma_gaussian_21%*% solve(sigma_gaussian_11)%*% w[1,]
  w[2,]<-rmvnorm(1,mu_star,sigma_star)
  mu_star<- cbind(sigma_gaussian_32,sigma_gaussian_31)%*% ginv(rbind(cbind(sigma_gaussian_22,sigma_gaussian_21),cbind(sigma_gaussian_12,sigma_gaussian_11)))%*% matrix(c(w[2,],w[1,]),ncol=1)
  w[3,]<-rmvnorm(1,mu_star,sigma_star)
  mu_star<- cbind(sigma_gaussian_43,sigma_gaussian_42,sigma_gaussian_41)%*% ginv(rbind(cbind(sigma_gaussian_33,sigma_gaussian_32,sigma_gaussian_31),cbind(sigma_gaussian_23,sigma_gaussian_22,sigma_gaussian_21),cbind(sigma_gaussian_13,sigma_gaussian_12,sigma_gaussian_11)))%*% matrix(c(w[3,],w[2,],w[1,]),ncol=1)
  w[4,]<-rmvnorm(1,mu_star,sigma_star)
  mu_star<- cbind(sigma_gaussian_54,sigma_gaussian_53,sigma_gaussian_52,sigma_gaussian_51)%*% ginv(rbind(cbind(sigma_gaussian_44,sigma_gaussian_43,sigma_gaussian_42,sigma_gaussian_41),cbind(sigma_gaussian_34,sigma_gaussian_33,sigma_gaussian_32,sigma_gaussian_31),cbind(sigma_gaussian_24,sigma_gaussian_23,sigma_gaussian_22,sigma_gaussian_21),cbind(sigma_gaussian_14,sigma_gaussian_13,sigma_gaussian_12,sigma_gaussian_11)))%*% matrix(c(w[4,],w[3,],w[2,],w[1,]),ncol=1)
  w[5,]<-rmvnorm(1,mu_star,sigma_star)
  # 步驟2: 將模擬好的第五列Gaussian random field往下拓展75列
  for(i in 1:75){
    mu_star<- cbind(sigma_gaussian_65,sigma_gaussian_64,sigma_gaussian_63,sigma_gaussian_62,sigma_gaussian_61)%*% ginv(rbind(cbind(sigma_gaussian_55,sigma_gaussian_54,sigma_gaussian_53,sigma_gaussian_52,sigma_gaussian_51),cbind(sigma_gaussian_45,sigma_gaussian_44,sigma_gaussian_43,sigma_gaussian_42,sigma_gaussian_41),cbind(sigma_gaussian_35,sigma_gaussian_34,sigma_gaussian_33,sigma_gaussian_32,sigma_gaussian_31),cbind(sigma_gaussian_25,sigma_gaussian_24,sigma_gaussian_23,sigma_gaussian_22,sigma_gaussian_21),cbind(sigma_gaussian_15,sigma_gaussian_14,sigma_gaussian_13,sigma_gaussian_12,sigma_gaussian_11))) %*% matrix(c(w[i+4,],w[i+3,],w[i+2,],w[i+1,],w[i,]),ncol=1)
    w[i+5,]<-rmvnorm(1,mu_star,sigma_star)
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
  rg6[[q]] <- x
}
heatmap(rg6[[1]], Colv = NA, Rowv = NA)
#3D plot
# Library
library(plotly)
# Plot
s <- plot_ly(z = rg6[[1]], type = "surface")
s 


#5.2驗證1000張Gamma的參數
z <- c()
sq_error <- c()
cube_error <- c()
N<-100
for (i in 1:N){
  z[i] <- rg6[[i]][10,10]
}
mu_hat <- mean(z); mu_hat

for (i in 1:N){
  sq_error[i] <- (rg6[[i]][10,10]-mu_hat)^2
}
sigma_z <- sqrt(sum(sq_error)/(N-1)); sigma_z

for (i in 1:N){
  cube_error[i] <- (rg6[[i]][10,10]-mu_hat)^3
}
skew <- N/(N-1)/(N-2)/sigma_z^3*sum(cube_error); skew

