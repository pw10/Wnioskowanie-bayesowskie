library(mlbench)

data("BostonHousing")

data0 <- BostonHousing[c(1:100), c("medv","rm", "crim", "lstat", "ptratio", "indus")]
data <- BostonHousing[c(101:500), c("medv","rm", "crim", "lstat", "ptratio", "indus")]
dataprog <- BostonHousing[c(501:506), c("medv","rm", "crim", "lstat", "ptratio", "indus")]


#dysponuj¹c licznym zbiorem na mniejszej próbie oszacowano MNK,
#na podstawie którego ustalono parametru rozk³ a priori
model0 <- lm(data = data0, medv ~ .)
summ0 <- summary(model0)

#wzory np. W5 7 slajd
#wektor oszacowan  parametrow modelu beta0 na podst wczesniejszego badania
beta0 <- as.vector(model0$coefficients)
#format: macierz
X0 <- matrix(data = c( rep(1, 100) , data0[,2], data0[,3], data0[,4], data0[,5],data0[,6] ), ncol = 6)
sigma0 <- solve(t(X0) %*% X0)
#parametry rozk³adu IG
alfa0 <- nrow(data0)-ncol(data0)
delta0 <- sum(model0$residuals^2)


#model MNK dla tych samych danych
modelMNK <- lm(data = data, medv ~ .)
betaMNK <- modelMNK$coefficients

#test normalnoœci reszt - > dane panelowe, wiêc nie
#badamy autokorelacji i heteroskedastycznoœci
shapiro.test(modelMNK$residuals)
#test odrzuca H0 o normalnoœci
#ale mo¿na zak³adaæ asymptotyczny rozk³ad normalny



  
#funkcja wylicaj¹ca par rozk³ a posteriori
calc_pars <- function(X, SIGMA, y, beta0, alfa0, n, delta0){
  
  #format: macierz
  X <- matrix(data = c( rep(1, n) , X[,1], X[,2], X[,3], X[,4], X[,5] ), ncol = 6)
  SIGMA <- as.matrix(SIGMA)
  y <- as.matrix(y)
  
  SIGMA1 <- solve( t(X) %*% X + solve(SIGMA))
  
  beta1 <- SIGMA1 %*% ( t(X) %*% y + solve(SIGMA) %*% beta0)
  
  alfa1 <- alfa0 + n
  delta1 <- delta0 + t(y) %*% y - t(beta1) %*% solve(SIGMA1) %*% beta1 + t(beta0)%*%solve(SIGMA)%*%beta0
  
  return(list(alfa1 = alfa1, delta1 = as.numeric(delta1), sigma1 = SIGMA1, beta1 = beta1))
}

pars <- calc_pars(X = data[,-1], SIGMA = sigma0, y = data[,1], beta0 = beta0, alfa0 = alfa0, n = nrow(data), delta0 = delta0)

pars$alfa1
pars$delta1
pars$sigma1
pars$beta1



#rozk³ad a posteriori parametru sigma^2
library(bayesAB)
plotInvGamma(shape = alfa0/2, scale = delta0/2)
plotInvGamma(shape = pars$alfa1/2, scale = pars$delta1/2)

pars$delta1/2/(pars$alfa1/2-1)

#obliczanie parametrów rozk³adów brzegowych bet
rozkl_brzeg_beta <- function(beta1, sigma1, alfa1, delta1){
  
  beta1 <- as.vector(beta1)
  res <- list()
  
  for (i in 1:length(beta1)){
    scale_par <- delta1/alfa1*sigma1[i,i]
    res <- c(res, list(c(ex = beta1[i], scale = scale_par, alfa = alfa1)))
    names(res)[i] <- paste0("beta", as.character(i))
  }
  
  return(res)
}

rozklady_brzegowe <- rozkl_brzeg_beta(beta1 = pars$beta1, alfa1 = pars$alfa1,
                                      delta1 = pars$delta1, sigma1 = pars$sigma1)


rozklady_brzegowe0 <- rozkl_brzeg_beta(beta1 = beta0, alfa1 = alfa0,
                                       delta1 = delta0, sigma1 = sigma0)


###wykresy 
{
  library(metRology)
  par(mfrow=c(3,2))
  
  #beta1
  dist_beta1_0 <- dt.scaled(df = rozklady_brzegowe0$beta1[[3]], mean = rozklady_brzegowe0$beta1[[1]], sd = rozklady_brzegowe0$beta1[[2]], x = seq(-170,150) )
  plot(dist_beta1_0, x = seq(-170,150), xlab = "", main=expression(paste("Rozk³ad a priori ", beta[1] )), type='l', lwd=3)
  dist_beta1 <- dt.scaled(df = rozklady_brzegowe$beta1[[3]], mean = rozklady_brzegowe$beta1[[1]], sd = rozklady_brzegowe$beta1[[2]], x = seq(-100,100) )
  plot(dist_beta1, x = seq(-100,100), xlab = "", main=expression(paste("Rozk³ad a posteriori ", beta[1] )), type='l', lwd=3)
  #beta2
  dist_beta2_0 <- dt.scaled(df = rozklady_brzegowe0$beta2[[3]], mean = rozklady_brzegowe0$beta2[[1]], sd = rozklady_brzegowe0$beta2[[2]], x = seq(5,9,0.01) )
  plot(dist_beta2_0, x = seq(5,9,0.01), xlab = "",main=expression(paste("Rozk³ad a priori ", beta[2] )), type='l', lwd=3)
  dist_beta2 <- dt.scaled(df = rozklady_brzegowe$beta2[[3]], mean = rozklady_brzegowe$beta2[[1]], sd = rozklady_brzegowe$beta2[[2]], x = seq(3,6,0.01) )
  plot(dist_beta2, x = seq(3,6,0.01), xlab = "", main=expression(paste("Rozk³ad a posteriori ", beta[2] )), type='l', lwd=3)
  #beta3
  dist_beta3_0<- dt.scaled(df = rozklady_brzegowe0$beta3[[3]], mean = rozklady_brzegowe0$beta3[[1]], sd = rozklady_brzegowe0$beta3[[2]], x = seq(-6,1,0.01) )
  plot(dist_beta3_0,  x = seq(-6,1,0.01), xlab = "",main=expression(paste("Rozk³ad a priori ", beta[3] )), type='l', lwd=3)
  dist_beta3 <- dt.scaled(df = rozklady_brzegowe$beta3[[3]], mean = rozklady_brzegowe$beta3[[1]], sd = rozklady_brzegowe$beta3[[2]], x = seq(-0.08,-0.06,0.00001) )
  plot(dist_beta3, x = seq(-0.08,-0.06,0.00001), xlab = "", main=expression(paste("Rozk³ad a posteriori ", beta[3] )), type='l', lwd=3)
  
  par(mfrow=c(3,2))
  #beta4
  dist_beta4_0 <- dt.scaled(df = rozklady_brzegowe0$beta4[[3]], mean = rozklady_brzegowe0$beta4[[1]], sd = rozklady_brzegowe0$beta4[[2]], x = seq(-.31,-.25,0.0001) )
  plot(dist_beta4_0, x = seq(-.31,-.25,0.0001), xlab = "", main=expression(paste("Rozk³ad a priori ", beta[4] )), type='l', lwd=3)
  dist_beta4 <- dt.scaled(df = rozklady_brzegowe$beta4[[3]], mean = rozklady_brzegowe$beta4[[1]], sd = rozklady_brzegowe$beta4[[2]], x = seq(-.57,-.54,0.0001) )
  plot(dist_beta4, x = seq(-.57,-.54,0.0001), xlab = "", main=expression(paste("Rozk³ad a posteriori ", beta[4] )), type='l', lwd=3)
  #beta5
  dist_beta5_0 <- dt.scaled(df = rozklady_brzegowe0$beta5[[3]], mean = rozklady_brzegowe0$beta5[[1]], sd = rozklady_brzegowe0$beta5[[2]], x = seq(-.3,0,0.001) )
  plot(dist_beta5_0, x = seq(-.3, 0 ,0.001), xlab = "", main=expression(paste("Rozk³ad a priori ", beta[5] )), type='l', lwd=3)
  dist_beta5 <- dt.scaled(df = rozklady_brzegowe$beta5[[3]], mean = rozklady_brzegowe$beta5[[1]], sd = rozklady_brzegowe$beta5[[2]], x = seq(-.95,-.8,0.001) )
  plot(dist_beta5, x = seq(-.95,-.8,0.001), xlab = "", main=expression(paste("Rozk³ad a posteriori ", beta[5] )), type='l', lwd=3)
  #beta6
  dist_beta6_0 <- dt.scaled(df = rozklady_brzegowe0$beta6[[3]], mean = rozklady_brzegowe0$beta6[[1]], sd = rozklady_brzegowe0$beta6[[2]], x = seq(-.2,-.12,0.0001) )
  plot(dist_beta6_0, x = seq(-.2, -.12, 0.0001), xlab = "", main=expression(paste("Rozk³ad a priori ", beta[6] )), type='l', lwd=3)
  dist_beta6 <- dt.scaled(df = rozklady_brzegowe$beta6[[3]], mean = rozklady_brzegowe$beta6[[1]], sd = rozklady_brzegowe$beta6[[2]], x = seq(0.01,0.04,0.0001) )
  plot(dist_beta6, x = seq(0.01,0.04,0.0001), xlab = "", main=expression(paste("Rozk³ad a posteriori ", beta[6] )), type='l', lwd=3)
}



#bayesowskie estymatory parametrów modelu - wartoœæ oczekiwana rozk³adu aposteriori
#czyli to ju¿ mamy
pars$beta1
#wyznaczyæ ich HPDI
library(HDInterval)

set.seed(1234)

hpdi <- list()

for( i in 1:length(rozklady_brzegowe)) {
  tmp <- hdi(rt.scaled(n = 1000,df = rozklady_brzegowe[[i]][[3]], mean = rozklady_brzegowe[[i]][[1]], sd = rozklady_brzegowe[[i]][[2]]))
  hpdi <- c(hpdi, list(c(tmp[1], tmp[2])))
  names(hpdi)[i] <- paste0("hpdi_beta", as.character(i))
}


#prognozowanie bayes

#W5 slajd10

alfa1 <- pars$alfa1

#format: macierz
data.matrix(dataprog[,-1])
Xt <- matrix(data = c( rep(1, nrow(dataprog)) , dataprog[,2], dataprog[,3], dataprog[,4], dataprog[,5], dataprog[,6] ), ncol = 6)
ex <- Xt %*% as.matrix(pars$beta1)

scale <- (pars$delta1 / alfa1) * (diag(1,ncol = nrow(dataprog), nrow =nrow(dataprog) ) +  Xt %*% as.matrix(pars$sigma1) %*% t(Xt) )



#rozk³ady yp
rozklady_brzegowe_y <- list()
  
for (i in 1:length(ex)){
    
  rozklady_brzegowe_y <- c(rozklady_brzegowe_y, list(c(ex = ex[i], scale = scale[i,i], alfa = alfa1)))
  names(rozklady_brzegowe_y)[i] <- paste0("yp", as.character(i))
}

#ANALOGICZNIE NARYSOWAC WYKRESY BRZEGOWYCH - JAK WCZESNIEJ 

#prognoza przedzia³owa bayes - hpdi 
hpdi_y <- list()

for( i in 1:length(rozklady_brzegowe_y)) {
  tmp <- hdi(rt.scaled(n = 1000,df = rozklady_brzegowe_y[[i]][[3]], mean = rozklady_brzegowe_y[[i]][[1]], sd = rozklady_brzegowe_y[[i]][[2]]))
  hpdi_y <- c(hpdi_y, list(c(tmp[1], tmp[2])))
  names(hpdi_y)[i] <- paste0("hpdi_y", as.character(i))
}


###restrykcje bayesowski sposob 
{
  
#indus i crim
model0R <- lm(data = data0, medv ~ rm + lstat + ptratio)
summ0R <- summary(model0R)

#wzory np. W5 7 slajd
#wektor oszacowan  parametrow modelu beta0 na podst wczesniejszego badania
beta0R <- as.vector(model0R$coefficients)
#format: macierz
X0R <- matrix(data = c( rep(1, 100) , data0[,2],  data0[,4], data0[,5]), ncol = 4)
sigma0R <- solve(t(X0R) %*% X0R)
#parametry rozk³adu IG
alfa0R <- nrow(data0)-ncol(data0)-2
delta0R <- sum(model0R$residuals^2)

calc_pars <- function(X, SIGMA, y, beta0, alfa0, n, delta0){
  
  #format: macierz
  X <- matrix(data = c( rep(1, n) , X[,1], X[,2], X[,3]), ncol = 4)
  SIGMA <- as.matrix(SIGMA)
  y <- as.matrix(y)
  
  SIGMA1 <- solve( t(X) %*% X + solve(SIGMA))
  
  beta1 <- SIGMA1 %*% ( t(X) %*% y + solve(SIGMA) %*% beta0)
  
  alfa1 <- alfa0 + n
  delta1 <- delta0 + t(y) %*% y - t(beta1) %*% solve(SIGMA1) %*% beta1 + t(beta0)%*%solve(SIGMA)%*%beta0
  
  return(list(alfa1 = alfa1, delta1 = as.numeric(delta1), sigma1 = SIGMA1, beta1 = beta1))
}

parsR <- calc_pars(X = data[,-1], SIGMA = sigma0R, y = data[,1], beta0 = beta0R, delta0 = delta0R, alfa0 = alfa0R, n = nrow(data))



pyM1_p1 <- sqrt(det(sigma0)/det(pars$sigma1))
pyM2_p1 <- sqrt(det(sigma0R)/det(parsR$sigma1))
cz1 <- pyM1_p1/pyM2_p1



cz2 <- ((delta0/delta0R)^46)*delta0*((parsR$delta1/pars$delta1)^246)/pars$delta1

cz3 <- 246/46

R12 <- cz1*cz2*cz3

log(cz1*cz3)

}

####test F
{
  #wyrzucamy indus i crim
  library(car)
  linearHypothesis(modelMNK, c("indus=0", "crim=0"))
  #P>0.05 - nie odrzucamy H0, ne sa laczie istotne
}

#ALGORYTM GIBSA
{
  #wzory np. W5 7 slajd
  #wektor oszacowan  parametrow modelu beta0 na podst wczesniejszego badania
  
  #za sime start bierzemy wartoœæ oczekiwan¹ z rozk³adu odwrotnego gamma - beta/(alfa-1)
  sigma_start <- (delta0/2)/(alfa0/2 - 1)
  losowanko <- 11111
  
  wyniki <- data.frame(1,2,3,4,5,6,7)
  colnames(wyniki) <- c(names(modelMNK$coefficients),'sigma')
  
  X <- matrix(data = c( rep(1, nrow(data)) ,data[,2], data[,3], data[,4], data[,5], data[,6]), ncol = 6)
  XtX <- t(X) %*% X
  Xty <- t(X) %*% as.matrix(data[,1])
  y <- data[,1]
  
  tmp_sigma_2 <- 1/sigma_start
  
  library(MASS)
  library(invgamma)
  
  for(i in 1:losowanko) {
    tmp_sigma_duza = solve(tmp_sigma_2 * XtX + solve(sigma0))
    tmp_beta <- tmp_sigma_duza %*% (tmp_sigma_2 * Xty + solve(sigma0)%*%beta0)
    
    wyniki[i,1:6] <- mvrnorm(n=1, mu=tmp_beta, Sigma=tmp_sigma_duza)
    
    
    tmp_delta <- delta0 + t(y - X%*% t(as.matrix(wyniki[i,1:6])) ) %*% (y-X%*% t(as.matrix(wyniki[i,1:6])) )
    wyniki[i,7] <- rinvgamma(n=1, shape = pars$alfa1/2,rate = tmp_delta/2)
    
    
    tmp_sigma_2 <- 1/wyniki[i,7]
  }
  
  wyniki <- wyniki[-(1:1111),]
  
  prognozy <- c()
  odchylenie <- c()
  for(i in 1:ncol(wyniki)) {
    prognozy[i] <- mean(wyniki[,i])
    odchylenie[i] <- sd(wyniki[,i])
  }
  round(prognozy,3)
  round(odchylenie, 3)
  
  
  par(mfrow=c(1,1))
  
  plot(density(wyniki$rm))
  plot(density(wyniki$crim))
  plot(density(wyniki$lstat))
  plot(density(wyniki$sigma))

}


###wykresy gibs
{
  par(mfrow=c(3,2))
  
  #beta1

  plot(dist_beta1, x = seq(-100,100), xlab = "", main=expression(paste("Rozk³ad a posteriori (analityczny)", beta[1] )), type='l', lwd=3)
  plot(density(wyniki$`(Intercept)`), xlab = "", main=expression(paste("Rozk³ad a posteriori (Gibbs) ", beta[1] )), lwd=3)

  #beta2

  plot(dist_beta2, x = seq(3,6,0.01), xlab = "", main=expression(paste("Rozk³ad a posteriori (analityczny)", beta[2] )), type='l', lwd=3)
  plot(density(wyniki$rm), xlab = "", main=expression(paste("Rozk³ad a posteriori (Gibbs) ", beta[2] )), lwd=3)
  
  #beta3
 
  plot(dist_beta3, x = seq(-0.08,-0.06,0.00001), xlab = "", main=expression(paste("Rozk³ad a posteriori (analityczny)", beta[3] )), type='l', lwd=3)
  plot(density(wyniki$crim), xlab = "", main=expression(paste("Rozk³ad a posteriori (Gibbs) ", beta[3] )), lwd=3)
  
  par(mfrow=c(3,2))
  #beta4

  plot(dist_beta4, x = seq(-.57,-.54,0.0001), xlab = "", main=expression(paste("Rozk³ad a posteriori (analityczny)", beta[4] )), type='l', lwd=3)
  plot(density(wyniki$lstat), xlab = "", main=expression(paste("Rozk³ad a posteriori (Gibbs) ", beta[4] )), lwd=3)
  
  #beta5

  plot(dist_beta5, x = seq(-.95,-.8,0.001), xlab = "", main=expression(paste("Rozk³ad a posteriori (analityczny)", beta[5] )), type='l', lwd=3)
  plot(density(wyniki$ptratio), xlab = "", main=expression(paste("Rozk³ad a posteriori (Gibbs) ", beta[5] )), lwd=3)
  
  #beta6

  plot(dist_beta6, x = seq(0.01,0.04,0.0001), xlab = "", main=expression(paste("Rozk³ad a posteriori (analityczny)", beta[6] )), type='l', lwd=3)
  plot(density(wyniki$indus), xlab = "", main=expression(paste("Rozk³ad a posteriori (Gibbs) ", beta[6] )), lwd=3)
  
  
  plotInvGamma(shape = pars$alfa1/2, scale = pars$delta1/2)
  plot(density(wyniki$sigma), xlab = "", main= "Rozk³ad a posteriori parametru sigma (Gibbs)", lwd = 3)
  
  
}

