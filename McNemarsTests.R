library(pROC)

linear <- read.csv("20170408_linear/IndexValues_2.csv")
densit <- read.csv("20170406_Density/IndexValues_0.5_GridSize_16.csv")

linear_lm1 <- glm(X0.1~Max.degree,data=linear)
linear_lm2 <- glm(X0.1~Connections.to.CA125+Connections.to.HE4, data=linear)
linear_lm3 <- glm(X0.1~Max.Page.Rank+Connections.to.CA125+Connections.to.HE4,data=linear)

density_lm1 <- glm(X0.1~Connections.to.CA125, data=densit)
density_lm2 <- glm(X0.1~Connections.to.CA125+Connections.to.HE4, data=densit)
density_lm3 <- glm(X0.1~Connections.to.CA125+Connections.to.HE4+Connections.to.MK, data=densit)

Results <- data.frame(CaseControl=densit$X0.1)

Results$L1 <- predict(linear_lm1,linear)
Results$L2 <- predict(linear_lm2,linear)
Results$L3 <- predict(linear_lm3,linear)
Results$D1 <- predict(density_lm1,densit)
Results$D2 <- predict(density_lm2,densit)
Results$D3 <- predict(density_lm3,densit)

L1 <- roc(CaseControl~L1, data=Results, plot=T)
L2 <- roc(CaseControl~L2, data=Results, plot=T)
L3 <- roc(CaseControl~L3, data=Results, plot=T)
D1 <- roc(CaseControl~D1, data=Results, plot=T)
D2 <- roc(CaseControl~D2, data=Results, plot=T)
D3 <- roc(CaseControl~D3, data=Results, plot=T)

threshold <- function(roc){
  temp <- which(roc$specificities>0.95)
  temp <- roc$sensitivities[min(temp)]
  temp <- which(roc$sensitivities==temp)
  temp <- roc$thresholds[max(temp)]
  return(temp)
}

L1_thresh <- threshold(L1)
L2_thresh <- threshold(L2)
L3_thresh <- threshold(L3)
D1_thresh <- threshold(D1)
D2_thresh <- threshold(D2)
D3_thresh <- threshold(D3)

#L1 vs D1
a <- length(with(Results, which(L1>L1_thresh & D1>D1_thresh)))
b <- length(with(Results, which(L1>L1_thresh & D1<D1_thresh)))
c <- length(with(Results, which(L1<L1_thresh & D1>D1_thresh)))
d <- length(with(Results, which(L1<L1_thresh & D1<D1_thresh)))

mat <- matrix(c(a,b,c,d),byrow=T, nrow=2)
mcnemar.test(mat)

#L2 vs D2
a <- length(with(Results, which(L2>L2_thresh & D2>D2_thresh)))
b <- length(with(Results, which(L2>L2_thresh & D2<D2_thresh)))
c <- length(with(Results, which(L2<L2_thresh & D2>D2_thresh)))
d <- length(with(Results, which(L2<L2_thresh & D2<D2_thresh)))

mat <- matrix(c(a,b,c,d),byrow=T, nrow=2)
mcnemar.test(mat)

#L3 vs D3
a <- length(with(Results, which(L3>L3_thresh & D3>D3_thresh)))
b <- length(with(Results, which(L3>L3_thresh & D3<D3_thresh)))
c <- length(with(Results, which(L3<L3_thresh & D3>D3_thresh)))
d <- length(with(Results, which(L3<L3_thresh & D3<D3_thresh)))

mat <- matrix(c(a,b,c,d),byrow=T, nrow=2)
mcnemar.test(mat)
