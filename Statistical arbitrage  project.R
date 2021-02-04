
load("/Users/auroracappadocian/Desktop/UCLA/statistical arbitrage/hw3/database.RData")
#### Problem Set 3 - Backtest
library(lubridate)
library(fBasics)
library(Jmisc)
library(matrixStats)

library(quadprog)

library(modopt.matlab)
library(Matrix)
library(dplyr)


#### Q1

n = length(allstocks) ## the number of stocks - 566
Tt = length(myday) ## the number of weekdays - 1504



#### Q2

myday_1 = as.Date(myday, format='%d-%B-%Y')
mon = month(myday_1)
mon_lag = c(NA, mon[1:Tt-1])
firstday = mon - mon_lag
index = which(firstday == 1 | firstday == -11 | is.na(firstday))

startIdx = index[13]
yrInterval = 12
shrink = rep(0, 60)
shrinkIdx = 1

## calculate arithmetic return
tri_return =  rbind(NA, diff(tri)) / tri 

for (i in 13:length(index)) {
  
  ## extract active day and stock
  active_now = which(isactivenow[index[i], ] == 1)
  past_1yr = tri_return[index[i - yrInterval]: index[i], active_now]
  past_1yr[is.na(past_1yr)] = 0 ## replace NA with 0
  t_1yr = nrow(past_1yr)
  n_1yr = ncol(past_1yr)
  
  ## calculate cov-var matrix
  S = var(past_1yr, na.rm=T)
  sigma2_bar = sum(diag(S)) / n_1yr
  shrinkage_target = sigma2_bar * diag(n_1yr)
  
  total_mse = 0
  for (j in 1:t_1yr) {
    Xt = past_1yr[j, ]
    total_mse = total_mse + sum((Xt %*% t(Xt) - S) ^2)
  }
  
  beta_hat = 1 - ((1 / ((t_1yr-1)*t_1yr)) * total_mse) / sum((S - shrinkage_target) ^2)
  shrink[shrinkIdx] = beta_hat
  shrinkIdx = shrinkIdx + 1
}



#### Q3

#### part 1 - alpha 1

## calculate industry dummy - rho
indus = c()
for (i in 1:n) {
  indus[i] = allstocks[[i]]$industrylist[2]$industry
}
indus = unique(indus)
rho = matrix(0, nrow = n, ncol = length(indus))

for (i in 1:n) {
  rho[i, which(indus == allstocks[[i]]$industrylist[2]$industry)] = 1
}


## calculate raw alpharev
w = 1/11 - 1/231 * seq(20, 0, -1)
alpharev_0 = matrix(0, nrow=Tt, ncol=n)

for (i in 21:length(myday)) {
  past_21d = tri_return[(i-20):i, ]
  past_21d[is.na(past_21d)] = 0 ## replace NA with 0
  alpharev_0[i, ] = -t(w) %*% past_21d
}


## calculate modified alpharev
alpharev = alpharev_0 %*% (diag(n) - rho %*% solve(t(rho) %*% rho) %*% t(rho))



#### part 2 - alpha 2

## calculate alpharec
w = 1/23 - 1/1035 * seq(44, 0, -1)
alpharec = matrix(0, nrow=Tt, ncol=n)

for (i in 45:length(myday)) {
  past_45d = rec[(i-44):i, ]
  past_45d[is.na(past_45d)] = 0 ## replace NA with 0
  alpharec[i, ] = t(w) %*% past_45d
}



#### part 3 - alpha 3

## calculate alphaval
alphaval = 1 / mtbv
alphaval[which(mtbv == 0)] = NA


#### part 4 - alpha 4

## calculate alphamom
#alphamom=rbind(matrix(0,nrow=231,ncol=n),tri[(252-21+1):length(myday),]/tri[1:(length(myday)-(252-21)),]-1)

alphamom=rbind(matrix(0,nrow=226,ncol=n),tri[227:length(myday),]/tri[1:(length(myday)-226),]-1)


demeanwin<-function(alphamean){
  for (i in 1:nrow(alphamean)){
    alphamean[i, ] = (alphamean[i, ] - mean(alphamean[i, ], na.rm = T)) / sd(alphamean[i, ], na.rm = T)
    alphamean[i, which(alphamean[i, ] < -3)] = -3
    alphamean[i, which(alphamean[i, ] > 3)] = 3
  }
  return(alphamean)
}

alpharev = demeanwin(alpharev)
alpharec = demeanwin(alpharec)
alphaval = demeanwin(alphaval)
alphamom = demeanwin(alphamom)

alphablend = demeanwin(0.5*alpharev+0.25*alpharec+0.15*alphaval+0.1*alphamom)



#### Q4
mu=1
lambda=0.1
w=matrix(0,nrow = 1, ncol = n)
trade=matrix(0,ncol=n, nrow=Tt)
back_weight = matrix(0,ncol=n,nrow = Tt)
daily_pnl = matrix(0, nrow = Tt, ncol=1)
tradesize = matrix(0, nrow = Tt, ncol=1)
booksize = matrix(0, nrow = Tt, ncol=1)
t0=246
r_star=300000
f_star=100000


retmat=rbind(rep(0,n),tri_return[2:nrow(tri_return),])
retmat[which(is.na(retmat))]=0
tcost[which(is.na(tcost))]=0




ctry = c()
for (i in 1:n) {
  ctry[i] = allstocks[[i]]$indexlist[2]$index
}
ctry = unique(ctry)[c(1:3,5:12)]
phi = matrix(0, nrow = n, ncol = length(ctry))

for (i in 1:n) {
  phi[i, which(ctry == allstocks[[i]]$indexlist[2]$index)] = 1
}


#calculate market beta
ew_mkt = rowMeans(retmat)
ew_mkt=as.matrix(ew_mkt)
dt = myday_1
mon = month(dt)
firDay = c(NA, mon[1:Tt-1])
firstDay = mon - firDay
firInd = which(firstDay != 0)
firInd = c(1,firInd)

yr = 12
mktbeta = matrix(0, nrow=length(firInd), ncol = n)
for (i in 13: length(firInd)) {
  pastyr = retmat[firInd[i-yr]:(firInd[i]-1), ]
  pastyr[which(is.na(pastyr))] = 0
  
  row = nrow(pastyr)
  col = ncol(pastyr)
  
  x = cbind(matrix(1, nrow=row, ncol=1), ew_mkt[firInd[i-yr]:(firInd[i]-1), ])
  beta = solve(t(x) %*% x) %*% t(x) %*% pastyr
  mktbeta[i, ] = beta[2, ]
  
}


#optimizer
w = matrix(0, nrow = n,ncol = 1)
for (i in t0:(Tt-1)){
  idx=which(isactivenow[i,]==1 & !is.na(alphablend[i-1,]))
  S=cov(retmat[(i-t0+1):(i-1),idx])
  target=mean(diag(S))*diag(1,length(idx))
  
  shrink_order=(year(myday_1[i])-1998)*12+month(myday_1[i])
  hat_beta=shrink[shrink_order]
  hat_sigma=(1-hat_beta)*target+hat_beta*S
  
  H=2*mu*as.matrix(rbind(cbind(hat_sigma, -hat_sigma),cbind(-hat_sigma,hat_sigma)))
  Hnew <- as.matrix(Matrix::nearPD(H)$mat)
  g=c((2*mu*hat_sigma%*%w[idx,1]-alphablend[i-1,idx]+lambda*tcost[i,idx]),
      (-2*mu*hat_sigma%*%w[idx,1]+alphablend[i-1,idx]+lambda*tcost[i,idx]))
  A=rbind(cbind(t(rho[idx,]),t(-rho[idx,])),cbind(t(-rho[idx,]),t(rho[idx,])),
          cbind(t(phi[idx,]),t(-phi[idx,])),cbind(t(-phi[idx,]),t(phi[idx,])))
  b=rbind(r_star*matrix(1,nrow = 40,ncol=1)-t(rho[idx,])%*%w[idx,1],
          r_star*matrix(1,nrow = 40,ncol=1)+t(rho[idx,])%*%w[idx,1],
          f_star*matrix(1,nrow = 11,ncol=1)-t(phi[idx,])%*%w[idx,1],
          f_star*matrix(1,nrow = 11,ncol=1)+t(phi[idx,])%*%w[idx,1])
  C = cbind(t(mktbeta[shrink_order,idx]), - t(mktbeta[shrink_order,idx]))
  d = - t(mktbeta[shrink_order,idx]) %*% w[idx,1]
  
  LB = rep(0,2 * length(idx))
  theta = pmin(volume[i,idx] * 0.01 * 1000, 150000, na.rm = T)
  pie = pmin(10 * theta, 0.025 * 50000000,na.rm = T)
  UB = c(pmax(0, pmin(t(theta), t(pie) - w[idx,1])), pmax(0, pmin(t(theta), t(pie) + w[idx,1])))
  
  solution = quadprog(Hnew, g, A, b, C, d, LB, UB)
  u=solution$x
  
  solution = quadprog(Hnew, g, A, b, C, d, LB, UB)
  u=solution$x
  if (is.na(u[1])) {
    UB = UB+0.1
    solution = quadprog(Hnew, g, A, b, C, d, LB, UB)
    u=solution$x
  }
  
  
  #[u, fval, exitflag, output]
  y = u[1:length(idx)]
  z = u[(length(idx)+1):length(u)]
  
  yy = rep(0,n)
  zz = rep(0,n)
  
  yy[idx] = y
  zz[idx] = z
  
  trade[i, ] = yy-zz
  tradesize[i] =sum(abs(trade[i, ]))*100
  
  w = matrix(t(w) * (1 + retmat[i, ]) + t(trade[i, ]),nrow = n,ncol = 1)
  back_weight[i+1, ] = t(w)
  booksize[i+1] = sum(abs(back_weight[i+1, ]))*100
  #daily_pnl[i] = sum(price[i-1,]*retmat[i, ] - price[i-1,]*tcost[i,]*(yy+zz))
  daily_pnl[i+1] = sum(back_weight[i+1, ] * retmat[i+1, ]) - sum(trade[i+1, ] * tcost[i+1, ])
  print(i)
}

pnl = cumsum(daily_pnl)
excess_return=daily_pnl/booksize
sharpe=mean(excess_return, na.rm =T)*sqrt(252)/sd(excess_return,na.rm =T)

#max  drawdown
dd=rep(0,Tt)
ddt=rep(0,Tt)
hwm_old=0
for (i in 257:Tt){
  hwm=max(pnl[1:i])
  if (hwm>hwm_old){
    hwm_old=hwm
  }
  else{
    ddt[i]=ddt[i-1]+1
  }
  dd[i]=hwm-pnl[i]
    
}
deepest_dd=max(dd)
longest_dd=max(ddt)

save(shrink,alpharev,alpharec,alphaval,alphamom,alphablend,
     lambda,mu,t0,trade,back_weight,pnl,booksize,tradesize,sharpe,longest_dd,deepest_dd,
     file='/Users/auroracappadocian/Desktop/UCLA/statistical arbitrage/hw3/data_hw3.RData')


