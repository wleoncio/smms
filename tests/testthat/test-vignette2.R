## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "##"
)

## ----message=FALSE------------------------------------------------------------
library(smms)
library(igraph) # For specifying the multi-state graph
library(msm) # To get the CAV dataset

dd = cav
dd = dd[!is.na(dd$pdiag),]

# Remove observations where the patient appears to go back to a previous state
# (assumed to be impossible):
id_wrong = unique(dd$PTNUM[which(dd$state!=dd$statemax)])
dd = dd[-which(dd$PTNUM %in% id_wrong),]

dd = dd[ ,-c(2, 5, 7, 9, 10)]
# rename relevant columns (necessary in current version):
colnames(dd)[1:2] <- c("patient","time")

# Covariates:
dd$dage_st = (dd$dage-mean(dd$dage))/sd(dd$dage)
dd$ihd = (dd$pdiag=="IHD")
colnames(dd)[1:2] <- c("patient","time")

X_data = aggregate(dd[,c("dage_st","ihd")],by=list(dd$patient),FUN=median)
#FUN does not matter
X_data = as.matrix(X_data[,2:3])

print(dd[1:11,])

## ----fig.height=2.5,fig.align = "center"--------------------------------------
# Specify the graph:
gg = graph_from_literal("1"--+"2"--+"3"--+"4", "1"--+"4",
                        "2"--+"4")
par(mar=c(1,1,1,1))
plot(gg,layout=layout_with_sugiyama(gg,layers=c(1,1,1,2))$layout,vertex.size=20)

## -----------------------------------------------------------------------------
S_01 = function(param, x, t){(1-pexp(t,exp(param[1]+param[2]*x[1]+param[3]*x[2])))}
S_12 = function(param, x, t){(1-pexp(t,exp(param[4]+param[5]*x[1]+param[6]*x[2])))}
S_23 = function(param, x, t){(1-pexp(t,exp(param[7]+param[8]*x[1]+param[9]*x[2])))}
S_03 = function(param, x, t){(1-pexp(t,exp(param[10]+param[11]*x[1]+param[12]*x[2])))}
S_13 = function(param, x, t){(1-pexp(t,exp(param[13]+param[14]*x[1]+param[15]*x[2])))}

f_01 = function(param, x, t){dexp(t,exp(param[1]+param[2]*x[1]+param[3]*x[2]))}
f_12 = function(param, x, t){dexp(t,exp(param[4]+param[5]*x[1]+param[6]*x[2]))}
f_23 = function(param, x, t){dexp(t,exp(param[7]+param[8]*x[1]+param[9]*x[2]))}
f_03 = function(param, x, t){dexp(t,exp(param[10]+param[11]*x[1]+param[12]*x[2]))}
f_13 = function(param, x, t){dexp(t,exp(param[13]+param[14]*x[1]+param[15]*x[2]))}

## -----------------------------------------------------------------------------
print(names_of_survival_density(gg))

## ----eval=F-------------------------------------------------------------------
#  startval <- c(-2.5,0.21,0.5,-1.1,-0.2,0.2,-1.2,-0.1,-0.5,
#                -3.1,0.5,0.3,-2.8,-1.1,1.1)
#
#  mlo <- smms(startval,dd,gg,X_data, mc_cores = 1,hessian_matrix=T)
#  #calculates hessian too, slow with 1 core

## ----echo=F-------------------------------------------------------------------
load("~/UiO/Projects/smms/vignettes/cav_expo_cov_optims")

## -----------------------------------------------------------------------------
# Compute AIC (higher values are better with this definition)
aic <- (-2*mlo$opt$objective)-2*length(mlo$opt$par) #-2851.2

# Look at estimates and 95% confidence intervals.
# On the -Inf to Inf scale:
print(round(est_ci(mlo$opt$par,mlo$hess),2))

# On the 0 to Inf scale (on the transition intensity scale):
round(exp(est_ci(mlo$opt$par,mlo$hess)),2)

## -----------------------------------------------------------------------------
tval <- seq(0.01,30,length=50)
p0_y0_ci <- occupancy_prob_ci_band("1",tval,mlo$opt$par,gg,xval=c(-1,0),mlo$hess)
p0_y1_ci <- occupancy_prob_ci_band("1",tval,mlo$opt$par,gg,xval=c(-1,1),mlo$hess)
p0_o0_ci <- occupancy_prob_ci_band("1",tval,mlo$opt$par,gg,xval=c(1,0),mlo$hess)
p0_o1_ci <- occupancy_prob_ci_band("1",tval,mlo$opt$par,gg,xval=c(1,1),mlo$hess)
p1_y0_ci <- occupancy_prob_ci_band("2",tval,mlo$opt$par,gg,xval=c(-1,0),mlo$hess)
p1_y1_ci <- occupancy_prob_ci_band("2",tval,mlo$opt$par,gg,xval=c(-1,1),mlo$hess)
p1_o0_ci <- occupancy_prob_ci_band("2",tval,mlo$opt$par,gg,xval=c(1,0),mlo$hess)
p1_o1_ci <- occupancy_prob_ci_band("2",tval,mlo$opt$par,gg,xval=c(1,1),mlo$hess)
p2_y0_ci <- occupancy_prob_ci_band("3",tval,mlo$opt$par,gg,xval=c(-1,0),mlo$hess)## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "##"
)

## ----message=FALSE------------------------------------------------------------
library(smms)
library(igraph) # For specifying the multi-state graph
library(msm) # To get the CAV dataset

dd = cav
dd = dd[!is.na(dd$pdiag),]

# Remove observations where the patient appears to go back to a previous state
# (assumed to be impossible):
id_wrong = unique(dd$PTNUM[which(dd$state!=dd$statemax)])
dd = dd[-which(dd$PTNUM %in% id_wrong),]

dd = dd[ ,-c(2, 5, 7, 9, 10)]
# rename relevant columns (necessary in current version):
colnames(dd)[1:2] <- c("patient","time")

# Covariates:
dd$dage_st = (dd$dage-mean(dd$dage))/sd(dd$dage)
dd$ihd = (dd$pdiag=="IHD")
colnames(dd)[1:2] <- c("patient","time")

X_data = aggregate(dd[,c("dage_st","ihd")],by=list(dd$patient),FUN=median)
#FUN does not matter
X_data = as.matrix(X_data[,2:3])

print(dd[1:11,])

## ----fig.height=2.5,fig.align = "center"--------------------------------------
# Specify the graph:
gg = graph_from_literal("1"--+"2"--+"3"--+"4", "1"--+"4",
                        "2"--+"4")
par(mar=c(1,1,1,1))
plot(gg,layout=layout_with_sugiyama(gg,layers=c(1,1,1,2))$layout,vertex.size=20)

## -----------------------------------------------------------------------------
S_01 = function(param, x, t){(1-pexp(t,exp(param[1]+param[2]*x[1]+param[3]*x[2])))}
S_12 = function(param, x, t){(1-pexp(t,exp(param[4]+param[5]*x[1]+param[6]*x[2])))}
S_23 = function(param, x, t){(1-pexp(t,exp(param[7]+param[8]*x[1]+param[9]*x[2])))}
S_03 = function(param, x, t){(1-pexp(t,exp(param[10]+param[11]*x[1]+param[12]*x[2])))}
S_13 = function(param, x, t){(1-pexp(t,exp(param[13]+param[14]*x[1]+param[15]*x[2])))}

f_01 = function(param, x, t){dexp(t,exp(param[1]+param[2]*x[1]+param[3]*x[2]))}
f_12 = function(param, x, t){dexp(t,exp(param[4]+param[5]*x[1]+param[6]*x[2]))}
f_23 = function(param, x, t){dexp(t,exp(param[7]+param[8]*x[1]+param[9]*x[2]))}
f_03 = function(param, x, t){dexp(t,exp(param[10]+param[11]*x[1]+param[12]*x[2]))}
f_13 = function(param, x, t){dexp(t,exp(param[13]+param[14]*x[1]+param[15]*x[2]))}

## -----------------------------------------------------------------------------
print(names_of_survival_density(gg))

## ----eval=F-------------------------------------------------------------------
#  startval <- c(-2.5,0.21,0.5,-1.1,-0.2,0.2,-1.2,-0.1,-0.5,
#                -3.1,0.5,0.3,-2.8,-1.1,1.1)
#
#  mlo <- smms(startval,dd,gg,X_data, mc_cores = 1,hessian_matrix=T)
#  #calculates hessian too, slow with 1 core

## ----echo=F-------------------------------------------------------------------
load("~/UiO/Projects/smms/vignettes/cav_expo_cov_optims")

## -----------------------------------------------------------------------------
# Compute AIC (higher values are better with this definition)
aic <- (-2*mlo$opt$objective)-2*length(mlo$opt$par) #-2851.2

# Look at estimates and 95% confidence intervals.
# On the -Inf to Inf scale:
print(round(est_ci(mlo$opt$par,mlo$hess),2))

# On the 0 to Inf scale (on the transition intensity scale):
round(exp(est_ci(mlo$opt$par,mlo$hess)),2)

## -----------------------------------------------------------------------------
tval <- seq(0.01,30,length=50)
p0_y0_ci <- occupancy_prob_ci_band("1",tval,mlo$opt$par,gg,xval=c(-1,0),mlo$hess)
p0_y1_ci <- occupancy_prob_ci_band("1",tval,mlo$opt$par,gg,xval=c(-1,1),mlo$hess)
p0_o0_ci <- occupancy_prob_ci_band("1",tval,mlo$opt$par,gg,xval=c(1,0),mlo$hess)
p0_o1_ci <- occupancy_prob_ci_band("1",tval,mlo$opt$par,gg,xval=c(1,1),mlo$hess)
p1_y0_ci <- occupancy_prob_ci_band("2",tval,mlo$opt$par,gg,xval=c(-1,0),mlo$hess)
p1_y1_ci <- occupancy_prob_ci_band("2",tval,mlo$opt$par,gg,xval=c(-1,1),mlo$hess)
p1_o0_ci <- occupancy_prob_ci_band("2",tval,mlo$opt$par,gg,xval=c(1,0),mlo$hess)
p1_o1_ci <- occupancy_prob_ci_band("2",tval,mlo$opt$par,gg,xval=c(1,1),mlo$hess)
p2_y0_ci <- occupancy_prob_ci_band("3",tval,mlo$opt$par,gg,xval=c(-1,0),mlo$hess)
p2_y1_ci <- occupancy_prob_ci_band("3",tval,mlo$opt$par,gg,xval=c(-1,1),mlo$hess)
p2_o0_ci <- occupancy_prob_ci_band("3",tval,mlo$opt$par,gg,xval=c(1,0),mlo$hess)
p2_o1_ci <- occupancy_prob_ci_band("3",tval,mlo$opt$par,gg,xval=c(1,1),mlo$hess)
p3_y0_ci <- occupancy_prob_ci_band("4",tval,mlo$opt$par,gg,xval=c(-1,0),mlo$hess)
p3_y1_ci <- occupancy_prob_ci_band("4",tval,mlo$opt$par,gg,xval=c(-1,1),mlo$hess)
p3_o0_ci <- occupancy_prob_ci_band("4",tval,mlo$opt$par,gg,xval=c(1,0),mlo$hess)
p3_o1_ci <- occupancy_prob_ci_band("4",tval,mlo$opt$par,gg,xval=c(1,1),mlo$hess)

## ----fig.width=8,fig.align = "center",fig.height=5----------------------------
# msm package  (to get non-parameteric prevalence curves)
oneway4.q <- rbind(c(0, 0.25, 0, 0.25), c(0, 0, 0.25, 0.25),c(0, 0, 0, 0.5),
                   c(0, 0, 0, 0))
rownames(oneway4.q) <- colnames(oneway4.q) <- c("Well", "Mild","Severe",
                                                "Death")
cav.msm <- msm(state ~ time, subject = patient, data = dd,qmatrix = oneway4.q,
               death = 4,method = "BFGS")

prev_np <- prevalence.msm(cav.msm,times=tval)

# Plot
par(mfrow=c(2,2))
par(bty="l")
par(mar=c(4,4,1,2))
par(cex=1)
plot(tval,prev_np$`Observed percentages`[,1],type="l",ylim=c(0,100),lwd=3,
     xlab=" ",col="dark grey",ylab="prevalence (%)",main="well")
polygon(c(tval, rev(tval)), c(p0_y0_ci$upper*100, rev(p0_y0_ci$lower*100)),
        col = adjustcolor("#0571b0",alpha.f=0.2), border=NA)
polygon(c(tval, rev(tval)), c(p0_y1_ci$upper*100, rev(p0_y1_ci$lower*100)),
        col = adjustcolor("#ca0020",alpha.f=0.2), border=NA)
polygon(c(tval, rev(tval)), c(p0_o0_ci$upper*100, rev(p0_o0_ci$lower*100)),
        col = adjustcolor("#0571b0",alpha.f=0.2), border=NA)
polygon(c(tval, rev(tval)), c(p0_o1_ci$upper*100, rev(p0_o1_ci$lower*100)),
        col = adjustcolor("#ca0020",alpha.f=0.2), border=NA)
lines(tval,p0_y0_ci$est*100,col="#0571b0",lwd=3)
lines(tval,p0_y1_ci$est*100,col="#ca0020",lwd=3)
lines(tval,p0_o0_ci$est*100,col="#0571b0",lwd=3,lty=2)
lines(tval,p0_o1_ci$est*100,col="#ca0020",lwd=3,lty=2)

plot(tval,prev_np$`Observed percentages`[,2],type="l",col="dark grey",
     ylim=c(0,100),lwd=3,xlab=" ",ylab="prevalence (%)",main="mild")
polygon(c(tval, rev(tval)), c(p1_y0_ci$upper*100, rev(p1_y0_ci$lower*100)),
        col = adjustcolor("#0571b0",alpha.f=0.2), border=NA)
polygon(c(tval, rev(tval)), c(p1_y1_ci$upper*100, rev(p1_y1_ci$lower*100)),
        col = adjustcolor("#ca0020",alpha.f=0.2), border=NA)
polygon(c(tval, rev(tval)), c(p1_o0_ci$upper*100, rev(p1_o0_ci$lower*100)),
        col = adjustcolor("#0571b0",alpha.f=0.2), border=NA)
polygon(c(tval, rev(tval)), c(p1_o1_ci$upper*100, rev(p1_o1_ci$lower*100)),
        col = adjustcolor("#ca0020",alpha.f=0.2), border=NA)
lines(tval,p1_y0_ci$est*100,col="#0571b0",lwd=3)
lines(tval,p1_y1_ci$est*100,col="#ca0020",lwd=3)
lines(tval,p1_o0_ci$est*100,col="#0571b0",lwd=3,lty=2)
lines(tval,p1_o1_ci$est*100,col="#ca0020",lwd=3,lty=2)


plot(tval,prev_np$`Observed percentages`[,3],type="l",col="dark grey",
     ylim=c(0,100),lwd=3,xlab="years after transplantation",
     ylab="prevalence (%)",main="severe")
polygon(c(tval, rev(tval)), c(p2_y0_ci$upper*100, rev(p2_y0_ci$lower*100)),
        col = adjustcolor("#0571b0",alpha.f=0.2), border=NA)
polygon(c(tval, rev(tval)), c(p2_y1_ci$upper*100, rev(p2_y1_ci$lower*100)),
        col = adjustcolor("#ca0020",alpha.f=0.2), border=NA)
polygon(c(tval, rev(tval)), c(p2_o0_ci$upper*100, rev(p2_o0_ci$lower*100)),
        col = adjustcolor("#0571b0",alpha.f=0.2), border=NA)
polygon(c(tval, rev(tval)), c(p2_o1_ci$upper*100, rev(p2_o1_ci$lower*100)),
        col = adjustcolor("#ca0020",alpha.f=0.2), border=NA)
lines(tval,p2_y0_ci$est*100,col="#0571b0",lwd=3)
lines(tval,p2_y1_ci$est*100,col="#ca0020",lwd=3)
lines(tval,p2_o0_ci$est*100,col="#0571b0",lwd=3,lty=2)
lines(tval,p2_o1_ci$est*100,col="#ca0020",lwd=3,lty=2)
legend("topright",legend=c("youger donor, no IHD","youger donor, IHD",
                           "older donor, no IHD","older donor, IHD"),
       col=c("#0571b0","#ca0020","#0571b0","#ca0020"),lwd=2,bty="n",
       lty=c(1,1,2,2),cex=0.7)

plot(tval,prev_np$`Observed percentages`[,4],type="l",col="dark grey",
     ylim=c(0,100),lwd=3,xlab="years after transplantation",
     ylab="prevalence (%)",main="death")
polygon(c(tval, rev(tval)), c(p3_y0_ci$upper*100, rev(p3_y0_ci$lower*100)),
        col = adjustcolor("#0571b0",alpha.f=0.2), border=NA)
polygon(c(tval, rev(tval)), c(p3_y1_ci$upper*100, rev(p3_y1_ci$lower*100)),
        col = adjustcolor("#ca0020",alpha.f=0.2), border=NA)
polygon(c(tval, rev(tval)), c(p3_o0_ci$upper*100, rev(p3_o0_ci$lower*100)),
        col = adjustcolor("#0571b0",alpha.f=0.2), border=NA)
polygon(c(tval, rev(tval)), c(p3_o1_ci$upper*100, rev(p3_o1_ci$lower*100)),
        col = adjustcolor("#ca0020",alpha.f=0.2), border=NA)
lines(tval,p3_y0_ci$est*100,col="#0571b0",lwd=3)
lines(tval,p3_y1_ci$est*100,col="#ca0020",lwd=3)
lines(tval,p3_o0_ci$est*100,col="#0571b0",lwd=3,lty=2)
lines(tval,p3_o1_ci$est*100,col="#ca0020",lwd=3,lty=2)


## ----message=F,error=F, fig.width=6,fig.align = "center",fig.height=4---------
tval <- seq(0.01,20,length=50)
Sy0 <- overall_survival_ci_band(tval,mlo$opt$par,gg,c(-1,0),mlo$hess)
Sy1 <- overall_survival_ci_band(tval,mlo$opt$par,gg,c(-1,1),mlo$hess)
So0 <- overall_survival_ci_band(tval,mlo$opt$par,gg,c(1,0),mlo$hess)
So1 <- overall_survival_ci_band(tval,mlo$opt$par,gg,c(1,1),mlo$hess)

par(mfrow=c(1,1))
par(bty="l")
par(mar=c(4,4,2,2))
par(cex=1)
plot.survfit.msm(cav.msm, col.surv="black",lwd.surv=2,
                 xlab="years after transplantation",
                 ylab="survival",main=" ",legend.pos=c(30,2),col=NULL,lwd=3)
# using msm package for the nonparametric estimate (for now)
polygon(c(tval, rev(tval)), c(Sy0$upper, rev(Sy0$lower)),
        col = adjustcolor("#0571b0",alpha.f=0.2), border=NA)
polygon(c(tval, rev(tval)), c(Sy1$upper, rev(Sy1$lower)),
        col = adjustcolor("#ca0020",alpha.f=0.2), border=NA)
polygon(c(tval, rev(tval)), c(So0$upper, rev(So0$lower)),
        col = adjustcolor("#0571b0",alpha.f=0.2), border=NA)
polygon(c(tval, rev(tval)), c(So1$upper, rev(So1$lower)),
        col = adjustcolor("#ca0020",alpha.f=0.2), border=NA)
lines(tval,Sy0$est,col="#0571b0",lwd=3)
lines(tval,Sy1$est,col="#ca0020",lwd=3)
lines(tval,So0$est,col="#0571b0",lwd=3,lty=2)
lines(tval,So1$est,col="#ca0020",lwd=3,lty=2)
legend("topright",legend=c("youger donor, no IHD","youger donor, IHD",
                           "older donor, no IHD","older donor, IHD"),
       col=c("#0571b0","#ca0020","#0571b0","#ca0020"),lwd=2,bty="n",
       lty=c(1,1,2,2),cex=0.7)


## -----------------------------------------------------------------------------
tval <- seq(10,40,length=50)
vt <- 10
t01_y0_ci <- transition_prob_ci_band("1-2",tval,vt,mlo$opt$par,gg,xval=c(-1,0),
                                     hessian=mlo$hess)
t01_y1_ci <- transition_prob_ci_band("1-2",tval,vt,mlo$opt$par,gg,xval=c(-1,1),
                                     hessian=mlo$hess)
t01_o0_ci <- transition_prob_ci_band("1-2",tval,vt,mlo$opt$par,gg,xval=c(1,0),hessian=mlo$hess)
t01_o1_ci <- transition_prob_ci_band("1-2",tval,vt,mlo$opt$par,gg,xval=c(1,1),hessian=mlo$hess)
t12_y0_ci <- transition_prob_ci_band("2-3",tval,vt,mlo$opt$par,gg,xval=c(-1,0),hessian=mlo$hess)
t12_y1_ci <- transition_prob_ci_band("2-3",tval,vt,mlo$opt$par,gg,xval=c(-1,1),hessian=mlo$hess)
t12_o0_ci <- transition_prob_ci_band("2-3",tval,vt,mlo$opt$par,gg,xval=c(1,0),hessian=mlo$hess)
t12_o1_ci <- transition_prob_ci_band("2-3",tval,vt,mlo$opt$par,gg,xval=c(1,1),hessian=mlo$hess)
t23_y0_ci <- transition_prob_ci_band("3-4",tval,vt,mlo$opt$par,gg,xval=c(-1,0),hessian=mlo$hess)
t23_y1_ci <- transition_prob_ci_band("3-4",tval,vt,mlo$opt$par,gg,xval=c(-1,1),hessian=mlo$hess)
t23_o0_ci <- transition_prob_ci_band("3-4",tval,vt,mlo$opt$par,gg,xval=c(1,0),hessian=mlo$hess)
t23_o1_ci <- transition_prob_ci_band("3-4",tval,vt,mlo$opt$par,gg,xval=c(1,1),hessian=mlo$hess)
t13_y0_ci <- transition_prob_ci_band("2-4",tval,vt,mlo$opt$par,gg,xval=c(-1,0),hessian=mlo$hess)
t13_y1_ci <- transition_prob_ci_band("2-4",tval,vt,mlo$opt$par,gg,xval=c(-1,1),hessian=mlo$hess)
t13_o0_ci <- transition_prob_ci_band("2-4",tval,vt,mlo$opt$par,gg,xval=c(1,0),hessian=mlo$hess)
t13_o1_ci <- transition_prob_ci_band("2-4",tval,vt,mlo$opt$par,gg,xval=c(1,1),hessian=mlo$hess)
t03_y0_ci <- transition_prob_ci_band("1-4",tval,vt,mlo$opt$par,gg,xval=c(-1,0),hessian=mlo$hess)
t03_y1_ci <- transition_prob_ci_band("1-4",tval,vt,mlo$opt$par,gg,xval=c(-1,1),hessian=mlo$hess)
t03_o0_ci <- transition_prob_ci_band("1-4",tval,vt,mlo$opt$par,gg,xval=c(1,0),hessian=mlo$hess)
t03_o1_ci <- transition_prob_ci_band("1-4",tval,vt,mlo$opt$par,gg,xval=c(1,1),hessian=mlo$hess)

## ----fig.width=8,fig.align = "center",fig.height=5----------------------------
par(mfrow=c(2,3))
par(bty="l")
par(mar=c(4,4,2,2))
par(cex=1)
plot(tval,t01_y0_ci$est,type="l",ylim=c(0,1),col="#0571b0",lwd=3,
     ylab="transition probability",xlab="time",main="well-mild")
polygon(c(tval, rev(tval)), c(t01_y0_ci$upper, rev(t01_y0_ci$lower)),
        col = adjustcolor("#0571b0",alpha.f=0.2), border=NA)
polygon(c(tval, rev(tval)), c(t01_y1_ci$upper, rev(t01_y1_ci$lower)),
        col = adjustcolor("#ca0020",alpha.f=0.2), border=NA)
polygon(c(tval, rev(tval)), c(t01_o0_ci$upper, rev(t01_o0_ci$lower)),
        col = adjustcolor("#0571b0",alpha.f=0.2), border=NA)
polygon(c(tval, rev(tval)), c(t01_o1_ci$upper, rev(t01_o1_ci$lower)),
        col = adjustcolor("#ca0020",alpha.f=0.2), border=NA)
lines(tval,t01_y1_ci$est,lwd=3,col="#ca0020")
lines(tval,t01_o0_ci$est,lwd=3,col="#0571b0",lty=2)
lines(tval,t01_o1_ci$est,lwd=3,col="#ca0020",lty=2)


plot(tval,t12_y0_ci$est,type="l",ylim=c(0,1),col="#0571b0",lwd=3,
     ylab="transition probability",xlab="time",main="mild-severe")
polygon(c(tval, rev(tval)), c(t12_y0_ci$upper, rev(t12_y0_ci$lower)),
        col = adjustcolor("#0571b0",alpha.f=0.2), border=NA)
polygon(c(tval, rev(tval)), c(t12_y1_ci$upper, rev(t12_y1_ci$lower)),
        col = adjustcolor("#ca0020",alpha.f=0.2), border=NA)
polygon(c(tval, rev(tval)), c(t12_o0_ci$upper, rev(t12_o0_ci$lower)),
        col = adjustcolor("#0571b0",alpha.f=0.2), border=NA)
polygon(c(tval, rev(tval)), c(t12_o1_ci$upper, rev(t12_o1_ci$lower)),
        col = adjustcolor("#ca0020",alpha.f=0.2), border=NA)
lines(tval,t12_y1_ci$est,lwd=3,col="#ca0020")
lines(tval,t12_o0_ci$est,lwd=3,col="#0571b0",lty=2)
lines(tval,t12_o1_ci$est,lwd=3,col="#ca0020",lty=2)


plot(tval,t23_y0_ci$est,type="l",ylim=c(0,1),col="#0571b0",lwd=3,
     ylab="transition probability",xlab="time",main="severe-death")
polygon(c(tval, rev(tval)), c(t23_y0_ci$upper, rev(t23_y0_ci$lower)),
        col = adjustcolor("#0571b0",alpha.f=0.2), border=NA)
polygon(c(tval, rev(tval)), c(t23_y1_ci$upper, rev(t23_y1_ci$lower)),
        col = adjustcolor("#ca0020",alpha.f=0.2), border=NA)
polygon(c(tval, rev(tval)), c(t23_o0_ci$upper, rev(t23_o0_ci$lower)),
        col = adjustcolor("#0571b0",alpha.f=0.2), border=NA)
polygon(c(tval, rev(tval)), c(t23_o1_ci$upper, rev(t23_o1_ci$lower)),
        col = adjustcolor("#ca0020",alpha.f=0.2), border=NA)
lines(tval,t23_y1_ci$est,lwd=3,col="#ca0020")
lines(tval,t23_o0_ci$est,lwd=3,col="#0571b0",lty=2)
lines(tval,t23_o1_ci$est,lwd=3,col="#ca0020",lty=2)

plot(tval,t13_y0_ci$est,type="l",ylim=c(0,1),col="#0571b0",lwd=3,
     ylab="transition probability",xlab="time",main="mild-death")
polygon(c(tval, rev(tval)), c(t13_y0_ci$upper, rev(t13_y0_ci$lower)),
        col = adjustcolor("#0571b0",alpha.f=0.2), border=NA)
polygon(c(tval, rev(tval)), c(t13_y1_ci$upper, rev(t13_y1_ci$lower)),
        col = adjustcolor("#ca0020",alpha.f=0.2), border=NA)
polygon(c(tval, rev(tval)), c(t13_o0_ci$upper, rev(t13_o0_ci$lower)),
        col = adjustcolor("#0571b0",alpha.f=0.2), border=NA)
polygon(c(tval, rev(tval)), c(t13_o1_ci$upper, rev(t13_o1_ci$lower)),
        col = adjustcolor("#ca0020",alpha.f=0.2), border=NA)
lines(tval,t13_y1_ci$est,lwd=3,col="#ca0020")
lines(tval,t13_o0_ci$est,lwd=3,col="#0571b0",lty=2)
lines(tval,t13_o1_ci$est,lwd=3,col="#ca0020",lty=2)

plot(tval,t03_y0_ci$est,type="l",ylim=c(0,1),col="#0571b0",lwd=3,
     ylab="transition probability",xlab="time",main="well-death")
polygon(c(tval, rev(tval)), c(t03_y0_ci$upper, rev(t03_y0_ci$lower)),
        col = adjustcolor("#0571b0",alpha.f=0.2), border=NA)
polygon(c(tval, rev(tval)), c(t03_y1_ci$upper, rev(t03_y1_ci$lower)),
        col = adjustcolor("#ca0020",alpha.f=0.2), border=NA)
polygon(c(tval, rev(tval)), c(t03_o0_ci$upper, rev(t03_o0_ci$lower)),
        col = adjustcolor("#0571b0",alpha.f=0.2), border=NA)
polygon(c(tval, rev(tval)), c(t03_o1_ci$upper, rev(t03_o1_ci$lower)),
        col = adjustcolor("#ca0020",alpha.f=0.2), border=NA)
lines(tval,t03_y1_ci$est,lwd=3,col="#ca0020")
lines(tval,t03_o0_ci$est,lwd=3,col="#0571b0",lty=2)
lines(tval,t03_o1_ci$est,lwd=3,col="#ca0020",lty=2)

plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = ''
     , xlab = '',ylim=c(0,1),xlim=c(0,40))
legend("right",legend=c("youger donor, no IHD","youger donor, IHD",
                        "older donor, no IHD","older donor, IHD"),
       col=c("#0571b0","#ca0020","#0571b0","#ca0020"),lwd=2,
       bty="n",lty=c(1,1,2,2),cex=0.7)
p2_y1_ci <- occupancy_prob_ci_band("3",tval,mlo$opt$par,gg,xval=c(-1,1),mlo$hess)
p2_o0_ci <- occupancy_prob_ci_band("3",tval,mlo$opt$par,gg,xval=c(1,0),mlo$hess)
p2_o1_ci <- occupancy_prob_ci_band("3",tval,mlo$opt$par,gg,xval=c(1,1),mlo$hess)
p3_y0_ci <- occupancy_prob_ci_band("4",tval,mlo$opt$par,gg,xval=c(-1,0),mlo$hess)
p3_y1_ci <- occupancy_prob_ci_band("4",tval,mlo$opt$par,gg,xval=c(-1,1),mlo$hess)
p3_o0_ci <- occupancy_prob_ci_band("4",tval,mlo$opt$par,gg,xval=c(1,0),mlo$hess)
p3_o1_ci <- occupancy_prob_ci_band("4",tval,mlo$opt$par,gg,xval=c(1,1),mlo$hess)

## ----fig.width=8,fig.align = "center",fig.height=5----------------------------
# msm package  (to get non-parameteric prevalence curves)
oneway4.q <- rbind(c(0, 0.25, 0, 0.25), c(0, 0, 0.25, 0.25),c(0, 0, 0, 0.5),
                   c(0, 0, 0, 0))
rownames(oneway4.q) <- colnames(oneway4.q) <- c("Well", "Mild","Severe",
                                                "Death")
cav.msm <- msm(state ~ time, subject = patient, data = dd,qmatrix = oneway4.q,
               death = 4,method = "BFGS")

prev_np <- prevalence.msm(cav.msm,times=tval)

# Plot
par(mfrow=c(2,2))
par(bty="l")
par(mar=c(4,4,1,2))
par(cex=1)
plot(tval,prev_np$`Observed percentages`[,1],type="l",ylim=c(0,100),lwd=3,
     xlab=" ",col="dark grey",ylab="prevalence (%)",main="well")
polygon(c(tval, rev(tval)), c(p0_y0_ci$upper*100, rev(p0_y0_ci$lower*100)),
        col = adjustcolor("#0571b0",alpha.f=0.2), border=NA)
polygon(c(tval, rev(tval)), c(p0_y1_ci$upper*100, rev(p0_y1_ci$lower*100)),
        col = adjustcolor("#ca0020",alpha.f=0.2), border=NA)
polygon(c(tval, rev(tval)), c(p0_o0_ci$upper*100, rev(p0_o0_ci$lower*100)),
        col = adjustcolor("#0571b0",alpha.f=0.2), border=NA)
polygon(c(tval, rev(tval)), c(p0_o1_ci$upper*100, rev(p0_o1_ci$lower*100)),
        col = adjustcolor("#ca0020",alpha.f=0.2), border=NA)
lines(tval,p0_y0_ci$est*100,col="#0571b0",lwd=3)
lines(tval,p0_y1_ci$est*100,col="#ca0020",lwd=3)
lines(tval,p0_o0_ci$est*100,col="#0571b0",lwd=3,lty=2)
lines(tval,p0_o1_ci$est*100,col="#ca0020",lwd=3,lty=2)

plot(tval,prev_np$`Observed percentages`[,2],type="l",col="dark grey",
     ylim=c(0,100),lwd=3,xlab=" ",ylab="prevalence (%)",main="mild")
polygon(c(tval, rev(tval)), c(p1_y0_ci$upper*100, rev(p1_y0_ci$lower*100)),
        col = adjustcolor("#0571b0",alpha.f=0.2), border=NA)
polygon(c(tval, rev(tval)), c(p1_y1_ci$upper*100, rev(p1_y1_ci$lower*100)),
        col = adjustcolor("#ca0020",alpha.f=0.2), border=NA)
polygon(c(tval, rev(tval)), c(p1_o0_ci$upper*100, rev(p1_o0_ci$lower*100)),
        col = adjustcolor("#0571b0",alpha.f=0.2), border=NA)
polygon(c(tval, rev(tval)), c(p1_o1_ci$upper*100, rev(p1_o1_ci$lower*100)),
        col = adjustcolor("#ca0020",alpha.f=0.2), border=NA)
lines(tval,p1_y0_ci$est*100,col="#0571b0",lwd=3)
lines(tval,p1_y1_ci$est*100,col="#ca0020",lwd=3)
lines(tval,p1_o0_ci$est*100,col="#0571b0",lwd=3,lty=2)
lines(tval,p1_o1_ci$est*100,col="#ca0020",lwd=3,lty=2)


plot(tval,prev_np$`Observed percentages`[,3],type="l",col="dark grey",
     ylim=c(0,100),lwd=3,xlab="years after transplantation",
     ylab="prevalence (%)",main="severe")
polygon(c(tval, rev(tval)), c(p2_y0_ci$upper*100, rev(p2_y0_ci$lower*100)),
        col = adjustcolor("#0571b0",alpha.f=0.2), border=NA)
polygon(c(tval, rev(tval)), c(p2_y1_ci$upper*100, rev(p2_y1_ci$lower*100)),
        col = adjustcolor("#ca0020",alpha.f=0.2), border=NA)
polygon(c(tval, rev(tval)), c(p2_o0_ci$upper*100, rev(p2_o0_ci$lower*100)),
        col = adjustcolor("#0571b0",alpha.f=0.2), border=NA)
polygon(c(tval, rev(tval)), c(p2_o1_ci$upper*100, rev(p2_o1_ci$lower*100)),
        col = adjustcolor("#ca0020",alpha.f=0.2), border=NA)
lines(tval,p2_y0_ci$est*100,col="#0571b0",lwd=3)
lines(tval,p2_y1_ci$est*100,col="#ca0020",lwd=3)
lines(tval,p2_o0_ci$est*100,col="#0571b0",lwd=3,lty=2)
lines(tval,p2_o1_ci$est*100,col="#ca0020",lwd=3,lty=2)
legend("topright",legend=c("youger donor, no IHD","youger donor, IHD",
                           "older donor, no IHD","older donor, IHD"),
       col=c("#0571b0","#ca0020","#0571b0","#ca0020"),lwd=2,bty="n",
       lty=c(1,1,2,2),cex=0.7)

plot(tval,prev_np$`Observed percentages`[,4],type="l",col="dark grey",
     ylim=c(0,100),lwd=3,xlab="years after transplantation",
     ylab="prevalence (%)",main="death")
polygon(c(tval, rev(tval)), c(p3_y0_ci$upper*100, rev(p3_y0_ci$lower*100)),
        col = adjustcolor("#0571b0",alpha.f=0.2), border=NA)
polygon(c(tval, rev(tval)), c(p3_y1_ci$upper*100, rev(p3_y1_ci$lower*100)),
        col = adjustcolor("#ca0020",alpha.f=0.2), border=NA)
polygon(c(tval, rev(tval)), c(p3_o0_ci$upper*100, rev(p3_o0_ci$lower*100)),
        col = adjustcolor("#0571b0",alpha.f=0.2), border=NA)
polygon(c(tval, rev(tval)), c(p3_o1_ci$upper*100, rev(p3_o1_ci$lower*100)),
        col = adjustcolor("#ca0020",alpha.f=0.2), border=NA)
lines(tval,p3_y0_ci$est*100,col="#0571b0",lwd=3)
lines(tval,p3_y1_ci$est*100,col="#ca0020",lwd=3)
lines(tval,p3_o0_ci$est*100,col="#0571b0",lwd=3,lty=2)
lines(tval,p3_o1_ci$est*100,col="#ca0020",lwd=3,lty=2)


## ----message=F,error=F, fig.width=6,fig.align = "center",fig.height=4---------
tval <- seq(0.01,20,length=50)
Sy0 <- overall_survival_ci_band(tval,mlo$opt$par,gg,c(-1,0),mlo$hess)
Sy1 <- overall_survival_ci_band(tval,mlo$opt$par,gg,c(-1,1),mlo$hess)
So0 <- overall_survival_ci_band(tval,mlo$opt$par,gg,c(1,0),mlo$hess)
So1 <- overall_survival_ci_band(tval,mlo$opt$par,gg,c(1,1),mlo$hess)

par(mfrow=c(1,1))
par(bty="l")
par(mar=c(4,4,2,2))
par(cex=1)
plot.survfit.msm(cav.msm, col.surv="black",lwd.surv=2,
                 xlab="years after transplantation",
                 ylab="survival",main=" ",legend.pos=c(30,2),col=NULL,lwd=3)
# using msm package for the nonparametric estimate (for now)
polygon(c(tval, rev(tval)), c(Sy0$upper, rev(Sy0$lower)),
        col = adjustcolor("#0571b0",alpha.f=0.2), border=NA)
polygon(c(tval, rev(tval)), c(Sy1$upper, rev(Sy1$lower)),
        col = adjustcolor("#ca0020",alpha.f=0.2), border=NA)
polygon(c(tval, rev(tval)), c(So0$upper, rev(So0$lower)),
        col = adjustcolor("#0571b0",alpha.f=0.2), border=NA)
polygon(c(tval, rev(tval)), c(So1$upper, rev(So1$lower)),
        col = adjustcolor("#ca0020",alpha.f=0.2), border=NA)
lines(tval,Sy0$est,col="#0571b0",lwd=3)
lines(tval,Sy1$est,col="#ca0020",lwd=3)
lines(tval,So0$est,col="#0571b0",lwd=3,lty=2)
lines(tval,So1$est,col="#ca0020",lwd=3,lty=2)
legend("topright",legend=c("youger donor, no IHD","youger donor, IHD",
                           "older donor, no IHD","older donor, IHD"),
       col=c("#0571b0","#ca0020","#0571b0","#ca0020"),lwd=2,bty="n",
       lty=c(1,1,2,2),cex=0.7)


## -----------------------------------------------------------------------------
tval <- seq(10,40,length=50)
vt <- 10
t01_y0_ci <- transition_prob_ci_band("1-2",tval,vt,mlo$opt$par,gg,xval=c(-1,0),
                                     hessian=mlo$hess)
t01_y1_ci <- transition_prob_ci_band("1-2",tval,vt,mlo$opt$par,gg,xval=c(-1,1),
                                     hessian=mlo$hess)
t01_o0_ci <- transition_prob_ci_band("1-2",tval,vt,mlo$opt$par,gg,xval=c(1,0),hessian=mlo$hess)
t01_o1_ci <- transition_prob_ci_band("1-2",tval,vt,mlo$opt$par,gg,xval=c(1,1),hessian=mlo$hess)
t12_y0_ci <- transition_prob_ci_band("2-3",tval,vt,mlo$opt$par,gg,xval=c(-1,0),hessian=mlo$hess)
t12_y1_ci <- transition_prob_ci_band("2-3",tval,vt,mlo$opt$par,gg,xval=c(-1,1),hessian=mlo$hess)
t12_o0_ci <- transition_prob_ci_band("2-3",tval,vt,mlo$opt$par,gg,xval=c(1,0),hessian=mlo$hess)
t12_o1_ci <- transition_prob_ci_band("2-3",tval,vt,mlo$opt$par,gg,xval=c(1,1),hessian=mlo$hess)
t23_y0_ci <- transition_prob_ci_band("3-4",tval,vt,mlo$opt$par,gg,xval=c(-1,0),hessian=mlo$hess)
t23_y1_ci <- transition_prob_ci_band("3-4",tval,vt,mlo$opt$par,gg,xval=c(-1,1),hessian=mlo$hess)
t23_o0_ci <- transition_prob_ci_band("3-4",tval,vt,mlo$opt$par,gg,xval=c(1,0),hessian=mlo$hess)
t23_o1_ci <- transition_prob_ci_band("3-4",tval,vt,mlo$opt$par,gg,xval=c(1,1),hessian=mlo$hess)
t13_y0_ci <- transition_prob_ci_band("2-4",tval,vt,mlo$opt$par,gg,xval=c(-1,0),hessian=mlo$hess)
t13_y1_ci <- transition_prob_ci_band("2-4",tval,vt,mlo$opt$par,gg,xval=c(-1,1),hessian=mlo$hess)
t13_o0_ci <- transition_prob_ci_band("2-4",tval,vt,mlo$opt$par,gg,xval=c(1,0),hessian=mlo$hess)
t13_o1_ci <- transition_prob_ci_band("2-4",tval,vt,mlo$opt$par,gg,xval=c(1,1),hessian=mlo$hess)
t03_y0_ci <- transition_prob_ci_band("1-4",tval,vt,mlo$opt$par,gg,xval=c(-1,0),hessian=mlo$hess)
t03_y1_ci <- transition_prob_ci_band("1-4",tval,vt,mlo$opt$par,gg,xval=c(-1,1),hessian=mlo$hess)
t03_o0_ci <- transition_prob_ci_band("1-4",tval,vt,mlo$opt$par,gg,xval=c(1,0),hessian=mlo$hess)
t03_o1_ci <- transition_prob_ci_band("1-4",tval,vt,mlo$opt$par,gg,xval=c(1,1),hessian=mlo$hess)

## ----fig.width=8,fig.align = "center",fig.height=5----------------------------
par(mfrow=c(2,3))
par(bty="l")
par(mar=c(4,4,2,2))
par(cex=1)
plot(tval,t01_y0_ci$est,type="l",ylim=c(0,1),col="#0571b0",lwd=3,
     ylab="transition probability",xlab="time",main="well-mild")
polygon(c(tval, rev(tval)), c(t01_y0_ci$upper, rev(t01_y0_ci$lower)),
        col = adjustcolor("#0571b0",alpha.f=0.2), border=NA)
polygon(c(tval, rev(tval)), c(t01_y1_ci$upper, rev(t01_y1_ci$lower)),
        col = adjustcolor("#ca0020",alpha.f=0.2), border=NA)
polygon(c(tval, rev(tval)), c(t01_o0_ci$upper, rev(t01_o0_ci$lower)),
        col = adjustcolor("#0571b0",alpha.f=0.2), border=NA)
polygon(c(tval, rev(tval)), c(t01_o1_ci$upper, rev(t01_o1_ci$lower)),
        col = adjustcolor("#ca0020",alpha.f=0.2), border=NA)
lines(tval,t01_y1_ci$est,lwd=3,col="#ca0020")
lines(tval,t01_o0_ci$est,lwd=3,col="#0571b0",lty=2)
lines(tval,t01_o1_ci$est,lwd=3,col="#ca0020",lty=2)


plot(tval,t12_y0_ci$est,type="l",ylim=c(0,1),col="#0571b0",lwd=3,
     ylab="transition probability",xlab="time",main="mild-severe")
polygon(c(tval, rev(tval)), c(t12_y0_ci$upper, rev(t12_y0_ci$lower)),
        col = adjustcolor("#0571b0",alpha.f=0.2), border=NA)
polygon(c(tval, rev(tval)), c(t12_y1_ci$upper, rev(t12_y1_ci$lower)),
        col = adjustcolor("#ca0020",alpha.f=0.2), border=NA)
polygon(c(tval, rev(tval)), c(t12_o0_ci$upper, rev(t12_o0_ci$lower)),
        col = adjustcolor("#0571b0",alpha.f=0.2), border=NA)
polygon(c(tval, rev(tval)), c(t12_o1_ci$upper, rev(t12_o1_ci$lower)),
        col = adjustcolor("#ca0020",alpha.f=0.2), border=NA)
lines(tval,t12_y1_ci$est,lwd=3,col="#ca0020")
lines(tval,t12_o0_ci$est,lwd=3,col="#0571b0",lty=2)
lines(tval,t12_o1_ci$est,lwd=3,col="#ca0020",lty=2)


plot(tval,t23_y0_ci$est,type="l",ylim=c(0,1),col="#0571b0",lwd=3,
     ylab="transition probability",xlab="time",main="severe-death")
polygon(c(tval, rev(tval)), c(t23_y0_ci$upper, rev(t23_y0_ci$lower)),
        col = adjustcolor("#0571b0",alpha.f=0.2), border=NA)
polygon(c(tval, rev(tval)), c(t23_y1_ci$upper, rev(t23_y1_ci$lower)),
        col = adjustcolor("#ca0020",alpha.f=0.2), border=NA)
polygon(c(tval, rev(tval)), c(t23_o0_ci$upper, rev(t23_o0_ci$lower)),
        col = adjustcolor("#0571b0",alpha.f=0.2), border=NA)
polygon(c(tval, rev(tval)), c(t23_o1_ci$upper, rev(t23_o1_ci$lower)),
        col = adjustcolor("#ca0020",alpha.f=0.2), border=NA)
lines(tval,t23_y1_ci$est,lwd=3,col="#ca0020")
lines(tval,t23_o0_ci$est,lwd=3,col="#0571b0",lty=2)
lines(tval,t23_o1_ci$est,lwd=3,col="#ca0020",lty=2)

plot(tval,t13_y0_ci$est,type="l",ylim=c(0,1),col="#0571b0",lwd=3,
     ylab="transition probability",xlab="time",main="mild-death")
polygon(c(tval, rev(tval)), c(t13_y0_ci$upper, rev(t13_y0_ci$lower)),
        col = adjustcolor("#0571b0",alpha.f=0.2), border=NA)
polygon(c(tval, rev(tval)), c(t13_y1_ci$upper, rev(t13_y1_ci$lower)),
        col = adjustcolor("#ca0020",alpha.f=0.2), border=NA)
polygon(c(tval, rev(tval)), c(t13_o0_ci$upper, rev(t13_o0_ci$lower)),
        col = adjustcolor("#0571b0",alpha.f=0.2), border=NA)
polygon(c(tval, rev(tval)), c(t13_o1_ci$upper, rev(t13_o1_ci$lower)),
        col = adjustcolor("#ca0020",alpha.f=0.2), border=NA)
lines(tval,t13_y1_ci$est,lwd=3,col="#ca0020")
lines(tval,t13_o0_ci$est,lwd=3,col="#0571b0",lty=2)
lines(tval,t13_o1_ci$est,lwd=3,col="#ca0020",lty=2)

plot(tval,t03_y0_ci$est,type="l",ylim=c(0,1),col="#0571b0",lwd=3,
     ylab="transition probability",xlab="time",main="well-death")
polygon(c(tval, rev(tval)), c(t03_y0_ci$upper, rev(t03_y0_ci$lower)),
        col = adjustcolor("#0571b0",alpha.f=0.2), border=NA)
polygon(c(tval, rev(tval)), c(t03_y1_ci$upper, rev(t03_y1_ci$lower)),
        col = adjustcolor("#ca0020",alpha.f=0.2), border=NA)
polygon(c(tval, rev(tval)), c(t03_o0_ci$upper, rev(t03_o0_ci$lower)),
        col = adjustcolor("#0571b0",alpha.f=0.2), border=NA)
polygon(c(tval, rev(tval)), c(t03_o1_ci$upper, rev(t03_o1_ci$lower)),
        col = adjustcolor("#ca0020",alpha.f=0.2), border=NA)
lines(tval,t03_y1_ci$est,lwd=3,col="#ca0020")
lines(tval,t03_o0_ci$est,lwd=3,col="#0571b0",lty=2)
lines(tval,t03_o1_ci$est,lwd=3,col="#ca0020",lty=2)

plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = ''
     , xlab = '',ylim=c(0,1),xlim=c(0,40))
legend("right",legend=c("youger donor, no IHD","youger donor, IHD",
                        "older donor, no IHD","older donor, IHD"),
       col=c("#0571b0","#ca0020","#0571b0","#ca0020"),lwd=2,
       bty="n",lty=c(1,1,2,2),cex=0.7)
