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
ddo = dd

# Change state names from 1,2,3,4 to well, mild, severe, death
tab = data.frame(state=1:4,name=c("well","mild","severe","death"))
dd$state = tab$name[match(dd$state,tab$state)]

print(dd[1:11,])

## ----fig.height=3,fig.align = "center"----------------------------------------
# Specify the graph:
gg = graph_from_literal("well"--+"mild"--+"severe"--+"death", "well"--+"death",
                        "mild"--+"death")
par(mar=c(1,1,1,1))
plot(gg,layout=layout_with_sugiyama(gg,layers=c(1,1,1,2))$layout,vertex.size=45)

## -----------------------------------------------------------------------------
f_01 = function(param, x, tt){dexp(tt,exp(param[1]))}
f_12 = function(param, x, tt){dexp(tt,exp(param[2]))}
f_23 = function(param, x, tt){dexp(tt,exp(param[3]))}
f_03 = function(param, x, tt){dexp(tt,exp(param[4]))}
f_13 = function(param, x, tt){dexp(tt,exp(param[5]))}

S_01 = function(param, x, tt){1-pexp(tt,exp(param[1]))}
S_12 = function(param, x, tt){1-pexp(tt,exp(param[2]))}
S_23 = function(param, x, tt){1-pexp(tt,exp(param[3]))}
S_03 = function(param, x, tt){1-pexp(tt,exp(param[4]))}
S_13 = function(param, x, tt){1-pexp(tt,exp(param[5]))}

## -----------------------------------------------------------------------------
print(names_of_survival_density(gg))

## ----eval=F-------------------------------------------------------------------
#  startval <- c(-2.5,-1.1,-1.2,-3.1,-2.8)
#
#  mlo <- smms(startval,dd,gg, mc_cores = 1, hessian_matrix = T)

## ----echo=F-------------------------------------------------------------------
load("~/UiO/Projects/smms/vignettes/cav_expo_noCov_optims")

## -----------------------------------------------------------------------------
# Compute AIC (higher values are better with this definition)
aic <- (-2*mlo$opt$objective)-2*length(mlo$opt$par) #-2887.1

# Look at estimates and 95% confidence intervals.
# On the -Inf to Inf scale:
print(round(est_ci(mlo$opt$par,mlo$hess),2))

# On the 0 to Inf scale (on the transition intensity scale):
round(exp(est_ci(mlo$opt$par,mlo$hess)),2)

## -----------------------------------------------------------------------------
tval <- seq(0.01,30,length=50)
# a sequence of time-points over which to compute the state occupancies
p0_ci <- occupancy_prob_ci_band("well",tval,mlo$opt$par,gg,hessian=mlo$hess)
#for computing the confidence bands, the hessian needs to be provided (but there
# exists a function computing only the occupancy probabilities too)
p1_ci <- occupancy_prob_ci_band("mild",tval,mlo$opt$par,gg,hessian=mlo$hess)
p2_ci <- occupancy_prob_ci_band("severe",tval,mlo$opt$par,gg,hessian=mlo$hess)
p3_ci <- occupancy_prob_ci_band("death",tval,mlo$opt$par,gg,hessian=mlo$hess)

## ----fig.width=6,fig.align = "center",fig.height=4----------------------------
# msm package (to get non-parameteric prevalence curves)
oneway4.q <- rbind(c(0, 0.25, 0, 0.25), c(0, 0, 0.25, 0.25),c(0, 0, 0, 0.5),
                   c(0, 0, 0, 0))
rownames(oneway4.q) <- colnames(oneway4.q) <- c("Well", "Mild","Severe",
                                                "Death")
cav.msm <- msm(state ~ time, subject = patient, data = ddo,qmatrix = oneway4.q,
               death = 4,method = "BFGS")

prev_np <- prevalence.msm(cav.msm,times=tval)

# Plot
par(bty="l")
par(mar=c(4,4,1,2))
par(cex=1)
plot(tval,prev_np$`Observed percentages`[,1],type="l",ylim=c(0,100),lwd=3,
     xlab="time",col="dark grey", ylab="prevalence (%)")
polygon(c(tval, rev(tval)), c(p0_ci$upper*100, rev(p0_ci$lower*100)),
        col = adjustcolor("#b2e2e2",alpha.f=0.2), border=NA)
lines(tval,p0_ci$est*100,col="#b2e2e2",lwd=3)

lines(tval,prev_np$`Observed percentages`[,2],col="dark grey",lwd=3)
polygon(c(tval, rev(tval)), c(p1_ci$upper*100, rev(p1_ci$lower*100)),
        col = adjustcolor("#66c2a4",alpha.f=0.2), border=NA)
lines(tval,p1_ci$est*100,col="#66c2a4",lwd=3)


lines(tval,prev_np$`Observed percentages`[,3],col="dark grey",lwd=3)
polygon(c(tval, rev(tval)), c(p2_ci$upper*100, rev(p2_ci$lower*100)),
        col = adjustcolor("#2ca25f",alpha.f=0.2), border=NA)
lines(tval,p2_ci$est*100,col="#2ca25f",lwd=3)

lines(tval,prev_np$`Observed percentages`[,4],col="dark grey",lwd=3)
polygon(c(tval, rev(tval)), c(p3_ci$upper*100, rev(p3_ci$lower*100)),
        col = adjustcolor("#006d2c",alpha.f=0.2), border=NA)
lines(tval,p3_ci$est*100,col="#006d2c",lwd=3)

legend("right",legend=c("non-parametric","well","mild","severe","death"),
       col=c("dark grey","#b2e2e2","#66c2a4","#2ca25f","#006d2c"),
       lwd=2,bty="n",cex=1)


## ----message=F,error=F, fig.width=6,fig.align = "center",fig.height=4---------
tval <- seq(0.01,30,length=50)
So <- overall_survival_ci_band(tval,mlo$opt$par,gg,hessian=mlo$hess)

par(bty="l")
par(mar=c(4,4,2,2))
par(cex=1)
plot.survfit.msm(cav.msm, col.surv="dark grey",lwd.surv=2,col=NULL,lwd=3,
                 xlab="years after transplantation",ylab="survival",main=" ",
                 legend.pos=c(30,2))
# using msm package for the nonparametric estimate (for now)
polygon(c(tval, rev(tval)), c(So$upper, rev(So$lower)),
        col = adjustcolor("#0571b0",alpha.f=0.2), border=NA)
lines(tval,So$est,col="#0571b0",lwd=3)

## -----------------------------------------------------------------------------
tval <- seq(10,40,length=50)
vt <- 10
t01_ci <- transition_prob_ci_band("well-mild",tval,vt,mlo$opt$par,gg,
                                  hessian=mlo$hess)
t12_ci <- transition_prob_ci_band("mild-severe",tval,vt,mlo$opt$par,gg,
                                  hessian=mlo$hess)
t23_ci <- transition_prob_ci_band("severe-death",tval,vt,mlo$opt$par,gg,
                                  hessian=mlo$hess)
t13_ci <- transition_prob_ci_band("mild-death",tval,vt,mlo$opt$par,gg,
                                  hessian=mlo$hess)
t03_ci <- transition_prob_ci_band("well-death",tval,vt,mlo$opt$par,gg,
                                  hessian=mlo$hess)

## ----fig.width=6,fig.align = "center",fig.height=4----------------------------
par(bty="l")
par(mar=c(4,4,2,2))
par(cex=1)
plot(tval,t01_ci$est,type="l",ylim=c(0,1),col="#5ab4ac",lwd=3,
     ylab="transition probability",xlab="time")
polygon(c(tval, rev(tval)), c(t01_ci$upper, rev(t01_ci$lower)),
        col = adjustcolor("#0571b0",alpha.f=0.2), border=NA)
polygon(c(tval, rev(tval)), c(t12_ci$upper, rev(t12_ci$lower)),
        col = adjustcolor("#5ab4ac",alpha.f=0.2), border=NA)
polygon(c(tval, rev(tval)), c(t23_ci$upper, rev(t23_ci$lower)),
        col = adjustcolor("#f6e8c3",alpha.f=0.2), border=NA)
polygon(c(tval, rev(tval)), c(t13_ci$upper, rev(t13_ci$lower)),
        col = adjustcolor("#d8b365",alpha.f=0.2), border=NA)
polygon(c(tval, rev(tval)), c(t03_ci$upper, rev(t03_ci$lower)),
        col = adjustcolor("#8c510a",alpha.f=0.2), border=NA)
lines(tval,t12_ci$est,lwd=3,col="#5ab4ac")
lines(tval,t23_ci$est,lwd=3,col="#f6e8c3")
lines(tval,t13_ci$est,lwd=3,col="#d8b365")
lines(tval,t03_ci$est,lwd=3,col="#8c510a")
legend(x=32,y=0.8,legend=c("well-mild","mild-severe","severe-death",
                           "mild-death","well-death"),
       col=c("#5ab4ac","#01665e","#f6e8c3","#d8b365","#8c510a"),
       lwd=2,bty="n",cex=0.8)


## ----eval=F-------------------------------------------------------------------
#  write_loglikelihood(gg,abs_exact = T)
