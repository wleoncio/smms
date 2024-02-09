dd <- msm::cav
dd <- dd[!is.na(dd$pdiag), ]

# Remove observations where the patient appears to go back to a previous state
# (assumed to be impossible):
id_wrong <- unique(dd$PTNUM[which(dd$state != dd$statemax)])
dd <- dd[-which(dd$PTNUM %in% id_wrong), ]

dd <- dd[, -c(2, 5, 7, 9, 10)]
# rename relevant columns (necessary in current version):
colnames(dd)[1:2] <- c("patient", "time")
ddo <- dd

# Change state names from 1,2,3,4 to well, mild, severe, death
tab <- data.frame(state = 1:4, name = c("well", "mild", "severe", "death"))
dd$state <- tab$name[match(dd$state, tab$state)]

print(dd[1:11, ])

## ----fig.height=3,fig.align = "center"----------------------------------------
# Specify the graph:
gg <- graph_from_literal(
  "well" - -+"mild" - -+"severe" - -+"death", "well" - -+"death",
  "mild" - -+"death"
)

## -----------------------------------------------------------------------------
f_01 <- function(param, x, tt) {
  dexp(tt, exp(param[1]))
}
f_12 <- function(param, x, tt) {
  dexp(tt, exp(param[2]))
}
f_23 <- function(param, x, tt) {
  dexp(tt, exp(param[3]))
}
f_03 <- function(param, x, tt) {
  dexp(tt, exp(param[4]))
}
f_13 <- function(param, x, tt) {
  dexp(tt, exp(param[5]))
}

S_01 <- function(param, x, tt) {
  1 - pexp(tt, exp(param[1]))
}
S_12 <- function(param, x, tt) {
  1 - pexp(tt, exp(param[2]))
}
S_23 <- function(param, x, tt) {
  1 - pexp(tt, exp(param[3]))
}
S_03 <- function(param, x, tt) {
  1 - pexp(tt, exp(param[4]))
}
S_13 <- function(param, x, tt) {
  1 - pexp(tt, exp(param[5]))
}

## -----------------------------------------------------------------------------
print(names_of_survival_density(gg))

## ----eval=F-------------------------------------------------------------------
startval <- c(-2.5,-1.1,-1.2,-3.1,-2.8)
mlo <- smms(startval,dd,gg, mc_cores = 1, hessian_matrix = T)

## -----------------------------------------------------------------------------
# Compute AIC (higher values are better with this definition)
aic <- (-2 * mlo$opt$objective) - 2 * length(mlo$opt$par) #-2887.1

# Look at estimates and 95% confidence intervals.
# On the -Inf to Inf scale:
print(round(est_ci(mlo$opt$par, mlo$hess), 2))

# On the 0 to Inf scale (on the transition intensity scale):
round(exp(est_ci(mlo$opt$par, mlo$hess)), 2)

## -----------------------------------------------------------------------------
tval <- seq(0.01, 30, length = 50)
# a sequence of time-points over which to compute the state occupancies
p0_ci <- occupancy_prob_ci_band("well", tval, mlo$opt$par, gg, hessian = mlo$hess)
# for computing the confidence bands, the hessian needs to be provided (but there
# exists a function computing only the occupancy probabilities too)
p1_ci <- occupancy_prob_ci_band("mild", tval, mlo$opt$par, gg, hessian = mlo$hess)
p2_ci <- occupancy_prob_ci_band("severe", tval, mlo$opt$par, gg, hessian = mlo$hess)
p3_ci <- occupancy_prob_ci_band("death", tval, mlo$opt$par, gg, hessian = mlo$hess)

## ----fig.width=6,fig.align = "center",fig.height=4----------------------------
# msm package (to get non-parameteric prevalence curves)
oneway4.q <- rbind(
  c(0, 0.25, 0, 0.25), c(0, 0, 0.25, 0.25), c(0, 0, 0, 0.5),
  c(0, 0, 0, 0)
)
rownames(oneway4.q) <- colnames(oneway4.q) <- c(
  "Well", "Mild", "Severe",
  "Death"
)
cav.msm <- msm(state ~ time,
  subject = patient, data = ddo, qmatrix = oneway4.q,
  death = 4, method = "BFGS"
)

prev_np <- prevalence.msm(cav.msm, times = tval)

## ----message=F,error=F, fig.width=6,fig.align = "center",fig.height=4---------
tval <- seq(0.01, 30, length = 50)
So <- overall_survival_ci_band(tval, mlo$opt$par, gg, hessian = mlo$hess)

## -----------------------------------------------------------------------------
tval <- seq(10, 40, length = 50)
vt <- 10
t01_ci <- transition_prob_ci_band("well-mild", tval, vt, mlo$opt$par, gg,
  hessian = mlo$hess
)
t12_ci <- transition_prob_ci_band("mild-severe", tval, vt, mlo$opt$par, gg,
  hessian = mlo$hess
)
t23_ci <- transition_prob_ci_band("severe-death", tval, vt, mlo$opt$par, gg,
  hessian = mlo$hess
)
t13_ci <- transition_prob_ci_band("mild-death", tval, vt, mlo$opt$par, gg,
  hessian = mlo$hess
)
t03_ci <- transition_prob_ci_band("well-death", tval, vt, mlo$opt$par, gg,
  hessian = mlo$hess
)

## ----eval=F-------------------------------------------------------------------
#  write_loglikelihood(gg,abs_exact = T)
