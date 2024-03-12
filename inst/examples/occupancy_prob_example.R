# Reduced msm::cav dataset
set.seed(114385)
cav_not_na <- msm::cav[!is.na(msm::cav$pdiag), ]
patient_sample <- sample(unique(cav_not_na$PTNUM), 7L)
dd <- cav_not_na[cav_not_na$PTNUM %in% patient_sample, ]

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

# Specify the graph:
gg <- igraph::graph_from_literal(
  "well" - -+"mild" - -+"severe" - -+"death", "well" - -+"death",
  "mild" - -+"death"
)

## -----------------------------------------------------------------------------
f_01 <- function(param, x, tt) dexp(tt, exp(param[1]))
f_12 <- function(param, x, tt) dexp(tt, exp(param[2]))
f_23 <- function(param, x, tt) dexp(tt, exp(param[3]))
f_03 <- function(param, x, tt) dexp(tt, exp(param[4]))
f_13 <- function(param, x, tt) dexp(tt, exp(param[5]))

S_01 <- function(param, x, tt) 1 - pexp(tt, exp(param[1]))
S_12 <- function(param, x, tt) 1 - pexp(tt, exp(param[2]))
S_23 <- function(param, x, tt) 1 - pexp(tt, exp(param[3]))
S_03 <- function(param, x, tt) 1 - pexp(tt, exp(param[4]))
S_13 <- function(param, x, tt) 1 - pexp(tt, exp(param[5]))

## -----------------------------------------------------------------------------
startval <- c(-2.5,-1.1,-1.2,-3.1,-2.8)
mlo <- smms(startval, dd, gg, mc_cores = 1L, hessian_matrix = TRUE)
tval <- seq(0.01, 30, length = 50)

occupancy_prob("well", tval, mlo$opt$par, gg)
occupancy_prob_ci_band("well", tval, mlo$opt$par, gg, hessian = mlo$hess)
