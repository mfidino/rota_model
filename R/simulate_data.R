library(nimble)

# simulate data for the rota two species model!

nsite <- 150
nvisit <- 4

# The intercept and slope terms for species occupancy,
#  as well as their interactions. Set up so they
#  co-occur more together when our environmental
#  gradient is positive.
spA_psi <- c(0.6, -1)
spB_psi <- c(-0.5, -0.75)
spAB_psi <- c(0.5, 1)

# The intercept and slope terms for species detection,
# as well as their interactions.
spA_rho <- c(0, 0.25)
spB_rho <- c(-0.7, 0)
spAB_rho <- c(0.5, 0.5)

# The covariate
set.seed(55)
x <- rnorm(nsite)

# make the design matrix (column of 1's for the intercept)
x <- cbind(1,x)

# calculate probability of each community state at each site
psi <- matrix(
  NA,
  ncol = 4, 
  nrow = nsite
)
colnames(psi) <- c("U", "A", "B", "AB")

# For no species present
U_linpred_psi <- 1 # exp(0)
A_linpred_psi <- exp(
  x %*% spA_psi
)
B_linpred_psi <- exp(
  x %*% spB_psi
)
AB_linpred_psi <- exp(
  x %*% spA_psi + x %*% spB_psi + x %*% spAB_psi
)
psi_denominator <- 
  U_linpred_psi + 
  A_linpred_psi + 
  B_linpred_psi +
  AB_linpred_psi
psi[,"U"] <- U_linpred_psi / psi_denominator
psi[,"A"] <- A_linpred_psi / psi_denominator
psi[,"B"] <- B_linpred_psi / psi_denominator
psi[,"AB"] <- AB_linpred_psi / psi_denominator

# ensure all rows sum to 1
all.equal(
  rowSums(psi),
  rep(1, nsite)
)

# simulate the community state
z_cat <- apply(
  psi,
  1,
  function(x){
    sample(
      x = c("U","A","B","AB"),
      size = 1,
      replace = TRUE,
      prob = x
    )
  }
)

# convert to numeric
z <- as.numeric(
  factor(
    z_cat,
    levels = colnames(psi)
  )
)
# now construct the observational data. This is a little
#  more complex as the you can't observe some states
#  given the true state.

rho <- array(
  0,
  dim = c(nsite,4,4)
)


# For no species present
U_linpred_rho <- 1 # exp(0)
A_linpred_rho <- exp(
  x %*% spA_rho
)
B_linpred_rho <- exp(
  x %*% spB_rho
)
AB_linpred_rho <- exp(
  x %*% spA_rho + x %*% spB_rho + x %*% spAB_rho
)
rho_denominator_A <- 
  U_linpred_rho +
  A_linpred_rho
rho_denominator_B <- 
  U_linpred_rho +
  B_linpred_rho
rho_denominator_AB <- 
  U_linpred_rho + 
  A_linpred_rho + 
  B_linpred_rho +
  AB_linpred_rho
# OS = observed state
# If true state is U then you detect it perfectly
rho[,1,1] <- 1 # OS = U
rho[,1,2] <- 0 # OS = A
rho[,1,3] <- 0 # OS = B
rho[,1,4] <- 0 # OS = AB
# If true state is A...
rho[,2,1] <- U_linpred_rho / rho_denominator_A # OS = U
rho[,2,2] <- A_linpred_rho / rho_denominator_A # OS = A
rho[,2,3] <- 0 # OS = B
rho[,2,4] <- 0 # OS = AB
# If true state is B...
rho[,3,1] <- U_linpred_rho / rho_denominator_B # OS = U
rho[,3,2] <- 0 # OS = A
rho[,3,3] <- B_linpred_rho / rho_denominator_B # OS = B
rho[,3,4] <- 0 # OS = AB
# If true is is AB...
rho[,4,1] <- U_linpred_rho / rho_denominator_AB # OS = U
rho[,4,2] <- A_linpred_rho / rho_denominator_AB # OS = A
rho[,4,3] <- B_linpred_rho / rho_denominator_AB # OS = B
rho[,4,4] <- AB_linpred_rho / rho_denominator_AB # OS = AB

# create the y matrix
y <- matrix(
  NA,
  ncol = nvisit,
  nrow = nsite
)
for(site in 1:nsite){
  y[site,] <-sample(
    x = 1:4,
    size = nvisit,
    replace = TRUE,
    prob = rho[site,z[site],]
  )
}

data_list <- list(
  y = y,
  x = x
)

constant_list <- list(
  nsite = nsite,
  nvisit = nvisit,
  nstate = ncol(psi),
  npar_psi = 2,
  npar_rho = 2
)

source("./nimble/rota_model.R")

# initial values function
inits <- function(){
  with(
    constant_list,
    {
      list(
        spA_psi = rnorm(npar_psi),
        spB_psi = rnorm(npar_psi),
        spAB_psi = rnorm(npar_psi),
        spA_rho = rnorm(npar_rho),
        spB_rho = rnorm(npar_rho),
        spAB_rho = rnorm(npar_rho),
        z = rep(nstate, nsite)
      )
    }
  )
}

longshot <- nimble::nimbleMCMC(
  code = rota_model,
  constants = constant_list,
  data = data_list,
  monitors = c(
    "spA_psi",
    "spB_psi",
    "spAB_psi",
    "spA_rho",
    "spB_rho",
    "spAB_rho"
  ),
  niter = 70000,
  nburnin = 20000,
  nchains = 2,
  thin = 2,
  inits = inits
)

my_summary <- MCMCvis::MCMCsummary(
  longshot,
  digits = 2
)
par_ests <- c(
  spAB_psi, 
  spAB_rho,
  spA_psi,
  spA_rho,
  spB_psi,
  spB_rho
)
jpeg("./plots/simulation_results.jpeg")
{
  par(mar = c(5,8,1,1))
  bbplot::blank(
    xlim = c(-3,3),
    ylim = c(0,13),
    bty = "l"
  )
  bbplot::axis_blank(
    side = 1
  )
  bbplot::axis_text(
    side = 1,
    line = 0.75
  )
  bbplot::axis_text(
    text = "Parameter estimate",
    side = 1,
    line = 2
  )
  bbplot::axis_blank(
    side = 2,
    at = seq(1,12),
    minor = FALSE
  )
  bbplot::axis_text(
    text = rev(row.names(my_summary)),
    side = 2,
    at = seq(1,12),
    las = 1,
    line = 1
  )
  for(i in 1:nrow(my_summary)){
    lines(
      x = my_summary[i,c(3,5)],
      y = rep(13 - i,2),
      lwd = 3
    )
  }
  points(
    x = my_summary[,4],
    y = rev(1:12),
    pch = 21,
    bg = "gray50",
    cex = 2.5
  )
  points(
    x = par_ests,
    y = rev(1:12),
    pch = 21,
    bg = "blue",
    cex = 2
  )
  legend(
    "topright",
    legend = c("Truth", "Estimate"),
    pch = 21,
    pt.cex = c(2,2.5),
    pt.bg = c("blue", "gray50"),
    bty = "n"
  )
}
dev.off()
