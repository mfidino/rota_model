library(nimble)
library(dplyr)
library(stringi)
library(MCMCvis)

# Read in the data
data <- read.csv(
  "./data/squirrel_coyote.csv"
)

dplyr::glimpse(data)


# figure out some general info about the data
nsite <- length(
  unique(
    data$Site
  )
)

ndays <- length(
  grep(
    "Day",
    colnames(
      data
    )
  )
)

# Determine community state at for each
#  secondary sampling period.

# This will store all of our detection 
#  non-detection data
com_state <- matrix(
  NA,
  ncol = ndays,
  nrow = nsite
)

# paste the species data together. The data has been
#  sorted by site and species already (alphabetically),
#  so it would go coyote data, then gray squirrel
#  data.
sp_combo <- data %>% 
  dplyr::group_by(Site) %>% 
  dplyr::summarise_at(
    dplyr::vars(
      dplyr::starts_with("Day")
    ),
    .funs = function(x) paste0(x, collapse = "-")
  )

# construct a map for what each combo means
combo_map <- data.frame(
  # what we have
  combo = c("NA-NA", "0-0", "1-0", "0-1", "1-1"),
  # what the model wants 
  value = c(NA, 1,2,3,4),
  # to remember what each one means
  meaning = c(
    "no sampling",
    "no species",
    "coyote",
    "gray squirrel",
    "coyote & gray squirrel"
  )
)

# use the combo_map and sp_combo to construct the 
#  detection matrix. Easiest way to do this is
#  with this stringi function for each column.
for(i in 1:ncol(com_state)){
  day_vec <- sp_combo[,paste0("Day_",i), drop = TRUE]
  day_vec <- stringi::stri_replace_all_fixed(
    day_vec,
    combo_map$combo,
    combo_map$value,
    vectorize_all = FALSE
  )
  # put into the detection matrix
  com_state[,i] <- as.numeric(day_vec)
}

# com_state is what is the data we will supply to the model,
#  so let's pull in our covariate data.
covs <- read.csv(
  "./data/site_covariates.csv"
)

# scale the covariates
cov_scale <- covs %>% 
  dplyr::summarise_if(
    is.numeric,
    scale
  )

# apply PCA to generate an urbanization score
(cov_pca <- prcomp(
  cov_scale
))

#                           PC1        PC2         PC3
# Impervious         -0.6360670  0.3094596 0.706861714
# Ndvi                0.6363000 -0.3078588 0.707350923
# Population_density -0.4365101 -0.8996987 0.001090312

# currently,these data represent a gradient of built environment when
#  negative to more vegetation when positive. Let's flip it so that
#  positive means more built environment. This is simple, just 
#  multiply the loadings (rotation) and the first principal 
#  component (x) by -1. 

cov_pca$rotation[,1] <- cov_pca$rotation[,1] * -1
cov_pca$x[,1] <- cov_pca$x[,1] * -1


# put together lists for analysis
data_list <- list(
  # detection /non-detection data
  y = com_state,
  # design matrix for occupancy and detection,
  #  currently assuming the same covariates
  x = cbind(1, cov_pca$x[,1])
)

constant_list <- list(
  nsite = nsite,
  nvisit = ncol(com_state),
  nstate = dplyr::n_distinct(data$Species)^2,
  npar_psi = ncol(data_list$x),
  npar_rho = ncol(data_list$x)
)


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

# load the model
source("./nimble/rota_model.R")

my_model <- nimble::nimbleMCMC(
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

saveRDS(
  my_model,
  "./nimble/mcmc_output.rds"
)

my_summary <- MCMCvis::MCMCsummary(
  my_model,
  digits = 2
)

#              mean    sd   2.5%   50%   97.5% Rhat n.eff
# spAB_psi[1] -0.19 0.630 -1.400 -0.20  1.1000    1  1402
# spAB_psi[2]  0.95 0.490  0.058  0.92  2.0000    1  1166
# spAB_rho[1]  0.36 0.250 -0.140  0.37  0.8200    1  4533
# spAB_rho[2] -0.37 0.190 -0.760 -0.36  0.0055    1  4743
# spA_psi[1]   0.57 0.560 -0.460  0.55  1.8000    1  1585
# spA_psi[2]  -0.70 0.380 -1.400 -0.71  0.0980    1  1770
# spA_rho[1]  -2.60 0.160 -2.900 -2.50 -2.2000    1  2628
# spA_rho[2]  -0.37 0.120 -0.610 -0.37 -0.1300    1  2622
# spB_psi[1]   1.20 0.490  0.250  1.20  2.2000    1  1667
# spB_psi[2]  -0.84 0.360 -1.600 -0.81 -0.2100    1  1285
# spB_rho[1]  -0.91 0.052 -1.000 -0.91 -0.8000    1 14993
# spB_rho[2]   0.26 0.035  0.190  0.26  0.3300    1 14632
