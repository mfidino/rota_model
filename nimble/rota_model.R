rota_model <- nimble::nimbleCode(
  {
    # latent state linear predictors
    # linear predictors. TS = True state
    # TS = U
    # The latent-state model
    for(site in 1:nsite){
      psi[site,1] <- 1 
      # TS = A
      psi[site,2] <- exp(
        inprod(
          x[site,1:npar_psi],
          spA_psi[1:npar_psi]
        )
      )
      # TS = B
      psi[site,3] <- exp(
        inprod(
          x[site,1:npar_psi],
          spB_psi[1:npar_psi]
        )
      )
      # TS = AB
      psi[site,4] <- exp(
        inprod(
          x[site,1:npar_psi],
          spA_psi[1:npar_psi]
        ) +
        inprod(
          x[site,1:npar_psi],
          spB_psi[1:npar_psi]
        ) +
        inprod(
          x[site,1:npar_psi],
          spAB_psi[1:npar_psi]
        )
      )
      z[site] ~ dcat(
        psi[site,1:nstate]
      )
    }
    # The observation model linear predictors.
    # OS = Observed state
    for(site in 1:nsite){
      # OS = U
      rho[site,1] <- 1
      # OS = A
      rho[site,2] <- exp(
        inprod(
          x[site,1:npar_rho],
          spA_rho[1:npar_rho]
        )
      )
      # OS = B
      rho[site,3] <- exp(
        inprod(
          x[site,1:npar_rho],
          spB_rho[1:npar_rho]
        )
      )
      # OS = AB
      rho[site,4] <- exp(
        inprod(
          x[site,1:npar_rho],
          spA_rho[1:npar_rho]
        ) +
        inprod(
          x[site,1:npar_rho],
          spB_rho[1:npar_rho]
        ) +
        inprod(
          x[site,1:npar_rho],
          spAB_rho[1:npar_rho]
        )
      )
    }
    # And use these to fill in the rho detection matrix (rdm)
    # TS = U
    rdm[1:nsite,1,1] <- rho[1:nsite,1] # ------- OS = U
    rdm[1:nsite,1,2] <- rep(0, nsite) # -------- OS = A
    rdm[1:nsite,1,3] <- rep(0, nsite) # -------- OS = B
    rdm[1:nsite,1,4] <- rep(0, nsite) # -------- OS = AB
    # TS = A
    rdm[1:nsite,2,1] <- rho[1:nsite,1] # ------- OS= U
    rdm[1:nsite,2,2] <- rho[1:nsite,2] # ------- OS = A
    rdm[1:nsite,2,3] <- rep(0, nsite) # -------- OS = B
    rdm[1:nsite,2,4] <- rep(0,nsite) # --------- OS = AB
    # TS = B
    rdm[1:nsite,3,1] <- rho[1:nsite,1] # ------- OS = U
    rdm[1:nsite,3,2] <- rep(0, nsite) # -------- OS = A
    rdm[1:nsite,3,3] <- rho[1:nsite,3] # ------- OS = B
    rdm[1:nsite,3,4] <- rep(0, nsite) # -------- OS = AB
    # TS = AB
    rdm[1:nsite,4,1] <- rho[1:nsite,1] # ------- OS = U
    rdm[1:nsite,4,2] <- rho[1:nsite,2] # ------- OS = A
    rdm[1:nsite,4,3] <- rho[1:nsite,3] # ------- OS = B
    rdm[1:nsite,4,4] <- rho[1:nsite,4] # ------- OS = AB
    # observational model. z indexes the correct
    #   part of rdm
    for(site in 1:nsite){
      for(visit in 1:nvisit){
        y[site,visit] ~ dcat(
          rdm[
            site,
            z[site],
            1:nstate
          ]
        )
      }
    }
    # priors
    for(psii in 1:npar_psi){
      spA_psi[psii] ~ dlogis(0,1)
      spB_psi[psii] ~ dlogis(0,1)
      spAB_psi[psii] ~ dlogis(0,1)
    }
    for(rhoi in 1:npar_rho){
      spA_rho[rhoi] ~ dlogis(0,1)
      spB_rho[rhoi] ~ dlogis(0,1)
      spAB_rho[rhoi] ~ dlogis(0,1)
    }
  }
)
