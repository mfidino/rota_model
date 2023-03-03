
# Example one: Species occur together as often as you would expect
#  by chance. No slope terms, second-order parameters set to 0.
a0 <- 1
b0 <- -0.5
c0 <- 0

# two, species, so four states
beta_psi <- rep(NA, 4)

# state 1: species not present
beta_psi[1] <- 1

# state 2: species A present
beta_psi[2] <- exp(a0)

# state 3: species B present
beta_psi[3] <- exp(b0)

# state 4: speecies A & B present
beta_psi[4] <- exp(a0 + b0 + c0)

# convert to probability
psi <- beta_psi / sum(beta_psi)

# see that they occur together as much as you would expect by chance

# calculate the marginal occupancy of species A,
#  which is the sum of all states that include species A
marginal_occ_a <- psi[2] + psi[4]

# calculate the marginal occupancy of species B,
#  which is the sum of all states that include species B
marginal_occ_b <- psi[3] + psi[4]

# the product of these two should be equivalent to psi[4]
c(
  "product of marginal occupancy" =  marginal_occ_a * marginal_occ_b,
  "Probability of state 4" = psi[4]
)
# they are the same!
# product of marginal occupancy 
# 0.2760043 
# Probability of state 4 
# 0.2760043 


# Example two: Species occur more together than you would expect 
#  by chance. No slope terms but a positive second-order parameter.
a0 <- 1
b0 <- -0.5
c0 <- 1.5

# two, species, so four states
beta_psi <- rep(NA, 4)

# state 1: species not present
beta_psi[1] <- 1

# state 2: species A present
beta_psi[2] <- exp(a0)

# state 3: species B present
beta_psi[3] <- exp(b0)

# state 4: speecies A & B present
beta_psi[4] <- exp(a0 + b0 + c0)

# convert to probability
psi <- beta_psi / sum(beta_psi)

# calculate the marginal occupancy of species A,
#  which is the sum of all states that include species A
marginal_occ_a <- psi[2] + psi[4]

# calculate the marginal occupancy of species B,
#  which is the sum of all states that include species B
marginal_occ_b <- psi[3] + psi[4]

# psi[4] should be greater than the product of the marginals
c(
  "product of marginal occupancy" =  marginal_occ_a * marginal_occ_b,
  "Probability of state 4" = psi[4]
)
# State 4 happens more than you would expect by chance
#  product of marginal occupancy 
#  0.5889609 
#  Probability of state 4 
#  0.6307955 

# However, look and see how common each community state is. They
#  occur together way more than they do alone!
round(psi,2)
# [1] 0.09 0.23 0.05 0.63

# Example three: Species occur more together than you would expect 
#  by chance, but that decreases along an environmental gradient.
a0 <- 1
a1 <- 0.7
b0 <- -0.5
b1 <- 1.2
c0 <- 1.5
c1 <- -2

# choose three covariate values for sake of example
x <- c(-1,0,1)

# two, species, so four states, but with 
#  three covariate values we need a matrix now.
beta_psi <- matrix(
  NA,
  ncol =4,
  nrow = length(x)
)

# fill in beta_psi with a loop
for(i in 1:nrow(beta_psi)){
  # state 1: species not present
  beta_psi[i,1] <- 1
  
  # state 2: species A present
  beta_psi[i,2] <- exp(a0 + a1 * x[i])
  
  # state 3: species B present
  beta_psi[i,3] <- exp(b0 + b1 * x[i])
  
  # state 4: speecies A & B present
  beta_psi[i,4] <- exp(a0 + a1 * x[i] + b0 + b1 * x[i] + c0 + c1 * x[i])
}


# convert to probability
psi <- sweep(
  beta_psi,
  1,
  rowSums(beta_psi),
  FUN = "/"
)

colnames(psi) <- c("none", "A", "B", "AB")

# look at probability of each state, we see that Pr(AB) goes down
#  as our environmental gradient increases (i.e., goes from
#  negative to positive).
round(psi,2)
#       none    A    B   AB
#  [1,] 0.09 0.13 0.02 0.76 # at x = -1
#  [2,] 0.09 0.23 0.05 0.63 # at x =  0
#  [3,] 0.07 0.36 0.13 0.44 # at x =  1
