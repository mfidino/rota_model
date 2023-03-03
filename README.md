# Fitting and interpretting the rota et al. (2016) occupancy model for interacting species

This repo is for this blog post here:


If you are here that means you are likely interested in running the code from that post. The R sub-folder pretty much has all you need.

1. `.R/fit_real_data.R`: Fits the occupancy model to coyote and gray squirrel camera trap data.
2. `./R/linear_predictor_examples.R`: Some examples I used in the blog post for how the first and second-order parameters work.
3. `./R/mcmc_utility.R`: Contains the `split_mcmc()` function I created to make predictions way easier with the MCMC output.
4. `./R/plot_results.R`: Plotting out the results from the model for both occupancy and detection.
5. `./R/simulate_data.R`: Not a part of the blog post, but a script that simualtes data and fits the model to those data.

The only other script you may want is in the `./nimble` sub-folder, which is `./nimble/rota_model.R`. It's the nimble code for the model I fitted to the coyote and gray squirrel data.
 
