  /* filename: well_titer_estimates.stan
   * author: Dylan H. Morris <dhmorris@princeton.edu>
   * 
   * description: estimate virus titers 
   * (log10 TCID50 / vol) directly 
   * from raw titration (well) data
   */

data {

  ////////////////////////////////
  // general data
  ////////////////////////////////
  int<lower = 1> n_total_datapoints;
  int<lower = 0, upper = n_total_datapoints> n_datapoints_used;
  int<lower = 1> n_titers;
  int<lower = 0, upper = 1> well_status[n_total_datapoints];
  int dilution[n_total_datapoints];
  
  int<lower = 1, upper = n_titers> titer_id[n_total_datapoints];


  ////////////////////////////////
  // parameters for priors
  ////////////////////////////////

  real titer_prior_mean;
  real<lower = 0> titer_prior_sd;

  ///////////////////////////////
  // flags
  ///////////////////////////////
  
  int<lower = 0, upper = 1> debug;
}


parameters{
  vector[n_titers] log10_titer;
}

model {

  // observation process: poisson single hit
  for (i_dat in 1:n_datapoints_used) {
    int i_titer = titer_id[i_dat];
    real ith_titer = log10_titer[i_titer];
    real dilute_dose = ith_titer + dilution[i_dat];
    real virions = log(2) * pow(10, dilute_dose);

    if(well_status[i_dat] == 0) {

      target += poisson_lpmf(0 | virions);
      
    } else if (well_status[i_dat] == 1) {

      target += poisson_lccdf(0 | virions);

    } else {
      // well status must be negative (0) or positive (1)
      reject("well_status data must be one or zero, found",
             well_status[i_dat]);
    }
  }

  // priors
  log10_titer ~ normal(titer_prior_mean,
                       titer_prior_sd);
}

generated quantities {
  vector[n_titers] sampled_titers_pred;

  for(t_id in 1:n_titers)
    sampled_titers_pred[t_id] = normal_rng(titer_prior_mean,
                                 titer_prior_sd);
}
