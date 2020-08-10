  /* filename: well_decay_model.stan
   * author: Dylan H. Morris <dhmorris@princeton.edu>
   * 
   * description: estimate half-lives 
   * (log10 TCID50 / vol) directly 
   * from raw titration (well) data
   */


data {

  ////////////////////////////////
  // general data
  ////////////////////////////////
  int<lower = 1> n_total_datapoints;
  int<lower = 0, upper = n_total_datapoints> n_datapoints_used;
  int<lower = 1> n_experiments;
  int<lower = 1> n_titers;
  int<lower = 0, upper = 1> well_status[n_total_datapoints];
  int dilution[n_total_datapoints];
  
  int<lower = 1> experiment_id[n_total_datapoints];

  int<lower = 1, upper = n_titers> titer_id[n_total_datapoints];
  int<lower = 1, upper = n_titers> titer_experiment_id[n_titers];
  vector<lower = 0>[n_titers] time;

  ////////////////////////////////
  // parameters for priors
  ////////////////////////////////

  real intercept_prior_mean;
  real<lower = 0> intercept_prior_sd;
  
  real log_hl_prior_mean;
  real<lower = 0> log_hl_prior_sd;

 
  real<lower = 0> mode_sd_intercept;
  real<lower = 0> sd_sd_intercept;

  real lower_lim_decay_rate;
  
  int<lower = 0, upper = 1> debug;

}

transformed data {
  real lld;
  lld = lower_lim_decay_rate;
}


parameters{
  vector[n_experiments] mean_intercept;
  vector[n_titers] error_intercept;
  vector[n_experiments] log_half_life;    
  vector<lower = 0>[n_experiments] sd_intercept;

}

transformed parameters {
  vector[n_titers] intercept;
  vector[n_titers] sampled_titer;


  vector<lower = lld>[n_experiments] decay_rate;

  decay_rate = log10(2) ./ exp(log_half_life);
  
  for(i_titer in 1:n_titers){
    int exp_id = titer_experiment_id[i_titer];
    real t = time[i_titer];
    
    real ith_predicted_titer = 0;
    real ith_intercept = mean_intercept[exp_id] +
      sd_intercept[exp_id] * error_intercept[i_titer];
    real decay = decay_rate[exp_id] * t;
    if(decay < 80){
      ith_predicted_titer = ith_intercept - decay;
    } else {
      ith_predicted_titer = ith_intercept - 80;
    }
    
    // save values
    intercept[i_titer] = ith_intercept;
    sampled_titer[i_titer] = ith_predicted_titer;
  }
}

model {

  // observation process: poisson single hit
  for (i_dat in 1:n_datapoints_used) {
    int i_titer = titer_id[i_dat];
    real dilute_dose = sampled_titer[i_titer] + dilution[i_dat];

    if(well_status[i_dat] == 0) {
      target += poisson_lpmf(0 | log(2) * pow(10, dilute_dose));

      // use debug to check log prob getting properly incremented
      if(debug) {
        print("time: ", time[i_titer]);
        print("dose: ", dilute_dose);
        print("lpmf: ", poisson_lpmf(0 | log(2) * pow(10, dilute_dose)));
      }
      
    } else if (well_status[i_dat] == 1) {
      // use debug to check log prob getting properly incremented
      if(debug) {
        print("lccdf: ", poisson_lccdf(0 | log(2) * pow(10, dilute_dose)));
      }

      target += poisson_lccdf(0 | log(2) * pow(10, dilute_dose));

      
    } else {
      // well status must be negative (0) or positive (1)
      reject("well_status data must be one or zero, found",
             well_status[i_dat]);
      
    }
  }

  // priors
  for(i_exp in 1:n_experiments)
    mean_intercept[i_exp] ~ normal(intercept_prior_mean,
                                   intercept_prior_sd);
  
  log_half_life ~ normal(log_hl_prior_mean,
                         log_hl_prior_sd);

  error_intercept ~ normal(0, 1);
  
  sd_intercept  ~ normal(mode_sd_intercept,
                         sd_sd_intercept);
}

generated quantities {
  vector[n_titers] intercepts_pred;
  vector[n_titers] sampled_titers_pred;

  for(i_titer in 1:n_titers){
    int i_exp = titer_experiment_id[i_titer];
    real t = time[i_titer];
    
    real ith_predicted_titer = 0;
    real ith_intercept = normal_rng(mean_intercept[i_exp],
                                    sd_intercept[i_exp]);
    
    ith_predicted_titer = ith_intercept - 
      decay_rate[i_exp] * t;
      
    // save values
    intercepts_pred[i_titer] = ith_intercept;
    sampled_titers_pred[i_titer] = ith_predicted_titer;
  }
}
