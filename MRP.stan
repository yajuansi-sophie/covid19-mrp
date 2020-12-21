data {
  int<lower = 0> N;  // number of tests in the sample 
  int<lower = 0, upper = 1> y[N];  // 1 if positive, 0 if negative
  int<lower = 1> J_time; //#weeks
  vector[N] male;  // -0.5 if female, 0.5 if male. 
  int<lower = 1, upper = 3> race[N];  // 1=white, 2=black, 3=other
  int<lower = 1, upper = 5> age[N];  // 1=0-17, 2=18-39, 3=40-64, 4=65-75, 5 = 75+
  int<lower = 1, upper = J_time> week[N]; //week
  int<lower = 1, upper = 2> county[N]; //county indcators
  int<lower = 1, upper = 5> sex_age[N]; //interactions between sex and age
  int<lower = 0> J_spec; //prior information for sensitivity and specifity
  int<lower = 0> y_spec [J_spec];
  int<lower = 0> n_spec [J_spec];
  int<lower = 0> J_sens;
  int<lower = 0> y_sens [J_sens];
  int<lower = 0> n_sens [J_sens];
  int<lower = 0> J;  // number of population cells, J = 2*3*5*2
  vector<lower = 0>[J] N_pop;  // hospital population sizes
  vector<lower = 0>[J] N_acs;  // community population sizes
  real intercept_prior_mean;
  real<lower = 0> intercept_prior_scale;
  real<lower = 0> coef_prior_scale;
  real<lower = 0> coef_prior_scale_time;
  real<lower = 0> logit_spec_prior_scale;
  real<lower = 0> logit_sens_prior_scale;
}
parameters {
  real mu_logit_spec;
  real mu_logit_sens;
  real<lower = 0> sigma_logit_spec;
  real<lower = 0> sigma_logit_sens;
  vector<offset = mu_logit_spec, multiplier = sigma_logit_spec>[J_spec] logit_spec;
  vector<offset = mu_logit_sens, multiplier = sigma_logit_sens>[J_sens] logit_sens;
  vector[2] b;  // intercept and coef for male
  real<lower = 0> sigma_race;
  real<lower = 0> sigma_age;
  real<lower = 0> sigma_time;
  real<lower = 0> sigma_county;
  real<lower = 0> sigma_sexage;
  vector<multiplier = sigma_race>[3] a_race;  
  vector<multiplier = sigma_age>[5] a_age; 
  vector<multiplier = sigma_county>[2] a_county;
  vector<multiplier = sigma_time>[J_time] a_time;
  vector<multiplier = sigma_sexage>[5] a_sexage;
}
transformed parameters { 
  vector[J_spec] spec; 
  vector[J_sens] sens;
  vector[N] a_sexage_trans;
  spec = inv_logit(logit_spec); 
  sens = inv_logit(logit_sens);
  for (i in 1:N)
  a_sexage_trans[i] = (male[i] + 0.5)*a_sexage[sex_age[i]];
}
model {
  vector[N] p;
  vector[N] p_sample;
  p = inv_logit(b[1] + b[2]*male  + a_race[race] + a_age[age] + a_sexage_trans + a_time[week] + a_county[county]); 
  p_sample = p*sens[1] + (1-p)*(1-spec[1]);
  y ~ bernoulli(p_sample);
  y_spec ~ binomial(n_spec, spec);
  y_sens ~ binomial(n_sens, sens);
  logit_spec ~ normal(mu_logit_spec, sigma_logit_spec);
  logit_sens ~ normal(mu_logit_sens, sigma_logit_sens);
  sigma_logit_spec ~ normal(0, logit_spec_prior_scale);
  sigma_logit_sens ~ normal(0, logit_sens_prior_scale);
  a_race ~ normal(0, sigma_race);
  a_age ~ normal(0, sigma_age);
  a_time ~ normal(0, sigma_time);
  a_county ~ normal(0, sigma_county);
  a_sexage ~ normal(0, sigma_sexage);
  b[1] + b[2] * mean(male) ~ normal(intercept_prior_mean, intercept_prior_scale);
  b[2] ~ normal(0, coef_prior_scale);
  sigma_race ~ normal(0, coef_prior_scale);
  sigma_age ~ normal(0, coef_prior_scale);
  sigma_time ~ normal(0, coef_prior_scale_time);
  sigma_county ~ normal(0, coef_prior_scale);
  sigma_sexage ~ normal(0, coef_prior_scale);
}
generated quantities {
  vector[J_time] p_avg_h;
  vector[J_time] p_avg_c;  
  vector[J] p_pop;  // population prevalence in the J poststratification cells
  int count;
  for (time in 1:J_time){
  count = 1;
  for (i_male in 0:1) {
    for (i_age in 1:5) {
      for (i_race in 1:3) {
          for (i_county in 1:2){
          p_pop[count] = inv_logit(b[1] + b[2] * (i_male - 0.5) + a_race[i_race] +                           a_age[i_age]+a_time[time]+
                          i_male*a_sexage[i_age] + a_county[i_county]);
          count += 1;
          }
        }
      }
    }
  p_avg_h[time] = sum(N_pop .* p_pop) / sum(N_pop);
  p_avg_c[time] = sum(N_acs .* p_pop) / sum(N_acs);
  }
}

