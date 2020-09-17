// Add constriants for data   
data {
  int<lower = 0> N;  // number of tests in the sample (3330 for Santa Clara)
  int<lower = 0, upper = 1> y[N];  // 1 if positive, 0 if negative
  int<lower = 1> J_time; //#weeks
  vector[N] male;  // -0.5 if female, 0.5 if male. Could add constraint back if we use 0/1 coding
  int<lower = 1, upper = 3> race[N];  // 1=white, 2=black, 3=other
  int<lower = 1, upper = 4> age[N];  // 1=0-17, 2=18-39, 3=40-64, 4=65+
  int<lower = 1, upper = J_time> week[N];
  real<lower = 1, upper = J_time> week_real[J_time]; //Introduced for Exponentiated quadratic covariance function calculation
  int<lower = 1, upper = 2> county[N];
  int<lower = 1, upper = 8> sex_age[N];
  //int<lower = 0> N_zip;  // number of zip codes (58 in this case)
  //int<lower = 1, upper = N_zip> zip[N];  // zip codes 1 through 58
  //vector[N_zip] x_zip;  // predictors at the zip code level
  int<lower = 0> J_spec;
  int<lower = 0> y_spec [J_spec];
  int<lower = 0> n_spec [J_spec];
  int<lower = 0> J_sens;
  int<lower = 0> y_sens [J_sens];
  int<lower = 0> n_sens [J_sens];
  int<lower = 0> J;  // number of population cells, J = 2*4*3*3
  vector<lower = 0>[J] N_pop;  // population sizes for poststratification
  vector<lower = 0>[J] N_acs;  // population sizes for poststratification
  real intercept_prior_mean;
  real<lower = 0> intercept_prior_scale;
  real<lower = 0> coef_prior_scale;
  //real<lower = 0> coef_prior_scale_time;
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
  //real<lower = 0> sigma_time;
  real<lower = 0> sigma_county;
  real<lower = 0> sigma_sexage;
  real<lower=0> rho;
  real<lower=0> alpha;
  real<lower=0> sigma;
  //real<lower = 0> sigma_zip;
  vector<multiplier = sigma_race>[3] a_race;  // varying intercepts for ethnicity
  vector<multiplier = sigma_age>[4] a_age;  // varying intercepts for age category
  //vector<multiplier = sigma_zip>[N_zip] a_zip;  // varying intercepts for zip code
  vector<multiplier = sigma_county>[2] a_county;
  //vector<multiplier = sigma_time>[J_time] a_time;
  vector[J_time] a_time;
  vector<multiplier = sigma_sexage>[8] a_sexage;
}
transformed parameters { 
  vector[J_spec] spec; 
  vector[J_sens] sens;
  vector[J_time] mu = rep_vector(0, J_time);  
  spec = inv_logit(logit_spec); 
  sens = inv_logit(logit_sens);
}
model {
  vector[N] p;
  vector[N] p_sample;
  matrix[J_time, J_time] L_K;
  matrix[J_time, J_time] K = cov_exp_quad(week_real, alpha, rho);
  real sq_sigma = square(sigma);
  p = inv_logit(b[1] + b[2]*male  + a_race[race] + a_age[age] +a_sexage[sex_age]+ a_time[week] + a_county[county]); 
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
  //a_time ~ normal(0, sigma_time);
  for (n in 1:J_time)
    K[n, n] = K[n, n] + sq_sigma;
  L_K = cholesky_decompose(K);
  rho ~ inv_gamma(5, 5);
  alpha ~ std_normal();
  sigma ~ std_normal();
  a_time ~ multi_normal(mu, L_K);
  a_county ~ normal(0, sigma_county);
  a_sexage ~ normal(0, sigma_sexage);
  //a_zip ~ normal(0, sigma_zip);
  // prior on centered intercept
  //b[1] + b[2] * mean(male) + b[3] * mean(x_zip[zip])
  //    ~ normal(intercept_prior_mean, intercept_prior_scale);
  b[1] + b[2] * mean(male)
      ~ normal(intercept_prior_mean, intercept_prior_scale);
  b[2] ~ normal(0, coef_prior_scale);
  sigma_race ~ normal(0, coef_prior_scale);
  sigma_age ~ normal(0, coef_prior_scale);
  //sigma_time ~ normal(0, coef_prior_scale_time);
  sigma_county ~ normal(0, coef_prior_scale);
  sigma_sexage ~ normal(0, coef_prior_scale);
  //sigma_zip ~ normal(0, coef_prior_scale);
  //b[3] ~ normal(0, coef_prior_scale / sd(x_zip[zip]));  // prior on scaled coefficient
}
generated quantities {
  vector[J_time] p_avg_h;
  vector[J_time] p_avg_c;  
  vector[J] p_pop;  // population prevalence in the J poststratification cells
  int count;
for (time in 1:J_time){
  count = 1;
  for (i_male in 0:1) {
    for (i_age in 1:4) {
      for (i_race in 1:3) {
          for (i_county in 1:2){
          p_pop[count] = inv_logit(b[1]
                                   + b[2] * (i_male - 0.5)
                                   + a_race[i_race]
                                   + a_age[i_age]+a_time[time]+a_sexage[i_male*2+i_age] + a_county[i_county]);
          count += 1;
      }
    }
  }
  }
  p_avg_h[time] = sum(N_pop .* p_pop) / sum(N_pop);
  p_avg_c[time] = sum(N_acs .* p_pop) / sum(N_acs);
}
  
}
// generated quantities {
//   real p_avg;
//   real p_male;
//   real p_female;
//   real p_0_17;
//   real p_18_39;
//   real p_40_64;
//   real p_65;
//   real p_white;
//   real p_black;
//   real p_other;
//   vector[J] p_pop;  // population prevalence in the J poststratification cells
//   // Further postratify
//   vector[J] sex_record; 
//   vector[J] age_record; 
//   vector[J] ethnicity_record;
//   int position_0_17;
//   int position_18_39; 
//   int position_40_64;
//   int position_65;
//   int position_white;
//   int position_black;
//   int position_other;
//   // Sex
//   vector[3*4] p_pop_male;
//   vector[3*4] p_pop_female;
//   vector[3*4] N_male;
//   vector[3*4] N_female;
//   // Age
//   vector[2*3] p_pop_0_17;
//   vector[2*3] p_pop_18_39;
//   vector[2*3] p_pop_40_64;
//   vector[2*3] p_pop_65;
//   vector[2*3] N_0_17;
//   vector[2*3] N_18_39;
//   vector[2*3] N_40_64;
//   vector[2*3] N_65;
//   // Ethnicity
//   vector[2*4] p_pop_white;
//   vector[2*4] p_pop_black;
//   vector[2*4] p_pop_other;
//   vector[2*4] N_white;
//   vector[2*4] N_black;
//   vector[2*4] N_other;
//   
//   int count;
//   count = 1;
//   for (i_male in 0:1) {
//     for (i_age in 1:4) {
//       for (i_eth in 1:3) {
//           p_pop[count] = inv_logit(b[1]
//                                    + b[2] * (i_male - 0.5)
//                                    + a_eth[i_eth]
//                                    + a_age[i_age]);
//           sex_record[count] = i_male;
//           age_record[count] = i_age;
//           ethnicity_record[count] = i_eth;
//           count += 1;
//       }
//     }
//   }
//   p_avg = sum(N_pop .* p_pop) / sum(N_pop);
//   // Sex
//   p_pop_female = p_pop[1:3*4];
//   p_pop_male = p_pop[1+3*4: J];
//   N_female = N_pop[1:3*4];
//   N_male = N_pop[1+3*4: J];
//   p_male = sum(N_male .* p_pop_male) / sum(N_male);
//   p_female = sum(N_female .* p_pop_female) / sum(N_female);
//   // Age
//   position_0_17 = 1;
//   position_18_39 = 1;
//   position_40_64 = 1;
//   position_65 = 1;
//   for (i in 1:J){
//     if(age_record[i] == 1){
//       p_pop_0_17[position_0_17] = p_pop[i];
//       N_0_17[position_0_17] = N_pop[i];
//       position_0_17 += 1;
//     }else if(age_record[i] == 2){
//       p_pop_18_39[position_18_39] = p_pop[i];
//       N_18_39[position_18_39] = N_pop[i];
//       position_18_39 += 1;
//     }else if(age_record[i] == 3){
//       p_pop_40_64[position_40_64] =  p_pop[i];
//       N_40_64[position_40_64] = N_pop[i];
//       position_40_64 += 1;
//     }else {
//       p_pop_65[position_65] = p_pop[i];
//       N_65[position_65] = N_pop[i];
//       position_65 += 1;
//     }
//   }
//   p_0_17 = sum(N_0_17 .* p_pop_0_17) / sum(N_0_17);
//   p_18_39 = sum(N_18_39 .* p_pop_18_39) / sum(N_18_39);
//   p_40_64 = sum(N_40_64 .* p_pop_40_64) / sum(N_40_64);
//   p_65 = sum(N_65 .*p_pop_65) / sum(N_65);
//   // Ethicity
//   position_white = 1;
//   position_black = 1;
//   position_other = 1;
//   for (i in 1:J){
//     if(ethnicity_record[i] == 1){
//       p_pop_white[position_white] = p_pop[i];
//       N_white[position_white] = N_pop[i];
//       position_white += 1;
//     }else if(ethnicity_record[i] == 2){
//       p_pop_black[position_black] =  p_pop[i];
//       N_black[position_black] = N_pop[i];
//       position_black += 1;
//     }else {
//       p_pop_other[position_other] = p_pop[i];
//       N_other[position_other] = N_pop[i];
//       position_other += 1;
//     }
//   }
//   p_white = sum(N_white .* p_pop_white) / sum(N_white);
//   p_black = sum(N_black .* p_pop_black) / sum(N_black);
//   p_other = sum(N_other .*p_pop_other) / sum(N_other);
//   
//   
//   
//   
//   
// }
