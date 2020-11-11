// Add constriants for data   
data {
  int<lower = 0> N;  // number of tests in the sample (3330 for Santa Clara)
  int<lower = 0, upper = 1> y[N];  // 1 if positive, 0 if negative
  int<lower = 1> J_time; //#weeks
  vector[N] male;  // -0.5 if female, 0.5 if male. Could add constraint back if we use 0/1 coding
  int<lower = 1, upper = 3> race[N];  // 1=white, 2=black, 3=other
  int<lower = 1, upper = 5> age[N];  // 1=0-17, 2=18-39, 3=40-64, 4=65-75, 5 = 75+
  int<lower = 1, upper = J_time> week[N];
  int<lower = 1, upper = 2> county[N];
  int<lower = 1, upper = 5> sex_age[N];
  //int<lower = 0> N_zip;  // number of zip codes (58 in this case)
  //int<lower = 1, upper = N_zip> zip[N];  // zip codes 1 through 58
  //vector[N_zip] x_zip;  // predictors at the zip code level
  int<lower = 0> J_spec;
  int<lower = 0> y_spec [J_spec];
  int<lower = 0> n_spec [J_spec];
  int<lower = 0> J_sens;
  int<lower = 0> y_sens [J_sens];
  int<lower = 0> n_sens [J_sens];
  int<lower = 0> J;  // number of population cells, J = 2*3*5*2
  vector<lower = 0>[J] N_pop;  // population sizes for poststratification
  vector<lower = 0>[J] N_acs;  // population sizes for poststratification
  real intercept_prior_mean;
  real<lower = 0> intercept_prior_scale;
  real<lower = 0> coef_prior_scale;
  real<lower = 0> coef_prior_scale_time;
  real<lower = 0> logit_spec_prior_scale;
  real<lower = 0> logit_sens_prior_scale;
  //inefficient code
  real offset;
  int n1;
  int n2;
  int n3;
  int n4;
  int n5;
  int n6;
  int n7;
  int n8;
  int n9;
  int n10;
  int n11;
  int n12;
  int n13;
  int n14;
  int n15;
  int n16;
  int n17;
  int n18;
  int n19;
  int n20;
  int n21;
  int n22;
  int n23;
  int n24;
  int n25;
  int n26;
  int n27;
  int n28;
  vector[n1] male_w1;
  int time_w1[n1];
  int<lower = 1, upper = 3> race_w1[n1];  
  int<lower = 1, upper = 5> age_w1[n1];  
  int<lower = 1, upper = 2> county_w1[n1];
  int<lower = 1, upper = 5> sex_age_w1[n1];
  vector[n2] male_w2;
  int time_w2[n2];
  int<lower = 1, upper = 3> race_w2[n2];  
  int<lower = 1, upper = 5> age_w2[n2];  
  int<lower = 1, upper = 2> county_w2[n2];
  int<lower = 1, upper = 5> sex_age_w2[n2];
  vector[n3] male_w3;
  int time_w3[n3];
  int<lower = 1, upper = 3> race_w3[n3];  
  int<lower = 1, upper = 5> age_w3[n3];  
  int<lower = 1, upper = 2> county_w3[n3];
  int<lower = 1, upper = 5> sex_age_w3[n3];
  vector[n4] male_w4;
  int time_w4[n4];
  int<lower = 1, upper = 3> race_w4[n4];  
  int<lower = 1, upper = 5> age_w4[n4];  
  int<lower = 1, upper = 2> county_w4[n4];
  int<lower = 1, upper = 5> sex_age_w4[n4];
  vector[n5] male_w5;
  int time_w5[n5];
  int<lower = 1, upper = 3> race_w5[n5];  
  int<lower = 1, upper = 5> age_w5[n5];  
  int<lower = 1, upper = 2> county_w5[n5];
  int<lower = 1, upper = 5> sex_age_w5[n5];
  vector[n6] male_w6;
  int time_w6[n6];
  int<lower = 1, upper = 3> race_w6[n6];  
  int<lower = 1, upper = 5> age_w6[n6];  
  int<lower = 1, upper = 2> county_w6[n6];
  int<lower = 1, upper = 5> sex_age_w6[n6];
  vector[n7] male_w7;
  int time_w7[n7];
  int<lower = 1, upper = 3> race_w7[n7];  
  int<lower = 1, upper = 5> age_w7[n7];  
  int<lower = 1, upper = 2> county_w7[n7];
  int<lower = 1, upper = 5> sex_age_w7[n7];
  vector[n8] male_w8;
  int time_w8[n8];
  int<lower = 1, upper = 3> race_w8[n8];  
  int<lower = 1, upper = 5> age_w8[n8];  
  int<lower = 1, upper = 2> county_w8[n8];
  int<lower = 1, upper = 5> sex_age_w8[n8];
  vector[n9] male_w9;
  int time_w9[n9];
  int<lower = 1, upper = 3> race_w9[n9];  
  int<lower = 1, upper = 5> age_w9[n9];  
  int<lower = 1, upper = 2> county_w9[n9];
  int<lower = 1, upper = 5> sex_age_w9[n9];
  vector[n10] male_w10;
  int time_w10[n10];
  int<lower = 1, upper = 3> race_w10[n10];  
  int<lower = 1, upper = 5> age_w10[n10];  
  int<lower = 1, upper = 2> county_w10[n10];
  int<lower = 1, upper = 5> sex_age_w10[n10];
  vector[n11] male_w11;
  int time_w11[n11];
  int<lower = 1, upper = 3> race_w11[n11];  
  int<lower = 1, upper = 5> age_w11[n11];  
  int<lower = 1, upper = 2> county_w11[n11];
  int<lower = 1, upper = 5> sex_age_w11[n11];
  vector[n12] male_w12;
  int time_w12[n12];
  int<lower = 1, upper = 3> race_w12[n12];  
  int<lower = 1, upper = 5> age_w12[n12];  
  int<lower = 1, upper = 2> county_w12[n12];
  int<lower = 1, upper = 5> sex_age_w12[n12];
  vector[n13] male_w13;
  int time_w13[n13];
  int<lower = 1, upper = 3> race_w13[n13];  
  int<lower = 1, upper = 5> age_w13[n13];  
  int<lower = 1, upper = 2> county_w13[n13];
  int<lower = 1, upper = 5> sex_age_w13[n13];
  vector[n14] male_w14;
  int time_w14[n14];
  int<lower = 1, upper = 3> race_w14[n14];  
  int<lower = 1, upper = 5> age_w14[n14];  
  int<lower = 1, upper = 2> county_w14[n14];
  int<lower = 1, upper = 5> sex_age_w14[n14];
  vector[n15] male_w15;
  int time_w15[n15];
  int<lower = 1, upper = 3> race_w15[n15];  
  int<lower = 1, upper = 5> age_w15[n15];  
  int<lower = 1, upper = 2> county_w15[n15];
  int<lower = 1, upper = 5> sex_age_w15[n15];
  vector[n16] male_w16;
  int time_w16[n16];
  int<lower = 1, upper = 3> race_w16[n16];  
  int<lower = 1, upper = 5> age_w16[n16];  
  int<lower = 1, upper = 2> county_w16[n16];
  int<lower = 1, upper = 5> sex_age_w16[n16];
  vector[n17] male_w17;
  int time_w17[n17];
  int<lower = 1, upper = 3> race_w17[n17];  
  int<lower = 1, upper = 5> age_w17[n17];  
  int<lower = 1, upper = 2> county_w17[n17];
  int<lower = 1, upper = 5> sex_age_w17[n17];
  vector[n18] male_w18;
  int time_w18[n18];
  int<lower = 1, upper = 3> race_w18[n18];  
  int<lower = 1, upper = 5> age_w18[n18];  
  int<lower = 1, upper = 2> county_w18[n18];
  int<lower = 1, upper = 5> sex_age_w18[n18];
  vector[n19] male_w19;
  int time_w19[n19];
  int<lower = 1, upper = 3> race_w19[n19];  
  int<lower = 1, upper = 5> age_w19[n19];  
  int<lower = 1, upper = 2> county_w19[n19];
  int<lower = 1, upper = 5> sex_age_w19[n19];
  vector[n20] male_w20;
  int time_w20[n20];
  int<lower = 1, upper = 3> race_w20[n20];  
  int<lower = 1, upper = 5> age_w20[n20];  
  int<lower = 1, upper = 2> county_w20[n20];
  int<lower = 1, upper = 5> sex_age_w20[n20];
  vector[n21] male_w21;
  int time_w21[n21];
  int<lower = 1, upper = 3> race_w21[n21];  
  int<lower = 1, upper = 5> age_w21[n21];  
  int<lower = 1, upper = 2> county_w21[n21];
  int<lower = 1, upper = 5> sex_age_w21[n21];
  vector[n22] male_w22;
  int time_w22[n22];
  int<lower = 1, upper = 3> race_w22[n22];  
  int<lower = 1, upper = 5> age_w22[n22];  
  int<lower = 1, upper = 2> county_w22[n22];
  int<lower = 1, upper = 5> sex_age_w22[n22];
  vector[n23] male_w23;
  int time_w23[n23];
  int<lower = 1, upper = 3> race_w23[n23];  
  int<lower = 1, upper = 5> age_w23[n23];  
  int<lower = 1, upper = 2> county_w23[n23];
  int<lower = 1, upper = 5> sex_age_w23[n23];
  vector[n24] male_w24;
  int time_w24[n24];
  int<lower = 1, upper = 3> race_w24[n24];  
  int<lower = 1, upper = 5> age_w24[n24];  
  int<lower = 1, upper = 2> county_w24[n24];
  int<lower = 1, upper = 5> sex_age_w24[n24];
  vector[n25] male_w25;
  int time_w25[n25];
  int<lower = 1, upper = 3> race_w25[n25];  
  int<lower = 1, upper = 5> age_w25[n25];  
  int<lower = 1, upper = 2> county_w25[n25];
  int<lower = 1, upper = 5> sex_age_w25[n25];
  vector[n26] male_w26;
  int time_w26[n26];
  int<lower = 1, upper = 3> race_w26[n26];  
  int<lower = 1, upper = 5> age_w26[n26];  
  int<lower = 1, upper = 2> county_w26[n26];
  int<lower = 1, upper = 5> sex_age_w26[n26];
  vector[n27] male_w27;
  int time_w27[n27];
  int<lower = 1, upper = 3> race_w27[n27];  
  int<lower = 1, upper = 5> age_w27[n27];  
  int<lower = 1, upper = 2> county_w27[n27];
  int<lower = 1, upper = 5> sex_age_w27[n27];
  vector[n28] male_w28;
  int time_w28[n28];
  int<lower = 1, upper = 3> race_w28[n28];  
  int<lower = 1, upper = 5> age_w28[n28];  
  int<lower = 1, upper = 2> county_w28[n28];
  int<lower = 1, upper = 5> sex_age_w28[n28];
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
  //real<lower = 0> sigma_zip;
  vector<multiplier = sigma_race>[3] a_race;  // varying intercepts for ethnicity
  vector<multiplier = sigma_age>[5] a_age;  // varying intercepts for age category
  //vector<multiplier = sigma_zip>[N_zip] a_zip;  // varying intercepts for zip code
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
  //a_zip ~ normal(0, sigma_zip);
  // prior on centered intercept
  //b[1] + b[2] * mean(male) + b[3] * mean(x_zip[zip])
  //    ~ normal(intercept_prior_mean, intercept_prior_scale);
  b[1] + b[2] * mean(male)
      ~ normal(intercept_prior_mean, intercept_prior_scale);
  b[2] ~ normal(0, coef_prior_scale);
  sigma_race ~ normal(0, coef_prior_scale);
  sigma_age ~ normal(0, coef_prior_scale);
  sigma_time ~ normal(0, coef_prior_scale_time);
  sigma_county ~ normal(0, coef_prior_scale);
  sigma_sexage ~ normal(0, coef_prior_scale);
  //sigma_zip ~ normal(0, coef_prior_scale);
  //b[3] ~ normal(0, coef_prior_scale / sd(x_zip[zip]));  // prior on scaled coefficient
}
generated quantities {
  vector[J_time] p_avg_h;
  vector[J_time] p_avg_c;  
  vector[J] p_pop;  // population prevalence in the J poststratification cells
  vector[J] p_pop_notime; // population prevalence in the J poststratification cells without considering time effect
  vector[J] p_pop_notime_up; 
  vector[J] p_pop_notime_low; 
  int count;
  real p_overall_h;
  real p_overall_c;
  real p_overall_h_up;
  real p_overall_h_low;
  real p_overall_c_up;
  real p_overall_c_low;
//posterior predictive check
  vector[n1] logit_p_w1 = b[1] + b[2]*male_w1  + a_race[race_w1] + a_age[age_w1] + (male_w1 + offset) .* a_sexage[sex_age_w1] + a_time[time_w1] + a_county[county_w1];
  int<lower = 0> y_rep_w1[n1] = bernoulli_logit_rng(logit_p_w1*sens[1] + (1-logit_p_w1)*(1-spec[1]));
  real<lower = 0> mean_y_rep_w1 = mean(to_vector(y_rep_w1));
  real<lower = 0> sd_y_rep_w1 = sd(to_vector(y_rep_w1));
  vector[n2] logit_p_w2 = b[1] + b[2]*male_w2  + a_race[race_w2] + a_age[age_w2] + (male_w2 + offset) .* a_sexage[sex_age_w2] + a_time[time_w2] + a_county[county_w2];
  int<lower = 0> y_rep_w2[n2] = bernoulli_logit_rng(logit_p_w2*sens[1] + (1-logit_p_w2)*(1-spec[1]));
  real<lower = 0> mean_y_rep_w2 = mean(to_vector(y_rep_w2));
  real<lower = 0> sd_y_rep_w2 = sd(to_vector(y_rep_w2));
  vector[n3] logit_p_w3 = b[1] + b[2]*male_w3  + a_race[race_w3] + a_age[age_w3] + (male_w3 + offset) .* a_sexage[sex_age_w3] + a_time[time_w3] + a_county[county_w3];
  int<lower = 0> y_rep_w3[n3] = bernoulli_logit_rng(logit_p_w3*sens[1] + (1-logit_p_w3)*(1-spec[1]));
  real<lower = 0> mean_y_rep_w3 = mean(to_vector(y_rep_w3));
  real<lower = 0> sd_y_rep_w3 = sd(to_vector(y_rep_w3));
  vector[n4] logit_p_w4 = b[1] + b[2]*male_w4  + a_race[race_w4] + a_age[age_w4] + (male_w4 + offset) .* a_sexage[sex_age_w4] + a_time[time_w4] + a_county[county_w4];
  int<lower = 0> y_rep_w4[n4] = bernoulli_logit_rng(logit_p_w4*sens[1] + (1-logit_p_w4)*(1-spec[1]));
  real<lower = 0> mean_y_rep_w4 = mean(to_vector(y_rep_w4));
  real<lower = 0> sd_y_rep_w4 = sd(to_vector(y_rep_w4));
  vector[n5] logit_p_w5 = b[1] + b[2]*male_w5  + a_race[race_w5] + a_age[age_w5] + (male_w5 + offset) .* a_sexage[sex_age_w5] + a_time[time_w5] + a_county[county_w5];
  int<lower = 0> y_rep_w5[n5] = bernoulli_logit_rng(logit_p_w5*sens[1] + (1-logit_p_w5)*(1-spec[1]));
  real<lower = 0> mean_y_rep_w5 = mean(to_vector(y_rep_w5));
  real<lower = 0> sd_y_rep_w5 = sd(to_vector(y_rep_w5));
  vector[n6] logit_p_w6 = b[1] + b[2]*male_w6 + a_race[race_w6] + a_age[age_w6] + (male_w6 + offset) .* a_sexage[sex_age_w6] + a_time[time_w6] + a_county[county_w6];
  int<lower = 0> y_rep_w6[n6] = bernoulli_logit_rng(logit_p_w6*sens[1] + (1-logit_p_w6)*(1-spec[1]));
  real<lower = 0> mean_y_rep_w6 = mean(to_vector(y_rep_w6));
  real<lower = 0> sd_y_rep_w6 = sd(to_vector(y_rep_w6));
  vector[n7] logit_p_w7 = b[1] + b[2]*male_w7 + a_race[race_w7] + a_age[age_w7] + (male_w7 + offset) .* a_sexage[sex_age_w7] + a_time[time_w7] + a_county[county_w7];
  int<lower = 0> y_rep_w7[n7] = bernoulli_logit_rng(logit_p_w7*sens[1] + (1-logit_p_w7)*(1-spec[1]));
  real<lower = 0> mean_y_rep_w7 = mean(to_vector(y_rep_w7));
  real<lower = 0> sd_y_rep_w7 = sd(to_vector(y_rep_w7));
  vector[n8] logit_p_w8 = b[1] + b[2]*male_w8  + a_race[race_w8] + a_age[age_w8] + (male_w8 + offset) .* a_sexage[sex_age_w8] + a_time[time_w8] + a_county[county_w8];
  int<lower = 0> y_rep_w8[n8] = bernoulli_logit_rng(logit_p_w8*sens[1] + (1-logit_p_w8)*(1-spec[1]));
  real<lower = 0> mean_y_rep_w8 = mean(to_vector(y_rep_w8));
  real<lower = 0> sd_y_rep_w8 = sd(to_vector(y_rep_w8));
  vector[n9] logit_p_w9 = b[1] + b[2]*male_w9  + a_race[race_w9] + a_age[age_w9] + (male_w9 + offset) .* a_sexage[sex_age_w9] + a_time[time_w9] + a_county[county_w9];
  int<lower = 0> y_rep_w9[n9] = bernoulli_logit_rng(logit_p_w9*sens[1] + (1-logit_p_w9)*(1-spec[1]));
  real<lower = 0> mean_y_rep_w9 = mean(to_vector(y_rep_w9));
  real<lower = 0> sd_y_rep_w9 = sd(to_vector(y_rep_w9));
  vector[n10] logit_p_w10 = b[1] + b[2]*male_w10  + a_race[race_w10] + a_age[age_w10] + (male_w10 + offset) .* a_sexage[sex_age_w10] + a_time[time_w10] + a_county[county_w10];
  int<lower = 0> y_rep_w10[n10] = bernoulli_logit_rng(logit_p_w10*sens[1] + (1-logit_p_w10)*(1-spec[1]));
  real<lower = 0> mean_y_rep_w10 = mean(to_vector(y_rep_w10));
  real<lower = 0> sd_y_rep_w10 = sd(to_vector(y_rep_w10));
  vector[n11] logit_p_w11 = b[1] + b[2]*male_w11  + a_race[race_w11] + a_age[age_w11] + (male_w11 + offset) .* a_sexage[sex_age_w11] + a_time[time_w11] + a_county[county_w11];
  int<lower = 0> y_rep_w11[n11] = bernoulli_logit_rng(logit_p_w11*sens[1] + (1-logit_p_w11)*(1-spec[1]));
  real<lower = 0> mean_y_rep_w11 = mean(to_vector(y_rep_w11));
  real<lower = 0> sd_y_rep_w11 = sd(to_vector(y_rep_w11));
  vector[n12] logit_p_w12 = b[1] + b[2]*male_w12  + a_race[race_w12] + a_age[age_w12] + (male_w12 + offset) .* a_sexage[sex_age_w12] + a_time[time_w12] + a_county[county_w12];
  int<lower = 0> y_rep_w12[n12] = bernoulli_logit_rng(logit_p_w12*sens[1] + (1-logit_p_w12)*(1-spec[1]));
  real<lower = 0> mean_y_rep_w12 = mean(to_vector(y_rep_w12));
  real<lower = 0> sd_y_rep_w12 = sd(to_vector(y_rep_w12));
  vector[n13] logit_p_w13 = b[1] + b[2]*male_w13  + a_race[race_w13] + a_age[age_w13] + (male_w13 + offset) .* a_sexage[sex_age_w13] + a_time[time_w13] + a_county[county_w13];
  int<lower = 0> y_rep_w13[n13] = bernoulli_logit_rng(logit_p_w13*sens[1] + (1-logit_p_w13)*(1-spec[1]));
  real<lower = 0> mean_y_rep_w13 = mean(to_vector(y_rep_w13));
  real<lower = 0> sd_y_rep_w13 = sd(to_vector(y_rep_w13));
  vector[n14] logit_p_w14 = b[1] + b[2]*male_w14  + a_race[race_w14] + a_age[age_w14] + (male_w14 + offset) .* a_sexage[sex_age_w14] + a_time[time_w14] + a_county[county_w14];
  int<lower = 0> y_rep_w14[n14] = bernoulli_logit_rng(logit_p_w14*sens[1] + (1-logit_p_w14)*(1-spec[1]));
  real<lower = 0> mean_y_rep_w14 = mean(to_vector(y_rep_w14));
  real<lower = 0> sd_y_rep_w14 = sd(to_vector(y_rep_w14));
  vector[n15] logit_p_w15 = b[1] + b[2]*male_w15  + a_race[race_w15] + a_age[age_w15] + (male_w15 + offset) .* a_sexage[sex_age_w15] + a_time[time_w15] + a_county[county_w15];
  int<lower = 0> y_rep_w15[n15] = bernoulli_logit_rng(logit_p_w15*sens[1] + (1-logit_p_w15)*(1-spec[1]));
  real<lower = 0> mean_y_rep_w15 = mean(to_vector(y_rep_w15));
  real<lower = 0> sd_y_rep_w15 = sd(to_vector(y_rep_w15));
  vector[n16] logit_p_w16 = b[1] + b[2]*male_w16  + a_race[race_w16] + a_age[age_w16] + (male_w16 + offset) .* a_sexage[sex_age_w16] + a_time[time_w16] + a_county[county_w16];
  int<lower = 0> y_rep_w16[n16] = bernoulli_logit_rng(logit_p_w16*sens[1] + (1-logit_p_w16)*(1-spec[1]));
  real<lower = 0> mean_y_rep_w16 = mean(to_vector(y_rep_w16));
  real<lower = 0> sd_y_rep_w16 = sd(to_vector(y_rep_w16));
  vector[n17] logit_p_w17 = b[1] + b[2]*male_w17  + a_race[race_w17] + a_age[age_w17] + (male_w17 + offset) .* a_sexage[sex_age_w17] + a_time[time_w17] + a_county[county_w17];
  int<lower = 0> y_rep_w17[n17] = bernoulli_logit_rng(logit_p_w17*sens[1] + (1-logit_p_w17)*(1-spec[1]));
  real<lower = 0> mean_y_rep_w17 = mean(to_vector(y_rep_w17));
  real<lower = 0> sd_y_rep_w17 = sd(to_vector(y_rep_w17));
  vector[n18] logit_p_w18 = b[1] + b[2]*male_w18  + a_race[race_w18] + a_age[age_w18] + (male_w18 + offset) .* a_sexage[sex_age_w18] + a_time[time_w18] + a_county[county_w18];
  int<lower = 0> y_rep_w18[n18] = bernoulli_logit_rng(logit_p_w18*sens[1] + (1-logit_p_w18)*(1-spec[1]));
  real<lower = 0> mean_y_rep_w18 = mean(to_vector(y_rep_w18));
  real<lower = 0> sd_y_rep_w18 = sd(to_vector(y_rep_w18));
  vector[n19] logit_p_w19 = b[1] + b[2]*male_w19  + a_race[race_w19] + a_age[age_w19] + (male_w19 + offset) .* a_sexage[sex_age_w19] + a_time[time_w19] + a_county[county_w19];
  int<lower = 0> y_rep_w19[n19] = bernoulli_logit_rng(logit_p_w19*sens[1] + (1-logit_p_w19)*(1-spec[1]));
  real<lower = 0> mean_y_rep_w19 = mean(to_vector(y_rep_w19));
  real<lower = 0> sd_y_rep_w19 = sd(to_vector(y_rep_w19));
  vector[n20] logit_p_w20 = b[1] + b[2]*male_w20  + a_race[race_w20] + a_age[age_w20] + (male_w20 + offset) .* a_sexage[sex_age_w20] + a_time[time_w20] + a_county[county_w20];
  int<lower = 0> y_rep_w20[n20] = bernoulli_logit_rng(logit_p_w20*sens[1] + (1-logit_p_w20)*(1-spec[1]));
  real<lower = 0> mean_y_rep_w20 = mean(to_vector(y_rep_w20));
  real<lower = 0> sd_y_rep_w20 = sd(to_vector(y_rep_w20));
  vector[n21] logit_p_w21 = b[1] + b[2]*male_w21  + a_race[race_w21] + a_age[age_w21] + (male_w21 + offset) .* a_sexage[sex_age_w21] + a_time[time_w21] + a_county[county_w21];
  int<lower = 0> y_rep_w21[n21] = bernoulli_logit_rng(logit_p_w21*sens[1] + (1-logit_p_w21)*(1-spec[1]));
  real<lower = 0> mean_y_rep_w21 = mean(to_vector(y_rep_w21));
  real<lower = 0> sd_y_rep_w21 = sd(to_vector(y_rep_w21));
  vector[n22] logit_p_w22 = b[1] + b[2]*male_w22  + a_race[race_w22] + a_age[age_w22] + (male_w22 + offset) .* a_sexage[sex_age_w22] + a_time[time_w22] + a_county[county_w22];
  int<lower = 0> y_rep_w22[n22] = bernoulli_logit_rng(logit_p_w22*sens[1] + (1-logit_p_w22)*(1-spec[1]));
  real<lower = 0> mean_y_rep_w22 = mean(to_vector(y_rep_w22));
  real<lower = 0> sd_y_rep_w22 = sd(to_vector(y_rep_w22));
  vector[n23] logit_p_w23 = b[1] + b[2]*male_w23  + a_race[race_w23] + a_age[age_w23] + (male_w23 + offset) .* a_sexage[sex_age_w23] + a_time[time_w23] + a_county[county_w23];
  int<lower = 0> y_rep_w23[n23] = bernoulli_logit_rng(logit_p_w23*sens[1] + (1-logit_p_w23)*(1-spec[1]));
  real<lower = 0> mean_y_rep_w23 = mean(to_vector(y_rep_w23));
  real<lower = 0> sd_y_rep_w23 = sd(to_vector(y_rep_w23));
  vector[n24] logit_p_w24 = b[1] + b[2]*male_w24  + a_race[race_w24] + a_age[age_w24] + (male_w24 + offset) .* a_sexage[sex_age_w24] + a_time[time_w24] + a_county[county_w24];
  int<lower = 0> y_rep_w24[n24] = bernoulli_logit_rng(logit_p_w24*sens[1] + (1-logit_p_w24)*(1-spec[1]));
  real<lower = 0> mean_y_rep_w24 = mean(to_vector(y_rep_w24));
  real<lower = 0> sd_y_rep_w24 = sd(to_vector(y_rep_w24));
  vector[n25] logit_p_w25 = b[1] + b[2]*male_w25  + a_race[race_w25] + a_age[age_w25] + (male_w25 + offset) .* a_sexage[sex_age_w25] + a_time[time_w25] + a_county[county_w25];
  int<lower = 0> y_rep_w25[n25] = bernoulli_logit_rng(logit_p_w25*sens[1] + (1-logit_p_w25)*(1-spec[1]));
  real<lower = 0> mean_y_rep_w25 = mean(to_vector(y_rep_w25));
  real<lower = 0> sd_y_rep_w25 = sd(to_vector(y_rep_w25));
  vector[n26] logit_p_w26 = b[1] + b[2]*male_w26  + a_race[race_w26] + a_age[age_w26] + (male_w26 + offset) .* a_sexage[sex_age_w26] + a_time[time_w26] + a_county[county_w26];
  int<lower = 0> y_rep_w26[n26] = bernoulli_logit_rng(logit_p_w26*sens[1] + (1-logit_p_w26)*(1-spec[1]));
  real<lower = 0> mean_y_rep_w26 = mean(to_vector(y_rep_w26));
  real<lower = 0> sd_y_rep_w26 = sd(to_vector(y_rep_w26));
  vector[n27] logit_p_w27 = b[1] + b[2]*male_w27  + a_race[race_w27] + a_age[age_w27] + (male_w27 + offset) .* a_sexage[sex_age_w27] + a_time[time_w27] + a_county[county_w27];
  int<lower = 0> y_rep_w27[n27] = bernoulli_logit_rng(logit_p_w27*sens[1] + (1-logit_p_w27)*(1-spec[1]));
  real<lower = 0> mean_y_rep_w27 = mean(to_vector(y_rep_w27));
  real<lower = 0> sd_y_rep_w27 = sd(to_vector(y_rep_w27));
  vector[n28] logit_p_w28 = b[1] + b[2]*male_w28  + a_race[race_w28] + a_age[age_w28] + (male_w28 + offset) .* a_sexage[sex_age_w28] + a_time[time_w28] + a_county[county_w28];
  int<lower = 0> y_rep_w28[n28] = bernoulli_logit_rng(logit_p_w28*sens[1] + (1-logit_p_w28)*(1-spec[1]));
  real<lower = 0> mean_y_rep_w28 = mean(to_vector(y_rep_w28));
  real<lower = 0> sd_y_rep_w28 = sd(to_vector(y_rep_w28));
//Marginal prevalence
  vector[J_time] p_male_h;
  vector[J_time] p_female_h;
  vector[J_time] p_0_17_h;
  vector[J_time] p_18_39_h;
  vector[J_time] p_40_64_h;
  vector[J_time] p_65_75_h;
  vector[J_time] p_75_h;
  vector[J_time] p_white_h;
  vector[J_time] p_black_h;
  vector[J_time] p_other_h;
  vector[J_time] p_male_c;
  vector[J_time] p_female_c;
  vector[J_time] p_0_17_c;
  vector[J_time] p_18_39_c;
  vector[J_time] p_40_64_c;
  vector[J_time] p_65_75_c;
  vector[J_time] p_75_c;
  vector[J_time] p_white_c;
  vector[J_time] p_black_c;
  vector[J_time] p_other_c;
// Further postratify
  vector[J] sex_record; 
  vector[J] age_record; 
  vector[J] ethnicity_record;
  int position_male;
  int position_female; 
  int position_0_17;
  int position_18_39; 
  int position_40_64;
  int position_65_75;
  int position_75;
  int position_white;
  int position_black;
  int position_other;
// Sex
  vector[2*3*5] p_pop_male;
  vector[2*3*5] p_pop_female;
  vector[2*3*5] N_male_h;
  vector[2*3*5] N_female_h;
  vector[2*3*5] N_male_c;
  vector[2*3*5] N_female_c;
// Age
  vector[2*2*3] p_pop_0_17;
  vector[2*2*3] p_pop_18_39;
  vector[2*2*3] p_pop_40_64;
  vector[2*2*3] p_pop_65_75;
  vector[2*2*3] p_pop_75;
  vector[2*2*3] N_0_17_h;
  vector[2*2*3] N_18_39_h;
  vector[2*2*3] N_40_64_h;
  vector[2*2*3] N_65_75_h;
  vector[2*2*3] N_75_h;
  vector[2*2*3] N_0_17_c;
  vector[2*2*3] N_18_39_c;
  vector[2*2*3] N_40_64_c;
  vector[2*2*3] N_65_75_c;
  vector[2*2*3] N_75_c;
// Ethnicity
  vector[2*2*5] p_pop_white;
  vector[2*2*5] p_pop_black;
  vector[2*2*5] p_pop_other;
  vector[2*2*5] N_white_h;
  vector[2*2*5] N_black_h;
  vector[2*2*5] N_other_h;
  vector[2*2*5] N_white_c;
  vector[2*2*5] N_black_c;
  vector[2*2*5] N_other_c;
for (time in 1:J_time){
  count = 1;
  for (i_male in 0:1) {
    for (i_age in 1:5) {
      for (i_race in 1:3) {
          for (i_county in 1:2){
          p_pop[count] = inv_logit(b[1]
                                   + b[2] * (i_male - 0.5)
                                   + a_race[i_race]
                                   + a_age[i_age]+a_time[time]+i_male*a_sexage[i_age] + a_county[i_county]);
          p_pop_notime[count] = inv_logit(b[1]
                                   + b[2] * (i_male - 0.5)
                                   + a_race[i_race]
                                   + a_age[i_age]+i_male*a_sexage[i_age] + a_county[i_county]);
          p_pop_notime_up[count] = inv_logit(b[1]
                                   + b[2] * (i_male - 0.5)
                                   + a_race[i_race]
                                   + a_age[i_age]+ sigma_time + i_male*a_sexage[i_age] + a_county[i_county]);  
          p_pop_notime_low[count] = inv_logit(b[1]
                                   + b[2] * (i_male - 0.5)
                                   + a_race[i_race]
                                   + a_age[i_age] - sigma_time + i_male*a_sexage[i_age] + a_county[i_county]);            
          sex_record[count] = i_male;
          age_record[count] = i_age;
          ethnicity_record[count] = i_race;
          count += 1;
      }
    }
  }
  }
  p_avg_h[time] = sum(N_pop .* p_pop) / sum(N_pop);
  p_avg_c[time] = sum(N_acs .* p_pop) / sum(N_acs);
// Sex  
  position_male = 1;
  position_female = 1;
  for (i in 1:J){
    if(sex_record[i] == 0){
      p_pop_male[position_male] = p_pop[i];
      N_male_h[position_male] = N_pop[i];
      N_male_c[position_male] = N_acs[i];
      position_male += 1;
     }else {
      p_pop_female[position_female] = p_pop[i];
      N_female_h[position_female] = N_pop[i];
      N_female_c[position_female] = N_acs[i];
      position_female += 1;
     }
   }
   p_male_h[time] = sum(N_male_h .* p_pop_male) / sum(N_male_h);
   p_female_h[time] = sum(N_female_h .* p_pop_female) / sum(N_female_h);
   p_male_c[time] = sum(N_male_c .* p_pop_male) / sum(N_male_c);
   p_female_c[time] = sum(N_female_c .* p_pop_female) / sum(N_female_c);
// Age
  position_0_17 = 1;
  position_18_39 = 1;
  position_40_64 = 1;
  position_65_75 = 1;
  position_75 = 1;
  for (i in 1:J){
     if(age_record[i] == 1){
       p_pop_0_17[position_0_17] = p_pop[i];
       N_0_17_h[position_0_17] = N_pop[i];
       N_0_17_c[position_0_17] = N_acs[i];
       position_0_17 += 1;
     }else if(age_record[i] == 2){
       p_pop_18_39[position_18_39] = p_pop[i];
       N_18_39_h[position_18_39] = N_pop[i];
       N_18_39_c[position_18_39] = N_acs[i];
       position_18_39 += 1;
     }else if(age_record[i] == 3){
       p_pop_40_64[position_40_64] =  p_pop[i];
       N_40_64_h[position_40_64] = N_pop[i];
       N_40_64_c[position_40_64] = N_acs[i];
       position_40_64 += 1;
     }else if(age_record[i] == 4){
       p_pop_65_75[position_65_75] = p_pop[i];
       N_65_75_h[position_65_75] = N_pop[i];
       N_65_75_c[position_65_75] = N_acs[i];
       position_65_75 += 1;
     }else {
       p_pop_75[position_75] = p_pop[i];
       N_75_h[position_75] = N_pop[i];
       N_75_c[position_75] = N_acs[i];
       position_75 += 1;
     }
   }
  p_0_17_h[time]= sum(N_0_17_h .* p_pop_0_17) / sum(N_0_17_h);
  p_18_39_h[time] = sum(N_18_39_h .* p_pop_18_39) / sum(N_18_39_h);
  p_40_64_h[time] = sum(N_40_64_h .* p_pop_40_64) / sum(N_40_64_h);
  p_65_75_h[time] = sum(N_65_75_h .*p_pop_65_75) / sum(N_65_75_h);
  p_75_h[time] = sum(N_75_h .*p_pop_75) / sum(N_75_h);
  p_0_17_c[time]= sum(N_0_17_c .* p_pop_0_17) / sum(N_0_17_c);
  p_18_39_c[time] = sum(N_18_39_c .* p_pop_18_39) / sum(N_18_39_c);
  p_40_64_c[time] = sum(N_40_64_c .* p_pop_40_64) / sum(N_40_64_c);
  p_65_75_c[time] = sum(N_65_75_c .*p_pop_65_75) / sum(N_65_75_c);
  p_75_c[time] = sum(N_75_c .*p_pop_75) / sum(N_75_c);
// Race
  position_white = 1;
  position_black = 1;
  position_other = 1;
  for (i in 1:J){
    if(ethnicity_record[i] == 1){
      p_pop_white[position_white] = p_pop[i];
      N_white_h[position_white] = N_pop[i];
      N_white_c[position_white] = N_acs[i];
      position_white += 1;
    }else if(ethnicity_record[i] == 2){
      p_pop_black[position_black] =  p_pop[i];
      N_black_h[position_black] = N_pop[i];
      N_black_c[position_black] = N_acs[i];
      position_black += 1;
    }else {
      p_pop_other[position_other] = p_pop[i];
      N_other_h[position_other] = N_pop[i];
      N_other_c[position_other] = N_acs[i];
      position_other += 1;
    }
  }
  p_white_h[time] = sum(N_white_h .* p_pop_white) / sum(N_white_h);
  p_black_h[time] = sum(N_black_h .* p_pop_black) / sum(N_black_h);
  p_other_h[time] = sum(N_other_h .* p_pop_other) / sum(N_other_h);
  p_white_c[time] = sum(N_white_c .* p_pop_white) / sum(N_white_c);
  p_black_c[time] = sum(N_black_c .* p_pop_black) / sum(N_black_c);
  p_other_c[time] = sum(N_other_c .* p_pop_other) / sum(N_other_c);
}
  p_overall_h = sum(N_pop .* p_pop_notime) / sum(N_pop);
  p_overall_c = sum(N_acs .* p_pop_notime) / sum(N_acs);
  p_overall_h_up = sum(N_pop .* p_pop_notime_up) / sum(N_pop);
  p_overall_c_up = sum(N_acs .* p_pop_notime_up) / sum(N_acs);
  p_overall_h_low = sum(N_pop .* p_pop_notime_low) / sum(N_pop);
  p_overall_c_low = sum(N_acs .* p_pop_notime_low) / sum(N_acs);
}
