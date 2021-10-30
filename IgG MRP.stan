functions { 
    real g(int i_age, int time, vector a_ageweek) {
    real output_1;
    if (i_age == 1) {
      output_1 = 0;
    }
    else {
       output_1 = (i_age - 1)*a_ageweek[(i_age-1)*time];
    }
    return output_1;}
}
data {
  int<lower = 0> N;  // number of tests in the sample 
  int<lower = 0, upper = 1> y[N];  // 1 if positive, 0 if negative
  int<lower = 1> J_time; //#weeks
  vector[N] male;  // -0.5 if female, 0.5 if male. 
  int<lower = 1, upper = 3> race[N];  // 1=white, 2=black, 3=other
  int<lower = 1, upper = 5> age[N];  // 1=0-17, 2=18-39, 3=40-64, 4=65-75, 5 = 75+
  int<lower = 1, upper = J_time> week[N]; //week
  int<lower = 1, upper = 2> county[N]; //county indcators
  int<lower = 1, upper = 4*J_time> age_week[N]; //interactions between age and week
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
  real<lower = 0> coef_prior_scale_interaction;
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
  real<lower = 0> sigma_ageweek;
  vector<multiplier = sigma_race>[3] a_race;  
  vector<multiplier = sigma_age>[5] a_age; 
  vector<multiplier = sigma_county>[2] a_county;
  vector<multiplier = sigma_time>[J_time] a_time;
  vector<multiplier = sigma_ageweek>[4*J_time] a_ageweek;
}
transformed parameters { 
  vector[J_spec] spec; 
  vector[J_sens] sens;
  vector[N] a_ageweek_trans;
  spec = inv_logit(logit_spec); 
  sens = inv_logit(logit_sens);
  for (i in 1:N){
  a_ageweek_trans[i] = (age[i] - 1)*a_ageweek[age_week[i]];
  }

}
model {
  vector[N] p;
  vector[N] p_sample;
  p = inv_logit(b[1] + b[2]*male  + a_race[race] + a_age[age] + a_time[week] + a_county[county] + a_ageweek_trans); 
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
  a_ageweek ~ normal(0, sigma_ageweek);
  b[1] + b[2] * mean(male) ~ normal(intercept_prior_mean, intercept_prior_scale);
  b[2] ~ normal(0, coef_prior_scale);
  sigma_race ~ normal(0, coef_prior_scale);
  sigma_age ~ normal(0, coef_prior_scale);
  sigma_time ~ normal(0, coef_prior_scale_time);
  sigma_county ~ normal(0, coef_prior_scale);
  sigma_ageweek ~ normal(0, coef_prior_scale_interaction);
}
generated quantities {
  vector[J_time] p_avg_h;
  vector[J_time] p_avg_c;  
  vector[J] p_pop;  // population prevalence in the J poststratification cells
  vector[J] p_pop_no_time;  // population prevalence in the J poststratification cells
  real p_avg_h_no_time;
  real p_avg_c_no_time;
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
  real p_male_h_no_time;
  real p_female_h_no_time;
  real p_0_17_h_no_time;
  real p_18_39_h_no_time;
  real p_40_64_h_no_time;
  real p_65_75_h_no_time;
  real p_75_h_no_time;
  real p_white_h_no_time;
  real p_black_h_no_time;
  real p_other_h_no_time;
  real p_male_c_no_time;
  real p_female_c_no_time;
  real p_0_17_c_no_time;
  real p_18_39_c_no_time;
  real p_40_64_c_no_time;
  real p_65_75_c_no_time;
  real p_75_c_no_time;
  real p_white_c_no_time;
  real p_black_c_no_time;
  real p_other_c_no_time;
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
  vector[2*3*5] p_pop_male_no_time;
  vector[2*3*5] p_pop_female_no_time;
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
  vector[2*2*3] p_pop_0_17_no_time;
  vector[2*2*3] p_pop_18_39_no_time;
  vector[2*2*3] p_pop_40_64_no_time;
  vector[2*2*3] p_pop_65_75_no_time;
  vector[2*2*3] p_pop_75_no_time;
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
  vector[2*2*5] p_pop_white_no_time;
  vector[2*2*5] p_pop_black_no_time;
  vector[2*2*5] p_pop_other_no_time;
  vector[2*2*5] N_white_h;
  vector[2*2*5] N_black_h;
  vector[2*2*5] N_other_h;
  vector[2*2*5] N_white_c;
  vector[2*2*5] N_black_c;
  vector[2*2*5] N_other_c;
  int count;
  for (time in 1:J_time){
  count = 1;
  for (i_male in 0:1) {
    for (i_age in 1:5) {
      for (i_race in 1:3) {
          for (i_county in 1:2){
          p_pop[count] = inv_logit(b[1] + b[2] * (i_male - 0.5) + a_race[i_race] + a_age[i_age]+a_time[time]+
                           a_county[i_county] + g(i_age, time, a_ageweek));
          p_pop_no_time[count] = inv_logit(b[1]
                                   + b[2] * (i_male - 0.5)
                                   + a_race[i_race]
                                   + a_age[i_age] + a_county[i_county]);  
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
      p_pop_male_no_time[position_male] = p_pop_no_time[i];
      N_male_h[position_male] = N_pop[i];
      N_male_c[position_male] = N_acs[i];
      position_male += 1;
     }else {
      p_pop_female[position_female] = p_pop[i];
      p_pop_female_no_time[position_female] = p_pop_no_time[i];
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
       p_pop_0_17_no_time[position_0_17] = p_pop_no_time[i];
       N_0_17_h[position_0_17] = N_pop[i];
       N_0_17_c[position_0_17] = N_acs[i];
       position_0_17 += 1;
     }else if(age_record[i] == 2){
       p_pop_18_39[position_18_39] = p_pop[i];
       p_pop_18_39_no_time[position_18_39] = p_pop_no_time[i];
       N_18_39_h[position_18_39] = N_pop[i];
       N_18_39_c[position_18_39] = N_acs[i];
       position_18_39 += 1;
     }else if(age_record[i] == 3){
       p_pop_40_64[position_40_64] =  p_pop[i];
       p_pop_40_64_no_time[position_40_64] =  p_pop_no_time[i];
       N_40_64_h[position_40_64] = N_pop[i];
       N_40_64_c[position_40_64] = N_acs[i];
       position_40_64 += 1;
     }else if(age_record[i] == 4){
       p_pop_65_75[position_65_75] = p_pop[i];
       p_pop_65_75_no_time[position_65_75] = p_pop_no_time[i];
       N_65_75_h[position_65_75] = N_pop[i];
       N_65_75_c[position_65_75] = N_acs[i];
       position_65_75 += 1;
     }else {
       p_pop_75[position_75] = p_pop[i];
       p_pop_75_no_time[position_75] = p_pop_no_time[i];
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
      p_pop_white_no_time[position_white] = p_pop_no_time[i];
      N_white_h[position_white] = N_pop[i];
      N_white_c[position_white] = N_acs[i];
      position_white += 1;
    }else if(ethnicity_record[i] == 2){
      p_pop_black[position_black] =  p_pop[i];
      p_pop_black_no_time[position_black] =  p_pop_no_time[i];
      N_black_h[position_black] = N_pop[i];
      N_black_c[position_black] = N_acs[i];
      position_black += 1;
    }else {
      p_pop_other[position_other] = p_pop[i];
      p_pop_other_no_time[position_other] = p_pop_no_time[i];
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
  p_avg_h_no_time = sum(N_pop .* p_pop_no_time) / sum(N_pop);
  p_avg_c_no_time = sum(N_acs .* p_pop_no_time) / sum(N_acs);
  p_male_h_no_time = sum(N_male_h .* p_pop_male_no_time) / sum(N_male_h);
  p_female_h_no_time = sum(N_female_h .* p_pop_female_no_time) / sum(N_female_h);
  p_male_c_no_time = sum(N_male_c .* p_pop_male_no_time) / sum(N_male_c);
  p_female_c_no_time = sum(N_female_c .* p_pop_female_no_time) / sum(N_female_c);
  p_0_17_h_no_time= sum(N_0_17_h .* p_pop_0_17_no_time) / sum(N_0_17_h);
  p_18_39_h_no_time = sum(N_18_39_h .* p_pop_18_39_no_time) / sum(N_18_39_h);
  p_40_64_h_no_time = sum(N_40_64_h .* p_pop_40_64_no_time) / sum(N_40_64_h);
  p_65_75_h_no_time = sum(N_65_75_h .*p_pop_65_75_no_time) / sum(N_65_75_h);
  p_75_h_no_time = sum(N_75_h .*p_pop_75_no_time) / sum(N_75_h);
  p_0_17_c_no_time = sum(N_0_17_c .* p_pop_0_17_no_time) / sum(N_0_17_c);
  p_18_39_c_no_time = sum(N_18_39_c .* p_pop_18_39_no_time) / sum(N_18_39_c);
  p_40_64_c_no_time = sum(N_40_64_c .* p_pop_40_64_no_time) / sum(N_40_64_c);
  p_65_75_c_no_time = sum(N_65_75_c .*p_pop_65_75_no_time) / sum(N_65_75_c);
  p_75_c_no_time = sum(N_75_c .*p_pop_75_no_time) / sum(N_75_c);
  p_white_h_no_time = sum(N_white_h .* p_pop_white_no_time) / sum(N_white_h);
  p_black_h_no_time = sum(N_black_h .* p_pop_black_no_time) / sum(N_black_h);
  p_other_h_no_time = sum(N_other_h .* p_pop_other_no_time) / sum(N_other_h);
  p_white_c_no_time = sum(N_white_c .* p_pop_white_no_time) / sum(N_white_c);
  p_black_c_no_time = sum(N_black_c .* p_pop_black_no_time) / sum(N_black_c);
  p_other_c_no_time = sum(N_other_c .* p_pop_other_no_time) / sum(N_other_c);
}

