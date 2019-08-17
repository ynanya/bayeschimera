

// The input data is a vector 'y' of length 'N'.
data {
  int N;
  int<lower=0> valid_num2_A_allele[N];
  int<lower=0> total_num_allele[N];
  int genotype[N];
  vector[3] theta;
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  real<lower=0, upper=1> chi; ##chimerism of donorz = chi
  simplex[3] a[N];  ## a[1]: Donor=AA, a[2]:Donor=AB, a[3]:Donor=BB

}

transformed parameters {
  vector<lower=0, upper=1>[N] estimated_freq_A_DonorAA;
  vector<lower=0, upper=1>[N] estimated_freq_A_DonorAB;
  vector<lower=0, upper=1>[N] estimated_freq_A_DonorBB;
  
  
  for(n in 1:N){

    if(genotype[n]==1){  //recipient is AA
      estimated_freq_A_DonorAA[n] = 0.99;
      estimated_freq_A_DonorAB[n] = 1-chi/2.0;
      estimated_freq_A_DonorBB[n] = 1-chi;
    }
    if(genotype[n]==2){  //recipient is AB
      estimated_freq_A_DonorAA[n] = (1+chi)/2.0;
      estimated_freq_A_DonorAB[n] = 1/2.0;
      estimated_freq_A_DonorBB[n] = (1-chi)/2.0;     
    }
    if(genotype[n]==3){  //recipient is BB
      estimated_freq_A_DonorAA[n] = chi;
      estimated_freq_A_DonorAB[n] = chi/2.0;
      estimated_freq_A_DonorBB[n] = 0.01;        
    }
  }

    
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  chi ~ beta(1,1);
  
  for(n in 1:N){
      a[n] ~ dirichlet(theta);
  }
  
  for (n in 1:N){
    
    vector[3] lp;
    lp[1] = log(a[n][1]) + binomial_lpmf(valid_num2_A_allele[n] | total_num_allele[n], estimated_freq_A_DonorAA[n]);
    lp[2] = log(a[n][2]) + binomial_lpmf(valid_num2_A_allele[n] | total_num_allele[n], estimated_freq_A_DonorAB[n]);
    lp[3] = log(a[n][3]) + binomial_lpmf(valid_num2_A_allele[n] | total_num_allele[n], estimated_freq_A_DonorBB[n]);
    target += log_sum_exp(lp);
  }
  

}

generated quantities{
  int donor_genotype[N];
  for(n in 1:N){
    if(a[n][1] == max(a[n])){
      donor_genotype[n] = 1;
    }
    if (a[n][2] == max(a[n])){
      donor_genotype[n] = 2;
    }
    if (a[n][3] == max(a[n])){
      donor_genotype[n] = 3;
    }
    
  }
  
}



