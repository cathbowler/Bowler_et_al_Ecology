data{
  int<lower = 1> N;
  int Fecundity[N];
  matrix[N,4] ModelMatrix;
  int<lower = 1>P;
  int Plot[N];
}
parameters{
  real epsilon[P];
  real<lower = 0> sigma;
  real<lower = 0> lambda;
  real alpha_intra;
  real alpha_InvForb;
  real alpha_NatForb;
  real alpha_unknown;
  
}
model{
  // create a vector of predictions
  vector[N] F_hat;
  vector[N] F_hat2;

  // set priors
  sigma ~ gamma(0.001, 0.001);
  epsilon ~ gamma(sigma, sigma);
  alpha_intra ~ normal(0, 1000);
  alpha_unknown ~ normal(0, 1000);
  alpha_InvForb ~ normal(0, 1000);
  alpha_NatForb ~ normal(0, 1000);
  lambda ~ gamma(0.001, 0.001);

  

  // implement the biological model
  for(i in 1:N){
    F_hat[i] = lambda * exp(alpha_intra * ModelMatrix[i,1] + alpha_InvForb * ModelMatrix[i,2] + alpha_NatForb * ModelMatrix[i,3] + alpha_unknown * ModelMatrix[i,4]); // + epsilon[Plot[i]];
   F_hat2[i] = F_hat[i]*epsilon[Plot[i]];
  }

  // calculate the likelihood
  Fecundity ~ poisson(F_hat2);
}
