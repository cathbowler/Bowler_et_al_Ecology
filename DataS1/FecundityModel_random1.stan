data{
  int<lower = 1> N;
  int Fecundity[N];
  matrix[N,5] ModelMatrix;
  int<lower = 1>P;
  int Plot[N];
}
parameters{
  real<lower = 0> epsilon[P];
  //real<lower = 0> sigma;
  real sigma_0;
  real<lower = 0> lambda;
  real alpha_intra;
  real alpha_InvGraminoid;
  real alpha_InvForb;
  real alpha_NatForb;
  real alpha_unknown;
}
transformed parameters{
  real<lower = 0> sigma;
  sigma = exp(sigma_0);
}
model{
  // create a vector of predictions
  vector[N] F_hat;
  vector[N] F_hat2;

  // set priors
  sigma_0 ~ normal(0, 1000);
  epsilon ~ gamma(sigma, sigma);
  alpha_intra ~ normal(0, 1000);
  alpha_InvGraminoid ~ normal(0, 1000);
  alpha_InvForb ~ normal(0, 1000);
  alpha_NatForb ~ normal(0, 1000);
  alpha_unknown ~ normal(0, 1000);
  lambda ~ gamma(0.001, 0.001);
  

  // implement the biological model
  for(i in 1:N){
    F_hat[i] = lambda * exp(alpha_intra * ModelMatrix[i,1] + alpha_InvGraminoid * ModelMatrix[i,2] + alpha_InvForb * ModelMatrix[i,3] + alpha_NatForb * ModelMatrix[i,4]);// + alpha_unknown * ModelMatrix[i,5]);
   F_hat2[i] = F_hat[i]*epsilon[Plot[i]];
  }

  // calculate the likelihood
  Fecundity ~ poisson(F_hat2);
}
