data {
  int<lower=0> P;
  int<lower=0> N;
  real<lower=-pi(), upper = pi()> y[N];
}

parameters {
  real alpha;
  real beta[P];
  real<lower=0> kappa;
}

model {
  alpha ~ normal(0, 1);
  for (j in 1:P){
    beta[j] ~ normal(0, 1);
  }
  kappa ~ inv_gamma(1, 1);
  for (t in (P + 1):N){
    real mu = alpha;
    for (j in 1:P){
      mu += beta[j] * tan(y[t - j] / 2); // tan takes values between -pi/2 and pi/2
    }
    y[t] ~ von_mises(atan(mu), kappa); // atan returns values between -pi and pi
  }
}
