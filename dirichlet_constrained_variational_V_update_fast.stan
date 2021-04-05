data {
  int<lower=0> J;   // number of data items
  int<lower=0> K;   // number of predictors
  matrix[J, K] X;   // predictor matrix
  matrix[J, K] X_squared;   // predictor matrix
  vector[J] y;      // outcome vector
  real intercept;
  vector[J] precisions; 
  real N;
  real<lower=0> alpha_zero;
}
parameters {
  simplex[K] beta;
}
model {
  beta ~ dirichlet(rep_vector(alpha_zero, K));
  target += N*precisions .* (y - (N*intercept) - 1.0) .* (X*beta);
  target += -precisions*square(N)*0.5 .* (X_squared*square(beta) + square(X*beta) - (square(X)*square(beta)));
}

