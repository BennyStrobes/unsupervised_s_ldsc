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
  for (j in 1:J) {
    target += N*precisions[j]*(y[j] - (N*intercept) - 1.0)*(X[j,:]*beta);
    target += -precisions[j]*square(N)*(X_squared[j,:]*square(beta))/2.0;
    target += -precisions[j]*square(N)*((X[j,:]*beta)*(X[j,:]*beta))/2.0;
    target += precisions[j]*square(N)*square(X[j,:])*square(beta)/2.0;
  }
}

