data {
  //int alpha;
  //int beta;
  int n;
  int d;
  //real vec_alpha[d];
  real X[d];
  
}
transformed data {
  //vector[d] Y;
  //Y[1] = 1;
  //Y[2] = 2;
}
parameters {
  real lambda;//[d];
  real Pi[d];
  //vector[n] Z;
}
model{

  //for(i in 1:n){
  //  target += log_mix(Pi, poisson_log_lpmf(X[i] | lambda[1]), poisson_log_lpmf(X[i] | lambda[2]));
  //}
  target+=  poisson(Pi | lambda);
}
