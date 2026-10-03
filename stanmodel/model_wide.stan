data {
  int N;
  real logr[N];
  int species[N];
  real n1[N];
  real n2[N];
  real n3[N];
  int Ngroup;
  int group[N];
}

parameters {
  real<lower=1> lambda[3,Ngroup];
  real<lower=0> alpha[3,3,Ngroup];
  real<lower=0> sigma;
}

transformed parameters{
  real expected[N];
  
  for (i in 1:N) {
    expected[i] = lambda[species[i], group[i]]/(1 + alpha[species[i],1,group[i]] * n1[i] + alpha[species[i],2,group[i]] * n2[i] + alpha[species[i],3,group[i]] * n3[i]);
  }
}

model {
  for (i in 1:3) {
    for (j in 1:Ngroup) {
      lambda[i,j] ~ normal(1, 10);
    }
  }
  
  for (i in 1:3) {
    for (j in 1:3) {
      for (k in 1:Ngroup) {
        alpha[i,j,k] ~ normal(0, 10);
      }
    }
  }
  
  sigma ~ normal(0, 1);
  
  for (i in 1:N) {
    logr[i] ~ normal(log(expected[i]), sigma);
  }
}

generated quantities {
  real<lower=0> rho[3,Ngroup];
  real<lower=0> fdiff[3,Ngroup];
  
  for (i in 1:Ngroup) {
    rho[1,i] = sqrt(alpha[1,2,i] * alpha[2,1,i]/(alpha[1,1,i] * alpha[2,2,i]));
    rho[2,i] = sqrt(alpha[1,3,i] * alpha[3,1,i]/(alpha[1,1,i] * alpha[3,3,i]));
    rho[3,i] = sqrt(alpha[2,3,i] * alpha[3,2,i]/(alpha[2,2,i] * alpha[3,3,i]));
    
    fdiff[1,i] = (lambda[2,i]-1)/(lambda[1,i]-1) * sqrt(alpha[1,2,i] * alpha[1,1,i]/(alpha[2,1,i]*alpha[2,2,i]));
    fdiff[2,i] = (lambda[3,i]-1)/(lambda[1,i]-1) * sqrt(alpha[1,3,i] * alpha[1,1,i]/(alpha[3,1,i]*alpha[3,3,i]));
    fdiff[3,i] = (lambda[3,i]-1)/(lambda[2,i]-1) * sqrt(alpha[2,3,i] * alpha[2,2,i]/(alpha[3,2,i]*alpha[3,3,i]));
  }
}
