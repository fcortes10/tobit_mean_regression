options(scipen = 999)
library(rstan)
library(data.table)
library(mvtnorm)
library(truncnorm)
library(kyotil)
library(bayesplot)
library(betareg)
library(loo)
library(rootSolve)

######################### MODELO CON COVARIABLES

data("ReadingSkills", package = "betareg")
data <- ReadingSkills
data$dyslexia <- ifelse(data$dyslexia == "no", 0, 1)
y <- data$accuracy + 0.01

x <- as.matrix(data[ , -1])
x <- cbind(intercept = 1, x)
x_gamma <- x

n <- length(y)
n_cov_gamma <- ncol(x_gamma)

modelo <- '

functions {

  real tobit_lpdf(real y, real mu, real sigma, real a, real b) {
    real prob;
    real alpha;
    real beta;
    
    alpha = (a-mu)/sigma;
    beta = (b-mu)/sigma;
    
    if(y >= b){
      prob = log(1-Phi_approx(beta));
    } else {
      prob = std_normal_lpdf(((y-mu)/sigma)) - log(sigma);
    }
  
    return(prob);
  
  }
  
  vector system(vector mu, vector theta, real[] x_r, int[] x_i){
    vector[1] z;
    
    real a;
    real b;
    
    real sigma;
    real gamma;
    
    real alpha;
    real beta;
    
    real X;
    real W;

    a = x_r[1];
    b = x_r[2];
    
    sigma = theta[1];
    gamma = theta[2];
    
    alpha = (a-mu[1])/sigma;
    beta = (b-mu[1])/sigma;
    
    X = exp(log_diff_exp(log(Phi_approx(beta)), log(Phi_approx(alpha))));
    W = exp(std_normal_lpdf(beta))*(1 - exp(std_normal_lpdf(alpha) - std_normal_lpdf(beta)));
    
    z[1] = (1-Phi_approx(beta)) + mu[1]*X - sigma*W - gamma;
    
    return(z);
  
  }

}

data {
  int n;
  int n_cov_gamma;
  int n_pars;
  vector[n] y;
  matrix[n, n_cov_gamma] x_gamma;
  real a;
  real b;
}

transformed data{
  vector[1] mu_guess;
  real x_r[2];
  int x_i[0];
  
  mu_guess[1] = 0.5;
  x_r = {a, b};
}

parameters {
  vector[n_cov_gamma] beta_gamma;
  real<lower = 0> sigma;
}

transformed parameters {
  vector[n] linpred_gamma;
  vector[n] gamma;
  
  vector[n] sigma_v;
  vector[n] d1;
  vector[1] mu[n];
  vector[n_pars] theta;
  
  //vector[n] theta_m[n_pars];
  
  linpred_gamma = x_gamma * beta_gamma;
  
  for(i in 1:n){
  
    gamma[i] = inv_logit(linpred_gamma[i]);
    sigma_v[i] = sigma;
    
    theta[1] = sigma_v[i];
    theta[2] = gamma[i];
    
    mu[i] = algebra_solver(system, mu_guess, theta, x_r, x_i,
                           1e-8, 1e-5, 1e4);
  
  }
  
}

model {

  // Priors
  
  for(i in 1:n_cov_gamma) beta_gamma[i] ~ normal(0, 5);
  
  sigma ~ cauchy(0, 1);
  
  for(i in 1:n){
    
    y[i] ~ tobit(mu[i, 1], sigma, a, b);
    
  }

}

'

rstan_options(auto_write = TRUE)
mydata <- list(n = n, n_cov_gamma = n_cov_gamma, n_pars = 2, y = y, 
               x_gamma = x_gamma, a = 0, b = 1)

fit <- stan(model_code = modelo,
            data = mydata, 
            iter = 2000, 
            chains = 1, 
            pars = c('beta_gamma', 'sigma'),
            include = TRUE,
            #init = list(list(beta_gamma = c(1.7, 0.2),
            #                 beta_alpha1 = c(-1.3, -0.1),
            #                 sigma = 0.2)),
            #init = list(list(beta_gamma = c(1, 3.17, -1.2, -0.57), 
            #                 beta_alpha1 = c(-0.21, 0.97, 1.47, -1.31),
            #                 sigma = 0.18)),
            #            list(gamma = 0.7, alpha1 = 0.3, sigma = 0.3)),
            control = list(adapt_delta = 0.8,
                           max_treedepth = 10))

summary(fit)$summary

plot(fit, pars = c('beta_gamma'))
plot(fit, pars = c('sigma'))
#plot(fit)

traceplot(fit, inc_warmup = TRUE, pars = c('beta_gamma'))
traceplot(fit, inc_warmup = TRUE, pars = c('sigma'))
#traceplot(fit, inc_warmup = TRUE)
#traceplot(fit, inc_warmup = FALSE)

stan_ac(fit, pars = c('beta_gamma'))
stan_ac(fit, pars = c('sigma'))
#stan_ac(fit)

pairs(fit, pars = c('beta_gamma', 'sigma'))
pairs(fit, pars = c('beta_gamma[2]', 'beta_alpha1[3]', 'sigma'))
#pairs(fit)

fwrite(as.data.table(extract(fit)), 'C:/Users/B35201/Desktop/Tesis/Tobit/STAN_tobit_dyslexia_mean.csv')
chain <- as.matrix(fread('C:/Users/B35201/Desktop/Tesis/Tobit/STAN_tobit_dyslexia_mean.csv'))

mtobit <- function(mu, sigma, a = 0, b = 1){
  
  alpha <- (a-mu)/sigma
  beta <- (b-mu)/sigma
  
  X <- exp(logDiffExp(pnorm(beta, log.p = TRUE), pnorm(alpha, log.p = TRUE)))
  W <- dnorm(beta)*(1 - exp(dnorm(alpha, log = TRUE) - dnorm(beta, log = TRUE)))
  
  m <- mu*X - sigma*W
  
  return(m)
}

g <- function(mu = 0, sigma = 1, valor = 0, a = 0, b = 1){
  g <- pnorm((b-mu)/sigma, lower.tail = FALSE)*b + mtobit(mu = mu, sigma = sigma, a = a, b = b) - valor
  return(g)
}

likelihood.i <- function(y, mu, sigma, a = 0, b = 1){
  if(y == 1){
    l <- log(pnorm((b-mu)/sigma, lower.tail = FALSE))
  } else {
    l <- dnorm(y, mu, sigma, log = TRUE)
  }
  return(l)
}

likelihood <- function(y, mu, sigma, a = 0, b = 1){
  n <- length(y)
  l <- numeric(n)
  for(i in 1:n){
    l[i] <- likelihood.i(y = y[i], mu = mu[i], sigma = sigma)
  }
  return(sum(l))
}

ll <- numeric(nrow(chain))
for(i in 1:nrow(chain)){
  
  print(i)
  
  gamma <- as.numeric(plogis(x_gamma %*% chain[i, 1:n_cov_gamma]))
  sigma <- chain[i, "sigma"]
  
  mu <- numeric(length(gamma))  
  for(j in 1:length(mu)){
    mu[j] <- multiroot(f = g, start = mean(y),
                       sigma = sigma, valor = gamma[j])$root
  }
  
  ll[i] <- likelihood(y = y, mu = mu, sigma = sigma)
  
}

gamma_mean <- as.numeric(plogis(x_gamma %*% colMeans(chain[ , 1:n_cov_gamma])))
sigma_mean <- mean(chain[ , "sigma"])

mu_mean <- numeric(length(gamma_mean))
for(j in 1:length(mu_mean)){
  mu_mean[j] <- multiroot(f = g, start = mean(y),
                          sigma = sigma_mean, valor = gamma_mean[j])$root
}

d.bar.theta <- -2*likelihood(y = y, mu = mu_mean, sigma = sigma_mean)
d.theta <- -2*ll

2*mean(d.theta)-d.bar.theta
#-7.636287

mean(d.theta)-d.bar.theta
#3.908139

waic(matrix(ll, ncol = 1))
#-8.4



c(mean = mean(chain[ , "beta_gamma.V1"]), quantile(chain[ , "beta_gamma.V1"], probs = c(0.025, 0.975)))
c(mean = mean(chain[ , "beta_gamma.V2"]), quantile(chain[ , "beta_gamma.V2"], probs = c(0.025, 0.975)))
c(mean = mean(chain[ , "beta_gamma.V3"]), quantile(chain[ , "beta_gamma.V3"], probs = c(0.025, 0.975)))
c(mean = mean(chain[ , "sigma"]), quantile(chain[ , "sigma"], probs = c(0.025, 0.975)))





######################### MODELO CON COVARIABLES MODELANDO SIGMA

data("ReadingSkills", package = "betareg")
data <- ReadingSkills
data$dyslexia <- ifelse(data$dyslexia == "no", 0, 1)
y <- data$accuracy + 0.01

x <- as.matrix(data[ , -1])
x <- cbind(intercept = 1, x, x[ , 1]*x[ , 2])
x_gamma <- x[ , c(1, 2)]
x_sigma <- x[ , c(1, 2)]#, 3, 4)]

n <- length(y)
n_cov_gamma <- ncol(x_gamma)
n_cov_sigma <- ncol(x_sigma)

modelo <- '

functions {

  real tobit_lpdf(real y, real mu, real sigma, real a, real b) {
    real prob;
    real alpha;
    real beta;
    
    alpha = (a-mu)/sigma;
    beta = (b-mu)/sigma;
    
    if(y >= b){
      prob = log(1-Phi_approx(beta));
    } else {
      prob = std_normal_lpdf(((y-mu)/sigma)) - log(sigma);
    }
  
    return(prob);
  
  }
  
  vector system(vector mu, vector theta, real[] x_r, int[] x_i){
    vector[1] z;
    
    real a;
    real b;
    
    real sigma;
    real gamma;
    
    real alpha;
    real beta;
    
    real X;
    real W;

    a = x_r[1];
    b = x_r[2];
    
    sigma = theta[1];
    gamma = theta[2];
    
    alpha = (a-mu[1])/sigma;
    beta = (b-mu[1])/sigma;
    
    X = exp(log_diff_exp(log(Phi_approx(beta)), log(Phi_approx(alpha))));
    W = exp(std_normal_lpdf(beta))*(1 - exp(std_normal_lpdf(alpha) - std_normal_lpdf(beta)));
    
    z[1] = (1-Phi_approx(beta)) + mu[1]*X - sigma*W - gamma;
    
    return(z);
  
  }

}

data {
  int n;
  int n_cov_gamma;
  int n_cov_sigma;
  int n_pars;
  vector[n] y;
  matrix[n, n_cov_gamma] x_gamma;
  matrix[n, n_cov_sigma] x_sigma;
  real a;
  real b;
}

transformed data{
  vector[1] mu_guess;
  real x_r[2];
  int x_i[0];
  
  mu_guess[1] = 0.5;
  x_r = {a, b};
}

parameters {
  vector[n_cov_gamma] beta_gamma;
  vector[n_cov_sigma] beta_sigma;
}

transformed parameters {
  vector[n] linpred_gamma;
  vector[n] gamma;
  
  vector[n] linpred_sigma;
  vector[n] sigma;
  
  vector[1] mu[n];
  vector[n_pars] theta;
  
  linpred_gamma = x_gamma * beta_gamma;
  linpred_sigma = x_sigma * beta_sigma;
  
  for(i in 1:n){
  
    gamma[i] = inv_logit(linpred_gamma[i]);
    sigma[i] = exp(linpred_sigma[i]);
    
    theta[1] = sigma[i];
    theta[2] = gamma[i];
    
    mu[i] = algebra_solver(system, mu_guess, theta, x_r, x_i,
                           1e-10, 1e-5, 1e5);
  
  }
  
}

model {

  // Priors
  
  for(i in 1:n_cov_gamma) beta_gamma[i] ~ normal(0, 5);
  for(i in 1:n_cov_sigma) beta_sigma[i] ~ normal(0, 5);
  
  for(i in 1:n){
    
    y[i] ~ tobit(mu[i, 1], sigma[i], a, b);
    
  }

}

'

rstan_options(auto_write = TRUE)
mydata <- list(n = n, n_cov_gamma = n_cov_gamma, n_cov_sigma = n_cov_sigma, n_pars = 2, 
               y = y, x_gamma = x_gamma, x_sigma = x_sigma, a = 0, b = 1)

fit <- stan(model_code = modelo,
            data = mydata, 
            iter = 4000, 
            chains = 1, 
            pars = c('beta_gamma', 'beta_sigma'),
            include = TRUE,
            control = list(adapt_delta = 0.99,
                           max_treedepth = 10))

summary(fit)$summary

plot(fit, pars = c('beta_gamma'))
plot(fit, pars = c('beta_sigma'))
#plot(fit)

traceplot(fit, inc_warmup = TRUE, pars = c('beta_gamma'))
traceplot(fit, inc_warmup = TRUE, pars = c('beta_sigma'))
#traceplot(fit, inc_warmup = TRUE)
#traceplot(fit, inc_warmup = FALSE)

stan_ac(fit, pars = c('beta_gamma'))
stan_ac(fit, pars = c('sigma'))
#stan_ac(fit)

pairs(fit, pars = c('beta_gamma', 'sigma'))
pairs(fit, pars = c('beta_gamma[2]', 'beta_alpha1[3]', 'sigma'))
#pairs(fit)

#fwrite(as.data.table(extract(fit)), 'C:/Users/B35201/Desktop/Tesis/Tobit/STAN_tobit_dyslexia_mean_var_sig_each.csv')
#write(as.data.table(extract(fit)), 'C:/Users/B35201/Desktop/Tesis/Tobit/STAN_tobit_dyslexia_mean_var_sig_any.csv')
#chain <- as.matrix(fread('C:/Users/B35201/Desktop/Tesis/Tobit/STAN_tobit_dyslexia_mean_var_sig_each.csv'))
chain <- as.matrix(fread('C:/Users/B35201/Desktop/Tesis/Tobit/STAN_tobit_dyslexia_mean_var_sig_any.csv'))

mtobit <- function(mu, sigma, a = 0, b = 1){
  
  alpha <- (a-mu)/sigma
  beta <- (b-mu)/sigma
  
  X <- exp(logDiffExp(pnorm(beta, log.p = TRUE), pnorm(alpha, log.p = TRUE)))
  W <- dnorm(beta)*(1 - exp(dnorm(alpha, log = TRUE) - dnorm(beta, log = TRUE)))
  
  m <- mu*X - sigma*W
  
  return(m)
}

g <- function(mu = 0, sigma = 1, valor = 0, a = 0, b = 1){
  g <- pnorm((b-mu)/sigma, lower.tail = FALSE)*b + mtobit(mu = mu, sigma = sigma, a = a, b = b) - valor
  return(g)
}

likelihood.i <- function(y, mu, sigma, a = 0, b = 1){
  if(y == 1){
    l <- log(pnorm((b-mu)/sigma, lower.tail = FALSE))
  } else {
    l <- dnorm(y, mu, sigma, log = TRUE)
  }
  return(l)
}

likelihood <- function(y, mu, sigma, a = 0, b = 1){
  n <- length(y)
  l <- numeric(n)
  for(i in 1:n){
    l[i] <- likelihood.i(y = y[i], mu = mu[i], sigma = sigma[i])
  }
  return(sum(l))
}

#######################################################################################
## SE TIENE QUE HABER CORRIDO LA FUNCION DE VEROSIMILITUD, NUMERO DE COVS y matriz X ##
#######################################################################################
#calcular DIC
dic_tobit <- function(y, chain){
  ll <- numeric(nrow(chain))
  for(i in 1:nrow(chain)){
    print(i)
    gamma <- as.numeric(plogis(x_gamma %*% chain[i, 1:n_cov_gamma]))
    sigma <- as.numeric(exp(x_sigma %*% chain[i, (n_cov_gamma+1):(n_cov_gamma+n_cov_sigma)]))
    
    mu <- numeric(length(gamma))  
    for(j in 1:length(mu)){
      mu[j] <- multiroot(f = g, start = mean(y), sigma = sigma[j], valor = gamma[j])$root
    }
    ll[i] <- likelihood(y = y, mu = mu, sigma = sigma)
  }
  gamma_mean <- as.numeric(plogis(x_gamma %*% colMeans(chain[ , 1:n_cov_gamma])))
  sigma_mean <- as.numeric(exp(x_sigma %*% colMeans(chain[ , (n_cov_gamma+1):(n_cov_gamma+n_cov_sigma)])))

  mu_mean <- numeric(length(gamma_mean))
  for(j in 1:length(mu_mean)){
    mu_mean[j] <- multiroot(f = g, start = mean(y), sigma = sigma_mean[j], valor = gamma_mean[j])$root
  }
  
  d.bar.theta <- -2*likelihood(y = y, mu = mu_mean, sigma = sigma_mean)
  d.theta <- -2*ll
  
  return(list(DIC = 2*mean(d.theta)-d.bar.theta,
              P_d = mean(d.theta)-d.bar.theta))
}

#calcular DIC_3
dic_3_tobit <- function(y, chain){
  ll <- numeric(nrow(chain))
  for(i in 1:nrow(chain)){
    print(i)
    gamma <- as.numeric(plogis(x_gamma %*% chain[i, 1:n_cov_gamma]))
    sigma <- as.numeric(exp(x_sigma %*% chain[i, (n_cov_gamma+1):(n_cov_gamma+n_cov_sigma)]))
    
    mu <- numeric(length(gamma))  
    for(j in 1:length(mu)){
      mu[j] <- multiroot(f = g, start = mean(y), sigma = sigma[j], valor = gamma[j])$root
    }
    ll[i] <- likelihood(y = y, mu = mu, sigma = sigma)
  }
  d.theta <- -2*ll
  
  py <- numeric(length(y))
  for(i in 1:length(y)){
    print(i)
    ll.ind <- numeric(nrow(chain))
    for(m in 1:nrow(chain)){
      gamma <- as.numeric(plogis(x_gamma[i, ] %*% chain[m, 1:n_cov_gamma]))
      sigma <- as.numeric(exp(x_sigma[i, ] %*% chain[m, (n_cov_gamma+1):(n_cov_gamma+n_cov_sigma)]))
      
      mu <- multiroot(f = g, start = mean(y), sigma = sigma, valor = gamma)$root
      
      ll.ind[m] <- likelihood.i(y = y[i], mu = mu, sigma = sigma)
    }
    py[i] <- sum(exp(ll.ind))/length(ll.ind)
  }
  
  return(list(DIC_3 = 2*mean(d.theta) + 2*sum(log(py)),
              P_d3 = mean(d.theta) + 2*sum(log(py))))
}

#calcular waic
waic_tobit <- function(y, chain){
  lppd <- numeric(length(y))
  pwaic <- numeric(length(y))
  for(i in 1:length(y)){
    print(i)
    ll.ind <- numeric(nrow(chain))
    for(m in 1:nrow(chain)){
      gamma <- as.numeric(plogis(x_gamma[i, ] %*% chain[m, 1:n_cov_gamma]))
      sigma <- as.numeric(exp(x_sigma[i, ] %*% chain[m, (n_cov_gamma+1):(n_cov_gamma+n_cov_sigma)]))
      
      mu <- multiroot(f = g, start = mean(y), sigma = sigma, valor = gamma)$root
      
      ll.ind[m] <- likelihood.i(y = y[i], mu = mu, sigma = sigma)
    }
    
    lppd[i] <- log(sum(exp(ll.ind))/length(ll.ind))
    pwaic[i] <- var(ll.ind)
  }
  return(list(waic = -2*(sum(lppd)-sum(pwaic)),
              P_waic = sum(pwaic)))
}

#any
dic_tobit(y, chain) #DIC -26.57135 #P_d 5.285454
dic_3_tobit(y, chain) #DIC_3 -27.76892 #P_d3 4.087887
waic_tobit(y, chain) #WAIC -26.92536 #P_waic 4.509667

#each
dic_tobit(y, chain) #DIC -28.16845 #P_d 3.849792
dic_3_tobit(y, chain) #DIC_3 -29.08533 #P_d3 2.932912
waic_tobit(y, chain) #WAIC -28.56784 #P_waic 3.19166

round(c(mean = mean(chain[ , "beta_gamma.V1"]), quantile(chain[ , "beta_gamma.V1"], probs = c(0.025, 0.975))), 3)
round(c(mean = mean(chain[ , "beta_gamma.V2"]), quantile(chain[ , "beta_gamma.V2"], probs = c(0.025, 0.975))), 3)
round(c(mean = mean(chain[ , "beta_gamma.V3"]), quantile(chain[ , "beta_gamma.V3"], probs = c(0.025, 0.975))), 3)
round(c(mean = mean(chain[ , "beta_sigma.V1"]), quantile(chain[ , "beta_sigma.V1"], probs = c(0.025, 0.975))), 3)
round(c(mean = mean(chain[ , "beta_sigma.V2"]), quantile(chain[ , "beta_sigma.V2"], probs = c(0.025, 0.975))), 3)
round(c(mean = mean(chain[ , "beta_sigma.V3"]), quantile(chain[ , "beta_sigma.V3"], probs = c(0.025, 0.975))), 3)
round(c(mean = mean(chain[ , "beta_sigma.V4"]), quantile(chain[ , "beta_sigma.V4"], probs = c(0.025, 0.975))), 4)
