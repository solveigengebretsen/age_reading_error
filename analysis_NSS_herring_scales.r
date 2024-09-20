library(xtable)
library(numDeriv)
# Read data
data_nss_scales = read.csv("SmartDots_Event_448_2022 NSS herring exchange (scales).csv", sep = ",")

# Only extract expert readers 
data_exp = data_nss_scales[data_nss_scales$expertise == 1, ]

# Remove missing ages (set to 0)
data_exp = data_exp[data_exp$age > 0, ]

# Define maximum and minimum age
min_age = 3
max_age = 12

### Compute modal age and weight for each fish ###
fish_id = unique(data_exp$FishID)
data_exp$modal_age_low = NA
data_exp$modal_age_closest = NA
data_exp$weight = 1

for (i in 1:length(fish_id)){
  subdat = data_exp[data_exp$FishID == fish_id[i], ]
  tbl = table(subdat$age)
  ages = names(tbl)
  mod = ages[which(tbl == max(tbl))]
  data_exp$modal_age_low[data_exp$FishID == fish_id[i]] = as.numeric(mod[1])
  n = length(mod)
  if(n>1){
    mean_age = mean(subdat$age)
    dist = abs(as.numeric(mod) - mean_age)
    min_dist = mod[which(dist == min(dist))]
    data_exp$modal_age_closest[data_exp$FishID == fish_id[i]] = as.numeric(min_dist[1])
    if(length(min_dist) > 1){
      nm = length(min_dist)
      data_exp$weight[data_exp$FishID == fish_id[i]]  = 1/nm
      for(j in 2:nm){
        extra = data_exp[data_exp$FishID == fish_id[i], ]
        extra$modal_age_closest = as.numeric(min_dist[j])
        extra$weight  = 1/nm
        data_exp = rbind(data_exp, extra)
      }
    }
  }
  else{
    data_exp$modal_age_closest[data_exp$FishID == fish_id[i]] = as.numeric(mod)
  }
}


# Empirical age-reading error matrix

mat = matrix(NA, ncol = 17, nrow = 17)
colnames(mat) = 3:19
rownames(mat) = 3:19

modes_all = unique(data_exp$modal_age_closest)
for (m in modes_all){
  subdat = data_exp$age[data_exp$modal_age_closest == m]
  subw = data_exp$weight[data_exp$modal_age_closest == m]
  freq = table(subdat, subw)
  ind = which(colnames(mat)%in%rownames(freq))
  tmp_vec = rep(0, 17)
  nw = length(unique(subw))
  nomin = 0
  denomin = 0
  for(i in 1:nw){
    nomin = nomin + freq[,i] * unique(subw)[i]
    denomin = denomin + sum(freq[,i])*unique(subw)[i]
  }
  tmp_vec[ind] = nomin/denomin
  mat[which(rownames(mat) == m), ] = tmp_vec
}

print(xtable(mat), type = "latex", file = "M1_NSS_scales.tex")


# Apply minimum and maximum age
data_exp$modal_age_closest[data_exp$modal_age_closest > max_age] = max_age
data_exp$modal_age_closest[data_exp$modal_age_closest < min_age] = min_age
data_exp$modal_age_low[data_exp$modal_age_low > max_age] = max_age
data_exp$modal_age_low[data_exp$modal_age_low < min_age] = min_age
data_exp$age[data_exp$age > max_age] = max_age
data_exp$age[data_exp$age < min_age] = min_age

### Full model ###

# Likelihood
lik_full = function(theta){
  beta0 = theta[1]
  beta1 = theta[2]
  alpha0 = theta[3]
  alpha1 = theta[4]
  phi = theta[5]
  true_age = data_exp$modal_age_closest # data_exp$modal_age_low ## if using the lowest mode as true age
  age = data_exp$age
  weight = data_exp$weight
  totLik = 0
  # See paper for definition of parametric matrix 
  for(i in 1:length(true_age)){
    paa = exp(beta0 + beta1 * true_age[i]) / (1 + exp(beta0 + beta1 * true_age[i]))
    qa = exp(alpha0 + alpha1 * true_age[i]) / (1 + exp(alpha0 + alpha1 * true_age[i]))
    if(true_age[i] == age[i]){
      if(true_age[i] == max_age){
        p = qa + paa - paa * qa
      } else if(true_age[i] == min_age){
        p = 1 - qa + paa * qa
      }else{
        p = paa
      }
    }
    else if(age[i] > true_age[i]){
      if(age[i] == max_age){
        p = (1 - paa) * qa * phi ^ (max_age - true_age[i] - 1)
      }else{
        p = (1 - paa) * qa * (1 - phi) * phi ^ (age[i] - 1 - true_age[i])
      }
    }
    else{
      d = true_age[i] - min_age
      p = (1 - paa) * (1 - qa) * phi ^ (true_age[i] - age[i] - 1) * (1 - phi) / (1 - phi ^ (d))
    }
    
    totLik = totLik + weight[i] * log(p)
  }
  return(totLik)
}

negLik_full = function(theta){
  -lik_full(theta)
}


fit_full = nlminb(start = c(2, 0, 0.5, 0, 0.2), negLik_full)


# Extract parameters, std. deviations and confidence intervals. 
params = fit_full$par
Sigma = solve(hessian(negLik_full, params))
sds = round(sqrt(diag(Sigma)), 7)
conf_int = cbind(fit_full$par - 1.96 * sqrt(diag(Sigma)), fit_full$par + 1.96 * sqrt(diag(Sigma)))

# Likelihood for full model
fit_full$objective


# Compute fitted matrix
fitted_matrix = matrix(0, nrow = 20, ncol = 20)
beta0 = fit_full$par[1]
beta1 = fit_full$par[2]
alpha0 = fit_full$par[3]
alpha1 = fit_full$par[4]
phi = fit_full$par[5]

for(i in min_age:max_age){
  paa = exp(beta0 + beta1 * i) / (1 + exp(beta0 + beta1*i))
  qa = exp(alpha0 + alpha1 * i) / (1 + exp(alpha0 + alpha1 * i))
  
  for(j in min_age:max_age){
    if(i==j){
      fitted_matrix[i, j] = paa
    }
    if(i == max_age & j == max_age){
      fitted_matrix[i, j] = qa + paa - qa * paa
    }
    if(i == min_age & j == min_age){
      fitted_matrix[i, j] = 1 - qa + paa * qa
    }
    if(j > i){
      if(j == max_age){
        fitted_matrix[i, j] = (1 - paa) * qa * phi ^ (max_age - i - 1)
      }else{
        fitted_matrix[i, j] = (1 - paa) * qa * (1 - phi) * phi ^ (j - 1 - i)
      }
      
    }
    if(i > j){
      d = i - min_age
      fitted_matrix[i, j] = (1 - paa) * (1 - qa) * phi ^ (i - j - 1)  * (1 - phi) / (1 - phi ^ (d))
    }
  }
}

fitted_matrix = fitted_matrix[rowSums(fitted_matrix) != 0, rowSums(fitted_matrix) != 0]
rownames(fitted_matrix) = min_age:max_age
colnames(fitted_matrix) = min_age:max_age
print(xtable(fitted_matrix), type = "latex", file = "M_full_scales_2022.tex")


### Model without age dependency ###

lik_0 = function(theta){
  beta0 = theta[1]
  beta1 = 0
  alpha0 = theta[2]
  alpha1 = 0
  phi = theta[3]
  true_age = data_exp$modal_age_closest # data_exp$modal_age_low ## if using the lowest mode as true age
  age = data_exp$age
  weight = data_exp$weight
  totLik = 0
  for(i in 1:length(true_age)){
    paa = exp(beta0 + beta1 * true_age[i]) / (1 + exp(beta0 + beta1 * true_age[i]))
    qa = exp(alpha0 + alpha1 * true_age[i]) / (1 + exp(alpha0 + alpha1 * true_age[i]))
    if(true_age[i] == age[i]){
      if(true_age[i] == max_age){
        p = qa + paa - paa * qa
      } else if(true_age[i] == min_age){
        p = 1 - qa + paa * qa
      }else{
        p = paa
      }
    }
    else if(age[i] > true_age[i]){
      if(age[i] == max_age){
        p = (1 - paa) * qa * phi ^ (max_age - true_age[i] - 1)
      }else{
        p = (1 - paa) * qa * (1 - phi) * phi ^ (age[i] - 1 - true_age[i])
      }
    }
    else{
      d = true_age[i] - min_age
      p = (1 - paa) * (1 - qa) * phi ^ (true_age[i] - age[i] - 1) * (1 - phi) / (1 - phi ^ (d))
    }
    totLik = totLik + weight[i] * log(p)
  }
  return(totLik)
}

negLik_0 = function(theta){
  -lik_0(theta)
}

fit_0 = nlminb(start = c(2, 0.5, 0.2), negLik_0, lower = c(-Inf, -Inf, 0.01), upper = c(Inf, Inf, 0.999))



# Likelihood for model without age dependency
fit_0$objective





### Model with age dependency on the diagonal ###

lik_1 = function(theta){
  beta0 = theta[1]
  beta1 = theta[2]
  alpha0 = theta[3]
  alpha1 = 0
  phi = theta[4]
  true_age = data_exp$modal_age_closest # data_exp$modal_age_low ## if using the lowest mode as true age
  age = data_exp$age
  weight = data_exp$weight
  totLik = 0
  for(i in 1:length(true_age)){
    paa = exp(beta0 + beta1 * true_age[i]) / (1 + exp(beta0 + beta1 * true_age[i]))
    qa = exp(alpha0 + alpha1 * true_age[i]) / (1 + exp(alpha0 + alpha1 * true_age[i]))
    if(true_age[i] == age[i]){
      if(true_age[i] == max_age){
        p = qa + paa - paa * qa
      } else if(true_age[i] == min_age){
        p = 1 - qa + paa * qa
      }else{
        p = paa
      }
    }
    else if(age[i] > true_age[i]){
      if(age[i] == max_age){
        p = (1 - paa) * qa * phi ^ (max_age - true_age[i] - 1)
      }else{
        p = (1 - paa) * qa * (1 - phi) * phi ^ (age[i] - 1 - true_age[i])
      }
    }
    else{
      d = true_age[i] - min_age
      p = (1 - paa) * (1 - qa) * phi ^ (true_age[i] - age[i] - 1) * (1 - phi) / (1 - phi ^ (d))
    }
    totLik = totLik + weight[i] * log(p)
  }
  return(totLik)
}

negLik_1 = function(theta){
  -lik_1(theta)
}

fit_1 = nlminb(start = c(2, 0, 0.5, 0.2), negLik_1, lower = c(-Inf, -Inf, -Inf, 0.01), upper = c(Inf, Inf, Inf, 0.999))

# Extract parameters, standard deviations, confidence intervals and likelihood
params = fit_1$par
Sigma = solve(hessian(negLik_1, params))
sds = round(sqrt(diag(Sigma)), 7)
conf_int = cbind(fit_1$par - 1.96 * sqrt(diag(Sigma)), fit_1$par + 1.96 * sqrt(diag(Sigma)))
likelihood = fit_1$objective





### Model without asymmetry ###

lik_3 = function(theta){
  beta0 = theta[1]
  beta1 = theta[2]
  alpha0 = 0
  alpha1 = 0
  phi = theta[3]
  true_age = data_exp$modal_age_closest # data_exp$modal_age_low ## if using the lowest mode as true age
  age = data_exp$age
  weight = data_exp$weight
  totLik = 0
  for(i in 1:length(true_age)){
    paa = exp(beta0 + beta1 * true_age[i]) / (1 + exp(beta0 + beta1 * true_age[i]))
    qa = exp(alpha0 + alpha1 * true_age[i]) / (1 + exp(alpha0 + alpha1 * true_age[i]))
    if(true_age[i] == age[i]){
      if(true_age[i] == max_age){
        p = qa + paa - paa * qa
      } else if(true_age[i] == min_age){
        p = 1 - qa + paa * qa
      }else{
        p = paa
      }
    }
    else if(age[i] > true_age[i]){
      if(age[i] == max_age){
        p = (1 - paa) * qa * phi ^ (max_age - true_age[i] - 1)
      }else{
        p = (1 - paa) * qa * (1 - phi) * phi ^ (age[i] - 1 - true_age[i])
      }
    }
    else{
      d = true_age[i] - min_age
      p = (1 - paa) * (1 - qa) * phi ^ (true_age[i] - age[i] - 1) * (1 - phi) / (1 - phi ^ (d))
    }
    totLik = totLik + weight[i] * log(p)
  }
  return(totLik)
}

negLik_3 = function(theta){
  -lik_3(theta)
}

fit_3 = nlminb(start = c(2, 0, 0.2), negLik_3, lower = c(-Inf, -Inf, 0.01), upper = c(Inf, Inf, 0.999))


# Likelihood for model without asymmetry
fit_3$objective

params = fit_3$par
Sigma = solve(hessian(negLik_3, params))
round(sqrt(diag(Sigma)), 7)
cbind(fit_3$par - 1.96 * sqrt(diag(Sigma)), fit_3$par + 1.96 * sqrt(diag(Sigma)))



fitted_matrix_b = matrix(0, nrow = 20, ncol = 20)

beta0 = fit_3$par[1]
beta1 = fit_3$par[2]
alpha0 = 0
alpha1 = 0
phi = fit_3$par[3]

for(i in min_age:max_age){
  paa = exp(beta0 + beta1 * i) / (1 + exp(beta0 + beta1 * i))
  qa = exp(alpha0 + alpha1 * i) / (1 + exp(alpha0 + alpha1 * i))
  
  for(j in min_age:max_age){
    if(i == j){
      fitted_matrix_b[i, j] = paa
    }
    if(i == max_age & j == max_age){
      fitted_matrix_b[i, j] = qa + paa - qa * paa
    }
    if(i == min_age & j == min_age){
      fitted_matrix_b[i, j] = 1 - qa + paa * qa
    }
    if(j > i){
      if(j == max_age){
        fitted_matrix_b[i, j] = (1 - paa) * qa * phi ^ (max_age - i - 1)
      }else{
        fitted_matrix_b[i, j] = (1 - paa) * qa * (1 - phi) * phi ^ (j - 1 - i)
      }
      
    }
    if(i > j){
      d = i - min_age
      fitted_matrix_b[i, j] = (1 - paa) * (1 - qa) * phi ^ (i - j - 1)  * (1 - phi) / (1 - phi ^ (d))
    }
  }
}

fitted_matrix_b = fitted_matrix_b[rowSums(fitted_matrix_b) != 0, rowSums(fitted_matrix_b) != 0]
rownames(fitted_matrix_b) = min_age:max_age
colnames(fitted_matrix_b) = min_age:max_age
print(xtable(fitted_matrix_b), type = "latex", file = "M_best_NSS_scales_2022_modal_mean.tex")

# Symmetric truth assumed in simulations
fitted_matrix_symmetric = fitted_matrix_b


### Model with quadratic age dependency on diagonal  ###

lik_4 = function(theta){
  beta0 = theta[1]
  beta1 = theta[2]
  alpha0 = theta[3]
  alpha1 = 0
  phi = theta[4]
  beta2 = theta[5]
  true_age = data_exp$modal_age_closest # data_exp$modal_age_low ## if using the lowest mode as true age
  age = data_exp$age
  weight = data_exp$weight
  totLik = 0
  for(i in 1:length(true_age)){
    paa = exp(beta0 + beta1 * true_age[i] + beta2 * true_age[i]^2) / (1 + exp(beta0 + beta1 * true_age[i]+ beta2 * true_age[i]^2))
    qa = exp(alpha0 + alpha1 * true_age[i]) / (1 + exp(alpha0 + alpha1 * true_age[i]))
    if(true_age[i] == age[i]){
      if(true_age[i] == max_age){
        p = qa + paa - paa * qa
      } else if(true_age[i] == min_age){
        p = 1 - qa + paa * qa
      }else{
        p = paa
      }
    }
    else if(age[i] > true_age[i]){
      if(age[i] == max_age){
        p = (1 - paa) * qa * phi ^ (max_age - true_age[i] - 1)
      }else{
        p = (1 - paa) * qa * (1 - phi) * phi ^ (age[i] - 1 - true_age[i])
      }
    }
    else{
      d = true_age[i] - min_age
      p = (1 - paa) * (1 - qa) * phi ^ (true_age[i] - age[i] - 1) * (1 - phi) / (1 - phi ^ (d))
    }
    totLik = totLik + weight[i] * log(p)
  }
  return(totLik)
}

negLik_4 = function(theta){
  -lik_4(theta)
}

fit_4 = nlminb(start = c(fit_1$par, 0), negLik_4)


fit_4$objective





## Simulation study, symmetric truth, fitting asymmetric model ##
lik_simulation = function(theta){
  beta0 = theta[1]
  beta1 = theta[2]
  alpha0 = theta[3]
  alpha1 = 0
  phi = theta[4]
  totLik = 0
  for(i in 1:length(true_age)){
    paa = exp(beta0 + beta1 * true_age[i]) / (1 + exp(beta0 + beta1*true_age[i]))
    qa = exp(alpha0 + alpha1 * true_age[i]) / (1 + exp(alpha0 + alpha1 * true_age[i]))
    if(true_age[i] == age[i]){
      if(true_age[i] == max_age){
        p = qa + paa - paa * qa
      } else if(true_age[i] == min_age){
        p = 1 - qa + paa * qa
      }else{
        p = paa
      }
    }
    else if(age[i] > true_age[i]){
      if(age[i] == max_age){
        p = (1 - paa) * qa * phi ^ (max_age - true_age[i] - 1)
      }else{
        p = (1 - paa) * qa * (1 - phi) * phi ^ (age[i] - 1 - true_age[i])
      }
    }
    else{
      d = true_age[i] - min_age
      p = (1 - paa) * (1 - qa) * phi ^ (true_age[i] - age[i] - 1) * (1 - phi) / (1 - phi ^ (d))
    }
    totLik = totLik + weight[i] * log(p)
  }
  return(totLik)
}

negLik_simulation = function(theta){
  -lik_simulation(theta)
}


### LOWEST MODAL AGE ###

fish_id = unique(data_exp$FishID)

B = 1000
mean_diff = c()
pars = c()
fitted_matrix_boot = matrix(0, nrow = 20, ncol = 20)
for (b in 1:B){
  true_ages = c()
  ages = c()
  df = data.frame()
  for (i in 1:length(fish_id)){
    subdata = data_exp[data_exp$FishID==fish_id[i], ]
    true_age = subdata$modal_age_low[1]
    ind = which(rownames(fitted_matrix_symmetric) == true_age)
    age_s = sample(colnames(fitted_matrix_symmetric), prob = fitted_matrix_symmetric[ind, ], size = dim(subdata)[1], replace = TRUE)
    true_ages = c(true_ages, true_age)
    age = min(as.numeric(names(table(age_s))[table(age_s)==max(table(age_s))]))
    ages = c(ages, age)
    df_sub = data.frame(true_age = age, age = age_s)
    df = rbind(df, df_sub)
  }
  mean_diff = c(mean_diff, mean(as.numeric(ages) != true_ages))
  true_age = as.numeric(df$true_age)
  age = as.numeric(df$age)
  df$weight = 1
  weight = as.numeric(df$weight)
  fit = nlminb(start = c(2, 0, 0, 0.2), negLik_simulation, lower = c(-Inf, -Inf, -Inf, 0.01), upper = c(Inf, Inf, Inf, 0.999))
  pars = rbind(pars, fit$par)

  beta0 = fit$par[1]
  beta1 = fit$par[2]
  alpha1 = 0
  alpha0 = fit$par[3]
  phi = fit$par[4]

  for(i in min_age:max_age){
    paa = exp(beta0 + beta1 * i) / (1 + exp(beta0 + beta1 * i))
    qa = exp(alpha0 + alpha1 * i) / (1 + exp(alpha0 + alpha1 * i))

    for(j in min_age:max_age){
      if(i==j & i!=min_age & i!=max_age){
        fitted_matrix_boot[i, j] = fitted_matrix_boot[i, j] + paa
      }
      if(i == max_age & j == max_age){
        fitted_matrix_boot[i, j] = fitted_matrix_boot[i, j] + qa + paa - qa * paa
      }
      if(i == min_age & j == min_age){
        fitted_matrix_boot[i, j] = fitted_matrix_boot[i, j] + 1 - qa + paa * qa
      }
      if(j > i){
        if(j == max_age){
          fitted_matrix_boot[i, j] = fitted_matrix_boot[i, j] + (1 - paa) * qa * phi ^ (max_age - i - 1)
        }else{
          fitted_matrix_boot[i, j] = fitted_matrix_boot[i, j] + (1 - paa) * qa * (1 - phi) * phi ^ (j - 1 - i)
        }

      }
      if(i > j){
        d = i - min_age
        fitted_matrix_boot[i, j] = fitted_matrix_boot[i, j] + (1 - paa) * (1 - qa) * phi ^ (i - j - 1)  * (1 - phi) / (1 - phi ^ (d))
      }
    }
  }
}

fitted_matrix_boot = fitted_matrix_boot[rowSums(fitted_matrix_boot) != 0, rowSums(fitted_matrix_boot) != 0] / B
rownames(fitted_matrix_boot) = min_age:max_age
colnames(fitted_matrix_boot) = min_age:max_age

print(xtable(fitted_matrix_boot), type = "latex", file = "M_boot_symmetric_truth_modallow.tex")

difference = fitted_matrix_boot - fitted_matrix_symmetric

print(xtable(round(difference, 5), digits = 3), type = "latex", file = "M_difference_symmetric_truth_modallow.tex")



### MODAL AGE CLOSEST TO MEAN ###
fish_id = unique(data_exp$FishID)

B = 500
mean_diff = c()
pars = c()
fitted_matrix_boot = matrix(0, nrow = 20, ncol = 20)
for (b in 1:B){
  true_ages = c()
  ages = c()
  df = data.frame()
  for (i in 1:length(fish_id)){
    subdata = data_exp[data_exp$FishID==fish_id[i], ]
    true_age = subdata$modal_age_low[1]
    ind = which(rownames(fitted_matrix_symmetric) == true_age)
    age_s = sample(colnames(fitted_matrix_symmetric), prob = fitted_matrix_symmetric[ind, ], size = dim(subdata)[1], replace = TRUE)
    true_ages = c(true_ages, true_age)
    tbl = table(age_s)
    ages_names = names(tbl)
    mod = ages_names[which(tbl == max(tbl))]
    n = length(mod)
    if(n>1){
      mean_age = mean(as.numeric(age_s))
      dist = abs(as.numeric(mod) - mean_age)
      min_dist = mod[which(dist == min(dist))]
      age = min_dist
      weight = 1/length(min_dist)
      for(j in age){
        df_sub = data.frame(true_age = j, age = age_s, weight = weight)
        df = rbind(df, df_sub)
      }
    }
    else{
      age = as.numeric(mod)
      weight = 1
      df_sub = data.frame(true_age = age, age = age_s, weight = weight)
      df = rbind(df, df_sub)
    }
  }
  true_age = as.numeric(df$true_age)
  age = as.numeric(df$age)
  weight = as.numeric(df$weight)
  fit = nlminb(start = c(2, 0, 0, 0.2), negLik_simulation, lower = c(-Inf, -Inf, -Inf, 0.01), upper = c(Inf, Inf, Inf, 0.999))
  pars = rbind(pars, fit$par)

  beta0 = fit$par[1]
  beta1 = fit$par[2]
  alpha0 = fit$par[3]
  alpha1 = 0
  phi = fit$par[4]

  for(i in min_age:max_age){
    paa = exp(beta0 + beta1 * i) / (1 + exp(beta0 + beta1*i))
    qa = exp(alpha0 + alpha1 * i) / (1 + exp(alpha0 + alpha1 * i))

    for(j in min_age:max_age){
      if(i==j & i != min_age & i != max_age){
        fitted_matrix_boot[i, j] = fitted_matrix_boot[i, j] + paa
      }
      if(i == max_age & j == max_age){
        fitted_matrix_boot[i, j] = fitted_matrix_boot[i, j] + qa + paa - qa * paa
      }
      if(i == min_age & j == min_age){
        fitted_matrix_boot[i, j] = fitted_matrix_boot[i, j] + 1 - qa + paa * qa
      }
      if(j > i){
        if(j == max_age){
          fitted_matrix_boot[i, j] = fitted_matrix_boot[i, j] + (1 - paa) * qa * phi ^ (max_age - i - 1)
        }else{
          fitted_matrix_boot[i, j] = fitted_matrix_boot[i, j] + (1 - paa) * qa * (1 - phi) * phi ^ (j - 1 - i)
        }

      }
      if(i > j){
        d = i - min_age
        fitted_matrix_boot[i, j] = fitted_matrix_boot[i, j] + (1 - paa) * (1 - qa) * phi ^ (i - j - 1)  * (1 - phi) / (1 - phi ^ (d))
      }
    }
  }
}

fitted_matrix_boot = fitted_matrix_boot[rowSums(fitted_matrix_boot) != 0, rowSums(fitted_matrix_boot) != 0] / B
rownames(fitted_matrix_boot) = min_age:max_age
colnames(fitted_matrix_boot) = min_age:max_age

print(xtable(fitted_matrix_boot), type = "latex", file = "M_boot_symmetric_truth_modalmean.tex")

difference = fitted_matrix_boot - fitted_matrix_symmetric

print(xtable(round(difference, 5), digits = 3), type = "latex", file = "M_difference_symmetric_truth_modalmean.tex")
