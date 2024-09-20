library(xtable)
library(numDeriv)

# Read data
data_mackerel = read.csv("SmartDots_Event_514_IMR mackerel age reading (internal).csv", sep = ",")

# Only extract expert readers 
data_exp = data_mackerel[data_mackerel$expertise == 1, ]

# Remove missing ages (set to 0)
data_exp = data_exp[data_exp$age > 0, ]

# Define maximum and minimum age
min_age = 1
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

mat = matrix(NA, ncol = 19, nrow = 19)
colnames(mat) = 1:19
rownames(mat) = 1:19

modes_all = unique(data_exp$modal_age_closest)
for (m in modes_all){
  subdat = data_exp$age[data_exp$modal_age_closest == m]
  subw = data_exp$weight[data_exp$modal_age_closest == m]
  freq = table(subdat, subw)
  ind = which(colnames(mat)%in%rownames(freq))
  tmp_vec = rep(0, 19)
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

print(xtable(mat), type = "latex", file = "M1_mackerel.tex")


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
print(xtable(fitted_matrix), type = "latex", file = "M_full_mackerel_2022.tex")


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



# Compute parametric matrix
fitted_matrix_b = matrix(0, nrow = 20, ncol = 20)

beta0 = fit_1$par[1]
beta1 = fit_1$par[2]
alpha0 = fit_1$par[3]
alpha1 = 0
phi = fit_1$par[4]

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
print(xtable(fitted_matrix_b), type = "latex", file = "M_best_mackerel_2022_modal_mean.tex")

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

fit_3 = nlminb(start = c(2, 0, 0.5, 0.2), negLik_3, lower = c(-Inf, -Inf, -Inf, 0.01), upper = c(Inf, Inf, Inf, 0.999))


# Likelihood for model without asymmetry
fit_3$objective


