  estBetaParams <- function(mu_weight, sigmasq_weight) {
    alpha_weight<- ((1 - mu_weight) / (sigmasq_weight) - 1 / mu_weight) * mu_weight ^ 2
    beta_weight <- alpha_weight * (1 / mu_weight - 1)
   return(params = list(alpha_weight = alpha_weight, beta_weight = beta_weight))
  }


