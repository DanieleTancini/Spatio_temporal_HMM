SS=1 

for(ss in 1:SS){
  #Packages
  library(ngspatial)
  library(igraph)
  library(LaplacesDemon)
  library(mvtnorm)
  library(bvartools)
  
  #Lattice generator
  generate_adjacency_matrix = function(n){
    adj_matrix = matrix(0, n^2, n^2)
    
    directions <- list(
      up = c(-1, 0),
      down = c(1, 0),
      left = c(0, -1),
      right = c(0, 1)
    )
    
    get_index = function(row, col){
      if (row < 1 || row > n || col < 1 || col > n){
        return(NULL)
      }
      return((row - 1) * n + col)
    }
    
    for (row in 1:n){
      for (col in 1:n){
        index = get_index(row, col)
        for (dir in directions) {
          neighbor_row = row + dir[1]
          neighbor_col = col + dir[2]
          neighbor_index = get_index(neighbor_row, neighbor_col)
          
          if (!is.null(neighbor_index)){
            adj_matrix[index, neighbor_index] = 1
          }
        }
      }
    }
    
    return(adj_matrix)
  }
  
  N = 3
  Com = adj_matrix = generate_adjacency_matrix(N)
  C = matrix(0, N^2, N^2)
  C[upper.tri(C)] = Com[upper.tri(Com)]
  print(adj_matrix)
  
  #Set Working Directory
  setwd("~/Desktop/Phd_proj/Repository lavori/Spatio-temporal/Simulation_study")
  
  #Functions
  source("U_sim_phase_transition_speed.R")
  source("Omega_sampler_general_model_speed.R")
  
  #Function MCMC ratio
  logbeta = function(Betan, Beta, k, j, u){
    I = as.numeric(table(factor(u[,1], 1:k))[j])
    (Betan[j] - Beta[j])*I
  }
  
  logbetas = function(Betan, Beta, k, j, u){
    I = as.numeric(table(factor(u[,-1], 1:k))[j])
    (Betan[j] - Beta[j])*I
  }
  
  loggamma = function(Gamman, Gamma, n, k, j, l, u, C){
    I = table(factor(matrix(u[,1], n,n),1:k),factor(t(matrix(u[,1],n,n)*t(C)),1:k))
    (Gamman[j,l] - Gamma[j,l])*I[j,l]
  }
  
  loggammas = function(Gamman, Gamma, n, Ti, k, j, l, u, C){
    I = matrix(0,k,k)
    for(t in 2:Ti){
      I = I + table(factor(matrix(u[,t], n,n),1:k),factor(t(matrix(u[,t],n,n)*t(C)),1:k))
    }
    (Gamman[j,l] - Gamma[j,l])*I[j,l]
  }
  
  logdelta = function(Deltan, Delta, Ti, k, j, l, u){
    I = matrix(0,k,k)
    for(t in 2:Ti){
      I = I + table(factor(u[,t-1],1:k),factor(u[,t],1:k))
    }
    (Deltan[j,l] - Delta[j,l])*I[j,l]
  }
  
  
  #Initial settings
  n = N^2
  Ti = 5
  k = 2
  U = 1:k
  
  #Parameters synthetic dataset
  Gamma_ = matrix(c(0,1,-1,0), 2,2)
  Delta_ = matrix(c(0,-1,-1,0), 2,2)
  Beta_ = c(2,0)
  Gammas_ = matrix(c(0,1,-1,0), 2,2)
  Betas_ = c(2,0)
  iniz_C = "no random"
  post_exp = "TRUE"
  
  #Generator with warning
  cond = k-1
  counter = 0
  set.seed(123456 + (ss-1))
  while(cond!=k){
    if(counter <= 10){
      u_real = U_sim(n, Ti, k, iniz_C = iniz_C, Beta = Beta_, Betas = Betas_, Delta = Delta_, Gamma = Gamma_, Gammas = Gammas_, iter = 20000, Com = Com, post_exp = post_exp) 
      cond = dim(table(u_real$mean_u))
    }else{
      u_real=NULL
      cond=k
    }
    counter = counter + 1
  }
  
  u_real$mean_u
  u_real$C
  
  #Simulated data observed
  d = 2
  mu_real = matrix(c(-3, -3, 3, 3), d, k)
  sigma_real = array(c(1,-0.5,-0.5,1,1,-0.5,-0.5,1), dim = c(d,d,k))
  y = matrix(0, d, n*Ti)
  pos_u = u_real$mean_u 
  for(i in 1:(n*Ti)) y[,i] = rmvnorm(1, mu_real[,as.vector(pos_u[i])], sigma_real[,,as.vector(pos_u[i])])
  yf = array(y, dim = c(d,n,Ti))
  
  #Matrix neigh
  C = u_real$C
  Com = u_real$Com
  
  #Prior
  log_prior_theta = function(x, mb, sb){
    dnorm(x, mb, sb, log = TRUE)
  }
  
  #Hyperparameter
  mb = 0
  sb = sqrt(1)
  sb_b = sqrt(1)
  
  mm = rep(0,d)
  Vm = diag(100, d)
  
  half_alpha = ceiling((d+1)/2) + 1
  alpha = half_alpha*2
  S = matrix(half_alpha,d,d)
  diag(S) = alpha
  
  #Iterations
  #R = 10000
  R = 1000 #short check
  
  #Iterations Gibbs in each exchange
  R_set = rep(5, R)
  
  #Constraint
  id_con = 2
  
  #Initial settings
  set.seed(123456 + (ss-1))
  beta = Beta = Betan = c(rnorm(k-1, 0, 1),0)
  betas = Betas = Betasn = c(rnorm(k-1, 0, 1),0)
  delta = Delta = Deltan = matrix(rnorm(k*k, 0, 1), k, k)
  gamma = Gamma = Gamman = matrix(rnorm(k*k, 0, 1), k, k)
  gammas = Gammas = Gammasn = matrix(rnorm(k*k, 0, 1), k, k)
  
  mu_iniz = matrix(rnorm(d*k, 0, 100), d, k)
  sigma_iniz = array(rinvwishart(alpha,S), dim = c(d, d, k))
  
  diag(delta) = 0
  diag(Delta) = 0
  diag(Deltan) = 0
  
  diag(gamma) = 0
  diag(Gamma) = 0
  diag(Gamman) = 0
  
  diag(gammas) = 0
  diag(Gammas) = 0
  diag(Gammasn) = 0
  
  #Initial u
  u = sample(U, n*Ti, replace = TRUE)
  u = matrix(u, n, Ti)
  U_ex = array(NA, dim = c(n,Ti,R))
  U_ex[,,1] = u
  
  #Cont
  mu = array(0, dim = c(d, k, R))
  sigma = array(0, dim = c(d, d, k, R))
  dens = matrix(0, k, n)
  Beta_es = matrix(0, k, R)
  Betas_es = matrix(0, k, R)
  Delta_es = array(0, dim = c(k,k,R))
  Gamma_es = array(0, dim = c(k,k,R))
  Gammas_es = array(0, dim = c(k,k,R))
  
  #Init
  mu[,,1] = mu_iniz
  sigma[,,,1] = sigma_iniz
  Beta_es[,1] = beta
  Betas_es[,1] = betas
  Gamma_es[,,1] = gamma
  Gammas_es[,,1] = gammas
  Delta_es[,,1] = delta
  
  #Full par
  solvesigmap1 = solve(Vm)
  solvesigmap1mup1 = solve(Vm)%*%mm
  
  #Acceptance
  accb = accbs = rep(0, k-1)
  accd = matrix(0, k, k)
  accg = accgs = matrix(0, k, k)
  
  #Sigma RW
  sb_rw = c(0.8)
  sbs_rw = c(0.8)
  sg_rw = matrix(c(0,0.8,0.8,0),2,2)
  sgs_rw = matrix(c(0,0.8,0.8,0),2,2)
  sd_rw = matrix(c(0,0.8,0.8,0),2,2)
  
  #Adaptive MCMC
  stop_ = 0.5*R
  lambda_b = lambda_bs = matrix(1,k-1,stop_)
  lambda_g = lambda_gs = lambda_d = array(1, dim = c(k,k,stop_))
  mu_b = matrix(beta, k, stop_)
  mu_bs = matrix(betas, k, stop_)
  sigma_b = sigma_bs = matrix(0,k-1,stop_)
  sigma_b[,1] = sb_rw
  sigma_bs[,1] = sbs_rw
  mu_g = array(gamma, dim = c(k,k,stop_))
  mu_gs = array(gammas, dim = c(k,k,stop_))
  mu_d = array(delta, dim = c(k,k,stop_))
  sigma_g = sigma_gs = sigma_d = array(0,dim = c(k,k,stop_))
  sigma_g[,,1] = sg_rw
  sigma_gs[,,1] = sgs_rw
  sigma_d[,,1] = sd_rw
  alph_opt = 0.4
  Cost = 0.01
  
  set.seed(123456 + (ss-1))
  for(r in 2:R){
    RR = R_set[r]
    ga = Cost/r
    for(v in 1:k){
      ns = ifelse(as.vector(u)==v,1,0)
      yn = y*matrix(rep(ns,d), nrow = d, byrow = TRUE)
      mu[,v,r] = rmvnorm(1,solve(solvesigmap1+sum(ns)*solve(sigma[,,v,r-1]))%*%(solvesigmap1mup1 + solve(sigma[,,v,r-1])%*%rowSums(yn)),solve((solvesigmap1+sum(ns)*solve(sigma[,,v,r-1]))))
      sum = colSums(yn)
      pos = ifelse(sum==0,1,0)
      posiz = which(pos==1,arr.ind = TRUE)
      ynn = as.matrix(yn[,-c(posiz)])
      if(dim(ynn)[2]==0){
        Sn = S
      }else{
        Sa = array((ynn-mu[,v,r]), dim = c(d,1,ncol(ynn)))
        St = matrix(rowSums(apply(Sa, 3, function(x) x%*%t(x))), nrow = d)
        Sn = S + St
      }
      an = alpha + sum(ns)
      sigma[,,v,r] = rinvwishart(an,Sn)
    }
    
    #Identif constraint
    norm_e = mu[id_con,,r]
    ind = order(norm_e)
    if(any(ind!=(1:k))){
      mu[,,r] =  mu[,ind,r]
      sigma[,,,r] = sigma[,,ind,r]
    }
    
    #Beta
    for(j in 1:(k-1)){
      Betan[j] = Beta[j] + rnorm(1, 0, sb_rw[j]) 
      
      omega = Omega_sim(n = n, Ti = Ti, k = k, u = u, Beta = Betan, Betas = Betas, Delta = Delta, Gamma = Gamma, Gammas = Gammas, iter = RR, Com = Com)
      omega = omega$u
      
      #Exchange beta
      ratio_b = log_prior_theta(Betan[j], mb, sb_b) - log_prior_theta(Beta[j], mb, sb_b) + logbeta(Betan =  Betan, Beta = Beta, k = k, j = j, u = u) + logbeta(Betan =  Beta, Beta = Betan, k = k, j = j, u = omega)
      
      alpha_b = min(1, exp(ratio_b))
      if(runif(1)<= alpha_b){
        Beta[j] = Betan[j]
        accb[j] = accb[j] +1
      }
      Betan = Beta
      Beta_es[j,r] = Beta[j]
      if(r <= stop_){
        lambda_b[j,r] = exp(log(lambda_b[j,r-1]) + ga*(alpha_b - alph_opt))
        mu_b[j,r] = mu_b[j,r-1] + ga*(Beta[j] - mu_b[j,r-1])
        sigma_b[j,r] = sigma_b[j,r-1] + ga*((Beta[j] - mu_b[j,r-1])^2 - sigma_b[j,r-1])
        sb_rw[j] = lambda_b[j,r]*sigma_b[j,r]
      }
    }
    
    
    #Betas
    for(j in 1:(k-1)){
      Betasn[j] = Betas[j] + rnorm(1, 0, sbs_rw[j]) 
      
      omega = Omega_sim(n = n, Ti = Ti, k = k, u = u, Beta = Beta, Betas = Betasn, Delta = Delta, Gamma = Gamma, Gammas = Gammas, iter = RR, Com = Com)
      omega = omega$u
      
      ratio_b = log_prior_theta(Betasn[j], mb, sb_b) - log_prior_theta(Betas[j], mb, sb_b) + logbetas(Betan =  Betasn, Beta = Betas, k = k, j = j, u = u) + logbetas(Betan =  Betas, Beta = Betasn, k = k, j = j, u = omega)
      
      alpha_b = min(1, exp(ratio_b))
      if(runif(1)<= alpha_b){
        Betas[j] = Betasn[j]
        accbs[j] = accbs[j] +1
      }
      Betasn = Betas
      Betas_es[j,r] = Betas[j]
      if(r <= stop_){
        lambda_bs[j,r] = exp(log(lambda_bs[j,r-1]) + ga*(alpha_b - alph_opt))
        mu_bs[j,r] = mu_bs[j,r-1] + ga*(Betas[j] - mu_bs[j,r-1])
        sigma_bs[j,r] = sigma_bs[j,r-1] + ga*((Betas[j] - mu_bs[j,r-1])^2 - sigma_bs[j,r-1])
        sbs_rw[j] = lambda_bs[j,r]*sigma_bs[j,r]
      }
    }
    
    #Gamma
    for(j in 1:(k-1)){
      for(l in (j+1):k){
        Gamman[j,l] = Gamma[j,l] + rnorm(1, 0, sg_rw[j,l])
        
        omega = Omega_sim(n = n, Ti = Ti, k = k, u = u, Beta = Beta, Betas = Betas, Delta = Delta, Gamma = Gamman, Gammas = Gammas, iter = RR, Com = Com)
        omega = omega$u
        
        ratio_g = log_prior_theta(Gamman[j,l], mb, sb) - log_prior_theta(Gamma[j,l], mb, sb) + loggamma(Gamman, Gamma, n, k, j, l, u = u, C = C) + loggamma(Gamma, Gamman, n, k, j, l, u = omega, C = C)
        
        alpha_g = min(1, exp(ratio_g))
        if(runif(1)<= alpha_g){
          Gamma[j,l] = Gamman[j,l]
          accg[j,l] = accg[j,l] +1
        }
        Gamman = Gamma
        Gamma_es[j,l,r] = Gamma[j,l]
        if(r <= stop_){
          lambda_g[j,l,r] = exp(log(lambda_g[j,l,r-1]) + ga*(alpha_g - alph_opt))
          mu_g[j,l,r] = mu_g[j,l,r-1] + ga*(Gamma[j,l] - mu_g[j,l,r-1])
          sigma_g[j,l,r] = sigma_g[j,l,r-1] + ga*((Gamma[j,l] - mu_g[j,l,r-1])^2 - sigma_g[j,l,r-1])
          sg_rw[j,l] = lambda_g[j,l,r]*sigma_g[j,l,r]
        }
      }
    }
    
    for(l in 1:(k-1)){
      for(j in (l+1):k){
        Gamman[j,l] = Gamma[j,l] + rnorm(1, 0, sg_rw[j,l])
        
        omega = Omega_sim(n = n, Ti = Ti, k = k, u = u, Beta = Beta, Betas = Betas, Delta = Delta, Gamma = Gamman, Gammas = Gammas, iter = RR, Com = Com)
        omega = omega$u
        
        ratio_g = log_prior_theta(Gamman[j,l], mb, sb) - log_prior_theta(Gamma[j,l], mb, sb) + loggamma(Gamman, Gamma, n, k, j, l, u = u, C = C) + loggamma(Gamma, Gamman, n, k, j, l, u = omega, C = C)
        
        alpha_g = min(1, exp(ratio_g))
        if(runif(1)<= alpha_g){
          Gamma[j,l] = Gamman[j,l]
          accg[j,l] = accg[j,l] +1
        }
        Gamman = Gamma
        Gamma_es[j,l,r] = Gamma[j,l]
        if(r <= stop_){
          lambda_g[j,l,r] = exp(log(lambda_g[j,l,r-1]) + ga*(alpha_g - alph_opt))
          mu_g[j,l,r] = mu_g[j,l,r-1] + ga*(Gamma[j,l] - mu_g[j,l,r-1])
          sigma_g[j,l,r] = sigma_g[j,l,r-1] + ga*((Gamma[j,l] - mu_g[j,l,r-1])^2 - sigma_g[j,l,r-1])
          sg_rw[j,l] = lambda_g[j,l,r]*sigma_g[j,l,r]
        }
      }
    }
    
    #Gammas
    for(j in 1:(k-1)){
      for(l in (j+1):k){
        Gammasn[j,l] = Gammas[j,l] + rnorm(1, 0, sgs_rw[j,l])
        
        omega = Omega_sim(n = n, Ti = Ti, k = k, u = u, Beta = Beta, Betas = Betas, Delta = Delta, Gamma = Gamma, Gammas = Gammasn, iter = RR, Com = Com)
        omega = omega$u
        
        ratio_g = log_prior_theta(Gammasn[j,l], mb, sb) - log_prior_theta(Gammas[j,l], mb, sb) + loggammas(Gammasn, Gammas, n, Ti, k, j, l, u = u, C = C) + loggammas(Gammas, Gammasn, n, Ti, k, j, l, u = omega, C = C)
        
        alpha_g = min(1, exp(ratio_g))
        if(runif(1)<= alpha_g){
          Gammas[j,l] = Gammasn[j,l]
          accgs[j,l] = accgs[j,l] +1
        }
        Gammasn = Gammas
        Gammas_es[j,l,r] = Gammas[j,l]
        if(r <= stop_){
          lambda_gs[j,l,r] = exp(log(lambda_gs[j,l,r-1]) + ga*(alpha_g - alph_opt))
          mu_gs[j,l,r] = mu_gs[j,l,r-1] + ga*(Gammas[j,l] - mu_gs[j,l,r-1])
          sigma_gs[j,l,r] = sigma_gs[j,l,r-1] + ga*((Gammas[j,l] - mu_gs[j,l,r-1])^2 - sigma_gs[j,l,r-1])
          sgs_rw[j,l] = lambda_gs[j,l,r]*sigma_gs[j,l,r]
        }
      }
    }
    
    for(l in 1:(k-1)){
      for(j in (l+1):k){
        Gammasn[j,l] = Gammas[j,l] + rnorm(1, 0, sgs_rw[j,l])
        
        omega = Omega_sim(n = n, Ti = Ti, k = k, u = u, Beta = Beta, Betas = Betas, Delta = Delta, Gamma = Gamma, Gammas = Gammasn, iter = RR, Com = Com)
        omega = omega$u
        
        ratio_g = log_prior_theta(Gammasn[j,l], mb, sb) - log_prior_theta(Gammas[j,l], mb, sb) + loggammas(Gammasn, Gammas, n, Ti, k, j, l, u = u, C = C) + loggammas(Gammas, Gammasn, n, Ti, k, j, l, u = omega, C = C)
        
        alpha_g = min(1, exp(ratio_g))
        if(runif(1)<= alpha_g){
          Gammas[j,l] = Gammasn[j,l]
          accgs[j,l] = accgs[j,l] +1
        }
        Gammasn = Gammas
        Gammas_es[j,l,r] = Gammas[j,l]
        if(r <= stop_){
          lambda_gs[j,l,r] = exp(log(lambda_gs[j,l,r-1]) + ga*(alpha_g - alph_opt))
          mu_gs[j,l,r] = mu_gs[j,l,r-1] + ga*(Gammas[j,l] - mu_gs[j,l,r-1])
          sigma_gs[j,l,r] = sigma_gs[j,l,r-1] + ga*((Gammas[j,l] - mu_gs[j,l,r-1])^2 - sigma_gs[j,l,r-1])
          sgs_rw[j,l] = lambda_gs[j,l,r]*sigma_gs[j,l,r]
        }
      }
    }
    
    #Delta
    for(j in 1:(k-1)){
      for(l in (j+1):k){
        Deltan[j,l] = Delta[j,l] + rnorm(1, 0, sd_rw[j,l])
        
        omega = Omega_sim(n = n, Ti = Ti, k = k, u = u, Beta = Beta, Betas = Betas, Delta = Deltan, Gamma = Gamma, Gammas = Gammas, iter = RR, Com = Com)
        omega = omega$u
        
        ratio_d = log_prior_theta(Deltan[j,l], mb, sb) - log_prior_theta(Delta[j,l], mb, sb) + logdelta(Deltan, Delta, Ti, k, j, l, u = u) + logdelta(Delta, Deltan, Ti, k, j, l, u = omega)
        
        alpha_d = min(1, exp(ratio_d))
        if(runif(1)<= alpha_d){
          Delta[j,l] = Deltan[j,l]
          accd[j,l] = accd[j,l] +1
        }
        Deltan = Delta 
        Delta_es[j,l,r] = Delta[j,l]
        if(r <= stop_){
          lambda_d[j,l,r] = exp(log(lambda_d[j,l,r-1]) + ga*(alpha_d - alph_opt))
          mu_d[j,l,r] = mu_d[j,l,r-1] + ga*(Delta[j,l] - mu_d[j,l,r-1])
          sigma_d[j,l,r] = sigma_d[j,l,r-1] + ga*((Delta[j,l] - mu_d[j,l,r-1])^2 - sigma_d[j,l,r-1])
          sd_rw[j,l] = lambda_d[j,l,r]*sigma_d[j,l,r]
        }
      }
    }
    
    for(l in 1:(k-1)){
      for(j in (l+1):k){
        Deltan[j,l] = Delta[j,l] + rnorm(1, 0, sd_rw[j,l])
        
        omega = Omega_sim(n = n, Ti = Ti, k = k, u = u, Beta = Beta, Betas = Betas, Delta = Deltan, Gamma = Gamma, Gammas = Gammas, iter = RR, Com = Com)
        omega = omega$u
        
        ratio_d = log_prior_theta(Deltan[j,l], mb, sb) - log_prior_theta(Delta[j,l], mb, sb) + logdelta(Deltan, Delta, Ti, k, j, l, u = u) + logdelta(Delta, Deltan, Ti, k, j, l, u = omega)
        
        alpha_d = min(1, exp(ratio_d))
        if(runif(1)<= alpha_d){
          Delta[j,l] = Deltan[j,l]
          accd[j,l] = accd[j,l] +1
        }
        Deltan = Delta
        Delta_es[j,l,r] = Delta[j,l]
        if(r <= stop_){
          lambda_d[j,l,r] = exp(log(lambda_d[j,l,r-1]) + ga*(alpha_d - alph_opt))
          mu_d[j,l,r] = mu_d[j,l,r-1] + ga*(Delta[j,l] - mu_d[j,l,r-1])
          sigma_d[j,l,r] = sigma_d[j,l,r-1] + ga*((Delta[j,l] - mu_d[j,l,r-1])^2 - sigma_d[j,l,r-1])
          sd_rw[j,l] = lambda_d[j,l,r]*sigma_d[j,l,r]
        }
      }
    }
    
    #Full U
    M = array(NA, dim= c(k,k,n))
    MG = array(NA, dim = c(k,k,n))
    res = matrix(NA, n, k)
    res2 = matrix(NA, n, k)
    res3 = matrix(NA, n, k)
    res4 = matrix(NA, n, k)
    fin = efin = efinb = efinc = eq_a = eq_b = eq_c = matrix(NA, n, k)
    matr = matrix(rep(c(1:k), rep(n,1)),k,n)
    
    for(t in 1:Ti){
      if(t==1){
        for(i in 1:k) dens[i,] = exp(loglik_normal(yf[,,t]-mu[,i,r],sigma[,,i,r]))
        Ma = C*t(matrix(rep(u[,1], n), n,n))
        for(l in 1:n){
          M[,,l] = table(matr, factor(matrix(rep(Ma[l,],k), k,n, byrow = TRUE),1:k))
          res2[l,] = Delta[1:k,u[l,2]]
          res3[l,] = Beta[1:k]
          MG[,,l] = M[,,l]*Gamma
          res[l,] = rowSums(MG[,,l])
          fin[l,] = res[l,]+res2[l,]+res3[l,]
          efin[l,] = exp(fin[l,])*dens[,l]
          eq_a[l,] = efin[l,]/sum(efin[l,])
          u[l,t] = sample(U, 1, prob = eq_a[l,], replace = TRUE)
        }
      }else if(t == Ti){
        for(i in 1:k) dens[i,] = exp(loglik_normal(yf[,,t]-mu[,i,r],sigma[,,i,r]))
        Ma = C*t(matrix(rep(u[,Ti], n), n,n))
        for(l in 1:n){
          M[,,l] = table(matr,factor(matrix(rep(Ma[l,],k), k,n, byrow = TRUE),1:k))
          res4[l,] = Delta[u[l,Ti-1],1:k]
          res3[l,] = Betas[1:k]
          MG[,,l] = M[,,l]*Gammas
          res[l,] = rowSums(MG[,,l])
          fin[l,] = res[l,]+res3[l,]+res4[l,]
          efinb[l,] = exp(fin[l,])*dens[,l]
          eq_b[l,] = efinb[l,]/sum(efinb[l,])
          u[l,t] = sample(U, 1, prob = eq_b[l,], replace = TRUE)
        }
      }else{
        for(i in 1:k) dens[i,] = exp(loglik_normal(yf[,,t]-mu[,i,r],sigma[,,i,r]))
        Ma = C*t(matrix(rep(u[,t], n), n,n))
        for(l in 1:n){
          M[,,l] = table(matr,factor(matrix(rep(Ma[l,],k), k,n, byrow = TRUE),1:k))
          res2[l,] = Delta[1:k,u[l,t+1]]
          res4[l,] = Delta[u[l,t-1],1:k]
          res3[l,] = Betas[1:k]
          MG[,,l] = M[,,l]*Gammas
          res[l,] = rowSums(MG[,,l])
          fin[l,] = res[l,]+res2[l,]+res3[l,]+res4[l,]
          efinc[l,] = exp(fin[l,])*dens[,l]
          eq_c[l,] = efinc[l,]/sum(efinc[l,])
          u[l,t] = sample(U, 1, prob = eq_c[l,], replace = TRUE)
        }
      }
    }
    U_ex[,,r] = u 
  }
  save.image(paste0("Sim_lattice_k2_n9_t5_sigmasqrt1_rv_r5_",ss,".Rdata"))
  rm(list = ls())
}
