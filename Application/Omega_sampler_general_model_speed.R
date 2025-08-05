#Omega simulator 

Omega_sim = function(n,Ti,k,u,Beta,Betas,Delta,Gamma,Gammas,iter,Com){
  R = iter
  
  #Initial Hidden states
  U = 1:k
  mU = matrix(U,k,n)
  
  Com = Com
  C = matrix(0,n,n)
  C[upper.tri(C, diag = FALSE)] = Com[upper.tri(Com, diag = FALSE)]
  
  M = array(NA, dim= c(k,k,n))
  MG = array(NA, dim = c(k,k,n))
  res = matrix(NA, n, k)
  res2 = matrix(NA, n, k)
  res3 = matrix(NA, n, k)
  res4 = matrix(NA, n, k)
  fin = matrix(NA, n, k)
  eq_a = eq_b = eq_c = matrix(NA, n, k)
  matr = matrix(rep(c(1:k), rep(n,1)),k,n)
  
  #Funzionante
  for(r in 2:R){
    for(t in 1:Ti){
      if(t==1){
        Ma = C*t(matrix(rep(u[,1], n), n,n))
        for(l in 1:n){
          M[,,l] = table(matr, factor(matrix(rep(Ma[l,],k), k,n, byrow = TRUE),1:k))
          res2[l,] = Delta[1:k,u[l,2]]
          res3[l,] = Beta[1:k]
          MG[,,l] = M[,,l]*Gamma
          res[l,] = rowSums(MG[,,l])
          fin[l,] = res[l,]+res2[l,]+res3[l,]
          eq_a[l,] = exp(fin[l,])/sum(exp(fin[l,]))
          u[l,t] = sample(U, 1, prob = eq_a[l,], replace = TRUE)
        }
      }else if(t == Ti){
        Ma = C*t(matrix(rep(u[,Ti], n), n,n))
        for(l in 1:n){
          M[,,l] = table(matr,factor(matrix(rep(Ma[l,],k), k,n, byrow = TRUE),1:k))
          res4[l,] = Delta[u[l,Ti-1],1:k]
          res3[l,] = Betas[1:k]
          MG[,,l] = M[,,l]*Gammas
          res[l,] = rowSums(MG[,,l])
          fin[l,] = res[l,]+res3[l,]+res4[l,]
          eq_b[l,] = exp(fin[l,])/sum(exp(fin[l,]))
          u[l,t] = sample(U, 1, prob = eq_b[l,], replace = TRUE)
        }
      }else{
        Ma = C*t(matrix(rep(u[,t], n), n,n))
        for(l in 1:n){
          M[,,l] = table(matr,factor(matrix(rep(Ma[l,],k), k,n, byrow = TRUE),1:k))
          res2[l,] = Delta[1:k,u[l,t+1]]
          res4[l,] = Delta[u[l,t-1],1:k]
          res3[l,] = Betas[1:k]
          MG[,,l] = M[,,l]*Gammas
          res[l,] = rowSums(MG[,,l])
          fin[l,] = res[l,]+res2[l,]+res3[l,]+res4[l,]
          eq_c[l,] = exp(fin[l,])/sum(exp(fin[l,]))
          u[l,t] = sample(U, 1, prob = eq_c[l,], replace = TRUE)
        }
      }
    }
  }
  res = list(u = u)
}  