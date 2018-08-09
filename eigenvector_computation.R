Sigma_w = diag(6)
Sigma_w[6,5] = Sigma_w[5,6] = Sigma_w[3,4] = Sigma_w[4,3] = Sigma_w[1,2] = Sigma_w[2,1] = -1
vecs = eigen(Sigma_w)$vectors
vals = eigen(Sigma_w)$values
z = rep(c(-1, 1), 3)
z = z / sum(z^2)
z
vecs %*% z
