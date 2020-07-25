# Rotation Matrix for two vectors in 3D Cartesian space

crossProduct <- function(ab, ac){
  
  ## vector cross product in Cartesian coordinates
  
  abci = ab[2] * ac[3] - ac[2] * ab[3]
  abcj = - (ab[1] * ac[3] - ac[1] * ab[3])
  abck = ab[1] * ac[2] - ac[1] * ab[2]
  return(c(abci, abcj, abck))
}

RotationMatrix <- function(A, B){

  # Rotate vector A around rotational axis unit vector k to new position B
  # In Cartesian coordinates
  # Based on Rodrigues' rotation formula

  # vector k = A x B
  k = crossProduct(A, B)
  
  print(k)
  # |k|
  kmod = sqrt(sum(k^2))
  print(kmod)
  
  # |A| & |B|
  Amod = sqrt(sum(A^2))
  Bmod = sqrt(sum(B^2))

  sin_theta = kmod/(Amod*Bmod)
  theta = asin(sin_theta)
  print(theta)
  
  khat <- k / (Amod*Bmod*sin_theta) # k rotational unit vector
  print(khat)
  
  I = crossprod(diag(length(khat)))
  K = rbind(c(0, -khat[3], khat[2]), 
            c(khat[3], 0, -khat[1]), 
            c(-khat[2], khat[1], 0))

  K2 = K %*% K

  R = I + sin_theta*K + (1-cos(theta))*K2
  print(R)
  return(R)
}

A = c(0,0,1)
B = c(0,1,0)
R1 <- RotationMatrix(A, B)
t(R1 %*% A) # this should return B

A = c(1,1,1) / sqrt(3)
B = c(1,0,0)
R1 <- RotationMatrix(A, B)
t(R1 %*% A) # this should return B


## rotation of north pole vector 45 degrees around y-axis
A = c(0,0,1)
B = c(-1,0,1) / sqrt(2)
R1 <- RotationMatrix(A, B)
t(R1 %*% A) # this should return B


## rotation of Mars North pole in J2000 frame to North pole in (lon, lat)
A = c(0.44615053, -0.40623628, 0.79744704) # Mars North Pole in J2000 (317.681*u.deg, 52.887*u.deg)
B = c(0,0,1) # north pole in (lon, lat)
R1 <- RotationMatrix(A, B)
t(R1 %*% A) # this should return B

# -0.4062363 -0.4461505  0.0000000
# [,1]       [,2]       [,3]
# [1,] 0.8892594  0.1008333 -0.4461505
# [2,] 0.1008333  0.9081876  0.4062363
# [3,] 0.4461505 -0.4062363  0.7974470

