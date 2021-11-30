set.seed(42069)
library(alphashape3d)
library(plot3D)
library(dirichletprocess)
library(mvnfast)
torus = rtorus(1500, 1, 3)
scatter3D(torus[,1],torus[,2],torus[,3],col = "blue", 
          theta=0, phi = 60, cex =0.5, pch=15, main = "original")


model = DirichletProcessMvnormal(torus)
model = Fit(model, 1500)

n = length(model$clusterParameters$mu)/3
if(n == 2){
  mus = rbind(model$clusterParameters$mu[,,1], model$clusterParameters$mu[,,2])
  sigmas = list(model$clusterParameters$sig[,,1], model$clusterParameters$sig[,,2])
  ws = model$weights 
  xSim = rmixn(1500, mu = mus, sigma = sigmas, w = ws)
} else {
  print("more than two clusters lad")
}
             
scatter3D(xSim[,1],xSim[,2],xSim[,3],col = "blue", 
          theta=0, phi = 60, cex =0.5, pch=15, main = "DPMG")
