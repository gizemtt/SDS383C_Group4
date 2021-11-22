# get spiral data
library(ider)
spiral = gendata(DataName = "Spiral", n = 750, noise =0.01, seed=1)
plot(spiral[,1], spiral[,2])

# get torus data
set.seed(42069)
library(alphashape3d)
library(plot3D)
torus = rtorus(1000, 1, 3)
scatter3D(torus[,1],torus[,2],torus[,3],col = "blue", 
          theta=0, phi = 60, cex =0.5, pch=15)
