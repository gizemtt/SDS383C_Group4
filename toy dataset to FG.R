# get spiral data
library(ider)
spiral = gendata(DataName = "Spiral", n = 750, noise =0.01, seed=1)
plot(spiral[,1], spiral[,2])

# get torus data
library(alphashape3d)
torus = rtorus(500, 3, 6)
scatter3D(torus[,1],torus[,2],torus[,3],col = "black", sphere.size=0.2, 
          theta=0,phi = 60)
