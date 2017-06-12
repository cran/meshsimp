#' Driver running the simplification process within R, reading the mesh from file.
#' 
#' @param n1,n2		Two different levels of simplification
#' @param wgeom		Weight for the geometric cost function
#' @param wdisp		Weight for the displacement cost function
#' @param wequi		Weight for the equidistribution cost function

fpath <- system.file("extdata", "pawn_2522.inp", package="meshsimp")
n1 <- 2000
n2 <- 1000
wgeom <- 1/3
wdisp <- 1/3
wequi <- 1/3

# Import the mesh
require(meshsimp)
mesh <- import.mesh.2.5D(file)

# Plot original mesh
plot.mesh.2.5D(mesh, main = sprintf("Original mesh, %i nodes", mesh$nnodes))

# Simplify the mesh, then plot
out1 <- simplify.mesh.2.5D(mesh, n1)
plot.mesh.2.5D(out1$mesh, out1$locations, main = sprintf("Simplified mesh, %i nodes", n1))

# Resume the simplification, then plot the final mesh
out2 <- simplify.mesh.2.5D(out1$mesh, n2, out1$locations)
plot.mesh.2.5D(out2$mesh, out2$locations, main = sprintf("Simplified mesh, %i nodes", n2))

