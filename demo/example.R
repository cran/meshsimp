#' Driver running the simplification process within R.
#' 
#' @param dataset	Name of dataset storing the mesh as an object of class mesh.2.5D (named 'mesh'):
#'					- pawn_2522: pawn geometry, 2522 nodes;
#'					- pawn_500: pawn geometry, 500 nodes;
#'					- pawn_300: pawn geometry, 300 nodes;
#'					- pawn_250: pawn geometry, 250 nodes.
#' @param n1,n2		Two different levels of simplification
#' @param wgeom		Weight for the geometric cost function
#' @param wdisp		Weight for the displacement cost function
#' @param wequi		Weight for the equidistribution cost function

dataset <- "pawn_2500"
n1 <- 400
n2 <- 300
wgeom <- 1/3
wdisp <- 1/3
wequi <- 1/3

# Load data
require(meshsimp)
if (dataset == "pawn_2522")
	data(pawn_2522)
if (dataset == "pawn_500")
	data(pawn_500)
if (dataset == "pawn_300")
	data(pawn_300)
if (dataset == "pawn_250")
	data(pawn_250)

# Plot original mesh
plot.mesh.2.5D(mesh, main = sprintf("Original mesh, %i nodes", mesh$nnodes))

# Simplify the mesh, then plot
out1 <- simplify.mesh.2.5D(mesh, n1)
plot.mesh.2.5D(out1$mesh, out1$locations, main = sprintf("Simplified mesh, %i nodes", n1))

# Resume the simplification, then plot the final mesh
out2 <- simplify.mesh.2.5D(out1$mesh, n2, out1$locations)
plot.mesh.2.5D(out2$mesh, out2$locations, main = sprintf("Simplified mesh, %i nodes", n2))



