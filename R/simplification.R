#'	Instantiate a \code{mesh.2.5D} object from file.
#' 
#'	@param	file	Absolute or relative path to the input mesh; the following file formats are supported:
#'					\itemize{
#'						\item AVS UCD ASCII (extension .inp);
#'						\item text files (extension .txt);
#'						\item Legacy VTK (extension .vtk).
#'					}
#'					Text files are assumed to be structured as the AVS UCD ASCII files (.inp extension).
#'	@description	This function reads a surface mesh from file, and returns an object of class \code{mesh.2.5D} holding the list of nodes and triangles in the mesh. The parsing of the file is carried out at the C++ level for the sake of efficiency.
#'	@seealso		\code{\link{simplify.mesh.2.5D}}, \code{\link{plot.mesh.2.5D}}
#'	@usage			import.mesh.2.5D(file)
#'	@return 		An object of class \code{mesh.2.5D}, provided with the following attributes:
#'					\itemize{
#'						\item \code{nnodes}: number of nodes in the mesh;
#'						\item \code{nodes}: \code{nnodes}-by-3 matrix collecting the coordinates of each vertex;
#'						\item \code{ntriangles}: number of triangles in the mesh;
#'						\item \code{triangles}: a \code{ntriangles}-by-3 (when \code{order} = 1) or \code{ntriangles}-by-6 (when \code{order} = 2) matrix. It specifies the triangles giving the row indices in \code{nodes} of the triangles vertices and (when \code{order} = 2) also of the triangles edges midpoints;
#'						\item \code{order}: either '1' or '2'. It specifies whether each mesh triangle should be represented by \eqn{3} nodes (the triangle vertices) or by \eqn{6} nodes (the triangle vertices and midpoints of the triangle edges). These are respectively used for linear (\code{order} = 1) and quadratic (\code{order} = 2) Finite Elements. Default is \code{order} = 1.
#'					}
#'	@export
#'	@examples
#'	## Import the mesh of a pawn
#'	fpath <- system.file("extdata", "pawn_250.inp", package="meshsimp")
#'	mesh <- import.mesh.2.5D(fpath)
#'	## Simplify the mesh down to 200 nodes; assume the components of the
#'	## edge cost functions are equally weighted and that the data locations
#'	## coincide with the vertices of the mesh
#'	out1 <- simplify.mesh.2.5D(mesh, 200)
#'	## Resume the simplification procedure, reducing the mesh down to 150 nodes
#'	out2 <- simplify.mesh.2.5D(out1$mesh, 150, out1$locations)

import.mesh.2.5D <- function(file)
{
	# Let C++ read the file and extract nodes and elements
	mod_meshsimp <- Module("mod_meshsimp", PACKAGE = "meshsimp")
	mesh <- mod_meshsimp$readMesh(file)
	
	# Convert from C++ 0-based indexing to R 1-based indexing
	mesh$triangles <- mesh$triangles + 1
	
	# Create a mesh.2.5D object
	out <- list(nnodes = mesh$nnodes, nodes = mesh$nodes, ntriangles = mesh$ntriangles, 
		triangles = mesh$triangles, order = 1)
	class(out) <- "mesh.2.5D"
	
	return(out)
}


#' Perform the simplification of a mesh, given as an object of class \code{mesh.2.5D}.
#' 
#' @param mesh 		A \code{mesh.2.5D} object, endowed with the following attributes:
#'					\itemize{
#'						\item \code{nnodes}: number of nodes in the mesh;
#'						\item \code{nodes}: \code{nnodes}-by-3 matrix collecting the coordinates of each vertex;
#'						\item \code{ntriangles}: number of triangles in the mesh;
#'						\item \code{triangles}: a \code{ntriangles}-by-3 (when \code{order} = 1) or \code{ntriangles}-by-6 (when \code{order} = 2) matrix. It specifies the triangles giving the row indices in \code{nodes} of the triangles vertices and (when \code{order} = 2) also of the triangles edges midpoints;
#'						\item \code{order}: either '1' or '2'. It specifies whether each mesh triangle should be represented by \eqn{3} nodes (the triangle vertices) or by \eqn{6} nodes (the triangle vertices and midpoints of the triangle edges). These are respectively used for linear (\code{order} = 1) and quadratic (\code{order} = 2) Finite Elements. Default is \code{order} = 1.
#'					}
#'	@param  target	Number of nodes which the mesh should feature at the end of the simplification process. In other terms, \code{target} is used as stopping criterium: the algorithm stops when the number of nodes in the mesh matches \code{target}.
#'	@param	loc 	#data-by-3 vector with data locations; default is NULL, i.e. data locations are assumed to coincide with the mesh vertices.
#'	@param	val		#data-by-1 vector with the observations associated with each data point; default is NULL.
#'	@param	wgeom	Weight for the geometric component of the edge cost function; default is 1/3. Note that the all weights should be positive and sum up to one.
#'	@param	wdisp 	Weight for the data displacement component of the edge cost function; default is 1/3. Note that the all weights should be positive and sum up to one.
#'	@param	wequi	Weight for the data equidistribution component of the edge cost function; default is 1/3. Note that the all weights should be positive and sum up to one.
#'	@param  file String specifying the path to the location where the simplified mesh will be stored; the following file formats are supported:
#'					\itemize{
#'						\item AVS UCD ASCII (extension .inp);
#'						\item text files (extension .txt);
#'						\item Legacy VTK (extension .vtk).
#'					}
#'					Text files are assumed to be structured as the AVS UCD ASCII files (.inp extension).
#'					If \code{outfile} is not provided, the mesh will not get saved to file but just returned as \code{mesh.2.5D} object.
#'	@description	Implementation of a mesh simplification strategy for a surface mesh with associated distributed data. The algorithm works by iteratively collapsing an edge into an internal point. The selection of the edge to contract at each iteration is driven by a cost functional, which measures the loss of geometrical accuracy and data information associated with the collapse. For a detailed description of the algorithm, please refer to: \cr
#'					Dassi, F., Ettinger, B., Perotto, S., Sangalli, L.M. (2015), \cr
#'					A mesh simplification strategy for a spatial regression analysis over the  cortical surface of the brain, \cr
#'					Applied Numerical Mathematics, Vol. 90, pp. 111-131. \cr
#'
#'					The whole (computing-intensive) procedure is carried out at the C++ level, thus ensuring high-performance. In detail, the function relies on the class \code{RcppSimplification} - a wrapper for the template class \code{simplification<Triangle, MeshType::DATA, DataGeo>} provided by the C++ library \code{meshsimplification}. Methods of class \code{RcppSimplification} are exposed to R through the API supplied by the \pkg{Rcpp} and \pkg{RcppEigen} packages.	
#'	@seealso		\code{\link{plot.mesh.2.5D}}, \code{\link{import.mesh.2.5D}}			
#'	@usage			simplify.mesh.2.5D(mesh,target,loc=NULL,val=NULL,wgeom=1/3,wdisp=1/3,wequi=1/3,file='')
#'	@return 		A list equipped with the following fields:
#'					\itemize{
#'						\item \code{mesh}: the simplified mesh as an instance of class \code{mesh.2.5D}; the order of the output mesh coincides with the order of the input mesh;
#'						\item \code{simplifier}: the object of class \code{RcppSimplification} used within the function to carry out the simplification procedure at the C++ level. This object is returned since eases (and speeds up) the extraction of useful information about the mesh, which are not directly made available by \code{mesh.2.5D} or which may rely on the connectivities of the mesh itself, as, e.g., the list of edges (see \code{\link{get.edges}});
#'						\item \code{locations}: #data-by-3 matrix holding the location of each data point over the simplified mesh;
#' 						\item \code{qoi}: vector listing the quantity of information associated with each triangle in the simplified mesh; see \code{\link{get.quantity.of.information}} for a rigorous definition of the quantity of information.
#'					}
#'	@export
#'	@examples
#'	## Import the mesh of a pawn 
#'  data(pawn_250)
#'	## Simplify the mesh down to 200 nodes; assume the components of the
#'	## edge cost functions are equally weighted and that the data locations
#'	## coincide with the vertices of the mesh
#'	out1 <- simplify.mesh.2.5D(mesh, 200)
#'	## Resume the simplification procedure, reducing the mesh down to 150 nodes
#'	out2 <- simplify.mesh.2.5D(out1$mesh, 150, out1$locations)

simplify.mesh.2.5D <- function(mesh, target, loc = NULL, val = NULL, wgeom = 1/3, wdisp = 1/3, wequi = 1/3, file = '')
{
	# Check whether the weights sum up to one
	s <- wgeom + wdisp + wequi
	tol <- 1e-15
	if (wgeom < 0 | wdisp < 0 | wequi < 0 | abs(s - 1) > tol)
		stop("The weights must be positive and sum up to one.") 
		
	# Check if mesh is a mesh.2.5D class object
	if (class(mesh) != "mesh.2.5D")
		stop("mesh must be an object of class mesh.2.5D.")
		
	# Load module
	mod_meshsimp <- Module("mod_meshsimp", PACKAGE = "meshsimp")
		
	# Extract vertices and triangles of the grid
	nodes <- mesh$nodes
	triangles <- mesh$triangles
	
	# Convert from R 1-based indexing to C++ 0-based indexing
	triangles <- triangles-1
		
	# If the mesh consists of quadratic elements, convert to linear elements
	if (mesh$order == 2)
	{
		res <- mod_meshsimp$getLinearFEMesh(nodes, triangles)
		nodes <- res$nodes
		triangles <- res$triangles
	}
		
	# Create an RcppSimplification object
	if (is.null(loc))
		simplifier <- new(mod_meshsimp$RcppSimplification, nodes, triangles, wgeom, wdisp, wequi)
	else if (is.null(val))
		simplifier <- new(mod_meshsimp$RcppSimplification, nodes, triangles, loc, wgeom, wdisp, wequi)
	else
		simplifier <- new(mod_meshsimp$RcppSimplification, nodes, triangles, loc, val, wgeom, wdisp, wequi)
		
	# Run the simplification
	simplifier$simplify(target, file)
	
	# Create the output mesh.2.5D object
	newmesh <- get.mesh.2.5D(simplifier, mesh$order)
	
	# Get list of data points
	loc <- get.data.locations(simplifier)
	
	# Get quantity of information for each triangle
	qoi <- get.quantity.of.information(simplifier)
	
	# Output
	return(list(mesh = newmesh, simplifier = simplifier, locations = loc, qoi = qoi))
}


#' Get the list of edges from an object of class \code{RcppSimplification}.
#'
#'	@param	x	An object of class \code{RcppSimplification}.
#'	@usage		get.edges(x)
#'	@seealso	\code{\link{simplify.mesh.2.5D}}
#'	@return		A #edges-by-2 matrix, where the \eqn{i}-th row stores the identifiers of the vertices at the end-points of the \eqn{i}-th edge.	
#'	@export

get.edges <- function(x)
{
	# Preliminary check
	if (class(x) != "Rcpp_RcppSimplification")
		stop("The input argument must be an object of class RcppSimplification.")
		
	# Run, converting from C++ 0-based indexing to R 1-based indexing
	out <- x$getEdges()
	out < out + 1
	return(out)
}


#'	Get the list of data locations from an object of class \code{RcppSimplification}.
#'
#'	@param	x	An object of class \code{RcppSimplification}.
#'	@usage		get.data.locations(x)
#'	@seealso	\code{\link{simplify.mesh.2.5D}}
#'	@return		A #data-by-3 matrix storing the coordinates of data locations.
#'	@export	

get.data.locations <- function(x)
{
	# Preliminary check
	if (class(x) != "Rcpp_RcppSimplification")
		stop("The input argument must be an object of class RcppSimplification.")
		
	# Run
	out <- x$getDataLocations()
	return(out)
}


#'	Get the quantity of information associated with each element of the triangulation, from an object of class \code{RcppSimplification}.
#'
#'	@param	x		An object of class \code{RcppSimplification}.
#'	@usage			get.quantity.of.information(x)
#'	@description	For each triangle \eqn{T} in the mesh, the associated quantity of information \eqn{N_T} is defined as: \deqn{N_T := n_f + \frac{1}{2}n_e + \frac{1}{\# \big( T_{v_1} \big)} n_1 + \frac{1}{\# \big( T_{v_2} \big)} n_2 + \frac{1}{\# \big( T_{v_3} \big)} n_3 \, ,} where \eqn{n_f} and \eqn{n_e} denote the number of data points associated with the face and the edges of the triangle \eqn{T}, respectively. For \eqn{j = 1, \, 2, \, 3}, \eqn{n_j} is the number of data points associated with the \eqn{j}-th vertex \eqn{v_j} of \eqn{T}, \eqn{T_{v_j}} is the patch of elements associated with \eqn{v_j} (i.e., the elements sharing \eqn{v_j}), and \eqn{\# \big( T_{v_j} \big)} denotes the cardinality of the patch \eqn{T_{v_j}}.
#'	@seealso		\code{\link{simplify.mesh.2.5D}}
#'	@return			A #elements-by-1 vector, storing the quantity of information for each element.
#'	@export

get.quantity.of.information <- function(x)
{
	# Preliminary check
	if (class(x) != "Rcpp_RcppSimplification")
		stop("The input argument must be an object of class RcppSimplification.")
		
	# Run
	out <- x$getQuantityOfInformation()
	return(out)
}


#'	Get surface mesh from an object of class \code{RcppSimplification}.
#'
#'	@param	x		An object of class \code{RcppSimplification}.
#'	@param  order	Either '1' or '2'. It specifies whether each mesh triangle should be represented by \eqn{3} nodes (the triangle vertices) or by \eqn{6} nodes (the triangle vertices and midpoints of the triangle edges). These are respectively used for linear (\code{order} = 1) and quadratic (\code{order} = 2) Finite Elements. Default is \code{order} = 1.
#'	@description	Extract the mesh from an object of class \code{RcppSimplification}. The mesh is returned as an instance of class \code{mesh.2.5D}. The order of the Finite Elements is compliant with the input parameter \code{order}.
#'	@usage			get.mesh.2.5D(x, order = 1)
#'	@seealso		\code{\link{simplify.mesh.2.5D}}
#'	@return			A \code{mesh.2.5D} object, endowed with the following attributes:
#'					\itemize{
#'						\item \code{nnodes}: number of nodes in the mesh;
#'						\item \code{nodes}: \code{nnodes}-by-3 matrix collecting the coordinates of each vertex;
#'						\item \code{ntriangles}: number of triangles in the mesh;
#'						\item \code{triangles}: a \code{ntriangles}-by-3 (when \code{order} = 1) or \code{ntriangles}-by-6 (when \code{order} = 2) matrix. It specifies the triangles giving the row indices in \code{nodes} of the triangles vertices and (when \code{order} = 2) also of the triangles edges midpoints;
#'						\item \code{order}: either '1' or '2'. It specifies whether each mesh triangle should be represented by \eqn{3} nodes (the triangle vertices) or by \eqn{6} nodes (the triangle vertices and midpoints of the triangle edges). These are respectively used for linear (\code{order} = 1) and quadratic (\code{order} = 2) Finite Elements. Default is \code{order} = 1.
#'					}

get.mesh.2.5D <- function(x, order = 1)
{
	# Preliminary check
	if (class(x) != "Rcpp_RcppSimplification")
		stop("The input argument must be an object of class RcppSimplification.")
		
	# Extract nodes and triangles, differentiating between
	# linear and quadratic Finite Elements
	if (order == 1)
	{
		nodes <- x$getNodes()
		triangles <- x$getElems()
	}
	else if (order == 2)
	{
		res <- x$getQuadraticFEMesh()
		nodes <- res$nodes
		triangles <- res$triangles
	}
	
	# Convert from C++ 0-based indexing to R 1-based indexing
	triangles <- triangles + 1
	
	# Create the mesh
	out <- list(nnodes = dim(nodes)[1], nodes = nodes, 
		ntriangles = dim(triangles)[1], triangles = triangles, order = order)
	class(out) <- "mesh.2.5D"
	return(out)
}


#'	Plot a mesh in a 3D perspective.
#'
#'	@param	x		An object of class \code{mesh.2.5D}, representing the mesh to plot.
#'	@param 	loc		#data-by-3 matrix collecting the locations of data points over the mesh.
#'					Default is NULL, i.e. the data points are assumed to coincide with the mesh vertices.
#'	@param	phi		Colatitude (in degrees) identifying the viewing direction; default is 40.
#'	@param	theta	Longitude (in degrees) identifying the viewing direction; default is 40.
#'	@param	...		Additional arguments passed to the plotting methods; these include:
#'					xlim, ylim, zlim, xlab, ylab, zlab, main, sub, r, d, scale, expand, box, axes, nticks, ticktype.
#'	@description	Plot a 2.5D surface mesh augmented with associated data locations. The package \pkg{plot3D} is used.
#'	@rdname			plot
#'	@method			plot mesh.2.5D
#'	@seealso		\code{\link{simplify.mesh.2.5D}}
#'	@export
#'	@examples
#'	## Import the mesh of a pawn
#'  data(pawn_250)
#'	## Plot the original mesh
#'	plot.mesh.2.5D(mesh, main = "Original mesh, 250 nodes")
#'	## Simplify the mesh down to 200 nodes; assume the components of the
#'	## edge cost functions are equally weighted and that the data locations
#'	## coincide with the vertices of the mesh
#'	out <- simplify.mesh.2.5D(mesh, 200)
#'	## Plot the simplified mesh
#'	plot.mesh.2.5D(out$mesh, loc = out$locations, main = "Simplified mesh, 200 nodes")

plot.mesh.2.5D <- function(x, loc = NULL, phi = 40, theta = 40, ...)
{
	# Preliminary check
	if (class(x) != "mesh.2.5D")
		stop("mesh must be an object of class mesh.2.5D.")
				
	# Extract vertices of each triangle 
	trs <- matrix(NA, nrow = 4*x$ntriangles-1, ncol = 3)
	for (i in 1:x$ntriangles)
	{
		trs[4*i-3,] = x$nodes[x$triangles[i,1],]
		trs[4*i-2,] = x$nodes[x$triangles[i,2],]
		trs[4*i-1,] = x$nodes[x$triangles[i,3],]
	}
	
	# Determine data locations
	if (is.null(loc))
		points = x$nodes
	else
		points = loc
	n = nrow(points)
	
    # Array of colors
	clr_trs <- matrix("blue", nrow = x$ntriangles, ncol = 1)
	clr_loc <- matrix("red", nrow = n, ncol = 1)
	
	# Plot
	polygon3D (x = trs[,1], y = trs[,2], z = trs[,3], ...,
      	colvar = NULL, phi = phi, theta = theta,
      	col = clr_trs, NAcol = "white", breaks = NULL,
      	facets = FALSE, bty = "b", pty = "s",
      	add = FALSE, plot = TRUE)
    points3D (x = points[,1], y = points[,2], z = points[,3], ...,
      	colvar = NULL, phi = phi, theta = theta,
      	col = clr_loc, breaks = NULL,
      	facets = FALSE, bty = "r", pty = "s",
      	add = TRUE, plot = TRUE)
}


#'	The mesh of a pawn.
#'
#'	@docType	data
#'	@name		mesh
'mesh'	
