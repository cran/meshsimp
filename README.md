The R package meshsimp
======================

The R package `meshsimp` offers a high-level access to the C++ library `meshsimplification`, implementing a simplification strategy for surface triangular meshes with associated distributed data. To expose the C++ code to R, the API provided by `Rcpp` and `RcppEigen` is used.

The tree structure of the package root folder is compliant with the standard. In detail:
- `src` contains C++ header and source files;
- `R` contains R source files;
- `demo` contains scripts examplifying the use of the package;
- `data` contains R datasets storing three-dimensional surface grids;
- `inst/extdata` contains three-dimensional surface grids in AVS UCD file format.

To install the package, first make sure that the package `devtools` is correctly installed on the current workstation. Then, from the package root folder type:

	R -e "library(devtools); install()" --silent 
	
Administrator privileges may be required.
