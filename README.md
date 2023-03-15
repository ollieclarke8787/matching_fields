# MatchingFields
MatchingFields is a package for Macaulay2. A matching field is a combinatorial object that parametrises certain toric degenerations of Grassmannians and partial flag varieties of type An. The main features of this package allow the user to: 
- construct a matching field,
- consruct a matching field polytope,
- check whether a matching field is coherent,
- construct a matching field ideal,
- check whether a matching field gives rise to a toric degeneration,
- compute the matroidal regular subdivision of the hypersimplex associated to a matching field for the Grassmannian (1).

The package is currently in development. Most of the core functionality is implemented. 
However, some functionality is currently undocumented.

## Notes:
1. A coherent matching field is naturally associated to a cone of the tropicalisation (of the open subset of the Grassmannian where the Pleucker coordinates are non-vanishing) and, therefore, is a subset of the corresponding Dressian, which (as a set) is the collection of weights that induce a regular matroidal subdivision of the hypersimplex (in the case of the Grassmannian).


# How to get the package

To get the package, you will need an up-to-date version of Macaulay2 (version 1.20+). 
The latest versions of Macaulay2 contain packages that are necessary for this package.
In particular, this package uses the packages: "Polyhedra", "Tropical", "Binomials", "SubalgebraBases", and "Matroids".
Use the following steps to get and install the package in Macaulay2:

1. Download the package
2. Open Macaulay2
3. Add the location of the downloaded package to the 'path':
path = join({"package_location_folder/"} , path)
4. Install the package:
installPackage "MatchingFields"

# Possible Issues

It is possible that some of the functions do not run because they require the package "FourTiTwo", 
which used for computing Groebner bases of toric ideals of matching fields. 
Some functions use it explicitly, so there are optional arguments such as in 'matchingFieldIdeal(... Strategy => "M2" ...) that can
be used to avoid running "FourTiTwo". However, other function may be running using the package implicitly. 
Please let me know if this is a problem.
