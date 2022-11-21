# CMZ Polytopes
The file 'CMZ_toric_degens.m2' contains scripts to produce matching field polytopes described in the paper:
"Toric degenerations of partial flag varieties and combinatorial mutations of matching field polytopes"
joint with Fatemeh Mohammadi and Francesca Zaffalon.
The paper is available on the arXiv (https://arxiv.org/abs/2206.13975).
In the paper, we consider matching field polytopes for the Grassmannians Gr(3,6) and Gr(3,7) and the full flag varieties Fl4 and Fl5.

The scripts in 'CMZ_toric_degens.m2' produce files that can be loaded directly into Polymake (https://polymake.org/doku.php), for a guide on scripting in Polymake see: https://polymake.org/doku.php/user_guide/howto/scripting
The Polymake files are already pre-computed. The following list of files should be present:
- fl4MFPolymakePolytopes
- fl5MFPolymakePolytopes
- gr36AllC1PolymakePolytopes
- gr36AllC2PolymakePolytopes
- gr36AllC3PolymakePolytopes
- gr36AllC4PolymakePolytopes
- gr36AllC5PolymakePolytopes
- gr37AllC1PolymakePolytopes
- gr37AllC2PolymakePolytopes
- gr37AllC3PolymakePolytopes
- gr37AllC4PolymakePolytopes
- gr37AllC5PolymakePolytopes
- gr37AllC6PolymakePolytopes

The files "gr36..." and "gr37..." contain the matching field polytopes for Gr(3,6) and Gr(3,7), respectively. In each file's name is another value "...CX..." for some 1 <= X <= 6. The value of X is the scaling coefficient 'c' from the aforementioned paper. The files "fl4..." and "fl5..." contain the matching field polytopes for the full flag varieties Fl4 and Fl5, respectively.

## Polytope identifiers

The polytopes are indexed by a number 'i' from 0 to (n! - 1) for some n depending on which Grassmannian or flag variety matching field polytopes are loaded. Each indexing number 'i' corresponds to the permutation that appears in the 'i'^th position in Macaulay2 when using the function permutations(n). For example, if n = 3, then the permutations, listed in order, are: {{0, 1, 2}, {0, 2, 1}, {1, 0, 2}, {1, 2, 0}, {2, 0, 1}, {2, 1, 0}}. In this case: index 0 is {0, 1, 2}, index 1 is {0, 2, 1}, and so on. Therefore, the polytope '$p0' in the file 'gr36AllC1PolymakePolytopes' refers to the polytope associated to the permutation {0, 1, 2, 3, 4, 5} and scaling coefficient c = 1. 

In the files that list the matching field polytopes for full flag varieties, the value of the scaling coefficient is the first digit of the identifier. The subsequent digits of the identifier are as above. So, for example, the polytope '$q20' in the file 'fl4MFPolymakePolytopes' is the matching field polytope associated to the permutation {0, 1, 2, 3} (i.e. the 0^th permutation) and scaling coefficient 2.

