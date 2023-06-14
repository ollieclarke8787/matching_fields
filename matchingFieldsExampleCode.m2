-- Section 3.1: Constructing matching fields
needsPackage "MatchingFields"

D = diagonalMatchingField(3, 6)
getTuples D

D' = diagonalMatchingField({1,2,3}, 6)
getTuples D'

L = matchingFieldFromPermutation(3, 6, {1,2,3,6,5,4})
getWeightMatrix L

L1 = grMatchingField matrix {{0,0,0,0,0,0}, {2,4,1,3,6,5}}
getTuples L1
getWeightMatrix L1

L2 = flMatchingField({1,2}, matrix {{0,0,0}, {3,1,2}})
getWeightMatrix L2
getTuples L2

L3 = grMatchingField(2, 4, {{1,2}, {1,3}, {4,1}, {2,3}, {4,2}, {3,4}})
isCoherent L3
getWeightMatrix L3

L4 = flMatchingField({1,2}, 3, {{{1}, {2}, {3}},  {{1,2}, {1,3}, {3,2}}})
isCoherent L4
getWeightMatrix L4


-- Section 3.2: Ideals and algebras of matching fields
restart
needsPackage "MatchingFields";
D = diagonalMatchingField(2, 4);
I = matchingFieldIdeal D

J = plueckerIdeal D
ideal leadTerm(1, J) == I
isToricDegeneration D

S = plueckerAlgebra D
transpose gens S
transpose subalgebraBasis S

M = matrix {{0,0,0,0,0,0},{18,3,15,6,9,12},{35,28,21,14,7,0}};
L = grMatchingField M;
T = plueckerAlgebra L;
numgens T
numgens sagbi T

-- Section 3.3: Other functions
restart
needsPackage "MatchingFields";
L = grMatchingField matrix {{0,0,0,0,0,0},{18,3,15,6,9,12},{35,28,21,14,7,0}};
P = matchingFieldPolytope L
vertices P
(volume P)*(dim P)!

Q = NOBody L
vertices Q
(volume Q)*(dim Q)!

(vertices Q)_{0 .. 19} == vertices P

L = grMatchingField matrix {{0,0,0,0,0},{1,3,2,5,4},{10,0,20,40,30}};
getWeightPluecker L
netList matroidSubdivision L

L = diagonalMatchingField(2, 6);
algebraicMatroid L
#algebraicMatroidBases L
netList (algebraicMatroidCircuits L)_{0 .. 6}

L = grMatchingField(3, 5, {{1,3,2}, {1,4,2}, {1,5,2}, {3,4,1}, {1,3,5}, {1,4,5}, {3,4,2}, {2,3,5}, {2,4,5}, {3,4,5}})
T = topeField L
isLinkage T
T2 = amalgamation(2, T)
getTuples T2
T23 = amalgamation(3, T2)
getTuples T23
