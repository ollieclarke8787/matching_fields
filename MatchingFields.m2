
-- Matching Fields package for macaulay2 by Oliver Clarke

newPackage(
    "MatchingFields",
    Version => "0.1",
    Date => "November 6, 2022",
    Authors => {
	{Name => "Oliver Clarke", Email => "oliver.clarke@ed.ac.uk", HomePage => "https://www.oliverclarkemath.com/"}
	},
    Headline => "Matching Fields in Macaulay2",
    Keywords => {"Grassmannian", "Flag Variety", "Polytopes", "Toric Degeneration", "SAGBI Basis"},
    DebuggingMode => false,
    PackageExports => {"Polyhedra", "Tropical", "Binomials", "SubalgebraBases", "Matroids"}
    )  

-- ###########
-- # Exports #
-- ###########

export {
    "GrMatchingField",
    "grMatchingField",
    "FlMatchingField",
    "flMatchingField",
    "getTuples",
    "getWeightMatrix",
    "getWeightPleucker",
    "getGrMatchingFields",
    "matchingFieldPolytope",
    "ExtraZeroRows",
    "diagonalMatchingField",
    "matchingFieldRingMap",
    "matchingFieldIdeal",
    "pleuckerIdeal",
    "pleuckerMap",
    "matchingFieldFromPermutation",
    "RowNum",
    "UsePrimePowers",
    "ScalingCoefficient",
    "PowerValue",
    "isToricDegeneration",
    "NOBody",
    "matroidSubdivision",
    "weightMatrixCone",
    "isCoherent",
    "linearSpanTropCone",
    "VerifyToricDegeneration",
    "algebraicMatroid",
    "algebraicMatroidBases",
    "algebraicMatroidCircuits"
    }



-- #############
-- # Main Code #
-- #############


-- Matching Field Types
GrMatchingField = new Type of HashTable;
-- Grassmannian Matching Fields MF have:
-- MF.n
-- MF.k
-- MF.tuples = List of k-subets of {1 .. n}  
-- MF.cache

FlMatchingField = new Type of HashTable;
-- Flag Matching Fields MF have:
-- MF.kList
-- MF.grMatchingFieldList
-- MF.cache

-- The cache table contains some of the following (or eventually computed):
-- weightMatrix
-- weightPleucker
-- mfPolytope 
-- matchingFieldIdeal
-- matchingFieldRingMap
-- ringP
-- ringX
-- X
-- mfSubring
-- mfNOBody
-- mfPleuckerIdeal

protect symbol tuples
protect symbol n
protect symbol k
protect symbol kList
protect symbol grMatchingFieldList

protect symbol weightMatrix
protect symbol weightPleucker
protect symbol mfPolytope

protect symbol ringP -- Polynomial ring in variables P_I, I in subsets(n, k)
protect symbol ringX -- Polynomial ring in variables x_(i,j), 1 <= i <= k, 1 <= j <= n
protect symbol X -- matroid with ringX variables
protect symbol mfRingMap
protect symbol pleuckerRingMap
protect symbol mfIdeal
protect symbol mfPleuckerIdeal
protect symbol mfSubring

protect symbol mfNOBody

protect symbol computedMatroidSubdivision
protect symbol computedWeightMatrixCone
protect symbol computedLinearSpanTropCone
protect symbol computedAlgebraicMatroid

---------------------------------------------
-- Matching Field constructor
grMatchingField = method()

-- MF from weight matrix
grMatchingField(Matrix) := M -> (
    -- uses min convention   
    -- returns matching f ield for weight matrix along with the weight for the Plucker variables
    -- NB we assume the MF is well defined from the weight matrix
    --    i.e. the minimum was uniquely attained for each plucker form
    Mk := numRows M;
    Mn := numColumns M;
    local subsetOrder;
    local subsetWeight;
    local weight;
    L := {};
    W := {};
    for I in subsets(Mn, Mk) do (
	subsetWeight = infinity;
	for ordering in permutations(I) do (
	    weight = sum for i from 0 to Mk-1 list M_(i, ordering_i);
	    if weight < subsetWeight then (
		subsetOrder = apply(ordering, i -> i + 1);
 		subsetWeight = weight;
 		);
	    );
	L = append(L, subsetOrder);
	W = append(W, subsetWeight);
	);
    MF := new GrMatchingField from {
	n => Mn, 
	k => Mk, 
	tuples => L, 
	cache => new CacheTable from {
	    weightMatrix => M,
	    weightPleucker => W
	    }
	}
    )

-- Matching Field from list of tuples
grMatchingField(ZZ, ZZ, List) := (Lk , Ln, L) -> (
    -- check user input:
    sortedTuples := sort (sort \ L);
    usualSubsets := subsets(1 .. Ln, Lk);
    if not (sortedTuples == sort usualSubsets) then (
	error("Unexpected tuples.");
	);
    -- put the subsets in the correct order:
    lookupPosition := new HashTable from for i from 0 to #usualSubsets - 1 list (usualSubsets_i => i);
    tupleList := sort(L, x -> lookupPosition#(sort(x)));
    new GrMatchingField from {
	n => Ln,
	k => Lk,
	tuples => tupleList,
	cache => new CacheTable from {}
	}
    )

-- Constructor for parital Flag Matching Fields
flMatchingField = method()
flMatchingField(List, Matrix) := (inputKList, inputWeightMatrix) -> (
    if numRows inputWeightMatrix < max inputKList then (
	error("expected a matrix with at least " | toString max inputKList | " rows");
	);
    grMatchingFields := for Mk in inputKList list grMatchingField(inputWeightMatrix^(toList(0 .. Mk - 1)));
    new FlMatchingField from {
	n => numColumns inputWeightMatrix, 
	kList => inputKList, 
	grMatchingFieldList => grMatchingFields, 
	cache => new CacheTable from {
	    weightMatrix => inputWeightMatrix,
	    weightPleucker => flatten for grMF in grMatchingFields list grMF.cache.weightPleucker
	    }
	}
    )

flMatchingField(Matrix) := inputWeightMatrix -> (
    flMatchingField(toList(1 .. numRows inputWeightMatrix), inputWeightMatrix)
    )

-- Flag matching field from list of tuples
-- tuple list should be a list of lists
flMatchingField(List, ZZ, List) := (LkList, Ln, L) -> (
    grMatchingFields := for kIndex from 0 to #LkList - 1 list (
	grMatchingField(LkList_kIndex, Ln, L_kIndex) 
	);
    new FlMatchingField from {
	n => Ln,
	kList => LkList,
	grMatchingFieldList => grMatchingFields,
	cache => new CacheTable from {}
	}
    )

net(GrMatchingField) := MF -> (
    "Grassmannan Matching Field for Gr(" | toString MF.k | ", " | toString MF.n | ")"  
    )

net(FlMatchingField) := MF -> (
    s := toString MF.kList;
    "Flag Matching Field for Fl(" | s_(1, #s - 2) | "; " | toString MF.n | ")"
    )

-----------------------
-- Setup weight vectors
-- unexported 
-- Called by weight-getters
setupWeights = method()
setupWeights(GrMatchingField) := MF -> (
    if not MF.cache.?weightMatrix then (
    	MF.cache.weightMatrix = computeWeightMatrix MF;
	);
    if not MF.cache.?weightPleucker then (
	MF.cache.weightPleucker = for tuple in MF.tuples list sum(for i from 0 to MF.k - 1 list MF.cache.weightMatrix_(i, tuple_i - 1));
	);
    )

setupWeights(FlMatchingField) := MF -> (
    if not MF.cache.?weightMatrix then (
    	MF.cache.weightMatrix = computeWeightMatrix MF;
	);
    if not MF.cache.?weightPleucker then (
	MF.cache.weightPleucker = flatten for grMF in MF.grMatchingFieldList list (
	    for tuple in grMF.tuples list sum(
		for i from 0 to grMF.k - 1 list MF.cache.weightMatrix_(i, tuple_i - 1)
		)
	    );
	);
    )


-- basic getters:
getTuples = method()
getTuples(GrMatchingField) := MF -> (
    MF.tuples
    )
getTuples(FlMatchingField) := MF -> (
    for grMF in MF.grMatchingFieldList list grMF.tuples
    )


getGrMatchingFields = method()
getGrMatchingFields(FlMatchingField) := MF -> (
    MF.grMatchingFieldList
    )


getWeightMatrix = method()
getWeightMatrix(GrMatchingField) := MF -> (
    if not MF.cache.?weightMatrix then (
    	setupWeights MF;
	);
    MF.cache.weightMatrix
    )

getWeightMatrix(FlMatchingField) := MF -> (
    if not MF.cache.?weightMatrix then (
	setupWeights MF;
	);
    MF.cache.weightMatrix
    )


getWeightPleucker = method()
getWeightPleucker(GrMatchingField) := MF -> (
    if not MF.cache.?weightPleucker then (
    	setupWeights MF;
	);
    MF.cache.weightPleucker
    )

getWeightPleucker(FlMatchingField) := MF -> (
    if not MF.cache.?weightPleucker then (
	setupWeights MF;
	);
    MF.cache.weightPleucker
    )

---------------------------
-- Matching Field polytope
-- The polytope with one vertex for each tuple of the matching field
-- given by the convex hull of the exponent vectors of the monomial map
--
-- Options:
-- ExtraZeroRows: adds this many rows of 0's to each vertex (thought of as a k by n matrix) 
--     	       	  of the polytope (used for constructing flag polytopes)
-- 
matchingFieldPolytope = method(
    Options => {
	ExtraZeroRows => 0 
	}
    )
matchingFieldPolytope(GrMatchingField) := opts -> MF -> (
    if MF.cache.?mfPolytope then (
	MF.cache.mfPolytope
	) else ( 
	-- construct a matching field polytope P_L from 
	--   L a matching field for Gr(k, n) Grassmannian 
    	points := {};
    	for I in MF.tuples do (
	    -- construct the point corresponding to I 
	    point := ();
	    for i in I do (
	    	point = point | (i - 1 : 0) | (1 : (1)) | (MF.n - i : 0);
	    	);
	    point = point | (MF.n * opts.ExtraZeroRows : 0);
	    point = toList point;
	    points = append(points, point);
	    ); 
    	points = transpose matrix points;
    	P := convexHull points;
    	if opts.ExtraZeroRows == 0 then MF.cache.mfPolytope = P;
	P    
	)
    )

matchingFieldPolytope(FlMatchingField) := opts -> MF -> (
    if opts.ExtraZeroRows == 0 and MF.cache.?mfPolytope then (
	MF.cache.mfPolytope
	) else (
	P := sum for grMF in MF.grMatchingFieldList list matchingFieldPolytope(grMF, 
	    ExtraZeroRows => (max MF.kList - grMF.k + opts.ExtraZeroRows)
	    );
	if opts.ExtraZeroRows == 0 then MF.cache.mfPolytope = P;
	P 
	)
    )

--------------------------
-- Diagonal matching field
-- Corresponds to the Gelfand-Tsetlin Cone of the Flag Variety / Grassmannian
diagonalMatchingField = method()
diagonalMatchingField(ZZ, ZZ) := (Lk, Ln) -> (
    M := matrix for i from 0 to Lk - 1 list for j from 0 to Ln - 1 list i*(Ln - j);
    grMatchingField(M)
    )

-- partial flag variety
diagonalMatchingField(List, ZZ) := (LkList, Ln) -> (
    M := matrix for i from 0 to max LkList - 1 list for j from 0 to Ln - 1 list i*(Ln - j);
    flMatchingField(LkList, M)
    )

-- full-flag variety 
-- by convention this is a (n-1) by (n) matrix
diagonalMatchingField(ZZ) := Ln -> (
    M := matrix for i from 0 to Ln - 2 list for j from 0 to Ln - 1 list i*(Ln - j);
    flMatchingField(M)
    )

------------------------
-- setting up the polynomials rings for the matching field
-- unexported method
-- The weight vector v / matrix m stored in the matching field uses min convention
-- so use weight vector {max v .. max v} - v for the M2 weight vector
--  and the same for each row of m
setupMatchingFieldRings = method()
setupMatchingFieldRings(GrMatchingField) := MF -> (
    local monomialOrder;
    if not MF.cache.?ringP then (
	p := symbol p;
	variables := for s in subsets(toList(1 .. MF.n), MF.k) list p_(if #s == 1 then s_0 else toSequence s);
	monomialOrder = (
	    maxVal := max (getWeightPleucker MF);
	    {Weights => for val in (getWeightPleucker MF) list maxVal - val} 
	    );
	MF.cache.ringP = QQ[variables, MonomialOrder => monomialOrder];
	);
    if not MF.cache.?ringX then (
	x := symbol x;
	monomialOrder = (
	    weights := for wRow in entries getWeightMatrix MF list (
		bigVal := toList(MF.n : (max wRow));
		bigVal - wRow
		); 
	    {Weights => flatten weights} 
	    );
	MF.cache.ringX = QQ[x_(1,1) .. x_(MF.k, MF.n), MonomialOrder => monomialOrder];
	);
    if not MF.cache.?X then (
	MF.cache.X = transpose genericMatrix(MF.cache.ringX, MF.n, MF.k);
	);
    )

setupMatchingFieldRings(FlMatchingField) := MF -> (
    local monomialOrder;
    if not MF.cache.?ringP then (
	p := symbol p;
	variables := flatten for Lk in MF.kList list for s in subsets(toList(1 .. MF.n), Lk) list p_(if #s == 1 then s_0 else toSequence s);
	monomialOrder = (
	    bigVals := flatten (
		currentIndex := 0;
		for grMF in MF.grMatchingFieldList list (
		    numberEntries := binomial(grMF.n, grMF.k);
		    currentIndex = currentIndex + numberEntries;
		    toList(numberEntries : max ((getWeightPleucker MF)_{currentIndex - numberEntries .. currentIndex - 1}))
		    )
	    	);
	    {Weights => (bigVals - (getWeightPleucker MF))} 
	    );
	MF.cache.ringP = QQ[variables, MonomialOrder => monomialOrder];
	);
    if not MF.cache.?ringX then (
	x := symbol x;
	monomialOrder = (
	    weights := for wRow in entries getWeightMatrix MF list (
		bigVal := toList(MF.n : (max wRow));
		bigVal - wRow
		); 
	    {Weights => flatten weights} 
	    );
	MF.cache.ringX = QQ[x_(1,1) .. x_(max MF.kList, MF.n), MonomialOrder => monomialOrder];
	);
    if not MF.cache.?X then (
	MF.cache.X = transpose genericMatrix(MF.cache.ringX, MF.n, max MF.kList);
	);
    )


-- gets the sign of a tuple
-- which is (-1) to power the number of descets modulo 2
-- 1 means even tuple, -1 means odd tuple
tupleSign = method()
tupleSign(List) := I -> (
    if #I <= 1 then 1 else (
    	(-1)^((sum for s in subsets(I, 2) list if s_0 > s_1 then 1 else 0) % 2)
	)
    )

-- matching field ring map: P_I -> x_(1,I_1) * x_(2, I_2) ... x_(k, I_k), for each tuple I
matchingFieldRingMap = method()
matchingFieldRingMap(GrMatchingField) := MF -> (
    setupMatchingFieldRings(MF);
    if not MF.cache.?mfRingMap then (
    	R := MF.cache.ringP;
    	S := MF.cache.ringX;
        MF.cache.mfRingMap = map(S, R, 
	    for tuple in MF.tuples list tupleSign(tuple) * (product for i from 0 to MF.k - 1 list (MF.cache.X)_(i, tuple_i - 1))
	    );
	);
    MF.cache.mfRingMap
    )

matchingFieldRingMap(FlMatchingField) := MF -> (
    setupMatchingFieldRings(MF);
    if not MF.cache.?mfRingMap then (
    	R := MF.cache.ringP;
    	S := MF.cache.ringX;
    	MF.cache.mfRingMap = map(S, R, 
	    flatten for grMF in MF.grMatchingFieldList list (
		for tuple in grMF.tuples list tupleSign(tuple) * (product for i from 0 to grMF.k - 1 list (MF.cache.X)_(i, tuple_i - 1))
	    	)
	    );
	);
    MF.cache.mfRingMap
    )	     

-- matching field ideal
matchingFieldIdeal = method()
matchingFieldIdeal(GrMatchingField) := MF -> (
    -- setting up MF rings is done by grMatchingFieldRingMap if necessary
    if not MF.cache.?mfIdeal then (
    	MF.cache.mfIdeal = kernel matchingFieldRingMap(MF)
    	);
    MF.cache.mfIdeal
    )

matchingFieldIdeal(FlMatchingField) := MF -> (
    -- setting up MF rings is done by grMatchingFieldRingMap if necessary
    if not MF.cache.?mfIdeal then (
    	MF.cache.mfIdeal = kernel matchingFieldRingMap(MF)
    	);
    MF.cache.mfIdeal
    )


-- Grassmannian ideal using the constructed pleucker variable ring
-- Sets the weight of the polynomial ring to be the MF pleucker weight
Grassmannian(GrMatchingField) := opts -> MF -> (
    if not MF.cache.?mfPleuckerIdeal then (
    	setupMatchingFieldRings(MF);
    	R := MF.cache.ringP;
    	MF.cache.mfPleuckerIdeal = Grassmannian(MF.k - 1, MF.n - 1, R);
	);
    MF.cache.mfPleuckerIdeal
    )


pleuckerIdeal = method()
pleuckerIdeal(GrMatchingField) := MF -> (
    Grassmannian(MF)
    )

pleuckerIdeal(FlMatchingField) := MF -> (
    if not MF.cache.?mfPleuckerIdeal then (
    	setupMatchingFieldRings(MF);
    	local i;
	local variableFromSubset;
	local generatorList;
	i = 0;
	variableFromSubset = new HashTable from flatten (
	    varsMatrix := vars MF.cache.ringP;
	    for grMF in MF.grMatchingFieldList list (
	    	for s in subsets(toList(1 .. grMF.n), grMF.k) list (
		    i = i + 1;
		    s => varsMatrix_(0, i-1)
		    )
	    	)
	    );
	--------------------------------
    	-- Grassmannian relations
	-- For example, see the Wiki-page on the Pleucker Embedding
	--
	-- TODO: remove the redundant generators
	-- E.g. for Gr(2,4) we get 4 copies of the same generator
	--
    	generatorList = flatten for grMF in MF.grMatchingFieldList list (
	    if grMF.k >= 2 and grMF.n - grMF.k >= 2 then (
	    	flatten for I in subsets(1 .. grMF.n, grMF.k - 1) list (
		    for J in subsets(1 .. grMF.n, grMF.k + 1) list (
			newGenerator := sum for jPosition from 0 to #J - 1 list (
			    j := J_jPosition;
			    if not member(j, I) then (
			    	IIndex := sort(I | {j});
			    	JIndex := delete(j, J);
			    	swapsToSortI := # for i in I list if i > j then i else continue;
				pI := variableFromSubset#IIndex;
			    	pJ := variableFromSubset#JIndex;
			    	(-1)^(jPosition + swapsToSortI)*pI*pJ
				) else continue
			    );
			if not zero(newGenerator) then (
			    newGenerator
			    ) else continue
			)
		    ) 
		) else continue
	    );
    	---------------------
	-- Incident Relations
	--
	generatorList = generatorList | flatten for grMFs in subsets(MF.grMatchingFieldList, 2) list (
	    grMF0 := grMFs_0;
	    grMF1 := grMFs_1;
	    flatten for I in subsets(1 .. grMF0.n, grMF0.k - 1) list (
		for J in subsets(1 .. grMF1.n, grMF1.k + 1) list (
		    sum for jPosition from 0 to #J - 1 list (
			j := J_jPosition;
			if not member(j, I) then (
			    IIndex := sort(I | {j});
			    JIndex := delete(j, J);
			    swapsToSortI := # for i in I list if i > j then i else continue;
			    pI := variableFromSubset#IIndex;
			    pJ := variableFromSubset#JIndex;
			    (-1)^(jPosition + swapsToSortI)*pI*pJ 
			    ) else continue
			)
		    )
		) 
	    );
	MF.cache.mfPleuckerIdeal = ideal(generatorList);
    	);    
    MF.cache.mfPleuckerIdeal
    )

----------------
-- Pleucker map is the determinantal map associated to Grassmannian / Flag variety
-- its kernel coincides with the pleucker ideal defined above 

pleuckerMap = method()
pleuckerMap(GrMatchingField) := MF -> (
    setupMatchingFieldRings(MF);
    if not MF.cache.?pleuckerRingMap then (
	R := MF.cache.ringP;
	S := MF.cache.ringX;
	matX := MF.cache.X;
	MF.cache.pleuckerRingMap = map(S, R, for s in subsets(MF.n, MF.k) list det(matX_s));
	);
    MF.cache.pleuckerRingMap
    )

pleuckerMap(FlMatchingField) := MF -> (
    setupMatchingFieldRings(MF);
    if not MF.cache.?pleuckerRingMap then (
	R := MF.cache.ringP;
	S := MF.cache.ringX;
	matX := MF.cache.X;
	MF.cache.pleuckerRingMap = map(S, R, flatten for grMF in MF.grMatchingFieldList list (
		for s in subsets(grMF.n, grMF.k) list det(matX_s^(toList(0 .. grMF.k - 1))))
		);
	);
    MF.cache.pleuckerRingMap
    )



-- matching field from permutation
-- Fix a permutation S, take a 'highly generic' weight matrix M
-- that induces the diagonal matching field 
-- Permute the 2nd row of M using S
-- See the paper: Clarke-Mohammadi-Zaffalon 2022
--
matchingFieldFromPermutation = method(
    Options => {
	RowNum => 2, -- which row to permute
	UsePrimePowers => false, -- Take N (in the definition of the weight matrix) to be a prime number
	ScalingCoefficient => 1, -- scale the permuted row by this coefficient, if > 2 then matrix may be non-generic unless prime power is true
	PowerValue => 0 -- Value of N in weight matrix, if supplied 0 then choose N to be n or nextPrime n depending on above options
	})
matchingFieldFromPermutation(ZZ, ZZ, List) := opts -> (Lk, Ln, S) -> (
    if # S != Ln or # set S < Ln then (
	error("expected a permutation of " | toString Ln | " distinct values");
	);
    if opts.ScalingCoefficient == 1 then (
	matchingFieldFromPermutationNoScaling(Lk, Ln, S, opts)
	) else (
	local N;
	local M;
	local W;
	if opts.PowerValue > 0 then (
	    N = opts.PowerValue;
	    ) else if opts.UsePrimePowers then (
	    N = nextPrime Ln;
	    ) else (
	    N = Ln
	    );
	if Lk == 1 then (
	    M = matrix {toList {Ln : 0}};
	    ) else (
    	    M = matrix {toList{Ln : 0}} || matrix for i from 1 to Lk - 1 list for j from 1 to Ln list (
	    	if i + 1 == opts.RowNum then (
	    	    (S_(j - 1))*opts.ScalingCoefficient*N^(i - 1)
	    	    ) else (
	    	    (Ln - j)*N^(i - 1)
	    	    )
	    	);
	    );
	grMatchingField M
	)
    );

-- The Flag matching field from permuting the second row
matchingFieldFromPermutation(List, ZZ, List) := opts -> (LkList, Ln, S) -> (
    if # S != Ln or # set S < Ln then (
	error("expected a permutation of " | toString Ln | " distinct values");
	);
    sortedLkList := sort LkList;
    grMatchingFields := for Lk in sortedLkList list matchingFieldFromPermutation(Lk, Ln, S, opts);
    lastGrMatchingField := grMatchingFields_(#sortedLkList - 1);
    new FlMatchingField from {
	n => Ln, 
	kList => sortedLkList, 
	grMatchingFieldList => grMatchingFields, 
	cache => new CacheTable from {
	    weightMatrix => lastGrMatchingField.cache.weightMatrix,
	    weightPleucker => flatten for grMF in grMatchingFields list grMF.cache.weightPleucker
	    }
	}
    );

-- matching field from permutation 
-- assume that there is no scaling coefficient so the tuples of the matching field can be written down quickly
-- unexported method (used by matchingFieldFromPermutation)
matchingFieldFromPermutationNoScaling = method(
    Options => {
	RowNum => 2, -- which row to permute
	UsePrimePowers => false, -- Take N (in the definition of the weight matrix) to be a prime number
	ScalingCoefficient => 1, -- scale the permuted row by this coefficient, if > 2 then matrix may be non-generic unless prime power is true
	PowerValue => 0 -- Value of N in weight matrix, if supplied 0 then choose N to be n or nextPrime n depending on above options
	})
matchingFieldFromPermutationNoScaling(ZZ, ZZ, List) := opts -> (Lk, Ln, S) -> (
    local IOrdered;
    local minIndex;
    local N;
    local M;
    local W;
    L := {};
    for I in subsets(Ln, Lk) do (
	if opts.RowNum <= Lk then (
	    -- find i in 0 .. rowNum-1 such that S_(I_i) is minimum
	    minIndex = 0;
	    for i from 1 to opts.RowNum - 1 do (
	    	if S_(I_i) < S_(I_minIndex) then (
		    minIndex = i;
		    );
	    	);
	    -- The elements I_0 .. I_(minIndex-1) are ordered in increasing order
	    IOrdered = for i from 0 to minIndex-1 list I_i + 1;
	    -- The next elements are I_(minIndex+1) .. I_(rowNum-1)
	    IOrdered = IOrdered | for i from minIndex+1 to opts.RowNum-1 list I_i + 1;
	    -- Then we get I_minIndex
	    IOrdered = append(IOrdered, I_minIndex + 1);
	    -- Then the rest of I in order
	    IOrdered = IOrdered | for i from opts.RowNum to Lk - 1 list I_i + 1;
	    L = append(L, IOrdered);
	    ) else (
	    L = append(L, apply(I, i -> i+1));
	    );
	);
    if opts.PowerValue > 0 then (
	N = opts.PowerValue;
	) else if opts.UsePrimePowers then (
	N = nextPrime Ln;
	) else (
	N = Ln
	);
    if Lk  == 1 then (
	M = matrix {toList {Ln : 0}};
	) else (
    	M = matrix {toList{Ln : 0}} || matrix for i from 1 to Lk - 1 list for j from 1 to Ln list (
	    if i + 1 == opts.RowNum then (
	    	(S_(j - 1))*N^(i - 1)
	    	) else (
	    	(Ln - j)*N^(i - 1)
	    	)
	    );
	);
    W = for I in L list (
	sum for i from 0 to Lk - 1 list M_(i, I_i - 1)
	);
    new GrMatchingField from {
	n => Ln, 
	k => Lk, 
	tuples => L, 
	cache => new CacheTable from {
	    weightMatrix => M,
	    weightPleucker => W
	    }
	}
    );



-- isToricDegeneration for a Matching Field
-- checks if the matching field ideal is equal to the initial ideal of the Grassmannian

isToricDegeneration = method ()
isToricDegeneration(GrMatchingField) := MF -> (
    (matchingFieldIdeal(MF) == ideal leadTerm(1, Grassmannian(MF)))
    )

isToricDegeneration(FlMatchingField) := MF -> (
    (matchingFieldIdeal(MF) == ideal leadTerm(1, pleuckerIdeal MF))
    )


-- subring of a matching field
-- the pleucker algebra inside inside a ring with term order
-- given by the weightMatrix
subring(GrMatchingField) := opts -> MF -> (
    if not MF.cache.?mfSubring then (
    	setupMatchingFieldRings(MF);
    	matX := MF.cache.X;
    	MF.cache.mfSubring = subring for s in subsets(MF.n, MF.k) list det(matX_s);
	);
    MF.cache.mfSubring
    )

subring(FlMatchingField) := opts -> MF -> (
    if not MF.cache.?mfSubring then (
    	setupMatchingFieldRings(MF);
    	matX := MF.cache.X;
    	MF.cache.mfSubring = subring flatten for grMF in MF.grMatchingFieldList list (
	    for s in subsets(grMF.n, grMF.k) list det(matX_s^(toList(0 .. grMF.k - 1)))
	    );
	);
    MF.cache.mfSubring
    )

-- Newton-Okounkov body for the Grassmannian from a matching field
NOBody = method()
NOBody(GrMatchingField) := MF -> (
    if not MF.cache.?mfNOBody then ( 
    	-- compute the initial algbera of the Pleucker algebra wrt the weight term order
    	initialAlgberaGens := first entries leadTerm subalgebraBasis subring MF; 
    	generatorExponents := apply(initialAlgberaGens, f -> (exponents(f))_0);
    	NOBodyVertices := apply(generatorExponents, v -> ((MF.k) / sum(v))*v); -- normalize the vertices
    	MF.cache.mfNOBody = convexHull transpose matrix NOBodyVertices;
	);
    MF.cache.mfNOBody
    )
-- TODO: NO-Body for partial flag varieties
-- >> add a grading matrix etc.

-----------------------
-- Regular Subdivision of a set of points
-- code is copied and modified from "Polyhedra" Package 
-- not exported
pointRegularSubdivision = method()
pointRegularSubdivision(Matrix, Matrix) := (points, weight) -> (
    -- Checking for input errors
    if numColumns weight != numColumns points or numRows weight != 1 then error("The weight must be a one row matrix with number of points many entries");
    P := convexHull(points || weight, matrix (toList(numRows points : {0}) | {{1}} ));
    F := select(faces (1,P), f -> #(f#1) == 0);
    V := vertices P;
    apply (F, f -> V_(f#0)^(toList(0 .. (numRows points - 1))))
    )
 
---------------------------------
-- matroidal subdivision from matching field
-- Take the Pleucker weight w of a matching field
-- Note that w lies in the Dressian
-- Compute the regular subdivision of the hypersimplex wrt w
matroidSubdivision = method()
matroidSubdivision(GrMatchingField) := MF -> (
    if not MF.cache.?computedMatroidSubdivision then (
    	SS := subsets(toList(1 .. MF.n), MF.k);
	hyperSimplex := sub(transpose matrix for s in SS list for i from 1 to MF.n list if member(i, s) then 1 else 0, QQ); -- sub to avoid entries in ZZ
	subdivisionPieces := pointRegularSubdivision(hyperSimplex, matrix {getWeightPleucker MF} );
	vertexLookup := new HashTable from for i from 0 to binomial(MF.n, MF.k) - 1 list hyperSimplex_{i} => SS_i;
	MF.cache.computedMatroidSubdivision = for piece in subdivisionPieces list for c from 0 to numColumns piece - 1 list vertexLookup#(piece_{c});
    	);
    MF.cache.computedMatroidSubdivision
    )


----------------------
-- weightMatrixCone
-- the cone whose interior points are weight matrices that induce the given matching field
--
-- TODO: check how to go between weightMatrixCone and the Tropicalisation
weightMatrixCone = method(
    Options => {
	ExtraZeroRows => 0 -- adds this many rows of 0 to each inequality, used for FlMatchingField cone
	}
    )
weightMatrixCone(GrMatchingField) := opts -> MF -> (
    if opts.ExtraZeroRows == 0 and MF.cache.?computedWeightMatrixCone then (
	MF.cache.computedWeightMatrixCone
	) else (
	-- form the matrix of inequalities A: such that the cone is Ax >= 0
	local inequalities;
	if MF.k > 1 then (
            inequalities = matrix ( 
	    	subsetList := subsets(1 .. MF.n, MF.k);
	    	flatten for i from 0 to binomial(MF.n, MF.k) - 1 list (
	    	    columnIndices := subsetList_i;
	    	    minimalTuple := MF.tuples_i;
		    for p in delete(minimalTuple, permutations columnIndices) list (
		    	-- the row vector the encodes: minimalTuple <= p
		    	-- E.g. if the minimal tuple is {1,2,3} then
		    	--      one of the inequalities is given by {1,2,3} <= {1,3,2}
		    	--      which cancels down further since 1 is in the same place
		    	for coord in (0,1) .. (MF.k - 1 + opts.ExtraZeroRows, MF.n) list (
		    	    sum {if coord_0 < MF.k and p_(coord_0) == coord_1 then 1 else 0, 
			    	if coord_0 < MF.k and minimalTuple_(coord_0) == coord_1 then -1 else 0}
		    	    )
		    	)
	    	    )
    	    	);
	    ) else (
	    inequalities = matrix {toList((MF.k + opts.ExtraZeroRows) * MF.n : 0)};
	    );
	C := coneFromHData(inequalities);   
	if opts.ExtraZeroRows == 0 then MF.cache.computedWeightMatrixCone = C;
        C
    	)
    )

weightMatrixCone(FlMatchingField) := opts -> MF -> (
    if opts.ExtraZeroRows == 0 and MF.cache.?computedWeightMatrixCone then (
	MF.cache.computedWeightMatrixCone
	) else (
    	kMax := max MF.kList;
    	weightMatrixConeList := for grMF in MF.grMatchingFieldList list (
	    weightMatrixCone(grMF, ExtraZeroRows => (kMax - grMF.k + opts.ExtraZeroRows))
	    );
	inequalityMatrix := facets (weightMatrixConeList_0);
	hyperplanesMatrix := hyperplanes (weightMatrixConeList_0);
	for i from 1 to #weightMatrixConeList - 1 do (
	    inequalityMatrix = inequalityMatrix || facets (weightMatrixConeList_i);
	    hyperplanesMatrix = hyperplanesMatrix || hyperplanes (weightMatrixConeList_i);
	    );
	C := coneFromHData(inequalityMatrix, hyperplanesMatrix);
	if opts.ExtraZeroRows == 0 then MF.cache. computedWeightMatrixCone = C;
	C
	)
    )
----------------------------------------------------------
-- isCoherent
-- check if a matching field is induced by a weight matrix
--
isCoherent = method()
isCoherent(GrMatchingField) := MF -> (
    if MF.cache.?weightMatrix then true else (
	C := weightMatrixCone MF;
	-- coherent iff C is full-dimensional
	(dim C) == (MF.k * MF.n)
	)
    )
isCoherent(FlMatchingField) := MF -> (
    if MF.cache.?weightMatrix then true else (
	C := weightMatrixCone MF;
	-- coherent iff C is full-dimensional
	(dim C) == ((max MF.kList)* MF.n)
	)
    )

------------------------
-- computeWeightMatrix
-- finds a weight matrix that induces the matching field
-- unexported (see getWeightMatrix)
computeWeightMatrix = method()
computeWeightMatrix(GrMatchingField) := MF -> (
    if not isCoherent MF then (
	error("expected a coherent matching field");
	);
    C := weightMatrixCone MF;
    CRays := rays C;
    -- construct an interior point of the cone
    weight := first entries transpose sum for c from 0 to numColumns CRays - 1 list CRays_{c};
    matrix for i from 0 to MF.k - 1 list weight_{i*MF.n .. (i+1)*MF.n - 1}
    )

computeWeightMatrix(FlMatchingField) := MF -> (
    if not isCoherent MF then (
	error("expected a coherent matching field");
	);
    C := weightMatrixCone MF;
    CRays := rays C;
    -- construct an interior point of the cone
    weight := first entries transpose sum for c from 0 to numColumns CRays - 1 list CRays_{c};
    matrix for i from 0 to (max MF.kList) - 1 list weight_{i*MF.n .. (i+1)*MF.n - 1}
    )

--------------------------------------------
-- compute linear span of the tropical cone
--
-- assume the initial ideal is toric
-- take the generators x^u - x^v of the initial ideal
-- form a matrix M with rows {.. 1_u .. -1_v ..}
-- kernel M is the linear span of the tropical cone 

removeZeroRows = method()
removeZeroRows(Matrix) := inputMatrix -> (
    nonZeroRowIndices := for i from 0 to numRows inputMatrix - 1 list if not zero inputMatrix^{i} then i else continue;
    inputMatrix^nonZeroRowIndices
    )

linearSpanTropCone = method(
    Options => {
	VerifyToricDegeneration => true
	}
    )
linearSpanTropCone(GrMatchingField) := opts -> MF -> (
    if not MF.cache.?computedLinearSpanTropCone then (
    	if opts.VerifyToricDegeneration and not isToricDegeneration MF then (
	    error("expected GrMatchingField that gives a toric degeneration.");
	    );
    	matchingFieldIdealExponents := exponents \ first entries gens matchingFieldIdeal MF; -- a list of pairs (since generators are binomials)
    	constraintMatrix := matrix apply(matchingFieldIdealExponents, exponentPair -> exponentPair_0 - exponentPair_1);
    	constraintMatrix = removeZeroRows reducedRowEchelonForm sub(constraintMatrix, QQ);
    	MF.cache.computedLinearSpanTropCone = ker constraintMatrix;
	);
    MF.cache.computedLinearSpanTropCone
    )


---------------------------
-- compute the algebraic matroid of the Matching field inside Grassmannian
-- this gives bases for the algebraic matroid of the grassmannian implied by the matching field
algebraicMatroid = method()
algebraicMatroid(GrMatchingField) := MF -> (
    if not MF.cache.?computedAlgebraicMatroid then (
	MF.cache.computedAlgebraicMatroid = matroid transpose gens linearSpanTropCone MF;
	);
    MF.cache.computedAlgebraicMatroid
    )

-- write down the bases of the algebraic matroid as subsets 
algebraicMatroidBases = method()
algebraicMatroidBases(GrMatchingField) := MF -> (
    SS := subsets(toList(1 .. MF.n), MF.k);
    for B in bases algebraicMatroid MF list (i -> SS_i) \ B
    )

-- write down the bases of the algebraic matroid as subsets 
algebraicMatroidCircuits = method()
algebraicMatroidCircuits(GrMatchingField) := MF -> (
    SS := subsets(toList(1 .. MF.n), MF.k);
    for B in circuits algebraicMatroid MF list (i -> SS_i) \ B
    )




-- #################
-- # Documentation #
-- #################

beginDocumentation()

doc ///
      Key
        MatchingFields
      Headline
        A package for working with matching fields for Grassmannians and partial flag varieties
      Description
        Text
	  A matching field $\Lambda$ for the Grassmannian Gr($k$, $n$), is a simple combinatorial object.
	  It may be thought of as a choice of initial term for each maximal minor of a generic $k \times n$ matrix 
	  of variables. For example, take $k = 2$ and $n = 4$. Let $X = (x_{i,j})$ be a generic $2 \times 4$ matrix of variables. 
	  Suppose that a matching field $\Lambda$ has tuples $\{12, 31, 14, 32, 24, 34\}$. This means that $\Lambda$
	  distinguishes the term $x_{1,1} x_{2,2}$ from the maximal minors on columns $1$ and $2$ of $X$: $x_{1,1} x_{2,2} - x_{1,2} x_{2,1}$.
	  Similarly for the terms $x_{1,3} x_{2,1}$, $x_{1,1} x_{2,4}$, and so on.
	  
	  If the terms of all maximal minors distinguished by a matching field are their initial terms with respect to a fixed weight matrix,
	  then we say that the matching field is coherent. Each such weight matrix induces a weight vector on the Pleucker coordinates of the 
	  Grassmannian. If the initial ideal of the Pleucker ideal of the Grassmannian with respect to this weight vector is a toric ideal,
	  i.e. a prime binomial ideal, then we say that the matching field gives rise to a toric degeneration of the Grassmannian.
	  By a result of Sturmfels (1996), a matching field gives rise to a toric degeneration if and only if the maximal minors of $X$ form
	  a subalgebra basis (or SAGBI basis) with respect to the order induced by the weight matrix.
	  
	  This concept naturally generalises to partial flag varieties under the Pleucker embedding.
	  
	  The MatchingFields package gives basic functions, to construct many of the well-studied examples of matching fields.
	  Given a matching field $L$, it is straight forward to check whether $L$ is coherent, what is a weight matrix that induces it,
	  and whether is gives rise to a toric degeneration. The package also produces polytopes associated to matching fields and Newton-Okounkov bodies.
        Example
	  L = grMatchingField(2, 4, {{1,2}, {3,1}, {1,4}, {3,2}, {2,4}, {3,4}})
	  isCoherent L
	  getWeightMatrix L
	  isToricDegeneration L
      SeeAlso
      Subnodes
///

doc ///
      Key
        pleuckerIdeal
	(pleuckerIdeal, FlMatchingField)
	(pleuckerIdeal, GrMatchingField)
      Headline
        The Pleucker ideal of a matching field
      Usage
      	I = pleuckerIdeal L
	I = pleuckerIdeal LL
      Inputs
        L: GrMatchingField 
	LL: FlMatchingField
      Outputs
        I: Ideal
	  The Pleucker ideal associated to the corresponding Grassmannian or partial flag variety
	  with the correct term order given by a weight that induced the matching field
      Description
        Text
	  The Pleucker ideal is the defining ideal of a partial flag variety embedded in a product of Grassmannians, where
	  each Grassmannian is embedded, by the Pleucker embedding, into a suitable projective space.
	  In the case of the Grassmannian Gr($k$, $n$), it is concretely given by kernel of the ring map 
	  $K[P_I : I \subseteq [n],\  |I| = k] \rightarrow K[x_{i,j} : i \in [k], \ j \in [n]]$ where $P_I$ is mapped 
	  to the $k \times k$ maximal minor of the matrix $(x_{i,j})$ whose columns are indexed by the set $I$.
	  It is well-known that this ideal has a Groebner basis consisting of homogeneous quadrics.
	  
	  The function @TO "pleuckerIdeal"@ takes a matching field, either for the Grassmannian or a partial flag variety
	  and outputs the Pleucker ideal for that Grassmannian or partial flag variety. The ambient polynomial ring that
	  contains this ideal is constructed to have the term order induced by the matching field.
	  
	Example
	  L = grMatchingField(2, 4, {{1,2}, {1,3}, {1,4}, {2,3}, {2,4}, {3,4}})
	  I = pleuckerIdeal L
	  (monoid ring I).Options.MonomialOrder
	  getWeightPleucker L
	
	Text
	  In the above example, the weights for the ambient ring are not the same as the Pleucker weights of the matching field.
	  This is because of the minimum-maximum convention problem. For compatibility with packages such as @TO "Tropical"@, we use
	  the minimum convention in @TO "MatchingFields"@ so the smallest weight with respect to the weight matrix that 
	  induces the matching field is the initial term of a Pleucker form. 
	  However, the monomial ordering given by @TO "Weights"@ uses the 
	  maximum convention, so the ambient ring has weights that are based on the negative of the induced Pleucker Weight.    
	  
	  Note that the given matching field must be coherent. If the matching field is not defined in terms of a weight
	  matrix, then the function will attempt to compute a weight matrix for the matching field. If the matching field is
	  not coherent then the function will produce an error.
        Example
	  L = grMatchingField(2, 4, {{1,2}, {1,3}, {4,1}, {2,3}, {2,4}, {3,4}})
	  isCoherent L
	  -- I = pleuckerIdeal L -- "error: expected a coherent matching field"
      
      SeeAlso
      Subnodes
///


doc ///
      Key
        matroidSubdivision
	(matroidSubdivision, GrMatchingField)
      Headline
        The matroid subdivision induced by the Pleucker weight of a coherent matching field
      Usage
        listOfBases = matroidSubdivision L
      Inputs
        L: GrMatchingField 
      Outputs
        listOfBases: List
	  Each element is a list of the vertices of a maximal cell of the matroid subdivision of the hypersimplex induced by 
	  the Pleucker weight of the matching field.
      Description
        Text
	  The hypersimplex $\Delta(k, n) \subseteq \RR^{n}$ is the convex hull of the characteristic vectors of all $k$-subsets
	  of $\{1, \dots, n\}$, and we label each vertex with with its corresponding subset. A regular subdivision of the vertices of $\Delta(k, n)$
	  is said to be matroidal if, for each maximal cell of the subdivision, the subsets labelling its vertices form the set of bases of a matroid.
	  The well-known result is: a point lies in the Dressian Dr($k$, $n$), the tropical prevariety of all $3$-term Pleucker relation in Gr($k$, $n$), if and only if
	  it induces a matroidal subdivision of the hypersimplex.
	Example
	  L = grMatchingField(2, 4, {{1,2}, {1,3}, {1,4}, {2,3}, {2,4}, {3,4}})
	  netList matroidSubdivision L -- an octahedron sliced into 2 pieces
      SeeAlso
      Subnodes
///

doc ///
      Key
        algebraicMatroid
	(algebraicMatroid, GrMatchingField)
      Headline
        The algebraic matroid of the cone that induces the matroid
      Usage
        M = algebraicMatroid L
      Inputs
        L: GrMatchingField 
      Outputs
        M: "matroid"
	  The algebraic matroid of the cone in Trop Gr($k$, $n$) that induces the matching field.
      Description
        Text
	  Let $V \subseteq \CC^n$ be an affine variety. 
	  The algebraic matroid of $V$ is a matroid whose independent sets $S \subseteq [n]$
	  are the subsets such that the projection from $V$ to the coordinates indexed by $S$
	  is a dominant morphism. Similarly, if $C \subseteq \RR^n$ is a polyhedral cone, then the algebraic matroid
	  of $C$ is the matroid whose independent sets $S \subseteq [n]$ are the subsets such that image of the
	  projection of $C$ onto the coordinates indexed by $S$ is full-dimensional.
	  
	  In the case of the affine cone of Grassmannian under the Pleucker embedding, 
	  there are a few different ways to compute its algebraic matroid. One way is to use its tropicalization.
	  The algebraic matroid of the Grassmannian is equal to the matroid whose bases are the union of all bases of the
	  algebraic matroid for all maximal cones of Trop Gr($k$, $n$).
	  
	  For each coherent matching field, we compute its cone in the tropicalization of the Grassmannian.
	  We compute the algebraic matroid of this cone. To view the bases of this matroid in terms of the $k$-subsets of $[n]$,
	  use the function @TO "algebraicMatroidBases"@.
	  
	Example
	  L = grMatchingField(2, 4, {{1,2}, {1,3}, {1,4}, {2,3}, {2,4}, {3,4}})
	  M = algebraicMatroid L
	  netList algebraicMatroidBases L
      SeeAlso
      Subnodes
///








-- #########
-- # Tests #
-- #########

TEST ///
L = grMatchingField matrix {
    {0,0,0,0}, 
    {1,3,2,4}};
tupleList = {{2,1},{3,1},{2,3},{4,1},{4,2},{4,3}};
assert(getTuples L == tupleList);
assert(getWeightPleucker L == {1, 1, 2, 1, 3, 2});
assert(isToricDegeneration L);
///

TEST ///
L = grMatchingField(2, 3, {{1,2}, {2,3}, {3,1}});
assert(isCoherent L == false);
///

TEST ///
L = diagonalMatchingField(2, 6);
assert(dim weightMatrixCone L == 12);
assert(numColumns rays weightMatrixCone L == 5);
assert(numColumns linealitySpace weightMatrixCone L == 7);
///

TEST ///
L = matchingFieldFromPermutation(2, 4, {2,1,4,3});
L' = matchingFieldFromPermutation(2, 4, {3,2,10,5});
assert(getTuples L == getTuples L');
///


end --
restart
loadPackage "MatchingFields"

L = flMatchingField({1,2}, 4, {{{1}, {2}, {3}, {4}}, {{1,2},{1,3},{3,2},{4,3},{4,1},{4,2}}})
isToricDegeneration L
isCoherent L
dim weightMatrixCone L


L = flMatchingField({1,2}, 4, {{{1}, {2}, {3}, {4}}, {{1,2},{3,1},{2,3},{4,3},{4,1},{4,2}}})
isToricDegeneration L

peek L.cache
W = computeWeightMatrix L
L' = grMatchingField W

C = weightMatrixCone L
linealitySpace C

D = diagonalMatchingField({1,2,3}, 6)
weightMatrixCone D

I = pleuckerIdeal D

D = diagonalMatchingField(3, 6)
D = diagonalMatchingField(2, 4)

C = weightMatrixCone D
rays C
linealitySpace C

peek D.cache
pleuckerMap D
vertices NOBody D
matroidSubdivision D

S = subring D

P = matchingFieldPolytope(D)
vertices P
(volume P) * (dim(P))!

ID = matchingFieldIdeal D
peek D.cache
I = Grassmannian(D)
J = ideal leadTerm(1, I)

L = matchingFieldFromPermutation(2, 8, {8,6,4,2,7,5,3,1}, RowNum => 2)
transpose gens matchingFieldIdeal L
transpose leadTerm(1, Grassmannian(L))
isToricDegeneration L


peek L
peek L.cache

vertices grNOBody L

hexMF = grMatchingField(matrix {
	{0,0,0,0,0,0},
	{15,0,12,3,6,9},
	{35,28,21,14,7,0}}
    )

isToricDegeneration hexMF
vertices grNOBody hexMF

matroidSubdivision hexMF

--------------------
D = diagonalMatchingField({1,2,3,4}, 6)
transpose gens matchingFieldIdeal D

hexMF = flMatchingField matrix {
	{0,0,0,0,0,0},
	{15,0,12,3,6,9},
	{35,28,21,14,7,0}}

I = matchingFieldIdeal hexMF
peek hexMF.cache
transpose gens matchingFieldIdeal hexMF


L = matchingFieldFromPermutation({1,3}, 5, {2,3,1,4,0})
isToricDegeneration L
peek L.cache
I = matchingFieldIdeal L
J = ker pleuckerMap L
leadTerm(1, J)
transpose leadTerm(1, J)
I == ideal leadTerm(1, J)

degree I
vertices matchingFieldPolytope L
dim matchingFieldPolytope L
L
volume matchingFieldPolytope L
(volume matchingFieldPolytope L) * (dim matchingFieldPolytope L)!

degree matchingFieldIdeal diagonalMatchingField({1,2,3}, 6)


for p in permutations(5) list (
    MF = matchingFieldFromPermutation({3}, 5, p, RowNum => 3);
    if not isToricDegeneration MF then (
        p
	) else (
	continue;
	)
    )



D = diagonalMatchingField({3},6)
I = pleuckerIdeal D
J = Grassmannian(2, 5, D.cache.ringP)
isSubset(J, I)

(gens J) % I

debug MatchingFields
