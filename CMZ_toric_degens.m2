needsPackage "MatchingFields"

polymakePolytope = (pointsMatrix, polytopeIndex) -> (
    -- produces polymake code to for the polytope whose vertices are the columns
    -- of the pointsMatrix
    points = entries transpose pointsMatrix;		    
    output = toString polytopeIndex |" = new Polytope(POINTS => [";
    for row in points do (
	output = output | "[1,"; --homogeneous coordinates
	for entry in row do (
	    output = output | toString entry | ",";
	    );
	output = output | "],";
	);
    output = output | "]);";
    output
    );

-- construct the matching fields arising from permutations for Gr(3,6)
for c from 1 to 5 do (
    f = openOut("./CMZ_Polytopes/gr36AllC" | toString c | "PolymakePolytopes");
    f << "use application \"polytope\";" << endl;
    i = 0; -- polytope number counter
    for s in permutations(6) do (
        L := matchingFieldFromPermutation(3, 6, s, UsePrimePowers => true, ScalingCoefficient => c);
	P := matchingFieldPolytope L;
	f << "my $p" | polymakePolytope(vertices P, i) << endl;
	f << "print \"P" | toString i |": \", $p" | toString i | "->F_VECTOR, \"\\n\";" << endl;
	i = i + 1;
	);
    f << close;
    )

--  construct polytopes for Gr(3,7)
for c from 1 to 6 do (
    f = openOut("./CMZ_Polytopes/gr37AllC" | toString c | "PolymakePolytopes");
    f << "use application \"polytope\";" << endl;
    i = 0; -- polytope number counter
    for s in permutations(7) do (
        L := matchingFieldFromPermutation(3, 7, s, UsePrimePowers => true, ScalingCoefficient => c);
	P := matchingFieldPolytope L;
	f << "my $p" | polymakePolytope(vertices P, i) << endl;
	f << "print \"P" | toString i |": \", $p" | toString i | "->F_VECTOR, \"\\n\";" << endl;
	i = i + 1;
	);
    f << close;
    )

