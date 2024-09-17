function alt_sum(a, b, theta)
	/*
	This function implements the algorithm of [Civino, Blondeau, Sala - Differential attacks: using alternative operations, Section 3]
	for an alternative operation over F_n. The inputs are:
	- a and b, two vectors in F_2^n to be summed;
	- theta is a list of n-d matrices, the matrices \Sigma_{e_i} in the notation of the paper, associated with the alternative
	  operation in the first n-d (non-weak) components; each theta[i] is a (n-d)*d matrix over F_2 that we then use to define the M_i.

	By definitio it holds a \circ e_i = aM_i + e_i, where M_i has theta[i] as upper right 
	block if i <= n-d, and is the identity otherwise. Then writing a = \Sum_i (a_i e_i) we have
	   a \circ b = \Sum(b \circ e_i)
	if the hamming weight of a is odd, and 
	   a \circ b = \Sum(b \circ e_i) + b
	if it is even. This function implements this logic and returns the sum a \circ b.
	*/

	// Recover sum parameters
	n := Ncols(a);
	d := n - #theta;

	V := VectorSpace(GF(2), n);
	ans := V!0;

	w := 0; // Hamming weight of a
	for i in [1..n] do
		if a[i] eq 1 then
			Mi := ScalarMatrix(GF(2), n, 1);
			// If i > n-d we need to add theta[i] to M_i
			if i le #theta then InsertBlock(~Mi, theta[i], 1, n-d+1); end if;
			ei := BasisElement(V, i);
			ans := ans + b*Mi + ei;
			w := 1 - w;
		end if;
	end for;

	// Check for the weight
	if w eq 0 then
		ans := ans + b;
	end if;
	return ans;
end function;

function cdot(a, b, theta)
	// \cdot function as defined in [Civino, Blondeau, Sala - Differential attacks: using alternative operations, Section 3.1.1]
	return a + b + alt_sum(a, b, theta);
end function;


function Uo(V, theta)
	// Computation of U_\circ for the sum defined by theta in the 
	// case d = n-3. Relies on [Blondeau, Civino, Sala, Remark 3.10]
	out := {};
	E := [x: x in V | Weight(x) eq 1];

	Include(~out,V!0);

	for u1 in V do
		for u2 in V do
			Include(~out, cdot(u1, u2, theta));
		end for;
	end for;

	return out;
end function;

function Ho(V, theta)
	// Given a sum \circ defined by theta returns the space H_\circ as
	// the intersection GL(+) \cap GL(\circ)
	
        Vset := {v : v in V};

	S := Sym(Vset);
	E := [x: x in V | Weight(x) eq 1];

	T := [map<V->V| x:->x+v> : v in E];
	tt := sub<S|[[t(v): v in Vset]:t in T]>;

	To := [map<V->V| x:->alt_sum(x, v, theta)> : v in E];
	tto := sub<S|[[t(v): v in Vset]:t in To]>;

	AGLt := Normalizer(S, tt);
	AGLo := Normalizer(S, tto);

	GLt := Stabilizer(AGLt, V!0);
	GLo := Stabilizer(AGLo, V!0);

	out := GLo meet GLt;
	return out;

end function;

function mat(lambda, V)
	// Given a linear map lambda over V returns
	// its representation as n x n matrix.
	// Sample usage: mat(Random(H0), V).
	n := Dimension(V);
	E := [x: x in V | Weight(x) eq 1];
	mlambda := Matrix(GF(2),n,n,[e^lambda:e in E]);
	return mlambda;
end function;
