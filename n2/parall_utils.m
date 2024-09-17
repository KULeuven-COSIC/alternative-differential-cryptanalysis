load "utils.m";

function splitv(v)
	// Given a vector v in (F_2)^n splits it in two vectors in (F_2)^{n/2}
	// n should be even
	n := Ncols(v) div 2;
	V := VectorSpace(GF(2), n);
	v1 := [];
	v2 := [];
	for i in [1 .. n] do
		Append(~v1, v[i]);
		Append(~v2, v[i+n]);
	end for;
	return V!v1, V!v2;
end function;

function mergev(v1, v2)
	// Given two vectors v1 and v2 returns the concatenation of v1 and v2
	n := Ncols(v1) + Ncols(v2);
	V := VectorSpace(GF(2), n);
	v := ElementToSequence(v1) cat ElementToSequence(v2);
	return V!v;
end function;

function parall_sum(a, b, theta)
	// Given two vectors a and b \in \F_2^(2n) applies a (2-)parallel sum
	// on the first and second half components and returns the result as a
	// single vector.
	a1, a2 := splitv(a);
	b1, b2 := splitv(b);
	c1 := alt_sum(a1, b1, theta);
	c2 := alt_sum(a2, b2, theta);
	c := mergev(c1, c2);
	return c;
end function;

function parall_cdot(a, b, theta)
	// \cdot function as defined in [Blodeau, Civino, Sala (2018) Section 3.1.1]
	return a + b + parall_sum(a, b, theta);
end function;

function parall_Uo(V, theta)
	// Parallel computation of U_\circ for the sum defined by theta in the 
	// case d = n-2. Relies on [Blondeu, Civino, Sala (2018) Remark 3.10]
	out := {};
	n := Dimension(V) div 2;
	E := [x: x in V | Weight(x) eq 1];

	Include(~out,V!0);

	u1 := parall_cdot(E[1], E[2], theta);
	u2 := parall_cdot(E[n+1], E[n+2], theta);

	Include(~out,u1);
	Include(~out,u2);
	Include(~out,u1+u2);
	return out;
end function;

function parall_Ho(V, theta)
	// Given a sum \circ defined by theta returns the space H_\circ as
	// the intersection GL(+) \cap GL(\circ)

	n := Dimension(V);
	Vset := {v : v in V};

	S := Sym(Vset);
	E := [x: x in V | Weight(x) eq 1];

	T := [map<V->V| x:-> x+v> : v in E];
	tt := sub<S|[[t(v): v in Vset]:t in T]>;

	To := [map<V->V| x:-> parall_sum(x, v, theta)> : v in E];
	tto := sub<S|[[t(v): v in Vset]:t in To]>;

	AGLt := Normalizer(S, tt);
	AGLo := Normalizer(S, tto);

	GLt := Stabilizer(AGLt, V!0);
	GLo := Stabilizer(AGLo, V!0);

	out := GLo meet GLt;
	return out;
end function;

function np2vec(a, b, n)
	// Given a pair of integer and the size of a single component
	// the concatenation of the corresponding vectors as a vector in (F_2)^(2n)
	va := n2vec(a, n);
	vb := n2vec(b, n);

	return mergev(va, vb);
end function;

function vec2np(v)
	// Given a vector splits it in a half and return the two resulting
	// vectors as integers
	v1, v2 := splitv(v);
	return vec2n(v1), vec2n(v2);
end function;
