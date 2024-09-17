// Compute the delta diff of optimal sboxes

load "utils.m";

function delta_diff(f, theta)
	// Function to compute the uniform differentiability with
	// respect to a sum \circ given by theta
	V := VectorSpace(GF(2), 4); // Working over V = \F_2^4
	Vs := [v: v in V];
	A := AssociativeArray();
	for i in [1..#Vs] do
		x := Vs[i]; 
		fx := f(x);
		for j in [i+1..#Vs] do
			y := Vs[j]; 
			fy := f(y);
			dt := alt_sum(x, y, theta);
			idt := vec2n(dt); // x \circ y = delta_in
			k := [idt, vec2n(alt_sum(fx, fy, theta))]; // [delta_in, delta_out]
			b, val := IsDefined(A, k);
			if b then
				val +:= 2;
				A[k] := val;
			else
				A[k] := 2;
			end if;
		end for;
	end for;
	return Max({A[x] : x in Keys(A)});
end function;

function get_sboxes(S, H0_mat)
	// Given an S-box and the space H0, compute all the equivalen S-boxes
	// with (potentially) different differential uniformity
	// S is a map, while H0_mat is H_o as a matrix subspace of GL(4,2)

	G := GL(4, 2);
	V := VectorSpace(GF(2), 4);
	
	// Compute the left equivalence classes
	gLeft := {};
	visited := {};
	for g in G do
		if g in H0_mat or g in visited then continue; end if;
		Include(~gLeft, g);

		for h in H0_mat do
			Include(~visited, h*g);
		end for;
	end for;

	// Compute the right equivalence classes
	gRight := {};
	visited := {};
	for g in G do
		if g in H0_mat or g in visited then continue; end if;
		Include(~gRight, g);

		for h in H0_mat do
			Include(~visited, g*h);
		end for;
	end for;

	// Generate all the sboxes
	sboxes := [];
	G1 := [];
	G2 := [];
	for g1 in gLeft do 
		for g2 in gRight do
			sbox := pmap<VectorSpace(GF(2),4) -> VectorSpace(GF(2),4)| [v -> (S(v*g1))*g2 : v in V] >;
			Append(~sboxes, sbox);
			Append(~G1, g1);
			Append(~G2, g2);
		end for;
	end for;
	// Return also G1 and G2 as the elements of GL such that
	// S(v * G1[i]) * G2[i] = sboxes[i]
	return sboxes, G1, G2;
end function;

// ---------------------------------------------------------
// Define the 16 optimal s-boxes as given in [Leander and Poschmann 2007, Table 6]
// Notice that the MAGMA notation is shifted, i.e. G0 will be G[1].
// ---------------------------------------------------------

G0 := [0, 1, 2,13, 4, 7,15, 6, 8, 11, 12, 9, 3, 14, 10, 5];
G1 := [0, 1, 2,13, 4, 7,15, 6, 8, 11, 14, 3, 5, 9, 10, 12];
G2 := [0, 1, 2,13, 4, 7,15, 6, 8, 11, 14, 3, 10, 12, 5, 9];
G3 := [0, 1, 2,13, 4, 7,15, 6, 8, 12, 5, 3, 10, 14, 11, 9];
G4 := [0, 1, 2,13, 4, 7,15, 6, 8, 12, 9, 11, 10, 14, 5, 3];
G5 := [0, 1, 2,13, 4, 7,15, 6, 8, 12, 11, 9, 10, 14, 3, 5];
G6 := [0, 1, 2,13, 4, 7,15, 6, 8, 12, 11, 9, 10, 14, 5, 3];
G7 := [0, 1, 2,13, 4, 7,15, 6, 8, 12, 14, 11, 10, 9, 3, 5];
G8 := [0, 1, 2,13, 4, 7,15, 6, 8, 14, 9, 5, 10, 11, 3, 12];
G9 := [0, 1, 2,13, 4, 7,15, 6, 8, 14, 11, 3, 5, 9, 10, 12];
G10:= [0, 1, 2,13, 4, 7,15, 6, 8, 14, 11, 5, 10, 9, 3, 12];
G11:= [0, 1, 2,13, 4, 7,15, 6, 8, 14, 11, 10, 5, 9, 12, 3];
G12:= [0, 1, 2,13, 4, 7,15, 6, 8, 14, 11, 10, 9, 3, 12, 5];
G13:= [0, 1, 2,13, 4, 7,15, 6, 8, 14, 12, 9, 5, 11, 10, 3];
G14:= [0, 1, 2,13, 4, 7,15, 6, 8, 14, 12, 11, 3, 9, 5, 10];
G15:= [0, 1, 2,13, 4, 7,15, 6, 8, 14, 12, 11, 9, 3, 10, 5];

G_n := [G0, G1, G2, G3, G4, G5, G6, G7, G8, G9, G10, G11, G12, G13, G14, G15];

V := VectorSpace(GF(2), 4);
G := [];

for Gi in G_n do 
	V_Gi := [n2vec(i, 4) : i in Gi];
	m := map< V -> V | x :-> V_Gi[vec2n(x) + 1]>;
	Append(~G, m);
end for;

// Sanity check: theta = 0 means that \circ is +
theta := [ZeroMatrix(GF(2), 2, 2), ZeroMatrix(GF(2), 2, 2)]; 
for i in [1..16] do
	dt_diff := delta_diff(G[i], theta);
	assert dt_diff eq 4; // Optimality
end for;

// Here we go with the actual sum defined by (1, 0)
V := VectorSpace(GF(2), 4);
c1 := Matrix(GF(2), 2, 2, [0, 0, 1, 0]);
c2 := Matrix(GF(2), 2, 2, [1, 0, 0, 0]);
theta := [c1, c2];
H0 := Ho(V, theta);

// We need the explicit elements of H0
H0_gen := Generators(H0);
E := [x: x in V | Weight(x) eq 1];
H0_mat := sub< GL(4,2) | [Matrix(GF(2),4, 4, [E[i]^h: i in [1..4]]): h in H0_gen]>;

for i in [1..16] do 
	printf "Computing S-box G[%o]...\n", i;
	sboxes, G1, G2 := get_sboxes(G[i], H0_mat);
	print "Generation done. Computing differentiability.";
	
	count_diff := [0,0,0,0,0,0,0,0];
	for ss in sboxes do 
		dd := delta_diff(ss, theta);
		count_diff[dd div 2] +:= 1; 
	end for;
	printf "Differentiability for G_%o: %o\n", i, count_diff;
end for;
