// Compute the delta diff of optimal sboxes load "utils.m";
load "AES.magma";
load "utils.m";

function delta_diff(f, theta)
	// Function to compute the uniform differentiability with
	// respect to a sum \circ given by theta
	V := VectorSpace(GF(2), 8);
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

V := VectorSpace(GF(2), 8);
c1 := Matrix(GF(2), 2, 6, [
    0, 0, 0, 0, 0, 0,
    1, 0, 0, 0, 0, 0
]);
c2 := Matrix(GF(2), 2, 6, [
    1, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0
]);
theta := [c1, c2];
dd := delta_diff(SBox, theta);
print dd;

