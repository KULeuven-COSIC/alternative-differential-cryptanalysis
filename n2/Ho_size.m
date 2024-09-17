// Here we compute the size of Ho for all possible operations
// defined in dimension n = 6 with d = n - 3 = 3 in relation with
// the dimension of the dimension of their error space

load "utils.m";

n := 6;

size2 := {}; // All possible sizes of Ho for sums with dim(Uo) = 2
size3 := {}; // All possible sizes of Ho for sums with dim(Uo) = 3

V3 := VectorSpace(GF(2), 3);
for b12 in V3 do
	for b13 in V3 do 
		for b23 in V3 do 
			
			// Check that b12, b13, b23 generate a valid operation
			S := sub<V3 | [b12, b13, b23]>;
			d := Dimension(S);
			if d le 1 then continue; end if;

			// Generate the 3 multiplication matrices
			M1 := ZeroMatrix(GF(2), 3, 3);
			InsertBlock(~M1, b12, 2, 1);
			InsertBlock(~M1, b13, 3, 1);

			M2 := ZeroMatrix(GF(2), 3, 3);
			InsertBlock(~M2, b12, 1, 1);
			InsertBlock(~M2, b23, 3, 1);

			M3 := ZeroMatrix(GF(2), 3, 3);
			InsertBlock(~M3, b13, 1, 1);
			InsertBlock(~M3, b23, 2, 1);

			// Compute Ho
			theta := [M1, M2, M3]; 
			H0 := Ho(theta, n);

			// Save the size
			if d eq 2 then Include(~size2, #H0);
			elif d eq 3 then Include(~size3, #H0);
			end if;
				

		end for;
	end for;
end for;

size2; // { 49152 } 
size3; // { 86016 }
