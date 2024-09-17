# Case n-3

For the code in this folder we assume `d = n-3`.

The file `utils.m` contains functions to operate with alternative sums. The structure is similar to the case `d = n-2`, but it is included to simplify the setup. An example of usage is given below.

```
load 'utils.m';

n := 6;

Z2 := GF(2);
V := VectorSpace(Z2, n);

E := [x: x in V | Weight(x) eq 1];

// Definition of alternative sum
c1 := Matrix(Z2, 3, 3, [
	0, 0, 0,
	1, 0, 1,
	1, 1, 1]);

c2 := Matrix(Z2, 3, 3, [
	1, 0, 1,
	0, 0, 0,
	0, 0, 1]);

c3 := Matrix(Z2, 3, 3, [
	1, 1, 1,
	0, 0, 1,
	0, 0, 0]);

theta := [c1, c2, c3]; 

// Space of errors
u0 := Uo(V, theta);
u0;

// Space of compatible maps
H0 := Ho(V, theta);
print mat(Random(H0), V);
```

## Size of Ho for `d = n-3`

The file `Ho_size.m` contains code to compute the size of the space of compatible maps `Ho` when `d = n-3`, depending on the dimension of the error space.
