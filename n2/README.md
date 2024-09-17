# Case n-2

For the code in this folder we generally assume `d = n-2`. However, most of the functions should work for any value of `d` and `n`.

The file `utils.m` contains functions to operate with alternative sums. An example of usage is given below.

```
load 'utils.m';
n := 6;

Z2 := GF(2);
V := VectorSpace(Z2, n);

E := [x: x in V | Weight(x) eq 1];

// Alternative sum defined by (1, 0)
c1 := Matrix(Z2, 2, 2, [0, 0, 1, 0]);
c2 := Matrix(Z2, 2, 2, [1, 0, 0, 0]);
theta := [c1, c2]; 

// Space of compatible maps
H0 := Ho(V, theta);

// Represents a map as matrix
print mat(Random(H0), V);

// Represent binary vectors as integers
a := n2vec(63, n);
b := n2vec(42, n);
c := alt_sum(a, b, theta);

print vec2n(c);

// Xor can be represented by the zero error
theta_xor := [ZeroMatrix(Z2, 2), ZeroMatrix(Z2, 2)];
c_xor := alt_sum(a, b, theta_xor);

print vec2n(c_xor);
```

## Parallel sums

The file `parall_utils.m` provides implementations for an alternative sum acting in parallel on two components. It builds on top of `utils.m`, which is loaded. An example of usage is given below.

```
load 'parall_utils.m';

n := 4; // Dimension of a single component (the whole space has dimension 2n)

Z2 := GF(2);
V := VectorSpace(Z2,n);
E := [x: x in V | Weight(x) eq 1];

// Alternative sum defined (on each component) by b = (0, 1)
c1 := Matrix(Z2, 2, 2, [0, 0, 1, 0]);
c2 := Matrix(Z2, 2, 2, [1, 0, 0, 0]);
theta := [c1, c2];

// Space of compatible maps
H0 := Ho(V, theta);

// Space of compatible maps
print mat(Random(H0), V);

// The sum of elements has now two components
a1 := 12;
a2 := 6;
a := np2vec(a1, a2, n);

b1 := 7;
b2 := 8;
b := np2vec(b1, b2, n);

c := parall_sum(a, b, theta);
c1, c2 := vec2np(c);
print c1, c2;

// We can compare again the result with xor
theta_xor := [ZeroMatrix(Z2, 2), ZeroMatrix(Z2, 2)];
c_xor := alt_sum(a, b, theta_xor);

print vec2np(c_xor);
```

## Analysis of optimal S-Boxes

The file `diff_opt.m` contains code to compute differentiability of optimal S-boxes with respect to a chosen sum. The computation for each S-box takes about 1 minute on a standard laptop.

