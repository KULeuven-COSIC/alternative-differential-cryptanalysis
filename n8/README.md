# Experiments on 8 bits sboxes

Running experiments on 8-bits sboxes with alternative operations.

To run the code just run `sage --python -O diff.py`. The code runs `magma` in
parallel to speed up the computation.For dimension `n = 8` and different `d`:
- prepare `template.m` with the correct input: load the file with the desired
  sbox and save it under the name `SBox`;
- generate a suitable number of operations; for `d = 7` and `d = 6` all
  operations can be generated with `gen_all_ops(n, d)`, for lower `d` it is
  better to use `gen_random_ops(n, d)` since the number of total operations
  blows up quickly
- then run the code sequentially or in parallel, uncommenting the corresponding
  lines
