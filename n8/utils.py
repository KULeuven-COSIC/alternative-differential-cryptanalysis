import itertools as itt

from sage.all import *

def id_matrix(n):
    """
    Couldn't find that as a sage builtin
    """
    return Matrix(GF(2), n).parent()(1)

class Column:
    def __init__(self, ls):
        """
        Helper class for a column of Theta, supporting additions and casting to
        operation matrix.
        Input:
        - ls: a list of n-d vectors in (F_2)^d
        """
        V = ls[0].parent()
        self.d = V.dimension()
        self.n = len(ls) + self.d
        self.V = V

        self.ls = ls
        self._M = None

    def is_zero(self):
        return all([vi == self.V(0) for vi in self.ls])

    def to_matrix(self):
        """
        Return the matrix M_i corresponding to self.
        """
        if self._M:
            return self._M

        n = self.n
        d = self.d
        M = id_matrix(n)

        for i in range(n-d):
            for j in range(d):
                M[i, n-d+j] = self.ls[i][j]

        self._M = M
        return M

    def to_magma(self):
        """
        Print the M_i in magma format
        """
        n, d = self.n, self.d
        out = f' := Matrix(GF(2), {n-d}, {d}, [\n'
        for j in range(n-d-1):
            out += '\t' + ', '.join([str(x) for x in self.ls[j]]) + ',\n'

        out += '\t' + ', '.join([str(x) for x in self.ls[n-d-1]]) + '\n'
        out += ']);\n\n'
        return out

    def __add__(self, other):
        if self == 0:
            return other
        if other == 0:
            return self

        assert self.d == other.d and self.n == other.n
        res = []
        for i in range(len(self.ls)):
            res.append(self.ls[i] + other.ls[i])
        return Column(res)

    def __radd__(self, other):
        return self.__add__(other)

    def __str__(self):
        s = f'Column for n={self.n}, d={self.d} with vectors {[i for i in self.ls]}'
        return s

    def __repr__(self):
        return self.__str__()

def _bi_to_theta(b_i, L):
    """
    Given a list of vectors in (F_2)^d, generate the matrix Theta defining the
    operation. Return Theta as a list of n-d columns.
    Input:
    - b_i: a list of (L^2-L)/2 vectors in (F_2)^d
    - L: n-d
    """
    b_i = list(b_i)

    V = b_i[0].parent()

    Theta = [ [V(0) for i in range(L)] for j in range(L) ]
    for i in range(L):
        for j in range(i+1, L):
            x = b_i.pop(0)
            Theta[i][j] = x
            Theta[j][i] = x

    assert b_i == [], "non empty list??"

    # Cast columns to Column object
    cols = []
    for j in range(L):
        ls = [Theta[i][j] for i in range(L)]
        cl = Column(ls)
        cols.append(cl)

    return cols

def gen_all_ops(n, d):
    """
    Generate all operations in dimension n with weak key space size d,
    according to [https://doi.org/10.1007/s10623-018-0516-z, Thm 3.4]
    Non scalable bruteforce approach.
    An operation is a list of n-d columns, each defining (the non trivial part
    of) a weak component.
    """
    assert n-d >= 2, "wrong input"

    V = VectorSpace(GF(2), d)
    L = n-d
    n_elt = (L**2 - L)//2

    out = []
    for b_i in itt.product(V, repeat=n_elt):
        theta = _bi_to_theta(b_i, L)

        # Check that the operation is valid: no F_2-linear combination of columns
        # of theta is zero
        is_valid = True
        for r in range(1, L+1):
            if not is_valid:
                break

            for ss in itt.combinations(theta, r):
                ccomb = sum(ss)
                if ccomb.is_zero():
                    is_valid = False
                    break
        if is_valid:
            out.append(theta)
    return out

def gen_random_ops(n, d):
    """
    Generate random operations in dimension n with weak key space size d,
    according to [https://doi.org/10.1007/s10623-018-0516-z, Thm 3.4]
    An operation is a list of n-d columns, each defining (the non trivial part
    of) a weak component.
    """
    assert n-d >= 2, "wrong input"

    V = VectorSpace(GF(2), d)
    L = n-d
    n_elt = (L**2 - L)//2

    while True:
        b_i = [V.random_element() for _ in range(n_elt)]
        theta = _bi_to_theta(b_i, L)

        # Check that the operation is valid: no F_2-linear combination of columns
        # of theta is zero
        M = Matrix(GF(2), 7)
        for i in range(len(theta)):
            c_i = theta[i].ls
            for j in range(len(theta)):
                M[j, i] = c_i[j][0]

        breakpoint()
        is_valid = True
        for r in range(1, L+1):
            if not is_valid:
                break

            for ss in itt.combinations(theta, r):
                ccomb = sum(ss)
                if ccomb.is_zero():
                    is_valid = False
                    break
        if is_valid:
            return theta

def alt_sum(a, b, theta=None):
    """
    Implementation from [https://doi.org/10.1007/s10623-018-0516-z, Sec 3] for
    an alternative operation over (F_2)^n.
    Input:
    - a, b: two vectors in (F_2)^n
    - theta: an alternative operation as a list of n-d Column objects
    """
    if not theta:
        return a+b

    V = a.parent()
    n = V.dimension()
    d = n - len(theta)

    out = V(0)

    w = 0 # Hamming weight parity of a

    for i in range(n):
        # Weak components
        if a[i] == 1:
            if i < n-d:
                # Weak components
                Mi = theta[i].to_matrix()
            else:
                # Regular components
                Mi = id_matrix(n)
            ei = V.basis()[i]
            out += b*Mi + ei
            w = 1 - w
    if w == 0:
        out += b
    return out

def theta_to_magma(theta):
    out = ''
    for i in range(len(theta)):
        out += f'c{i+1} ' + theta[i].to_magma()
    out += 'theta := [' + ', '.join([f'c{i+1}' for i in range(len(theta))]) + '];\n'
    return out


def vec2int(v):
    """
    Maps (F_2)^n -> N, lsb on the right
    """
    v = list(v)[::-1]
    out = 0
    for i in range(len(v)):
        out += 2**i * ZZ(v[i])
    return out

def int2vec(a, n):
    """
    Maps N -> (F_2)^n, lsb on the right
    """
    V = VectorSpace(GF(2), n)
    b = bin(a)[2:].zfill(n)
    return V([int(x) for x in b])

if __name__ == "__main__":
    # Testing the sum
    assert len(gen_all_ops(5, 2)) == 42
    assert len(gen_all_ops(6, 4)) == 15

    # Theta to magma
    # ops = gen_all_ops(8, 6)
    # print(theta_to_magma(ops[0]))

    # Generate all ops for given n and d
    # n = 3
    # ops = gen_all_ops(3, 1)
    # theta = ops[0]

    # for _ in range(10):
    #     a, b = [randint(0, 7) for _ in range(2)]
    #     print(f'{a = } {b = }')
    #     a, b = [int2vec(i, n) for i in [a, b]]

    #     v1 = vec2int(a+b)
    #     v2 = vec2int(alt_sum(a, b, theta))
    #     print(f'-> {v1 = } {v2 = }')

    # Generate random ops
    ops = [gen_random_ops(8, 5) for _ in range(20)]



