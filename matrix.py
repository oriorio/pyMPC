"""
Matrix - A skinny & improved version of https://github.com/Ismael-VC/pymatrix
VanderMatrix - a Vandermonde matrix, supporting fast inverse calculation
"""

from utils import inner_product, copy
from serializer import *
from itertools import izip, product
from operator import add, sub, mul

__all__ = ['Matrix', 'MatrixError']

class MatrixError(Exception):
    """Invalid operation attempted on a Matrix object"""
    pass

class Matrix(object):
    """A generic matrix object. Elements could be of any type supporting arithmetic operations.
    Supports various initialization forms:
        Matrix()
            Initialize a (0 x 0) matrix.
        Matrix(n_rows, n_cols, fill=0)
            Initialize a (@n_rows x @n_cols) matrix filled with @fill.
        Matrix(n)
            Initialize a (@n x @n) square matrix filled optional keyword argument @fill.
        Matrix.matrix([ [row0 values], [row1 vlaues], ... ])
            A 2d array where each inner array is of the same length.
        Matrix.matrix([diagonal values])
            A 1d array defines a diagonal matrix with values given set as the diagonal.
            Equivalent to Matrix.diag method.
        Matrix.matrix(n)
            An integer @n will return a (@n x @n) identity matrix.
            Equivalent to Matrix.identity method.
    """

    def __init__(self, n_rows=0, n_cols=None, fill=0):
        if n_cols is None:
            n_cols = n_rows
        self._nrows = n_rows
        self._ncols = n_cols
        self._grid = [[copy.deepcopy(fill) for i in xrange(n_cols)] for j in xrange(n_rows)]
        self._one = None
    
    def is_zero(self, proof=None):
        """Returns True iff matrix is a zero matrix"""
        if self.find() == (-1,-1):
            return True
        return False
    
    def one(self):
        """Returns matrix element's `1`"""
        if self._one is not None:
            return self._one
        
        row, col = self.find()
        if (row, col) == (-1, -1):
            raise MatrixError("Cannot determine `one` on zero matrix")
        x = self._grid[row][col]
        self._one = x / x
        return self._one
    
    def __str__(self):
        if self.ncols() == 0 or self.nrows() == 0:
            return ""
        maxlen = max(len(repr(e)) for e in self)
        return '\n'.join(
            ' '.join(repr(e).rjust(maxlen) for e in row) for row in self.rows()
        )

    def __repr__(self):
        return '<{} {}x{}>'.\
        format(self.__class__.__name__, self._nrows, self._ncols) \
        + "\n" + str(self)

    def __deepcopy__(self, memo=None):
        if memo is None:
            memo = {}
        return self.map(lambda e: copy.deepcopy(e, memo))
        
    def __getitem__(self, key):
        if isinstance(key, tuple):
            _assert(len(key) == 2, "Invalid indexing %r" % (key,))
            i, j = key
            return self._grid[i][j]
        else:
            return self._grid[key]

    def __contains__(self, v):
        for e in self:
            if e == v:
                return True
        return False

    def __neg__(self):
        return self.map(lambda e: -copy.deepcopy(e))

    def __pos__(self):
        return self.copy()

    def __eq__(self, other):
        return self._grid == other._grid

    def __ne__(self, other):
        return not self.__eq__(other)

    def __add__(self, other):
        _assert(isinstance(other, Matrix), "Cannot add %s to a Matrix" % type(other))
        _assert((self._nrows, self._ncols) == (other._nrows, other._ncols), "Cannot add matrices of different sizes")
        return self._op_elementwise(add, other)

    def __sub__(self, other):
        _assert(isinstance(other, Matrix), "Cannot subtract %s to a Matrix" % type(other))
        _assert((self._nrows, self._ncols) == (other._nrows, other._ncols), "Cannot subtract matrices of different sizes")
        return self._op_elementwise(sub, other)

    def __mul__(self, other):
        if isinstance(other, Matrix):
            _assert(self.ncols() == other.nrows(), "Incompatible sizes for multiplication")
            m = [[inner_product(self.row(i), col) for col in other.cols()] for i in xrange(self.nrows())]
            return Matrix.from_list(m)
        else:
            return self.map(lambda e: e * other)

    def __rmul__(self, other):
        return self * other

    def __pow__(self, power):
        _assert(isinstance(power, int), "Invalid power %r" % power)
        _assert(self.is_square(), "Cannot calculate powers of non-square matrix")
        
        y = abs(power)
        if y == 1:
            return Matrix.identity(self.nrows(), one=self.one())
        
        m = self.copy()
        if power < 0:
            m = m.inv()
        for i in xrange(y - 1):
            m = m * self
        return m

    def _op_elementwise(self, f, other):
        """Returns a new matrix where [m_ij] = [f(self_ij, other_ij)]"""
        _assert((self.nrows(), self.ncols()) == (other.nrows(), other.ncols()), "Incompatible sizes for elementwise %s" % f.__name__)
        m = (
            map(lambda x : f(*x), izip(row1, row2))
                for row1, row2 in izip(self.rows(), other.rows())
        )
        return Matrix.from_list(m)
        
    def nrows(self):
        """Get number of rows"""
        return self._nrows
    
    def ncols(self):
        """Get number of columns"""
        return self._ncols

    def __iter__(self):
        for row in self.rows():
            for x in row:
                yield x
    
    def row(self, n):
        """Returns an iterator over the specified row"""
        row = self[n]
        for x in row:
            yield x

    def col(self, n):
        """Returns an iterator over the specified column"""
        for row in xrange(self._nrows):
            yield self[row][n]

    def rows(self):
        """Returns a row iterator for each row in the matrix"""
        for row in xrange(self._nrows):
            yield self.row(row)

    def cols(self):
        """Returns a column iterator for each column in the matrix"""
        for col in xrange(self._ncols):
            yield self.col(col)
    
    def elements(self):
        """Iterator returning the tuple (row, col, element)"""
        for row in xrange(self._nrows):
            for col in xrange(self._ncols):
                yield row, col, self[row][col]
    
    def find(self, val=None, start=None, transpose=False):
        """Find first occurnce of @x in the matrix. If x not supplied, stops at
        any non-zero element. Returns a tuple (row, col) or (-1, -1) if not found.
        If @start supplied as a (x,y) - search will be on rows >= x, cols >= y.
        if @transpose is set True, search is performed column by column instead of by rows.
        """
        if val is None:
            check = lambda e, val : bool(e)
        else:
            check = lambda e, val : (e == val)
        
        start_x, start_y = start if (start is not None) else (0,0)
        if transpose is False:
            row_col = product(xrange(start_x, self.nrows()), xrange(start_y, self.ncols()))
        else:
            row_col = ((r,c) for c,r in product(xrange(start_y, self.ncols()), xrange(start_x, self.nrows())))
        
        for r,c in row_col:
                if check(self[r][c], val):
                    return (r,c)
        return (-1, -1)

    def rowvec(self, n):
        """Returns the specified row as a new row vector"""
        return Matrix.from_list([self.row(n)])

    def colvec(self, n):
        """Returns the specified column as a new column vector"""
        return Matrix.from_list([self.col(n)])

    def append_row(self, new_row):
        """Append a new row to the matrix"""
        return self.set_row(self.nrows(), new_row)
    
    def append_col(self, new_col):
        """Append a new column to the matrix"""
        return self.set_col(self.ncols(), new_col)
    
    def delete_row(self, n):
        """Delete row number @n of the matrix"""
        try:
            self._grid.pop(n)
        except:
            raise MatrixError("Invalid row index %r" % n)
        self._nrows -= 1
        if self._nrows == 0:
            self._ncols = 0
    
    def delete_col(self, n):
        """Delete column number @n of the matrix"""
        try:
            for i in xrange(self.nrows()):
                self._grid[i].pop(n)
        except:
            raise MatrixError("Invalid row index %r" % n)
        self._ncols -= 1
        if self._ncols == 0:
            self._nrows = 0
        
    def set_row(self, n, new_row):
        """Append a new / override an existing row"""
        new_row = list(new_row)
        _assert(len(new_row) == self.ncols(), "Invalid size for new row %d != %d" % (len(new_row), self.ncols()))
        _assert(0 <= n <= self.nrows(), "Invalid row number %d!" % n)
        row = [copy.deepcopy(x) for x in new_row]
        if n == self.nrows():
            self._nrows += 1
            self._grid.append(row)
        else:
            self._grid[n] = row
    
    def set_col(self, n, new_col):
        """Append a new / override an existing column"""
        new_col = list(new_col)
        _assert(len(new_col) == self.nrows(), "Invalid size for new col %d != %d" % (len(new_col), self.nrows()))
        _assert(0 <= n <= self.ncols(), "Invalid column number %d!" % n)
        
        if n == self.ncols():
            self._ncols += 1
            for i, e in enumerate(new_col):
                self._grid[i].append(copy.deepcopy(e))
        else:
            for i, e in enumerate(new_col):
                self._grid[i][n] = copy.deepcopy(e)
    
    def rowop_multiply(self, row, m):
        """In-place row operation. Multiplies the row @row by @m"""
        self._grid[row] = [m * x for x in self.row(row)]

    def rowop_swap(self, r1, r2):
        """In-place row operation. Interchanges the two specified rows"""
        if r1 != r2:
            self._grid[r1], self._grid[r2] = self._grid[r2], self._grid[r1]

    def rowop_add(self, r1, m, r2):
        """In-place row operation. Adds @m times row @r2 to row @r1"""
        if m:
            self._grid[r1] = [x + (m * y) for x,y in izip(self.row(r1), self.row(r2))]

    def copy(self):
        """Returns a copy of the matrix"""
        return self.__deepcopy__()

    def transpose(self):
        """Returns the transpose of the matrix as a new object"""
        m = Matrix(self._ncols, self._nrows)
        for row, col, element in self.elements():
            m[col][row] = element
        return m

    def t(self):
        """Shorthand for transpose method"""
        return self.transpose()

    def map(self, func):
        """Forms a new matrix by applying `func` to each element"""
        m = Matrix(self._nrows, self._ncols)
        for row, col, e in self.elements():
            m[row][col] = func(e)
        return m

    def is_square(self):
        """True iff the matrix is square"""
        return self._nrows == self._ncols
    
    def ref(self, mirror=None):
        """Transform to the row echelon form of the matrix using the forward phase of the 
        Gaussian elimination algorithm. Leading coefficients are normalized to `1`.
        If a @mirror matrix is supplied, the same sequence of row operations will be applied to it.
        Note that both of matrices are altered in-place!
        """
        # Start with the top row and work downwards.
        min_col = 0
        for top_row in xrange(self.nrows()):

            # Find a non-zero element in the leftmost column that is not all zeros
            # If such a column doesn't exist, were done
            row, col = self.find(start=(top_row, min_col), transpose=True)
            if (row, col) == (-1,-1):
                break
            #else:
            min_col = col

            # Bring the non-zero entry to the top of this column
            self.rowop_swap(top_row, row)
            if mirror:
                mirror.rowop_swap(top_row, row)

            # Normalize this entry '1'
            multiplier = self[top_row][col] ** (-1)
            self.rowop_multiply(top_row, multiplier)
            if mirror:
                mirror.rowop_multiply(top_row, multiplier)

            # Make all entries below this leading '1' zero
            for row in xrange(top_row + 1, self.nrows()):
                multiplier = -self[row][col]
                self.rowop_add(row, multiplier, top_row)
                if mirror:
                    mirror.rowop_add(row, multiplier, top_row)
    
    def rref(self, mirror=None):
        """Transform to the reduced row echelon form of the matrix using the 
        Gaussian elimination algorithm. If a @mirror matrix is supplied, the 
        same sequence of row operations will be applied to it. Note that both 
        matrices are altered in-place!
        """
        # Run the forward phase of the algorithm to determine the row echelon form.
        self.ref(mirror)

        # The backward phase of the algorithm. For each row, starting at the bottom
        # and working up, find the column containing the leading 1 and zero all the
        # entries above it
        for last_row in xrange(self.nrows() - 1, 0, -1):
            for col in xrange(self.ncols()):
                if self[last_row][col]: # Found leading coefficient?
                    for row in xrange(last_row):
                        multiplier = -self[row][col]
                        self.rowop_add(row, multiplier, last_row)
                        if mirror:
                            mirror.rowop_add(row, multiplier, last_row)
                    break # Move to next row
    
    def rank(self):
        """Returns the rank of the matrix.
        Counts the amount of zero lines in the row echoelon form of the matrix.
        """
        if self.is_zero():
            return 0
        
        # Find first zero row of the row echoelon form
        ref = self.copy()
        ref.ref()
        cur_row = 0
        for i in xrange(self.ncols()):
            if ref[cur_row][i]:
                cur_row += 1
        return cur_row
    
    def inv(self):
        """Returns the inverse matrix if it exists or raise MatrixError"""
        _assert(self.is_square(), "Non-square matrix doesn't have an inverse")
        orig = self.copy()
        inverse = Matrix.identity(self.nrows(), self.one())
        orig.rref(inverse)
        if orig != Matrix.identity(self.nrows(), self.one()):
            raise MatrixError("Matrix is not invertible")
        return inverse
    
    def __serialize__(self):
        grid = [[serialize(x) for x in row] for row in self.rows()]
        return self._nrows, self._ncols, grid
    
    @classmethod
    def __unserialize__(cls, serialized):
        nrows, ncols, grid_s = serialized
        self = cls()
        self._nrows, self._ncols = nrows, ncols
        self._grid = [[unserialize(x) for x in row] for row in grid_s]
        return self

    @classmethod
    def from_list(cls, lst):
        """Instantiates a new matrix object from a list of iterables / constants.
        On the latter case - a diagonal matrix is created.
        """
        l = list(lst)
        if l == []:
            return cls(0,0)
        
        if hasattr(l[0], '__iter__'):
            l = map(list, l)
            nrows, ncols = len(l), len(l[0])
            m = cls(nrows, ncols)
            for i in xrange(nrows):
                _assert(len(l[i]) == ncols, "Invalid grid for from_list")
                for j, e in enumerate(l[i]):
                    m[i][j] = copy.deepcopy(e)
            return m
        else:
            return cls.diag(l)
    
    @classmethod
    def diag(cls, l):
        """Instantiates a new len(l) x len(l) diagonal matrix with the values of l"""
        n = len(l)
        if n == 0:
            return cls(0,0)
        
        m = cls(0, n)
        zero = l[0] - l[0]
        zerow = [zero] * n
        for i in xrange(n):
            zerow[i] = l[i]
            m.append_row(zerow)
            zerow[i] = zero
        return m

    @classmethod
    def identity(cls, n, one=1):
        """Instantiates a new n x n identity matrix"""
        return cls.diag([one] * n)

    @classmethod
    def matrix(cls, *args):
        """Convenience function for instantiating Matrix objects"""
        if len(args) == 0:
            return cls(0,0)
        elif isinstance(args[0], int):
            return cls.identity(args[0])
        if isinstance(args[0], (list, tuple)):
            return cls.from_list(*args)
        else:
            raise NotImplementedError


class VanderMatrix(Matrix):
    """ A Vandermonde matrix.
    Each row is 1, x, x^2, ..., x^n - where n is a given value (or default - square matrix).
    Supports quick inverse calculation.
    """
    def __init__(self, xs, ncols=None):
        """
        @xs - the values defining each row in the matrix. Can be of any type supporting arithmetics
        @ncols - number of columns of the matrix. If None - a square matrix generated (set to len(xs))
        """
        if ncols is None:
            ncols = len(xs)
        _assert(len(xs) > 0 and ncols > 0, "VanderMatrix cannot be empty")
        
        Matrix.__init__(self, 0, ncols)
        self._xs = [copy.deepcopy(x) for x in xs]
        for x in xs:
            self.append_row([x**i for i in xrange(ncols)])
        self._one = copy.deepcopy(self[0,0])
    
    def inv(self):
        """Returns the inverse VanderMatrix.
        Quick calculation according to:
        https://arxiv.org/pdf/1211.1566.pdf
        """
        _assert(self.nrows() == self.ncols(), "Cannot invert a non-square matrix")
        
        n = self.nrows()
        one = self._one
        zero = one - one
        
        # Calc U matrix:
        #   every row has i zeros, then a growing product of differences with xs[i]
        U = Matrix(0, n)
        for i,x in enumerate(self._xs):
            # New row, initialize low part of triangle
            row = [zero] * i
            
            # Next, every non-zero element is a longer product. Build iteratively:
            diffs = ((x - self._xs[j]) for j in xrange(i))
            cur = reduce(mul, diffs, one).inv()
            row.append(cur)
            for j in xrange(i + 1, n):
                cur *= (x - self._xs[j]).inv()
                row.append(cur)
            
            # Append new row
            U.append_row(row)
        
        
        # Calculate L matrix
        L = Matrix(0, n)
        for i in xrange(n):
            row = []
            for j in xrange(i):
                cur = prev[j-1] - (self._xs[i-1] * prev[j])
                row.append(cur)
            row.append(one)
            row.extend([zero] * (n - i - 1))
            L.append_row(row)
            prev = row
        
        return (U * L).t()


def _assert(cond, msg):
    if not cond:
        raise MatrixError(msg)
