import unittest
from gf import GF, Extension, Zmod, FFE
from utils import powers, factor
from utils import random

class TestGF(unittest.TestCase):

    def test_generator(self):
        mul_ints = lambda ints : reduce(int.__mul__, ints, 1)
        test_cases = ( (3,1,6), (5,1,3), (5,1,4), (17,1,2), (97,1,2), (3,3,6), (3,2,6), (5,2,4) )
        for p, t, n in test_cases:
            K = GF(p, t)
            L = Extension(K, n // t)
            g = L.generator()
            self.assertEqual(len(set(powers(g, len(L)-1))), len(L)-1)
    
    def test_all_generators(self):
        mul_ints = lambda ints : reduce(int.__mul__, ints, 1)
        test_cases = ( (3,1,6), (17,1,2), (3,3,6), (5,2,4) )
        rand_checks = 5
        for p, t, n in test_cases:
            if t > 1:
                K = GF(p, t)
            else:
                K = Zmod(p)
            L = Extension(K, n // t)
            all_gens = list(L.all_generators())
            factors = factor(len(L) - 1)
            self.assertEqual(len(all_gens), (len(L) / mul_ints(factors)) * mul_ints(x-1 for x in factors))
            
            for i in xrange(rand_checks):
                rand_g = random.choice(all_gens)
                self.assertEqual(len(set(powers(rand_g, len(L)-1))), len(L)-1)
    
    def test_trace_pows(self):
        test_cases = ( (3,1,6), (5,1,3), (5,1,4), (17,1,2), (97,1,2), (3,3,6), (3,2,6), (5,2,4) )
        n_tests = 10
        for p, t, n in test_cases:
            K = GF(p, t)
            L = Extension(K, n // t)
            q = len(K)
            for i in xrange(n_tests):
                x = L.rand()
                self.assertEqual(L.trace(x), L.trace(x**q))
    
    def test_trace_counts(self):
        test_cases = ( (3,1,6), (3,3,6), (3,2,6), (5,1,3), (5,2,4), (17,1,2) )
        for p, t, n in test_cases:
            if t > 1:
                K = GF(p, t)
            else:
                K = Zmod(p)
            L = Extension(K, n // t)
            kgen, lgen = K.generator(), L.generator()
            traces = [L.trace(x) for x in powers(lgen, len(L)-1)]
            counts = [traces.count(x) for x in powers(kgen, len(L)-1)]
            self.assertEqual(counts[0], p**(n - t))
            self.assertEqual(len(set(counts)), 1)
    
    def test_trace_linearity(self):
        test_cases = ( (3,6), (5,3), (5,4), (17,2) )
        n_tests = 10
        for p, n in test_cases:
            f = Zmod(p)
            L = Extension(f, n)
            tr = L.trace
            for i in xrange(n_tests):
                x, y = L.rand(), L.rand()
                a, b = f.rand(), f.rand()
                a_l, b_l = map(L.extend, (a,b))
                self.assertEqual(tr(x + y), tr(x) + tr(y))
                self.assertEqual(tr(a_l*x + b_l*y), a*tr(x) + b*tr(y))

if __name__ == '__main__':
    unittest.main()
