import unittest
from transformation import Truncinator, Resampler
from polynomial import Polynomial
from gf import Zmod
from utils import random

class TestTruncinator(unittest.TestCase):
    
    def test_rand_truncinator(self):
        number_of_tests = 5
        for i in xrange(number_of_tests):
            d = random.randint(2, 10)
            d2 = random.randint(d + 1 , min(3 * d, 120))
            n = random.randint(d2, d2 + 3)
            p = random.choice([257, 443, 701, 997])
            F = Zmod(p)
            
            one_to_p = range(1, p)
            random.shuffle(one_to_p)
            pts = map(F, one_to_p[:n])
            
            coef2 = [F.rand() for i in xrange(d2)]
            coef1 = coef2[:d]
            
            # Evaluate the polynomial at all points
            vals1 = []
            vals2 = []
            for pt in pts:
                v1 = sum(((c * (pt**i)) for i,c in enumerate(coef1)), F.zero())
                v2 = sum(((c * (pt**i)) for i,c in enumerate(coef2)), F.zero())
                vals1.append(v1)
                vals2.append(v2)
            
            t = Truncinator(pts, d2, d, F)
            myvals = t.reduce(vals2)
            
            self.assertEqual(myvals, vals1, \
                             "n=%d, d=%d, p=%d, coef:\n%r \npts:\n%r" % (n, d, p, coef2, pts) \
                            )

    def test_rand_resampler(self):
        number_of_tests = 5
        for i in xrange(number_of_tests):
            p = random.choice([257, 443, 701, 997])
            F = Zmod(p)
            d = random.randint(2, min(80, p-1))
            n_in = random.randint(d + 1 , min(3 * d, 120))
            n_out = random.randint(1, min(p - n_in, 120))
            
            zero_to_p = range(p)
            random.shuffle(zero_to_p)
            srcs = map(F, zero_to_p[:n_in])
            dsts = map(F, zero_to_p[-n_out:])
            
            poly = Polynomial.rand(F, d)
            
            r = Resampler(srcs, dsts, F)
            myvals = r.apply(map(poly, srcs))
            realvals = map(poly, dsts)
            
            self.assertEqual(myvals, realvals, \
                             "n_in=%d, n_out=%d, d=%d, p=%d, poly:\n%r \npts:\n%r\n%r" % \
                                (n_in, n_out, d, p, poly, srcs, dsts) \
                            )

if __name__ == '__main__':
    unittest.main()
