import unittest
from secretsharing import NullSSS, ShamirSSS
from gf import Zmod, Extension, GF
from utils import random

class TestSSS(unittest.TestCase):

    def test_nullSSS_sharing(self):
        n = 10
        t = 3
        ss = NullSSS(n, t)
        secret = Zmod(11).rand()
        shares = ss.share(secret)
        self.assertEqual(len(shares), n)
        self.assertEqual(shares[0], secret)
    
    def test_nullSSS_reconstruction(self):
        n = 101
        t = 3
        ss = NullSSS(n, t)
        secret = Zmod(97).rand()
        shares = ss.share(secret)
        recon = ss.reconstruct(ss.preprocess(shares))
        self.assertEqual(secret, recon)
        
    def test_shamirSSS_Zmod(self):
        number_of_tests = 10
        for i in xrange(number_of_tests):
            n = random.randint(2, 256)
            t = random.randint(0, n-1)
            p = random.choice([257, 443, 701, 997])
        
            F = Zmod(p)
            secret = F.rand()
            
            one_to_p = range(1, p)
            random.shuffle(one_to_p)
            xs = map(F, one_to_p[:n])
            pos = F(one_to_p[-1])
            
            s = ShamirSSS(n, t, F, xs)
            shares = s.share(secret, pos=pos)
            res = s.reconstruct(s.preprocess(shares), pos=pos)

            points = zip(map(repr, xs), map(repr, shares))
            self.assertEqual(res, secret, "n=%d, t=%d, p=%d, s=%r, pos=%r, pts:\n%r" % (n, t, p, secret, pos, points))

    def test_shamirSSS_Extension(self):
        number_of_tests = 10
        for i in xrange(number_of_tests):
            p, m1, m2 = random.choice([(7,1,3), (97,1,2), (3,3,6), (2,4,8)])
            ext = m2 // m1
            F = Extension(GF(p, m1), ext)
            g = F.generator()
            
            n = random.randint(2, min(len(F)-2, 50))
            t = random.randint(0, n-1)
        
            ind = range(len(F)-1)
            random.shuffle(ind)
            xs = g.multi_pows(ind[:n])
            pos = g**ind[-1]
            
            secret = F.rand()
            ss = ShamirSSS(n, t, F, xs)
            shares = ss.share(secret, pos=pos)
            res = ss.reconstruct(ss.preprocess(shares), pos=pos)

            points = zip(map(repr, xs), map(repr, shares))
            self.assertEqual(res, secret, "n=%d, t=%d, p=%d, s=%r, pos=%r, pts:\n%r" % (n, t, p, secret, pos, points))
            
if __name__ == '__main__':
    unittest.main()
