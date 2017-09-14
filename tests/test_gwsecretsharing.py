import unittest
from secretsharing import GWSSS
from gw_params import GWParams
from gf import FFE, GF, Extension
from polynomial import Polynomial
from utils import powers
from utils import random

class TestGWSSS(unittest.TestCase):

    def test_manual(self):
        test_cases = ( (7,1,3), (97,1,2), (3,3,6), (2,4,8) )
        max_n = 20
        rand_checks = 2
        for p, m1, m2 in test_cases:
            ext = m2 // m1
            f = GF(p, m1)
            gf = Extension(f, ext)
            n = random.randint(2, min(max_n, len(gf) - 2))
            t = (random.randint(0, n-2) + random.randint(0, n-2)) / 2 #We want (0, n-2) but more in the middle
            g = gf.generator()
            pts = list(powers(g, n))
            pos = g.inv()
            dual, mus, mus_basis, _ = GWParams(gf, pts, t).generate(pos, tries_limit=30)

            # Check that f(pos) = sum(  tr(z_i*f(pt_i)) * dual_i )
            # For a few random f's of degree t
            for test_n in xrange(rand_checks):
                pol = Polynomial.rand(gf, t)
                evals = pol.evaluate(pts)
                eval_pos = pol(pos)
                res = gf.zero()
                for i in xrange(gf.get_extension()):
                    a = sum((gf.trace(mu * evals[pid]) for (pid, mu) in mus[i]), f.zero())
                    res += gf.extend(a) * dual[i]
                self.assertEqual(res, eval_pos, "%r: %r != %r" % ((p,m1,m2,n,t), res, eval_pos))
    
    def test_GWSSSS(self):
        test_cases = ( (7,1,3), (97,1,2), (3,3,6), (2,4,8) )
        max_n = 20
        zero_pos_checks = 1
        rand_checks = 1
        lost_node_checks = 1
        precalculate_params = True
        for p, m1, m2 in test_cases:
            ext = m2 // m1
            f = GF(p, m1)
            gf = Extension(f, ext)
            n = random.randint(3, min(max_n, len(gf) - rand_checks - 1))
            t = (random.randint(1, n-2) + random.randint(1, n-2)) / 2 #We want (1, n-2) but more in the middle
            inds = range(len(gf)-1)
            random.shuffle(inds)
            g = gf.generator()
            pts = g.multi_pows(inds[:n])
            positions = []
            if zero_pos_checks:
                positions += [gf.zero()] * zero_pos_checks
            if rand_checks:
                positions += g.multi_pows(inds[-rand_checks:])
            if lost_node_checks:
                pts_cpy = [p for p in pts]
                random.shuffle(pts_cpy)
                positions += pts_cpy[:lost_node_checks]
            param_generator = GWParams(gf, pts, t)
            params = None
            if precalculate_params:
                params = {pos : param_generator.generate(pos, tries_limit=30) for pos in set(positions)}
            ss = GWSSS(n, t, gf, pts, gf.zero(), params)
            for pos in positions:
                secret = gf.rand()
                shares = ss.share(secret, pos=pos)
                traces = [ss.preprocess(s, i, pos) for i, s in enumerate(shares)]
                res = ss.reconstruct(traces, pos)
                self.assertEqual(res, secret, "%r: %r != %r\n%r" % ((p,m1,m2,n,t), res, secret, pos))
    

if __name__ == '__main__':
    unittest.main()
