"""
An example run of the whole program.
Made by Ori Chen
"""

from master import Master
from gf import Zmod, Extension
from circuit import AdditionGate, MultiplicationGate, ScalarMultGate, NullGate, Circuit
from transformation import Truncinator, Resampler
from gw_params import GWParams
from config import Config
from timeit import default_timer as timer
from utils import log, powers

def main():
    global init_time, run_time
    init_time = timer()
    log.info("Main start")
    
    # Finite field of some arbitrary prime power
    p = 97                # Arbitrary prime
    F = Zmod(p)           # Finite field of this prime
    F = Extension(F, 2)   # Extension field of this field
    #F = Extension(F, 2)   # Extension field of this field
    
    # Generator of the field
    g = F.generator()
    log.info("Generator generated")
    
    example = 2
    if example == 1:
        # Simple example
        c, circuit_emulator = GF_circuit(F)
        n = 2
    elif example == 2:
        c, circuit_emulator = simple_circuit(F)
        n = 6
    elif example == 3:
        # Benchmark example - lots of multiplication gates, arbitrary number of parties
        n = 30
        gates = 50
        c, circuit_emulator = mul_circuit(n, gates, F)
    else:
        return (None, None)
    
    log.info("Circuit generated")

    # Define party x-axis values for the secret sharing schemes
    pts = g.multi_pows(range(1,n+1))
    
    # General secret sharing parameters
    # Secret position must be zero for BGW!
    secret_pos = F.zero()
    threshold, threshold2 = Config.THRESHOLD(n), Config.THRESHOLD2(n)
    
    ss2_args = (F, pts, secret_pos)
    
    # Generate parameters for regular / GW15 reconstruction of output
    if Config.USE_GW_SSS1:
        params_factory1 = GWParams(F, pts, threshold)
        recon_params1 = {secret_pos : params_factory1.generate(secret_pos, tries_limit=Config.GW_PARAM_GEN_LIMIT)}
        ss_args = (F, pts, secret_pos, recon_params1)
        log.info("Set-up: GW15 reconstruction parameters generated (low deg)")
    else:
        ss_args = (F, pts, secret_pos)
    
    if Config.USE_GW_SSS2:
        params_factory2 = GWParams(F, pts, threshold2)
        recon_params2 = {secret_pos : params_factory2.generate(secret_pos, tries_limit=Config.GW_PARAM_GEN_LIMIT)}
        ss2_args = (F, pts, secret_pos, recon_params2)
        log.info("Set-up: GW15 reconstruction parameters generated (high deg)")
    else:
        ss2_args = (F, pts, secret_pos)

    # Define and calculate Truncinator for BGW multiplication if needed
    if Config.MULT_METHOD == Config.MULT_METHODS.BGW:
        mult_transformation = Truncinator(pts, threshold2, threshold, F)
        log.info("Set-up: Truncinator generated")
    
    # Define and calculate Resampler for DIK multiplication if needed
    if Config.MULT_METHOD == Config.MULT_METHODS.DIK:
        n_honest = n - threshold
        assert n + n_honest < len(F), "Field not large enough for resampler"
        dsts = g.multi_pows(-i for i in xrange(n_honest))
        mult_transformation = Resampler(pts, dsts, F)
        log.info("Set-up: Resampler generated")
    
    # Create Master with all the above
    m = Master(n, c, mult_transformation, ss_args=ss_args, ss2_args=ss2_args)

    log.info("Master & Node init done")
    
    # Initialize all Nodes: set circuit and truncinator & tell them to choose a random input and share their inputs
    m.init(F)
    log.info("Preprocessing done")
    
    # Main MPC - evaluate the whole circuit + output reconstruction
    m.run()
    log.info("MPC done")
    
    return circuit_emulator, F

def mul_circuit(n, gates, F):
    """Returns a circuit and a function for calculating:
    (input_0 * input_1 * ... * input_n) * (input_0 * input_1 * ... * input_n) * ...
    The result is again multiplied by input_0.
    Total of multiplications is set by @gates
    """
    c = Circuit(n)
    prev_layer = ["INPUT" + str(i) for i in xrange(n)]
    last_gate = prev_layer[0]
    layers, extra = divmod(gates, n)
    for k in xrange(layers + 1):
        cur = n
        if k == layers:
            cur = extra
        for i in xrange(cur):
            gname = "%d_%d" % (k,i)
            gt = MultiplicationGate(gname)
            c.add_gate(gt, [last_gate, "INPUT" + str(i)])
            last_gate = gname
        prev_layer = ["%d_%d" % (k,j) for j in xrange(n)]
    c.add_gate(NullGate(label="out"), [last_gate], True)

    def general_emulator(gates, n, *inputs):
        res = inputs[0]
        for i in xrange(gates):
            res = res * inputs[i % n]
        return res
    circuit_emulator = lambda *inputs : general_emulator(gates, n, *inputs)
    
    return c, circuit_emulator

def simple_circuit(F):
    """An exmaple circuit for 6 parties over @F:
    C(x,y,z,a,b,c) = c*((2*x + y*a + z) + 3*(2*x * b*z))
    Where `2` is 1+1 and `3` is 1+1+1 in the field.
    """
    c = Circuit(6)
    v1 = F.one()
    v2 = v1 + v1
    v3 = v2 + v1
    fg = AdditionGate("F")
    gg = AdditionGate("G")
    hg = ScalarMultGate(v2, "H")
    jg = MultiplicationGate("J")
    kg = ScalarMultGate(v3, "K")
    lg = MultiplicationGate("L")
    mg = MultiplicationGate("M")
    ng = MultiplicationGate("N")
    c.add_gate(mg, ["INPUT2", "INPUT4"])
    c.add_gate(hg, ["INPUT0"])
    c.add_gate(lg, ["INPUT1", "INPUT3"])
    c.add_gate(jg, ["M", "H"])
    c.add_gate(kg, ["J"])
    c.add_gate(fg, ["H","L", "INPUT2"])
    c.add_gate(gg, ["F", "K"])
    c.add_gate(ng, ["G", "INPUT5"], True)
    
    circuit_emulator = lambda x,y,z,a,b,c : (c*((v2*x + y*a + z) + v3*(v2*x * b*z)))
    
    return c, circuit_emulator

def GF_circuit(F):
    """An exmaple circuit for 2 parties over @F:
    C(x,y) = (x + y) + 2*x
    Where `2` is 1+1 in the field.
    """
    c = Circuit(2)
    two = F.one() + F.one()
    fg = AdditionGate("F")
    gg = AdditionGate("G")
    hg = ScalarMultGate(two, "H")
    c.add_gate(hg, ["INPUT0"])
    c.add_gate(fg, ["INPUT0", "INPUT1"])
    c.add_gate(gg, ["F", "H"], True)
        
    circuit_emulator = lambda x,y : (x + y) + (two*x)
    return c, circuit_emulator

if __name__ == '__main__':
    start = timer()
    circuit_emulator, F = main()
    print "Total %.2f seconds" % (timer() - start)
    
"""
For testing the result for an Extension, calculate c(ext, s) where @ext is the extension size and s is a list of the coeffcients of inputs:
from polynomial import *
f = F.get_subfield()
tuplize = lambda s, n : [map(int, s[i:i+n]) for i in xrange(0, len(s), n)]
c = lambda n, s : circuit_emulator(*map(lambda x : F(Polynomial(f, map(f, x))), tuplize(s,n)))

For Zmod:
circuit_emulator(*map(F,[...]))
"""