"""
Module to hold all configuration of the program and MPC
"""

from utils import enum

class Config(object):

    # Are we enabling debugging (breakpoints & interactive shell)
    DEBUGGING = False
    DEBUGGING = DEBUGGING and __debug__ # DO NOT ALTER LINE - turn off debugging if asserts are off
    
    MULT_METHODS = enum("DIK", "BGW")
    MULT_METHOD = MULT_METHODS.DIK
    
    # TODO Future feature - party goes down
    PARTY_RESTORE = False

    # MPC secret sharing threshold for n parties, and for secret sharing during multiplication
    THRESHOLD = staticmethod( lambda n : ((n-1) / 2) - int(Config.PARTY_RESTORE) )
    THRESHOLD2 = staticmethod( lambda n : Config.THRESHOLD(n) * 2 )
    
    # GW15 reconstruction scheme ?
    USE_GW_SSS1 = True
    USE_GW_SSS2 = False
    # When generating parameters for GW15 reconstruction - size of quality check subgroup
    # First value is the fraction size (5.0 -> 20%, 100.0 -> 1% etc')
    # second is minimum number of parties in such a subgroup
    GW_PARAMS_SUBGROUP_SIZE = (5.0, 3)
    
    # Amount of schemes to test before choosing in main.py
    GW_PARAM_GEN_LIMIT = 60

    # Serialization compression parameters
    COMPRESS_SERIALIZATION = True
    COMPRESS_LIM = 4000
    
