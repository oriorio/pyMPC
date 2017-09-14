import os, sys
import logging as log
import traceback
import operator
import copy
import random
from itertools import repeat, izip

log.basicConfig(level=log.INFO)

def enum(*sequential, **named):
    """Create an enum. Can be made out of a list of strings or a dictionary mapping to values.
    Access to constants under enum E are done by E.VAL.
    Reverse mapping (value -> constant name) is made by E.reverse_mapping[val]
    """
    enums = dict(zip(sequential, range(len(sequential))), **named)
    reverse = dict((value, key) for key, value in enums.iteritems())
    enums['reverse_mapping'] = reverse
    return type('Enum', (), enums)

def factor(n):
    """Return a sorted list of n's prime factors"""
    res = []
    d = 2
    while d*d <= n:
        if (n % d) == 0:
            res.append(d)
            n //= d
            while (n % d) == 0:
                n //= d
        d = (d + 1) | 1 #Next odd number
    if n > 1:
       res.append(n)
    return res

def powers(v, n=None):
    """Returns a iterator of @v**0 .. @v**(@n-1) or endlessly if @n not supplied"""
    if n is None:
        it = repeat(None)
    else:
        it = xrange(1,n)
    x = v**0
    yield x
    for _ in it:
        x *= v
        yield x

def inner_product(vec1, vec2):
    """Performs inner product of two iterables (stops after shortest iterable is exhausted)"""
    return reduce(operator.add, (x*y for x,y in izip(vec1, vec2)))
