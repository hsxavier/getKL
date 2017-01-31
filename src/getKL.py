#!/usr/bin/env python 

import healpy as hp
from sympy.physics.wigner import gaunt
import numpy as np
import sys

def ThreeYlmIntegral(l1, m1, l2, m2, l3, m3):
    if (m3%2==0):
        fac1 = 1
    else:
        fac1 =-1
    return fac1 * gaunt(l1,l2,l3,m1,m2,-m3).evalf()


# Get input:

l1 = int(sys.argv[1])
l2 = int(sys.argv[2])
l3 = int(sys.argv[3])
m1 = int(sys.argv[4])
m2 = int(sys.argv[5])
m3 = int(sys.argv[6])

print gaunt(l1,l2,l3,m1,m2,m3)
#print ThreeYlmIntegral(l1,m1,l2,m2,l3,m3)
