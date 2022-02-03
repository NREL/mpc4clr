"""
Distribution grid load tester.
"""

import pytest
import numpy as np
import os
import math

import julia
from julia.api import Julia

jl = Julia(compiled_modules=False)  
jl = julia.Julia()

current_dir = os.getcwd()
os.chdir(os.path.dirname(current_dir) + "/mpc_three_phase_balanced_system")
    
jl.include("data_handler_threephase.jl")
jl.include("lr_threephase_opf_mpc.jl")

NETWORK = jl.load_case_data(datafile="13buscase")

abs_tol = 1e-05


class TestLoad:

    def test_num_load(self):

        bus = NETWORK[0]

        num_load = 0

        for idx in range(1, len(bus)+1):
            pl = bus[idx].d_P
            ql = bus[idx].d_Q
            if pl != 0.0 or ql != 0.0:
                num_load += 1

        assert num_load == 9

    def test_simple_load(self):

        bus = NETWORK[0]

        spec = np.array([5, 0.002874791974852071, 0.0016757581360946747, 0.8639363184095202], dtype=object)

        assert bus[5].index == spec[0]
        assert math.isclose(bus[5].d_P, spec[1], abs_tol=abs_tol)
        assert math.isclose(bus[5].d_Q, spec[2], abs_tol=abs_tol)
        assert math.isclose(bus[5].cosphi, spec[3], abs_tol=abs_tol)

          
# pytest test_load.py